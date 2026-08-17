# =============================================================================
# 03_GSEA_metabolism_and_figures_v2.R
#
# Publication-oriented visualization of cross-cohort pre-ranked GSEA results.
#
# Design principles
#   - Each figure is exported separately.
#   - No descriptive plot titles or subtitles are placed inside the panels.
#   - Panel letters are retained for manuscript assembly.
#   - Pathway labels use sentence case.
#   - Adult-enriched = blue.
#   - Pediatric-enriched = yellow.
#   - Positive NES = Pediatric-enriched.
#   - Negative NES = Adult-enriched.
#
# Main outputs
#   A. Global cross-cohort NES heatmap
#   B. Metabolism-focused NES heatmap
#   C. Shared leading-edge genes
#   D. Representative pediatric-enriched GSEA curves
#   E. Representative adult-enriched GSEA curves
#
# Extended GSEA curves are exported separately for supplementary use.
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project.
# =============================================================================

SCRIPT_BUILD <- "GSEA_FIGURES_2026-08-17_v2"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

required_packages <- c(
  "fgsea",
  "msigdbr",
  "tidyverse",
  "here",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(tidyverse)
  library(here)
  library(patchwork)
})


# =============================================================================
# 01. Configuration
# =============================================================================

ALPHA <- 0.05

ADULT_COLOR <- "#3E76BC"
PEDIATRIC_COLOR <- "#FCCE24"

dataset_order <- c(
  "GSE172274",
  "GSE179277",
  "GSE231409"
)

gsea_long_file <- here::here(
  "results",
  "gsea",
  "GSEA_all_cohorts_long.csv"
)

concordance_file <- here::here(
  "results",
  "gsea",
  "cross_cohort",
  "tables",
  "GSEA_pathway_concordance_complete.csv"
)

leading_edge_file <- here::here(
  "results",
  "gsea",
  "cross_cohort",
  "tables",
  "GSEA_leading_edge_overlap.csv"
)

rank_files <- c(
  GSE172274 = here::here(
    "results",
    "deseq2",
    "GSE172274",
    "primary",
    "GSE172274_GSEA_ranked_genes.csv"
  ),
  GSE179277 = here::here(
    "results",
    "deseq2",
    "GSE179277",
    "primary",
    "GSE179277_GSEA_ranked_genes.csv"
  ),
  GSE231409 = here::here(
    "results",
    "deseq2",
    "GSE231409",
    "primary",
    "GSE231409_GSEA_ranked_genes.csv"
  )
)

output_dir <- here::here(
  "results",
  "gsea",
  "figures_final"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

curves_dir <- file.path(
  output_dir,
  "supplementary_curves"
)

diagnostics_dir <- file.path(
  output_dir,
  "diagnostics"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  curves_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 02. Validate inputs
# =============================================================================

required_files <- c(
  gsea_long_file,
  concordance_file,
  leading_edge_file,
  rank_files
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required files were not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# =============================================================================
# 03. Read results
# =============================================================================

gsea_long <- read.csv(
  gsea_long_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

concordance <- read.csv(
  concordance_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

leading_edge <- read.csv(
  leading_edge_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gsea_long <- gsea_long %>%
  mutate(
    Dataset = as.character(Dataset),
    Collection = as.character(Collection),
    pathway = as.character(pathway),
    NES = as.numeric(NES),
    pval = as.numeric(pval),
    padj = as.numeric(padj)
  )

concordance <- concordance %>%
  mutate(
    Collection = as.character(Collection),
    pathway = as.character(pathway),
    N_cohorts_tested = as.integer(N_cohorts_tested),
    N_FDR05 = as.integer(N_FDR05),
    Mean_NES = as.numeric(Mean_NES),
    DirectionConsistent = as.logical(DirectionConsistent)
  )

leading_edge <- leading_edge %>%
  mutate(
    Collection = as.character(Collection),
    pathway = as.character(pathway),
    N_LE_shared_all3 = as.integer(N_LE_shared_all3)
  )


# =============================================================================
# 04. Clean pathway labels
# =============================================================================

# Remove database prefixes, convert underscores to spaces, and use sentence case.
# This keeps pathway names visually clean and avoids all-uppercase labels.

concordance <- concordance %>%
  mutate(
    PathwayLabel = pathway,
    PathwayLabel = sub(
      "^HALLMARK_",
      "",
      PathwayLabel
    ),
    PathwayLabel = sub(
      "^GOBP_",
      "",
      PathwayLabel
    ),
    PathwayLabel = sub(
      "^REACTOME_",
      "",
      PathwayLabel
    ),
    PathwayLabel = sub(
      "^KEGG_",
      "",
      PathwayLabel
    ),
    PathwayLabel = gsub(
      "_",
      " ",
      PathwayLabel
    ),
    PathwayLabel = stringr::str_to_sentence(
      tolower(
        PathwayLabel
      )
    ),
    EnrichmentGroup = case_when(
      Mean_NES > 0 ~ "Pediatric-enriched",
      Mean_NES < 0 ~ "Adult-enriched",
      TRUE ~ "Neutral"
    )
  )


# =============================================================================
# 05. Global representative pathway selection
# =============================================================================

# Only pathways:
#   - tested in all 3 cohorts,
#   - direction-consistent across all 3 cohorts,
#   - significant in at least 2 of 3 cohorts
# are considered for the primary heatmap.
#
# Representative terms are selected to reduce redundancy.

global_theme_patterns <- c(
  "ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
  "HUMORAL_IMMUNE_RESPONSE$",
  "ORGANIC_ANION_TRANSPORT$",
  "BILE_ACID_METABOL",
  "FATTY_ACID_METABOL",
  "LIPID_CATABOLIC_PROCESS",
  "XENOBIOTIC_METABOL",
  "STEROID_METABOL",
  "ARACHIDONATE_METABOL",
  "CYTOCHROME_P450",
  "CHOLESTEROL",
  "OXIDATIVE_PHOSPHORYLATION",
  "GLYCOLYSIS",
  "PEROXISOME",
  "REACTIVE_OXYGEN_SPECIES",
  "INTERFERON_ALPHA_RESPONSE",
  "INTERFERON_GAMMA_RESPONSE",
  "TNF"
)

global_candidates <- concordance %>%
  filter(
    N_cohorts_tested == 3,
    DirectionConsistent,
    N_FDR05 >= 2
  )

selected_global <- data.frame()

for (pattern_i in global_theme_patterns) {

  candidate_i <- global_candidates %>%
    filter(
      grepl(
        pattern_i,
        pathway,
        ignore.case = TRUE
      )
    ) %>%
    arrange(
      desc(N_FDR05),
      desc(abs(Mean_NES))
    ) %>%
    slice_head(
      n = 1
    )

  if (nrow(candidate_i) > 0) {
    selected_global <- bind_rows(
      selected_global,
      candidate_i
    )
  }
}

selected_global <- selected_global %>%
  distinct(
    Collection,
    pathway,
    .keep_all = TRUE
  )

if (nrow(selected_global) < 15) {

  supplemental_global <- global_candidates %>%
    anti_join(
      selected_global %>%
        dplyr::select(
          Collection,
          pathway
        ),
      by = c(
        "Collection",
        "pathway"
      )
    ) %>%
    arrange(
      desc(N_FDR05),
      desc(abs(Mean_NES))
    ) %>%
    slice_head(
      n = 15 - nrow(selected_global)
    )

  selected_global <- bind_rows(
    selected_global,
    supplemental_global
  )
}

selected_global <- selected_global %>%
  arrange(
    EnrichmentGroup,
    Mean_NES
  )

write.csv(
  selected_global,
  file.path(
    tables_dir,
    "Selected_global_concordant_pathways.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 06. FIGURE A — Global cross-cohort NES heatmap
# =============================================================================

figure_A_data <- gsea_long %>%
  semi_join(
    selected_global %>%
      dplyr::select(
        Collection,
        pathway
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  left_join(
    selected_global %>%
      dplyr::select(
        Collection,
        pathway,
        PathwayLabel,
        EnrichmentGroup
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = dataset_order
    ),
    Significance = case_when(
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    EnrichmentGroup = factor(
      EnrichmentGroup,
      levels = c(
        "Adult-enriched",
        "Pediatric-enriched"
      )
    )
  )

figure_A_order <- selected_global %>%
  arrange(
    EnrichmentGroup,
    Mean_NES
  ) %>%
  pull(
    PathwayLabel
  ) %>%
  unique()

figure_A_data$PathwayLabel <- factor(
  figure_A_data$PathwayLabel,
  levels = figure_A_order
)

figure_A <- ggplot(
  figure_A_data,
  aes(
    x = Dataset,
    y = PathwayLabel,
    fill = NES
  )
) +
  geom_tile(
    linewidth = 0.45,
    colour = "white"
  ) +
  geom_text(
    aes(
      label = Significance
    ),
    size = 4.2,
    fontface = "bold"
  ) +
  facet_grid(
    EnrichmentGroup ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = ADULT_COLOR,
    mid = "white",
    high = PEDIATRIC_COLOR,
    midpoint = 0,
    name = "NES"
  ) +
  labs(
    tag = "A",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_A_GSEA_global_concordance_heatmap.tiff"
  ),
  plot = figure_A,
  width = 8.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_A_GSEA_global_concordance_heatmap.pdf"
  ),
  plot = figure_A,
  width = 8.5,
  height = 8.5,
  units = "in"
)


# =============================================================================
# 07. Metabolism-focused pathway selection
# =============================================================================

metabolic_regex <- paste(
  c(
    "METABOL",
    "OXIDATIVE_PHOSPHORYLATION",
    "GLYCOLYSIS",
    "PEROXISOME",
    "CHOLESTEROL",
    "FATTY_ACID",
    "BILE_ACID",
    "LIPID",
    "STEROID",
    "XENOBIOTIC",
    "ARACHIDONATE",
    "CYTOCHROME_P450",
    "TRICARBOXYLIC",
    "TCA_CYCLE",
    "AMINO_ACID",
    "NUCLEOTIDE"
  ),
  collapse = "|"
)

metabolic_candidates <- concordance %>%
  filter(
    N_cohorts_tested == 3,
    DirectionConsistent,
    N_FDR05 >= 2,
    grepl(
      metabolic_regex,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    desc(N_FDR05),
    desc(abs(Mean_NES))
  )

write.csv(
  metabolic_candidates,
  file.path(
    tables_dir,
    "All_concordant_metabolic_pathways.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 08. Reduce metabolic redundancy
# =============================================================================

metabolic_theme_patterns <- c(
  "BILE_ACID_METABOL",
  "FATTY_ACID_METABOL",
  "LIPID_CATABOLIC_PROCESS",
  "XENOBIOTIC_METABOL",
  "STEROID_METABOL",
  "ARACHIDONATE_METABOL",
  "CYTOCHROME_P450",
  "CHOLESTEROL",
  "OXIDATIVE_PHOSPHORYLATION",
  "GLYCOLYSIS",
  "PEROXISOME",
  "TRICARBOXYLIC|TCA_CYCLE",
  "AMINO_ACID_METABOL",
  "NUCLEOTIDE_METABOL"
)

selected_metabolic <- data.frame()

for (pattern_i in metabolic_theme_patterns) {

  candidate_i <- metabolic_candidates %>%
    filter(
      grepl(
        pattern_i,
        pathway,
        ignore.case = TRUE
      )
    ) %>%
    arrange(
      desc(N_FDR05),
      desc(abs(Mean_NES))
    ) %>%
    slice_head(
      n = 1
    )

  if (nrow(candidate_i) > 0) {
    selected_metabolic <- bind_rows(
      selected_metabolic,
      candidate_i
    )
  }
}

selected_metabolic <- selected_metabolic %>%
  distinct(
    Collection,
    pathway,
    .keep_all = TRUE
  ) %>%
  arrange(
    Mean_NES
  )

write.csv(
  selected_metabolic,
  file.path(
    tables_dir,
    "Selected_metabolic_pathways_for_figure.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 09. FIGURE B — Metabolism-focused NES heatmap
# =============================================================================

figure_B_data <- gsea_long %>%
  semi_join(
    selected_metabolic %>%
      dplyr::select(
        Collection,
        pathway
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  left_join(
    selected_metabolic %>%
      dplyr::select(
        Collection,
        pathway,
        PathwayLabel,
        EnrichmentGroup
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = dataset_order
    ),
    Significance = case_when(
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

figure_B_order <- selected_metabolic %>%
  arrange(
    Mean_NES
  ) %>%
  pull(
    PathwayLabel
  ) %>%
  unique()

figure_B_data$PathwayLabel <- factor(
  figure_B_data$PathwayLabel,
  levels = figure_B_order
)

figure_B <- ggplot(
  figure_B_data,
  aes(
    x = Dataset,
    y = PathwayLabel,
    fill = NES
  )
) +
  geom_tile(
    linewidth = 0.45,
    colour = "white"
  ) +
  geom_text(
    aes(
      label = Significance
    ),
    size = 4.2,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = ADULT_COLOR,
    mid = "white",
    high = PEDIATRIC_COLOR,
    midpoint = 0,
    name = "NES"
  ) +
  labs(
    tag = "B",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.ticks = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_B_GSEA_metabolism_heatmap.tiff"
  ),
  plot = figure_B,
  width = 8.5,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_B_GSEA_metabolism_heatmap.pdf"
  ),
  plot = figure_B,
  width = 8.5,
  height = 7,
  units = "in"
)


# =============================================================================
# 10. Select pathways for leading-edge analysis and GSEA curves
# =============================================================================

curve_theme_patterns <- c(
  "HUMORAL_IMMUNE_RESPONSE$",
  "ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
  "ORGANIC_ANION_TRANSPORT$",
  "BILE_ACID_METABOL",
  "FATTY_ACID_METABOL",
  "XENOBIOTIC_METABOL"
)

curve_pathways_table <- data.frame()

for (pattern_i in curve_theme_patterns) {

  candidate_i <- concordance %>%
    filter(
      N_cohorts_tested == 3,
      DirectionConsistent,
      N_FDR05 >= 2,
      grepl(
        pattern_i,
        pathway,
        ignore.case = TRUE
      )
    ) %>%
    arrange(
      desc(N_FDR05),
      desc(abs(Mean_NES))
    ) %>%
    slice_head(
      n = 1
    )

  if (nrow(candidate_i) > 0) {
    curve_pathways_table <- bind_rows(
      curve_pathways_table,
      candidate_i
    )
  }
}

curve_pathways_table <- curve_pathways_table %>%
  distinct(
    Collection,
    pathway,
    .keep_all = TRUE
  )

write.csv(
  curve_pathways_table,
  file.path(
    tables_dir,
    "Selected_pathways_for_GSEA_curves.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. FIGURE C — Shared leading-edge genes
# =============================================================================

leading_edge_selected <- leading_edge %>%
  semi_join(
    curve_pathways_table %>%
      dplyr::select(
        Collection,
        pathway
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  left_join(
    concordance %>%
      dplyr::select(
        Collection,
        pathway,
        PathwayLabel,
        Mean_NES,
        EnrichmentGroup
      ),
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  arrange(
    N_LE_shared_all3
  )

write.csv(
  leading_edge_selected,
  file.path(
    tables_dir,
    "Leading_edge_selected_pathways.csv"
  ),
  row.names = FALSE
)

leading_edge_selected$PathwayLabel <- factor(
  leading_edge_selected$PathwayLabel,
  levels = leading_edge_selected$PathwayLabel
)

figure_C <- ggplot(
  leading_edge_selected,
  aes(
    x = N_LE_shared_all3,
    y = PathwayLabel,
    colour = EnrichmentGroup
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = N_LE_shared_all3,
      y = PathwayLabel,
      yend = PathwayLabel
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 4
  ) +
  scale_colour_manual(
    values = c(
      "Adult-enriched" = ADULT_COLOR,
      "Pediatric-enriched" = PEDIATRIC_COLOR
    )
  ) +
  labs(
    tag = "C",
    x = "Leading-edge genes shared by all three cohorts",
    y = NULL,
    colour = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 9
    ),
    axis.text.y = element_text(
      size = 9
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_C_GSEA_shared_leading_edge.tiff"
  ),
  plot = figure_C,
  width = 8.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_C_GSEA_shared_leading_edge.pdf"
  ),
  plot = figure_C,
  width = 8.5,
  height = 5,
  units = "in"
)


# =============================================================================
# 12. Retrieve gene sets for enrichment curves
# =============================================================================

msig_human <- msigdbr::msigdbr(
  species = "Homo sapiens"
)

selected_curve_gene_sets <- msig_human %>%
  filter(
    gs_name %in% curve_pathways_table$pathway
  ) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()

pathway_gene_lists <- split(
  selected_curve_gene_sets$gene_symbol,
  selected_curve_gene_sets$gs_name
)


# =============================================================================
# 13. Read complete cohort rankings
# =============================================================================

rank_172274 <- read.csv(
  rank_files["GSE172274"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rank_179277 <- read.csv(
  rank_files["GSE179277"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rank_231409 <- read.csv(
  rank_files["GSE231409"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stats_172274 <- rank_172274$ranking_statistic
names(stats_172274) <- rank_172274$SYMBOL
stats_172274 <- sort(
  stats_172274,
  decreasing = TRUE
)

stats_179277 <- rank_179277$ranking_statistic
names(stats_179277) <- rank_179277$SYMBOL
stats_179277 <- sort(
  stats_179277,
  decreasing = TRUE
)

stats_231409 <- rank_231409$ranking_statistic
names(stats_231409) <- rank_231409$SYMBOL
stats_231409 <- sort(
  stats_231409,
  decreasing = TRUE
)


# =============================================================================
# 14. Select one representative Pediatric and one representative Adult pathway
# =============================================================================

representative_pediatric <- curve_pathways_table %>%
  filter(
    Mean_NES > 0
  ) %>%
  arrange(
    desc(N_FDR05),
    desc(abs(Mean_NES))
  ) %>%
  slice_head(
    n = 1
  )

representative_adult <- curve_pathways_table %>%
  filter(
    Mean_NES < 0
  ) %>%
  arrange(
    desc(N_FDR05),
    desc(abs(Mean_NES))
  ) %>%
  slice_head(
    n = 1
  )

write.csv(
  bind_rows(
    representative_pediatric,
    representative_adult
  ),
  file.path(
    tables_dir,
    "Representative_pathways_for_main_curves.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. FIGURE D — Representative Pediatric-enriched GSEA curves
# =============================================================================

if (nrow(representative_pediatric) == 1) {

  pathway_D <- representative_pediatric$pathway[1]
  genes_D <- pathway_gene_lists[[pathway_D]]

  results_D <- gsea_long %>%
    filter(
      pathway == pathway_D
    )

  pretty_D <- representative_pediatric$PathwayLabel[1]

  D_172274 <- results_D %>%
    filter(
      Dataset == "GSE172274"
    )

  D_179277 <- results_D %>%
    filter(
      Dataset == "GSE179277"
    )

  D_231409 <- results_D %>%
    filter(
      Dataset == "GSE231409"
    )

  plot_D1 <- fgsea::plotEnrichment(
    genes_D,
    stats_172274
  ) +
    labs(
      tag = "GSE172274",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(D_172274$NES[1], 2),
        "   FDR = ",
        format(
          D_172274$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = PEDIATRIC_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  plot_D2 <- fgsea::plotEnrichment(
    genes_D,
    stats_179277
  ) +
    labs(
      tag = "GSE179277",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(D_179277$NES[1], 2),
        "   FDR = ",
        format(
          D_179277$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = PEDIATRIC_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  plot_D3 <- fgsea::plotEnrichment(
    genes_D,
    stats_231409
  ) +
    labs(
      tag = "GSE231409",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(D_231409$NES[1], 2),
        "   FDR = ",
        format(
          D_231409$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = PEDIATRIC_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  figure_D <- (
    plot_D1 |
      plot_D2 |
      plot_D3
  ) +
    plot_annotation(
      tag_levels = NULL
    )

  figure_D <- figure_D +
    plot_annotation(
      theme = theme(
        plot.margin = margin(
          10,
          10,
          10,
          10
        )
      )
    )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Figure_D_GSEA_pediatric_",
        gsub(
          "[^A-Za-z0-9_.-]",
          "_",
          pathway_D
        ),
        ".tiff"
      )
    ),
    plot = figure_D,
    width = 10.5,
    height = 3.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Figure_D_GSEA_pediatric_",
        gsub(
          "[^A-Za-z0-9_.-]",
          "_",
          pathway_D
        ),
        ".pdf"
      )
    ),
    plot = figure_D,
    width = 10.5,
    height = 3.3,
    units = "in"
  )

  writeLines(
    pretty_D,
    file.path(
      tables_dir,
      "Figure_D_pathway_name.txt"
    )
  )
}


# =============================================================================
# 16. FIGURE E — Representative Adult-enriched GSEA curves
# =============================================================================

if (nrow(representative_adult) == 1) {

  pathway_E <- representative_adult$pathway[1]
  genes_E <- pathway_gene_lists[[pathway_E]]

  results_E <- gsea_long %>%
    filter(
      pathway == pathway_E
    )

  pretty_E <- representative_adult$PathwayLabel[1]

  E_172274 <- results_E %>%
    filter(
      Dataset == "GSE172274"
    )

  E_179277 <- results_E %>%
    filter(
      Dataset == "GSE179277"
    )

  E_231409 <- results_E %>%
    filter(
      Dataset == "GSE231409"
    )

  plot_E1 <- fgsea::plotEnrichment(
    genes_E,
    stats_172274
  ) +
    labs(
      tag = "GSE172274",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(E_172274$NES[1], 2),
        "   FDR = ",
        format(
          E_172274$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = ADULT_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  plot_E2 <- fgsea::plotEnrichment(
    genes_E,
    stats_179277
  ) +
    labs(
      tag = "GSE179277",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(E_179277$NES[1], 2),
        "   FDR = ",
        format(
          E_179277$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = ADULT_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  plot_E3 <- fgsea::plotEnrichment(
    genes_E,
    stats_231409
  ) +
    labs(
      tag = "GSE231409",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(E_231409$NES[1], 2),
        "   FDR = ",
        format(
          E_231409$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = ADULT_COLOR,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  figure_E <- (
    plot_E1 |
      plot_E2 |
      plot_E3
  ) +
    plot_annotation(
      tag_levels = NULL
    )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Figure_E_GSEA_adult_",
        gsub(
          "[^A-Za-z0-9_.-]",
          "_",
          pathway_E
        ),
        ".tiff"
      )
    ),
    plot = figure_E,
    width = 10.5,
    height = 3.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Figure_E_GSEA_adult_",
        gsub(
          "[^A-Za-z0-9_.-]",
          "_",
          pathway_E
        ),
        ".pdf"
      )
    ),
    plot = figure_E,
    width = 10.5,
    height = 3.3,
    units = "in"
  )

  writeLines(
    pretty_E,
    file.path(
      tables_dir,
      "Figure_E_pathway_name.txt"
    )
  )
}


# =============================================================================
# 17. Export all selected enrichment curves separately for supplementary use
# =============================================================================

for (pathway_i in curve_pathways_table$pathway) {

  genes_i <- pathway_gene_lists[[pathway_i]]

  if (is.null(genes_i)) {
    next
  }

  result_i <- gsea_long %>%
    filter(
      pathway == pathway_i
    )

  result_172274 <- result_i %>%
    filter(
      Dataset == "GSE172274"
    )

  result_179277 <- result_i %>%
    filter(
      Dataset == "GSE179277"
    )

  result_231409 <- result_i %>%
    filter(
      Dataset == "GSE231409"
    )

  direction_color <- if (
    mean(
      result_i$NES,
      na.rm = TRUE
    ) < 0
  ) {
    ADULT_COLOR
  } else {
    PEDIATRIC_COLOR
  }

  supplementary_1 <- fgsea::plotEnrichment(
    genes_i,
    stats_172274
  ) +
    labs(
      tag = "GSE172274",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(result_172274$NES[1], 2),
        "   FDR = ",
        format(
          result_172274$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = direction_color,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  supplementary_2 <- fgsea::plotEnrichment(
    genes_i,
    stats_179277
  ) +
    labs(
      tag = "GSE179277",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(result_179277$NES[1], 2),
        "   FDR = ",
        format(
          result_179277$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = direction_color,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  supplementary_3 <- fgsea::plotEnrichment(
    genes_i,
    stats_231409
  ) +
    labs(
      tag = "GSE231409",
      x = "Rank in ordered gene list",
      y = "Enrichment score",
      caption = paste0(
        "NES = ",
        round(result_231409$NES[1], 2),
        "   FDR = ",
        format(
          result_231409$padj[1],
          digits = 2,
          scientific = TRUE
        )
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        colour = direction_color,
        size = 10
      ),
      plot.caption = element_text(
        size = 7.5,
        hjust = 0.5
      )
    )

  supplementary_panel <- (
    supplementary_1 |
      supplementary_2 |
      supplementary_3
  )

  safe_name <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    pathway_i
  )

  ggsave(
    filename = file.path(
      curves_dir,
      paste0(
        "GSEA_",
        safe_name,
        ".tiff"
      )
    ),
    plot = supplementary_panel,
    width = 10.5,
    height = 3.3,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      curves_dir,
      paste0(
        "GSEA_",
        safe_name,
        ".pdf"
      )
    ),
    plot = supplementary_panel,
    width = 10.5,
    height = 3.3,
    units = "in"
  )
}


# =============================================================================
# 18. Export leading-edge gene table
# =============================================================================

leading_edge_gene_table <- leading_edge_selected %>%
  dplyr::select(
    Collection,
    pathway,
    PathwayLabel,
    EnrichmentGroup,
    Mean_NES,
    N_LE_GSE172274,
    N_LE_GSE179277,
    N_LE_GSE231409,
    N_LE_shared_all3,
    Shared_LE_all3
  )

write.csv(
  leading_edge_gene_table,
  file.path(
    tables_dir,
    "Table_shared_leading_edge_genes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 19. Figure summary
# =============================================================================

summary_text <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  "",
  "Figures are exported separately and contain no descriptive titles/subtitles.",
  "",
  "Color convention:",
  paste(
    "  Adult-enriched:",
    ADULT_COLOR
  ),
  paste(
    "  Pediatric-enriched:",
    PEDIATRIC_COLOR
  ),
  "",
  "Main outputs:",
  "  Figure_A_GSEA_global_concordance_heatmap.tiff",
  "  Figure_B_GSEA_metabolism_heatmap.tiff",
  "  Figure_C_GSEA_shared_leading_edge.tiff",
  "  Figure_D_GSEA_pediatric_<PATHWAY>.tiff",
  "  Figure_E_GSEA_adult_<PATHWAY>.tiff",
  "",
  paste(
    "Global representative pathways:",
    nrow(selected_global)
  ),
  paste(
    "Metabolic representative pathways:",
    nrow(selected_metabolic)
  ),
  paste(
    "Pathways available for supplementary curves:",
    nrow(curve_pathways_table)
  ),
  "",
  "Pathway labels are shown in sentence case.",
  "Heatmap asterisks represent cohort-level GSEA FDR thresholds:",
  "  *  FDR < 0.05",
  "  ** FDR < 0.01"
)

writeLines(
  summary_text,
  file.path(
    output_dir,
    "GSEA_figure_summary.txt"
  )
)

cat(
  paste(
    summary_text,
    collapse = "\n"
  ),
  "\n"
)


# =============================================================================
# 20. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_GSEA_figures.txt"
  )
)


# =============================================================================
# 21. Final message
# =============================================================================

message("")
message("============================================================")
message("GSEA figures generated successfully.")
message("============================================================")
message("Figures saved to: ", output_dir)
message("============================================================")
