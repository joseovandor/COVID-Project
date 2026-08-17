# =============================================================================
# 05_previous_signature_concordance_v2.R
#
# Reviewer-facing concordance analysis of the two previously reported gene sets
# using ONLY results from the final harmonized meta-analysis pipeline.
#
# Historical sets
#   - Original manuscript: 17 genes
#   - First revised manuscript: 19 genes
#
# Main reviewer-facing figures
#   Figure A. Forest plot of the original 17-gene set
#   Figure B. Heatmap of cohort-specific effects for the original 17-gene set
#   Figure C. Forest plot of the previous 19-gene set
#   Figure D. Heatmap of cohort-specific effects for the previous 19-gene set
#
# Interpretation
#   These figures are not intended to "rescue" previously reported genes.
#   They show how those historical candidates behave when re-evaluated using
#   the final pipeline: cohort-specific direction, pooled random-effects
#   estimate, 95% CI, final FDR, and between-study heterogeneity.
#
# Positive log2FC = Pediatric-enriched
# Negative log2FC = Adult-enriched
# =============================================================================

SCRIPT_BUILD <- "META_SIGNATURE_CONCORDANCE_2026-08-17_v3"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

required_packages <- c(
  "tidyverse",
  "here"
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
  library(tidyverse)
  library(here)
})


# =============================================================================
# 01. Configuration
# =============================================================================

ALPHA <- 0.05

meta_all_file <- here::here(
  "results",
  "meta_analysis",
  "primary",
  "meta_analysis_REM_all_genes.csv"
)

cohort_files <- c(
  GSE172274 = here::here(
    "results", "deseq2", "GSE172274", "primary",
    "GSE172274_meta_input.csv"
  ),
  GSE179277 = here::here(
    "results", "deseq2", "GSE179277", "primary",
    "GSE179277_meta_input.csv"
  ),
  GSE231409 = here::here(
    "results", "deseq2", "GSE231409", "primary",
    "GSE231409_meta_input.csv"
  )
)

output_dir <- here::here(
  "results",
  "meta_analysis",
  "signature_concordance"
)

figures_dir <- file.path(
  output_dir,
  "figures"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 02. Historical gene lists exactly as previously reported
# =============================================================================

original_17 <- c(
  "SAMD9L",
  "FGD2",
  "NUPR1",
  "STARD4",
  "ATP10A",
  "SAMD9",
  "EPSTI1",
  "IFIH1",
  "NT5C3A",
  "MT2A",
  "IL18BP",
  "IFITM3",
  "TOR1B",
  "CXCL14",
  "IL13",
  "SIX6",
  "HSPB6"
)

previous_19 <- c(
  "GJD2",
  "DCAF8L2",
  "PLCZ1",
  "SPACA3",
  "EDDM3B",
  "CT55",
  "C3orf85",
  "TCF23",
  "TBR1",
  "OR10G4",
  "CDC14C",
  "OR5H14",
  "SCN7A",
  "SATL1",
  "OR5A1",
  "DDX39B",
  "CD83",
  "CXCL8",
  "IL3RA"
)


# Historical enrichment direction as reported in each submitted manuscript.
# This is used only for visual grouping of the previously reported candidates.
# Statistical direction in the FINAL analysis is always read from the new
# cohort-level and pooled estimates.

original_17_pediatric <- c(
  "SAMD9L",
  "FGD2",
  "NUPR1",
  "STARD4",
  "ATP10A",
  "SAMD9",
  "EPSTI1",
  "IFIH1",
  "NT5C3A",
  "MT2A",
  "IL18BP"
)

original_17_adult <- c(
  "IFITM3",
  "TOR1B",
  "CXCL14",
  "IL13",
  "SIX6",
  "HSPB6"
)

previous_19_adult <- c(
  "GJD2",
  "DCAF8L2",
  "PLCZ1",
  "SPACA3",
  "EDDM3B",
  "CT55",
  "C3orf85",
  "TCF23",
  "TBR1",
  "OR10G4",
  "CDC14C",
  "OR5H14",
  "SCN7A",
  "SATL1",
  "OR5A1"
)

previous_19_pediatric <- c(
  "DDX39B",
  "CD83",
  "CXCL8",
  "IL3RA"
)

ADULT_COLOR <- "#3E76BC"
PEDIATRIC_COLOR <- "#FCCE24"


# =============================================================================
# 03. Validate input files
# =============================================================================

required_files <- c(
  meta_all_file,
  cohort_files
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# =============================================================================
# 04. Read final meta-analysis and cohort-level results
# =============================================================================

meta_all <- read.csv(
  meta_all_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gse172274 <- read.csv(
  cohort_files["GSE172274"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gse179277 <- read.csv(
  cohort_files["GSE179277"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gse231409 <- read.csv(
  cohort_files["GSE231409"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# =============================================================================
# 05. Keep only columns required for reviewer-facing comparison
# =============================================================================

gse172274_small <- gse172274 %>%
  dplyr::select(
    SYMBOL,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  rename(
    GSE172274_log2FC = log2FoldChange,
    GSE172274_SE = lfcSE,
    GSE172274_pvalue = pvalue,
    GSE172274_padj = padj
  )

gse179277_small <- gse179277 %>%
  dplyr::select(
    SYMBOL,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  rename(
    GSE179277_log2FC = log2FoldChange,
    GSE179277_SE = lfcSE,
    GSE179277_pvalue = pvalue,
    GSE179277_padj = padj
  )

gse231409_small <- gse231409 %>%
  dplyr::select(
    SYMBOL,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  rename(
    GSE231409_log2FC = log2FoldChange,
    GSE231409_SE = lfcSE,
    GSE231409_pvalue = pvalue,
    GSE231409_padj = padj
  )

meta_small <- meta_all %>%
  dplyr::select(
    SYMBOL,
    k,
    RE_log2FC,
    RE_SE,
    RE_CI_lower,
    RE_CI_upper,
    RE_pvalue,
    RE_FDR,
    FE_log2FC,
    FE_pvalue,
    FE_FDR,
    tau2,
    I2,
    Q_pvalue,
    direction_pattern,
    direction_consistent
  )


# =============================================================================
# 06. Build one final-pipeline table containing all historical genes
# =============================================================================

historical_genes <- union(
  original_17,
  previous_19
)

historical_results <- data.frame(
  SYMBOL = historical_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(
    gse172274_small,
    by = "SYMBOL"
  ) %>%
  left_join(
    gse179277_small,
    by = "SYMBOL"
  ) %>%
  left_join(
    gse231409_small,
    by = "SYMBOL"
  ) %>%
  left_join(
    meta_small,
    by = "SYMBOL"
  ) %>%
  mutate(
    Historical_set = case_when(
      SYMBOL %in% original_17 ~ "Original 17-gene set",
      SYMBOL %in% previous_19 ~ "Previous 19-gene set",
      TRUE ~ "Other"
    ),
    
    HistoricalEnrichment = case_when(
      SYMBOL %in% original_17_adult ~ "Adult-enriched",
      SYMBOL %in% original_17_pediatric ~ "Pediatric-enriched",
      SYMBOL %in% previous_19_adult ~ "Adult-enriched",
      SYMBOL %in% previous_19_pediatric ~ "Pediatric-enriched",
      TRUE ~ "Unclassified"
    ),
    
    Final_significant = (
      !is.na(RE_FDR) &
        RE_FDR < ALPHA
    ),
    
    Final_status = case_when(
      is.na(k) ~ "Not meta-analysis eligible",
      RE_FDR < ALPHA &
        direction_consistent %in% TRUE ~
        "FDR < 0.05; direction consistent",
      RE_FDR < ALPHA &
        direction_consistent %in% FALSE ~
        "FDR < 0.05; direction discordant",
      RE_FDR >= ALPHA &
        direction_consistent %in% FALSE ~
        "FDR >= 0.05; direction discordant",
      RE_FDR >= ALPHA ~
        "FDR >= 0.05",
      TRUE ~ "Not retained"
    )
  )

write.csv(
  historical_results,
  file.path(
    tables_dir,
    "Supplementary_Table_previous_signatures_final_pipeline.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 07. Split final-pipeline results by historical list
# =============================================================================

original_17_results <- historical_results %>%
  filter(
    SYMBOL %in% original_17
  )

previous_19_results <- historical_results %>%
  filter(
    SYMBOL %in% previous_19
  )

original_17_results$SYMBOL <- factor(
  original_17_results$SYMBOL,
  levels = original_17
)

previous_19_results$SYMBOL <- factor(
  previous_19_results$SYMBOL,
  levels = previous_19
)

original_17_results <- original_17_results %>%
  arrange(
    SYMBOL
  )

previous_19_results <- previous_19_results %>%
  arrange(
    SYMBOL
  )

original_17_results$SYMBOL <- as.character(
  original_17_results$SYMBOL
)

previous_19_results$SYMBOL <- as.character(
  previous_19_results$SYMBOL
)

write.csv(
  original_17_results,
  file.path(
    tables_dir,
    "Original_17_final_pipeline_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  previous_19_results,
  file.path(
    tables_dir,
    "Previous_19_final_pipeline_results.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 08. Forest plot data: original 17-gene set
# =============================================================================

forest_original_17 <- original_17_results %>%
  mutate(
    GeneLabel = ifelse(
      !is.na(RE_FDR) & RE_FDR < ALPHA,
      paste0(SYMBOL, " *"),
      SYMBOL
    ),
    
    HistoricalEnrichment = factor(
      HistoricalEnrichment,
      levels = c(
        "Adult-enriched",
        "Pediatric-enriched"
      )
    )
  )

# Keep ALL 17 historical genes in the plot. Genes without a valid pooled
# estimate remain represented in the data table; the forest geometry is drawn
# only where the final REM estimate and CI are available.
forest_original_17$GeneLabel <- factor(
  forest_original_17$GeneLabel,
  levels = rev(
    c(
      ifelse(
        original_17_adult %in%
          forest_original_17$SYMBOL[
            !is.na(forest_original_17$RE_FDR) &
              forest_original_17$RE_FDR < ALPHA
          ],
        paste0(original_17_adult, " *"),
        original_17_adult
      ),
      ifelse(
        original_17_pediatric %in%
          forest_original_17$SYMBOL[
            !is.na(forest_original_17$RE_FDR) &
              forest_original_17$RE_FDR < ALPHA
          ],
        paste0(original_17_pediatric, " *"),
        original_17_pediatric
      )
    )
  )
)


# =============================================================================
# 09. Forest plot: original 17-gene set
# =============================================================================

p_forest_original_17 <- ggplot(
  forest_original_17,
  aes(
    x = RE_log2FC,
    y = GeneLabel,
    colour = HistoricalEnrichment
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = RE_CI_lower,
      xmax = RE_CI_upper
    ),
    height = 0.18,
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    size = 3.2,
    na.rm = TRUE
  ) +
  facet_grid(
    HistoricalEnrichment ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_colour_manual(
    values = c(
      "Adult-enriched" = ADULT_COLOR,
      "Pediatric-enriched" = PEDIATRIC_COLOR
    )
  ) +
  labs(
    title = "Original 17-gene set re-evaluated with the final meta-analysis",
    subtitle = "* Final random-effects FDR < 0.05",
    x = "Pooled log2 fold change (Pediatric vs Adult)",
    y = NULL,
    colour = "Previously reported group"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold"
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    figures_dir,
    "Forest_original_17_final_meta_analysis.tiff"
  ),
  p_forest_original_17,
  width = 9,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    figures_dir,
    "Forest_original_17_final_meta_analysis.pdf"
  ),
  p_forest_original_17,
  width = 9,
  height = 8
)


# =============================================================================
# 10. Forest plot data: previous 19-gene set
# =============================================================================

forest_previous_19 <- previous_19_results %>%
  mutate(
    GeneLabel = ifelse(
      !is.na(RE_FDR) & RE_FDR < ALPHA,
      paste0(SYMBOL, " *"),
      SYMBOL
    ),
    
    HistoricalEnrichment = factor(
      HistoricalEnrichment,
      levels = c(
        "Adult-enriched",
        "Pediatric-enriched"
      )
    )
  )

forest_previous_19$GeneLabel <- factor(
  forest_previous_19$GeneLabel,
  levels = rev(
    c(
      ifelse(
        previous_19_adult %in%
          forest_previous_19$SYMBOL[
            !is.na(forest_previous_19$RE_FDR) &
              forest_previous_19$RE_FDR < ALPHA
          ],
        paste0(previous_19_adult, " *"),
        previous_19_adult
      ),
      ifelse(
        previous_19_pediatric %in%
          forest_previous_19$SYMBOL[
            !is.na(forest_previous_19$RE_FDR) &
              forest_previous_19$RE_FDR < ALPHA
          ],
        paste0(previous_19_pediatric, " *"),
        previous_19_pediatric
      )
    )
  )
)


# =============================================================================
# 11. Forest plot: previous 19-gene set
# =============================================================================

p_forest_previous_19 <- ggplot(
  forest_previous_19,
  aes(
    x = RE_log2FC,
    y = GeneLabel,
    colour = HistoricalEnrichment
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = RE_CI_lower,
      xmax = RE_CI_upper
    ),
    height = 0.18,
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    size = 3.2,
    na.rm = TRUE
  ) +
  facet_grid(
    HistoricalEnrichment ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_colour_manual(
    values = c(
      "Adult-enriched" = ADULT_COLOR,
      "Pediatric-enriched" = PEDIATRIC_COLOR
    )
  ) +
  labs(
    title = "Previous 19-gene set re-evaluated with the final meta-analysis",
    subtitle = "* Final random-effects FDR < 0.05",
    x = "Pooled log2 fold change (Pediatric vs Adult)",
    y = NULL,
    colour = "Previously reported group"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold"
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    figures_dir,
    "Forest_previous_19_final_meta_analysis.tiff"
  ),
  p_forest_previous_19,
  width = 9,
  height = 8.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    figures_dir,
    "Forest_previous_19_final_meta_analysis.pdf"
  ),
  p_forest_previous_19,
  width = 9,
  height = 8.5
)


# =============================================================================
# 12. Heatmap preparation: original 17-gene set
# =============================================================================

heatmap_original_17 <- original_17_results %>%
  mutate(
    GeneLabel = ifelse(
      !is.na(RE_FDR) & RE_FDR < ALPHA,
      paste0(SYMBOL, " *"),
      SYMBOL
    )
  ) %>%
  dplyr::select(
    SYMBOL,
    GeneLabel,
    HistoricalEnrichment,
    GSE172274_log2FC,
    GSE179277_log2FC,
    GSE231409_log2FC,
    RE_log2FC
  ) %>%
  pivot_longer(
    cols = c(
      GSE172274_log2FC,
      GSE179277_log2FC,
      GSE231409_log2FC,
      RE_log2FC
    ),
    names_to = "Dataset",
    values_to = "log2FC"
  ) %>%
  mutate(
    Dataset = recode(
      Dataset,
      "GSE172274_log2FC" = "GSE172274",
      "GSE179277_log2FC" = "GSE179277",
      "GSE231409_log2FC" = "GSE231409",
      "RE_log2FC" = "Pooled REM"
    ),
    HistoricalEnrichment = factor(
      HistoricalEnrichment,
      levels = c(
        "Adult-enriched",
        "Pediatric-enriched"
      )
    )
  )

heatmap_original_17$Dataset <- factor(
  heatmap_original_17$Dataset,
  levels = c(
    "GSE172274",
    "GSE179277",
    "GSE231409",
    "Pooled REM"
  )
)

heatmap_original_17$GeneLabel <- factor(
  heatmap_original_17$GeneLabel,
  levels = rev(
    unique(
      forest_original_17$GeneLabel
    )
  )
)


# =============================================================================
# 13. Heatmap: original 17-gene set
# =============================================================================

p_heatmap_original_17 <- ggplot(
  heatmap_original_17,
  aes(
    x = Dataset,
    y = GeneLabel,
    fill = log2FC
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  facet_grid(
    HistoricalEnrichment ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = ADULT_COLOR,
    mid = "white",
    high = PEDIATRIC_COLOR,
    midpoint = 0,
    name = "log2FC
Pediatric vs Adult",
    na.value = "grey90"
  ) +
  labs(
    title = "Cohort-level effect concordance of the original 17-gene set",
    subtitle = "* Final random-effects FDR < 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold"
    )
  )

ggsave(
  file.path(
    figures_dir,
    "Heatmap_original_17_cohort_concordance.tiff"
  ),
  p_heatmap_original_17,
  width = 9,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    figures_dir,
    "Heatmap_original_17_cohort_concordance.pdf"
  ),
  p_heatmap_original_17,
  width = 9,
  height = 8
)


# =============================================================================
# 14. Heatmap preparation: previous 19-gene set
# =============================================================================

heatmap_previous_19 <- previous_19_results %>%
  mutate(
    GeneLabel = ifelse(
      !is.na(RE_FDR) & RE_FDR < ALPHA,
      paste0(SYMBOL, " *"),
      SYMBOL
    )
  ) %>%
  dplyr::select(
    SYMBOL,
    GeneLabel,
    HistoricalEnrichment,
    GSE172274_log2FC,
    GSE179277_log2FC,
    GSE231409_log2FC,
    RE_log2FC
  ) %>%
  pivot_longer(
    cols = c(
      GSE172274_log2FC,
      GSE179277_log2FC,
      GSE231409_log2FC,
      RE_log2FC
    ),
    names_to = "Dataset",
    values_to = "log2FC"
  ) %>%
  mutate(
    Dataset = recode(
      Dataset,
      "GSE172274_log2FC" = "GSE172274",
      "GSE179277_log2FC" = "GSE179277",
      "GSE231409_log2FC" = "GSE231409",
      "RE_log2FC" = "Pooled REM"
    ),
    HistoricalEnrichment = factor(
      HistoricalEnrichment,
      levels = c(
        "Adult-enriched",
        "Pediatric-enriched"
      )
    )
  )

heatmap_previous_19$Dataset <- factor(
  heatmap_previous_19$Dataset,
  levels = c(
    "GSE172274",
    "GSE179277",
    "GSE231409",
    "Pooled REM"
  )
)

heatmap_previous_19$GeneLabel <- factor(
  heatmap_previous_19$GeneLabel,
  levels = rev(
    unique(
      forest_previous_19$GeneLabel
    )
  )
)


# =============================================================================
# 15. Heatmap: previous 19-gene set
# =============================================================================

p_heatmap_previous_19 <- ggplot(
  heatmap_previous_19,
  aes(
    x = Dataset,
    y = GeneLabel,
    fill = log2FC
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  facet_grid(
    HistoricalEnrichment ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = ADULT_COLOR,
    mid = "white",
    high = PEDIATRIC_COLOR,
    midpoint = 0,
    name = "log2FC
Pediatric vs Adult",
    na.value = "grey90"
  ) +
  labs(
    title = "Cohort-level effect concordance of the previous 19-gene set",
    subtitle = "* Final random-effects FDR < 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold"
    )
  )

ggsave(
  file.path(
    figures_dir,
    "Heatmap_previous_19_cohort_concordance.tiff"
  ),
  p_heatmap_previous_19,
  width = 9,
  height = 8.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    figures_dir,
    "Heatmap_previous_19_cohort_concordance.pdf"
  ),
  p_heatmap_previous_19,
  width = 9,
  height = 8.5
)


# =============================================================================
# 16. Reviewer-oriented compact table
# =============================================================================

reviewer_table <- historical_results %>%
  dplyr::select(
    Historical_set,
    HistoricalEnrichment,
    SYMBOL,
    GSE172274_log2FC,
    GSE172274_SE,
    GSE172274_padj,
    GSE179277_log2FC,
    GSE179277_SE,
    GSE179277_padj,
    GSE231409_log2FC,
    GSE231409_SE,
    GSE231409_padj,
    k,
    RE_log2FC,
    RE_CI_lower,
    RE_CI_upper,
    RE_pvalue,
    RE_FDR,
    FE_log2FC,
    FE_FDR,
    tau2,
    I2,
    Q_pvalue,
    direction_pattern,
    direction_consistent,
    Final_status
  )

write.csv(
  reviewer_table,
  file.path(
    tables_dir,
    "Reviewer_Table_historical_genes_final_pipeline.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 17. Summary
# =============================================================================

summary_text <- c(
  paste("Script build:", SCRIPT_BUILD),
  "",
  "Reviewer-facing historical-signature re-evaluation",
  "",
  paste(
    "Original 17 genes retained at final RE FDR < 0.05:",
    sum(
      original_17_results$RE_FDR < ALPHA,
      na.rm = TRUE
    ),
    "/",
    length(original_17)
  ),
  paste(
    "Previous 19 genes retained at final RE FDR < 0.05:",
    sum(
      previous_19_results$RE_FDR < ALPHA,
      na.rm = TRUE
    ),
    "/",
    length(previous_19)
  ),
  "",
  "Figures:",
  "  Forest_original_17_final_meta_analysis.tiff",
  "  Heatmap_original_17_cohort_concordance.tiff",
  "  Forest_previous_19_final_meta_analysis.tiff",
  "  Heatmap_previous_19_cohort_concordance.tiff",
  "",
  "All previously reported genes are displayed regardless of current significance.",
  "Adult-enriched and pediatric-enriched historical groups are shown separately.",
  "Forest plot points and confidence-interval lines are colored by the originally reported group.",
  "An asterisk marks genes with final random-effects FDR < 0.05.",
  "The heatmaps show cohort-specific effect directions and the pooled REM effect.",
  "Genes are not excluded from these plots because they fail FDR; the purpose is",
  "to transparently show why previously reported candidates are or are not",
  "supported by the final harmonized analysis."
)

writeLines(
  summary_text,
  file.path(
    output_dir,
    "reviewer_signature_concordance_summary.txt"
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
# 18. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo_previous_signature_concordance.txt"
  )
)


# =============================================================================
# 19. Final message
# =============================================================================

message("")
message("============================================================")
message("Historical signature concordance figures completed.")
message("============================================================")
message("Figures saved to: ", figures_dir)
message("Tables saved to: ", tables_dir)
message("============================================================")
