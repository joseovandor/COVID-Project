# =============================================================================
# 01_GSEA_by_cohort.R
#
# Pre-ranked GSEA performed independently in each RNA-seq cohort.
#
# Reviewer-facing objective
#   The reviewer requested a pathway-level cross-cohort analysis using COMPLETE
#   ranked gene lists rather than over-representation analysis of a small set
#   of significant genes. This script performs the first stage of that analysis:
#   identical pre-ranked GSEA within each cohort.
#
# Cohorts
#   - GSE172274
#   - GSE179277
#   - GSE231409
#
# Ranking metric
#   - DESeq2 Wald statistic exported by each final cohort-level pipeline.
#   - Positive statistic = Pediatric-enriched.
#   - Negative statistic = Adult-enriched.
#   - No p-value, FDR, direction, or effect-size filtering is applied before GSEA.
#
# Gene-set collections
#   - MSigDB Hallmark
#   - MSigDB Reactome
#   - MSigDB KEGG (all currently available KEGG subcollections)
#   - MSigDB Gene Ontology Biological Process
#
# GSEA engine
#   - fgsea::fgseaMultilevel()
#
# Primary outputs
#   results/gsea/by_cohort/<DATASET>/
#       GSEA_Hallmark.csv
#       GSEA_Reactome.csv
#       GSEA_KEGG.csv
#       GSEA_GO_BP.csv
#       GSEA_all_collections.csv
#       GSEA_significant_FDR05.csv
#       GSEA_summary.txt
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project.
# =============================================================================

SCRIPT_BUILD <- "GSEA_BY_COHORT_2026-08-17_v1"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

required_packages <- c(
  "fgsea",
  "msigdbr",
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
    "\n\nInstall CRAN packages with:\n",
    "install.packages(c('msigdbr', 'tidyverse', 'here'))\n\n",
    "Install fgsea with Bioconductor using:\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install('fgsea')"
  )
}

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(tidyverse)
  library(here)
})


# =============================================================================
# 01. Configuration
# =============================================================================

ALPHA <- 0.05

# Conventional gene-set-size limits for pathway-level GSEA.
MIN_GENESET_SIZE <- 15
MAX_GENESET_SIZE <- 500

# One process is used for deterministic, Windows-safe execution.
FGSEA_NPROC <- 1

dataset_ids <- c(
  "GSE172274",
  "GSE179277",
  "GSE231409"
)

input_files <- c(
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

results_root <- here::here(
  "results",
  "gsea"
)

by_cohort_dir <- file.path(
  results_root,
  "by_cohort"
)

diagnostics_dir <- file.path(
  results_root,
  "diagnostics"
)

dir.create(
  by_cohort_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 02. Validate cohort-ranked input files
# =============================================================================

missing_input_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "The following GSEA-ranked input files were not found:\n",
    paste(
      paste0(
        names(missing_input_files),
        ": ",
        missing_input_files
      ),
      collapse = "\n"
    ),
    "\n\nRun the final DESeq2 cohort scripts first."
  )
}

message("All three cohort-level GSEA ranked files were found.")


# =============================================================================
# 03. Retrieve MSigDB human gene sets
# =============================================================================

# Retrieve the human MSigDB collections once so that exactly the same gene-set
# definitions are used for all three cohorts.

message("Retrieving MSigDB human gene sets with msigdbr...")

msig_human <- msigdbr::msigdbr(
  species = "Homo sapiens"
)

required_msig_columns <- c(
  "gs_name",
  "gene_symbol",
  "gs_collection",
  "gs_subcollection"
)

missing_msig_columns <- setdiff(
  required_msig_columns,
  colnames(msig_human)
)

if (length(missing_msig_columns) > 0) {
  stop(
    "The installed msigdbr version does not provide expected columns: ",
    paste(missing_msig_columns, collapse = ", ")
  )
}

# Document all collections available in the installed MSigDB release.
msig_collection_inventory <- msig_human %>%
  distinct(
    gs_collection,
    gs_subcollection
  ) %>%
  arrange(
    gs_collection,
    gs_subcollection
  )

write.csv(
  msig_collection_inventory,
  file.path(
    diagnostics_dir,
    "MSigDB_collection_inventory.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 04. Define pathway collections
# =============================================================================

# Hallmark
msig_hallmark <- msig_human %>%
  filter(
    gs_collection == "H"
  ) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()

# Reactome
msig_reactome <- msig_human %>%
  filter(
    gs_collection == "C2",
    gs_subcollection == "CP:REACTOME"
  ) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()

# KEGG
# MSigDB can expose more than one KEGG subcollection depending on release.
# All subcollections beginning with "CP:KEGG" are included and pathway names
# are de-duplicated before GSEA.
msig_kegg <- msig_human %>%
  filter(
    gs_collection == "C2",
    grepl(
      "^CP:KEGG",
      gs_subcollection
    )
  ) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()

# Gene Ontology: Biological Process
msig_go_bp <- msig_human %>%
  filter(
    gs_collection == "C5",
    gs_subcollection == "GO:BP"
  ) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()


# =============================================================================
# 05. Validate pathway collections
# =============================================================================

collection_sizes <- data.frame(
  Collection = c(
    "Hallmark",
    "Reactome",
    "KEGG",
    "GO_BP"
  ),
  N_pathways = c(
    length(unique(msig_hallmark$gs_name)),
    length(unique(msig_reactome$gs_name)),
    length(unique(msig_kegg$gs_name)),
    length(unique(msig_go_bp$gs_name))
  ),
  N_unique_genes = c(
    length(unique(msig_hallmark$gene_symbol)),
    length(unique(msig_reactome$gene_symbol)),
    length(unique(msig_kegg$gene_symbol)),
    length(unique(msig_go_bp$gene_symbol))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  collection_sizes,
  file.path(
    diagnostics_dir,
    "GSEA_gene_set_collection_sizes.csv"
  ),
  row.names = FALSE
)

print(collection_sizes)

if (any(collection_sizes$N_pathways == 0)) {
  stop(
    "At least one requested MSigDB collection contains zero pathways. ",
    "Inspect results/gsea/diagnostics/MSigDB_collection_inventory.csv ",
    "and verify the installed msigdbr/MSigDB release."
  )
}


# =============================================================================
# 06. Convert pathway tables to fgsea pathway lists
# =============================================================================

pathways_hallmark <- split(
  msig_hallmark$gene_symbol,
  msig_hallmark$gs_name
)

pathways_reactome <- split(
  msig_reactome$gene_symbol,
  msig_reactome$gs_name
)

pathways_kegg <- split(
  msig_kegg$gene_symbol,
  msig_kegg$gs_name
)

pathways_go_bp <- split(
  msig_go_bp$gene_symbol,
  msig_go_bp$gs_name
)


# =============================================================================
# 07. Preallocate combined result containers
# =============================================================================

all_gsea_results <- list()

input_diagnostics <- data.frame(
  Dataset = character(),
  N_rows_input = integer(),
  N_unique_symbols = integer(),
  N_missing_symbol = integer(),
  N_duplicate_symbol = integer(),
  N_invalid_statistic = integer(),
  Minimum_statistic = numeric(),
  Maximum_statistic = numeric(),
  stringsAsFactors = FALSE
)


# =============================================================================
# 08. Run identical pre-ranked GSEA separately in each cohort
# =============================================================================

for (dataset_id in dataset_ids) {

  message("")
  message("============================================================")
  message("Running GSEA for: ", dataset_id)
  message("============================================================")

  dataset_output_dir <- file.path(
    by_cohort_dir,
    dataset_id
  )

  dir.create(
    dataset_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08A. Read complete ranked gene list
  # ---------------------------------------------------------------------------

  ranked_table <- read.csv(
    input_files[dataset_id],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_rank_columns <- c(
    "Dataset",
    "SYMBOL",
    "ranking_statistic"
  )

  missing_rank_columns <- setdiff(
    required_rank_columns,
    colnames(ranked_table)
  )

  if (length(missing_rank_columns) > 0) {
    stop(
      dataset_id,
      " GSEA ranked file is missing required columns: ",
      paste(
        missing_rank_columns,
        collapse = ", "
      )
    )
  }

  ranked_table <- ranked_table %>%
    mutate(
      Dataset = as.character(Dataset),
      SYMBOL = as.character(SYMBOL),
      ranking_statistic = as.numeric(ranking_statistic)
    )


  # ---------------------------------------------------------------------------
  # 08B. Validate ranking
  # ---------------------------------------------------------------------------

  n_missing_symbol <- sum(
    is.na(ranked_table$SYMBOL) |
      ranked_table$SYMBOL == ""
  )

  n_duplicate_symbol <- sum(
    duplicated(ranked_table$SYMBOL)
  )

  n_invalid_statistic <- sum(
    is.na(ranked_table$ranking_statistic) |
      !is.finite(ranked_table$ranking_statistic)
  )

  if (n_missing_symbol > 0) {
    stop(
      dataset_id,
      " contains missing SYMBOL values in its GSEA ranking."
    )
  }

  if (n_duplicate_symbol > 0) {
    stop(
      dataset_id,
      " contains duplicated SYMBOL values in its GSEA ranking. ",
      "Resolve SYMBOL duplication in the cohort-level DESeq2 pipeline rather ",
      "than resolving it after ranking."
    )
  }

  if (n_invalid_statistic > 0) {
    stop(
      dataset_id,
      " contains invalid ranking statistics."
    )
  }

  if (
    any(
      ranked_table$Dataset != dataset_id,
      na.rm = TRUE
    )
  ) {
    stop(
      "Unexpected Dataset values were found in ",
      basename(
        input_files[dataset_id]
      )
    )
  }

  input_diagnostics <- bind_rows(
    input_diagnostics,
    data.frame(
      Dataset = dataset_id,
      N_rows_input = nrow(ranked_table),
      N_unique_symbols = length(
        unique(
          ranked_table$SYMBOL
        )
      ),
      N_missing_symbol = n_missing_symbol,
      N_duplicate_symbol = n_duplicate_symbol,
      N_invalid_statistic = n_invalid_statistic,
      Minimum_statistic = min(
        ranked_table$ranking_statistic
      ),
      Maximum_statistic = max(
        ranked_table$ranking_statistic
      ),
      stringsAsFactors = FALSE
    )
  )


  # ---------------------------------------------------------------------------
  # 08C. Create named pre-ranked statistic vector
  # ---------------------------------------------------------------------------

  gene_stats <- ranked_table$ranking_statistic

  names(gene_stats) <- ranked_table$SYMBOL

  gene_stats <- sort(
    gene_stats,
    decreasing = TRUE
  )

  # The final cohort pipeline defines:
  #   positive Wald statistic -> Pediatric-enriched
  #   negative Wald statistic -> Adult-enriched


  # ---------------------------------------------------------------------------
  # 08D. Hallmark GSEA
  # ---------------------------------------------------------------------------

  set.seed(172274)

  gsea_hallmark <- fgsea::fgseaMultilevel(
    pathways = pathways_hallmark,
    stats = gene_stats,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    eps = 0,
    scoreType = "std",
    nproc = FGSEA_NPROC
  )

  gsea_hallmark <- as.data.frame(
    gsea_hallmark
  )

  gsea_hallmark$Collection <- "Hallmark"
  gsea_hallmark$Dataset <- dataset_id

  if ("leadingEdge" %in% colnames(gsea_hallmark)) {
    gsea_hallmark$leadingEdge <- vapply(
      gsea_hallmark$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  gsea_hallmark <- gsea_hallmark %>%
    mutate(
      EnrichmentDirection = case_when(
        NES > 0 ~ "Pediatric",
        NES < 0 ~ "Adult",
        TRUE ~ "Zero"
      )
    ) %>%
    arrange(
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_hallmark,
    file.path(
      dataset_output_dir,
      "GSEA_Hallmark.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08E. Reactome GSEA
  # ---------------------------------------------------------------------------

  set.seed(179277)

  gsea_reactome <- fgsea::fgseaMultilevel(
    pathways = pathways_reactome,
    stats = gene_stats,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    eps = 0,
    scoreType = "std",
    nproc = FGSEA_NPROC
  )

  gsea_reactome <- as.data.frame(
    gsea_reactome
  )

  gsea_reactome$Collection <- "Reactome"
  gsea_reactome$Dataset <- dataset_id

  if ("leadingEdge" %in% colnames(gsea_reactome)) {
    gsea_reactome$leadingEdge <- vapply(
      gsea_reactome$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  gsea_reactome <- gsea_reactome %>%
    mutate(
      EnrichmentDirection = case_when(
        NES > 0 ~ "Pediatric",
        NES < 0 ~ "Adult",
        TRUE ~ "Zero"
      )
    ) %>%
    arrange(
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_reactome,
    file.path(
      dataset_output_dir,
      "GSEA_Reactome.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08F. KEGG GSEA
  # ---------------------------------------------------------------------------

  set.seed(231409)

  gsea_kegg <- fgsea::fgseaMultilevel(
    pathways = pathways_kegg,
    stats = gene_stats,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    eps = 0,
    scoreType = "std",
    nproc = FGSEA_NPROC
  )

  gsea_kegg <- as.data.frame(
    gsea_kegg
  )

  gsea_kegg$Collection <- "KEGG"
  gsea_kegg$Dataset <- dataset_id

  if ("leadingEdge" %in% colnames(gsea_kegg)) {
    gsea_kegg$leadingEdge <- vapply(
      gsea_kegg$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  gsea_kegg <- gsea_kegg %>%
    mutate(
      EnrichmentDirection = case_when(
        NES > 0 ~ "Pediatric",
        NES < 0 ~ "Adult",
        TRUE ~ "Zero"
      )
    ) %>%
    arrange(
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_kegg,
    file.path(
      dataset_output_dir,
      "GSEA_KEGG.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08G. GO Biological Process GSEA
  # ---------------------------------------------------------------------------

  set.seed(20260817)

  gsea_go_bp <- fgsea::fgseaMultilevel(
    pathways = pathways_go_bp,
    stats = gene_stats,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    eps = 0,
    scoreType = "std",
    nproc = FGSEA_NPROC
  )

  gsea_go_bp <- as.data.frame(
    gsea_go_bp
  )

  gsea_go_bp$Collection <- "GO_BP"
  gsea_go_bp$Dataset <- dataset_id

  if ("leadingEdge" %in% colnames(gsea_go_bp)) {
    gsea_go_bp$leadingEdge <- vapply(
      gsea_go_bp$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  gsea_go_bp <- gsea_go_bp %>%
    mutate(
      EnrichmentDirection = case_when(
        NES > 0 ~ "Pediatric",
        NES < 0 ~ "Adult",
        TRUE ~ "Zero"
      )
    ) %>%
    arrange(
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_go_bp,
    file.path(
      dataset_output_dir,
      "GSEA_GO_BP.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08H. Combine all collections for this cohort
  # ---------------------------------------------------------------------------

  gsea_all_dataset <- bind_rows(
    gsea_hallmark,
    gsea_reactome,
    gsea_kegg,
    gsea_go_bp
  ) %>%
    relocate(
      Dataset,
      Collection,
      pathway,
      NES,
      pval,
      padj,
      EnrichmentDirection
    ) %>%
    arrange(
      Collection,
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_all_dataset,
    file.path(
      dataset_output_dir,
      "GSEA_all_collections.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08I. Significant pathways for this cohort
  # ---------------------------------------------------------------------------

  gsea_significant_dataset <- gsea_all_dataset %>%
    filter(
      !is.na(padj),
      padj < ALPHA
    ) %>%
    arrange(
      Collection,
      padj,
      desc(
        abs(
          NES
        )
      )
    )

  write.csv(
    gsea_significant_dataset,
    file.path(
      dataset_output_dir,
      "GSEA_significant_FDR05.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08J. Collection-level summary for this cohort
  # ---------------------------------------------------------------------------

  gsea_dataset_summary <- gsea_all_dataset %>%
    group_by(
      Collection
    ) %>%
    summarise(
      N_tested_pathways = n(),
      N_FDR05 = sum(
        padj < ALPHA,
        na.rm = TRUE
      ),
      N_FDR05_Pediatric = sum(
        padj < ALPHA &
          NES > 0,
        na.rm = TRUE
      ),
      N_FDR05_Adult = sum(
        padj < ALPHA &
          NES < 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    gsea_dataset_summary,
    file.path(
      dataset_output_dir,
      "GSEA_collection_summary.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 08K. Dataset text summary
  # ---------------------------------------------------------------------------

  dataset_summary_text <- c(
    paste(
      "Script build:",
      SCRIPT_BUILD
    ),
    paste(
      "Dataset:",
      dataset_id
    ),
    "",
    "GSEA design:",
    "  Analysis: pre-ranked GSEA",
    "  Ranking metric: complete DESeq2 Wald statistic",
    "  Positive statistic: Pediatric-enriched",
    "  Negative statistic: Adult-enriched",
    "  Pre-GSEA significance filtering: none",
    paste(
      "  Ranked genes:",
      length(gene_stats)
    ),
    paste(
      "  Minimum gene-set size:",
      MIN_GENESET_SIZE
    ),
    paste(
      "  Maximum gene-set size:",
      MAX_GENESET_SIZE
    ),
    "  Multiple-testing correction: fgsea Benjamini-Hochberg padj within each collection/cohort run",
    "",
    "Pathway collections:",
    paste(
      "  Hallmark:",
      nrow(gsea_hallmark),
      "tested pathways"
    ),
    paste(
      "  Reactome:",
      nrow(gsea_reactome),
      "tested pathways"
    ),
    paste(
      "  KEGG:",
      nrow(gsea_kegg),
      "tested pathways"
    ),
    paste(
      "  GO Biological Process:",
      nrow(gsea_go_bp),
      "tested pathways"
    ),
    "",
    paste(
      "Total FDR < 0.05 pathways across exported collection-level analyses:",
      nrow(gsea_significant_dataset)
    ),
    paste(
      "FDR < 0.05 Pediatric-enriched pathways:",
      sum(
        gsea_significant_dataset$NES > 0,
        na.rm = TRUE
      )
    ),
    paste(
      "FDR < 0.05 Adult-enriched pathways:",
      sum(
        gsea_significant_dataset$NES < 0,
        na.rm = TRUE
      )
    )
  )

  writeLines(
    dataset_summary_text,
    file.path(
      dataset_output_dir,
      "GSEA_summary.txt"
    )
  )

  cat(
    paste(
      dataset_summary_text,
      collapse = "\n"
    ),
    "\n"
  )

  all_gsea_results[[dataset_id]] <- gsea_all_dataset
}


# =============================================================================
# 09. Export input diagnostics
# =============================================================================

write.csv(
  input_diagnostics,
  file.path(
    diagnostics_dir,
    "GSEA_ranked_input_diagnostics.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 10. Combine cohort-level GSEA results
# =============================================================================

gsea_all_cohorts <- bind_rows(
  all_gsea_results
) %>%
  arrange(
    Collection,
    pathway,
    Dataset
  )

write.csv(
  gsea_all_cohorts,
  file.path(
    results_root,
    "GSEA_all_cohorts_long.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. Export significant pathways across all cohort-specific runs
# =============================================================================

gsea_all_cohorts_significant <- gsea_all_cohorts %>%
  filter(
    !is.na(padj),
    padj < ALPHA
  ) %>%
  arrange(
    Collection,
    pathway,
    Dataset
  )

write.csv(
  gsea_all_cohorts_significant,
  file.path(
    results_root,
    "GSEA_all_cohorts_significant_FDR05.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 12. Export preliminary cross-cohort pathway coverage
# =============================================================================

# This is NOT yet the final pathway concordance analysis. It simply identifies
# which pathways were tested in one, two, or all three cohorts so that script 02
# can perform the formal cross-cohort concordance step.

pathway_coverage <- gsea_all_cohorts %>%
  distinct(
    Collection,
    pathway,
    Dataset
  ) %>%
  count(
    Collection,
    pathway,
    name = "N_cohorts_tested"
  ) %>%
  arrange(
    Collection,
    desc(
      N_cohorts_tested
    ),
    pathway
  )

write.csv(
  pathway_coverage,
  file.path(
    results_root,
    "GSEA_pathway_coverage_across_cohorts.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 13. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_GSEA_by_cohort.txt"
  )
)

package_versions <- c(
  paste(
    "R version:",
    R.version.string
  ),
  paste(
    "fgsea version:",
    as.character(
      packageVersion(
        "fgsea"
      )
    )
  ),
  paste(
    "msigdbr version:",
    as.character(
      packageVersion(
        "msigdbr"
      )
    )
  ),
  paste(
    "tidyverse version:",
    as.character(
      packageVersion(
        "tidyverse"
      )
    )
  ),
  paste(
    "here version:",
    as.character(
      packageVersion(
        "here"
      )
    )
  )
)

if ("db_version" %in% colnames(msig_human)) {

  msigdb_versions <- unique(
    msig_human$db_version
  )

  package_versions <- c(
    package_versions,
    paste(
      "MSigDB version(s):",
      paste(
        msigdb_versions,
        collapse = ", "
      )
    )
  )
}

writeLines(
  package_versions,
  file.path(
    diagnostics_dir,
    "package_versions_GSEA.txt"
  )
)


# =============================================================================
# 14. Final summary
# =============================================================================

final_summary <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  "",
  "Cohort-specific pre-ranked GSEA completed.",
  "",
  "Datasets:",
  paste0(
    "  ",
    dataset_ids
  ),
  "",
  "Ranking:",
  "  Complete final-pipeline DESeq2 Wald statistics",
  "  Positive NES = Pediatric-enriched",
  "  Negative NES = Adult-enriched",
  "",
  "Collections:",
  "  Hallmark",
  "  Reactome",
  "  KEGG",
  "  GO Biological Process",
  "",
  paste(
    "Total cohort-pathway results:",
    nrow(gsea_all_cohorts)
  ),
  paste(
    "Total cohort-pathway results with FDR < 0.05:",
    nrow(gsea_all_cohorts_significant)
  ),
  "",
  "Next analysis:",
  "  02_GSEA_cross_cohort_concordance.R",
  "  This will compare NES, FDR, direction, and leading-edge genes across cohorts."
)

writeLines(
  final_summary,
  file.path(
    results_root,
    "GSEA_by_cohort_summary.txt"
  )
)

cat(
  paste(
    final_summary,
    collapse = "\n"
  ),
  "\n"
)

message("")
message("============================================================")
message("Cohort-specific pre-ranked GSEA completed successfully.")
message("============================================================")
message("Results saved to: ", results_root)
message("============================================================")
