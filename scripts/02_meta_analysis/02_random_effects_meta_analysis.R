# =============================================================================
# 02_random_effects_meta_analysis.R
#
# Random-effects meta-analysis of cohort-level differential-expression results
#
# Input
#   results/meta_analysis/harmonization/meta_input_combined_long.csv
#
# Datasets
#   - GSE172274
#   - GSE179277
#   - GSE231409
#
# Primary model
#   - Random-effects meta-analysis
#   - REML estimator for between-study variance (tau^2)
#   - Effect measure: unshrunken DESeq2 log2FoldChange
#   - Standard error: DESeq2 lfcSE
#   - Positive pooled effect = higher expression in Pediatric samples
#
# Sensitivity model
#   - Fixed-effect inverse-variance model
#
# Multiple testing
#   - Benjamini-Hochberg FDR calculated across all genes with k >= 2
#   - Separate k=3 and k=2 tables are exported for reporting, but the primary
#     FDR is not recalculated within those subsets.
#
# Heterogeneity
#   - tau^2
#   - I^2
#   - Cochran's Q
#   - Q-test p-value
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project. Paths are resolved
#   with here::here(); no machine-specific working directory is required.
#
# The script is intentionally linear and does not define custom functions.
# =============================================================================

SCRIPT_BUILD <- "META_REM_2026-08-17_v1"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

required_packages <- c(
  "metafor",
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
  library(metafor)
  library(tidyverse)
  library(here)
})


# =============================================================================
# 01. Configuration
# =============================================================================

ALPHA <- 0.05

project_root <- here::here()

harmonization_dir <- here::here(
  "results",
  "meta_analysis",
  "harmonization"
)

input_file <- file.path(
  harmonization_dir,
  "meta_input_combined_long.csv"
)

results_root <- here::here(
  "results",
  "meta_analysis"
)

primary_dir <- file.path(
  results_root,
  "primary"
)

diagnostics_dir <- file.path(
  primary_dir,
  "diagnostics"
)

dir.create(
  primary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", project_root)
message("Input file: ", input_file)
message("Primary meta-analysis directory: ", primary_dir)


# =============================================================================
# 02. Validate harmonized input
# =============================================================================

if (!file.exists(input_file)) {
  stop(
    "Harmonized meta-analysis input not found:\n",
    input_file,
    "\nRun 01_harmonize_meta_inputs.R first."
  )
}

meta_long <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "Dataset",
  "SYMBOL",
  "log2FoldChange",
  "lfcSE"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_long)
)

if (length(missing_columns) > 0) {
  stop(
    "Harmonized meta-analysis input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

meta_long <- meta_long %>%
  mutate(
    Dataset = as.character(Dataset),
    SYMBOL = as.character(SYMBOL),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE)
  )

invalid_rows <- (
  is.na(meta_long$SYMBOL) |
    meta_long$SYMBOL == "" |
    is.na(meta_long$log2FoldChange) |
    !is.finite(meta_long$log2FoldChange) |
    is.na(meta_long$lfcSE) |
    !is.finite(meta_long$lfcSE) |
    meta_long$lfcSE <= 0
)

if (any(invalid_rows)) {
  stop(
    "Invalid rows were found in meta_input_combined_long.csv. ",
    "Run and inspect the harmonization diagnostics before meta-analysis."
  )
}

if (anyDuplicated(meta_long[, c("SYMBOL", "Dataset")]) > 0) {
  stop(
    "Duplicated SYMBOL-Dataset combinations were found in the harmonized input."
  )
}


# =============================================================================
# 03. Gene coverage
# =============================================================================

gene_coverage <- meta_long %>%
  distinct(
    SYMBOL,
    Dataset
  ) %>%
  count(
    SYMBOL,
    name = "k"
  )

meta_eligible_genes <- gene_coverage %>%
  filter(
    k >= 2
  ) %>%
  arrange(
    SYMBOL
  )

genes_for_meta <- meta_eligible_genes$SYMBOL

write.csv(
  meta_eligible_genes,
  file.path(
    diagnostics_dir,
    "genes_eligible_for_meta_analysis.csv"
  ),
  row.names = FALSE
)

coverage_summary <- meta_eligible_genes %>%
  count(
    k,
    name = "N_genes"
  ) %>%
  arrange(
    k
  )

write.csv(
  coverage_summary,
  file.path(
    diagnostics_dir,
    "meta_analysis_gene_coverage_summary.csv"
  ),
  row.names = FALSE
)

print(coverage_summary)

if (length(genes_for_meta) == 0) {
  stop(
    "No genes are available in at least two datasets."
  )
}


# =============================================================================
# 04. Preallocate result table
# =============================================================================

meta_results <- data.frame(
  SYMBOL = genes_for_meta,
  k = NA_integer_,

  GSE172274_log2FC = NA_real_,
  GSE172274_SE = NA_real_,

  GSE179277_log2FC = NA_real_,
  GSE179277_SE = NA_real_,

  GSE231409_log2FC = NA_real_,
  GSE231409_SE = NA_real_,

  RE_log2FC = NA_real_,
  RE_SE = NA_real_,
  RE_CI_lower = NA_real_,
  RE_CI_upper = NA_real_,
  RE_z = NA_real_,
  RE_pvalue = NA_real_,

  FE_log2FC = NA_real_,
  FE_SE = NA_real_,
  FE_CI_lower = NA_real_,
  FE_CI_upper = NA_real_,
  FE_z = NA_real_,
  FE_pvalue = NA_real_,

  tau2 = NA_real_,
  I2 = NA_real_,
  Q = NA_real_,
  Q_df = NA_real_,
  Q_pvalue = NA_real_,

  direction_pattern = NA_character_,
  direction_consistent = NA,

  stringsAsFactors = FALSE
)


# =============================================================================
# 05. Gene-by-gene random-effects meta-analysis
# =============================================================================

for (i in seq_along(genes_for_meta)) {

  current_gene <- genes_for_meta[i]

  gene_data <- meta_long %>%
    filter(
      SYMBOL == current_gene
    ) %>%
    arrange(
      Dataset
    )

  k_gene <- nrow(gene_data)

  meta_results$k[i] <- k_gene


  # ---------------------------------------------------------------------------
  # Store cohort-level estimates
  # ---------------------------------------------------------------------------

  row_172274 <- gene_data[
    gene_data$Dataset == "GSE172274",
    ,
    drop = FALSE
  ]

  if (nrow(row_172274) == 1) {
    meta_results$GSE172274_log2FC[i] <- row_172274$log2FoldChange
    meta_results$GSE172274_SE[i] <- row_172274$lfcSE
  }

  row_179277 <- gene_data[
    gene_data$Dataset == "GSE179277",
    ,
    drop = FALSE
  ]

  if (nrow(row_179277) == 1) {
    meta_results$GSE179277_log2FC[i] <- row_179277$log2FoldChange
    meta_results$GSE179277_SE[i] <- row_179277$lfcSE
  }

  row_231409 <- gene_data[
    gene_data$Dataset == "GSE231409",
    ,
    drop = FALSE
  ]

  if (nrow(row_231409) == 1) {
    meta_results$GSE231409_log2FC[i] <- row_231409$log2FoldChange
    meta_results$GSE231409_SE[i] <- row_231409$lfcSE
  }


  # ---------------------------------------------------------------------------
  # Direction pattern
  # ---------------------------------------------------------------------------

  direction_172274 <- ifelse(
    is.na(meta_results$GSE172274_log2FC[i]),
    "NA",
    ifelse(
      meta_results$GSE172274_log2FC[i] > 0,
      "+",
      ifelse(
        meta_results$GSE172274_log2FC[i] < 0,
        "-",
        "0"
      )
    )
  )

  direction_179277 <- ifelse(
    is.na(meta_results$GSE179277_log2FC[i]),
    "NA",
    ifelse(
      meta_results$GSE179277_log2FC[i] > 0,
      "+",
      ifelse(
        meta_results$GSE179277_log2FC[i] < 0,
        "-",
        "0"
      )
    )
  )

  direction_231409 <- ifelse(
    is.na(meta_results$GSE231409_log2FC[i]),
    "NA",
    ifelse(
      meta_results$GSE231409_log2FC[i] > 0,
      "+",
      ifelse(
        meta_results$GSE231409_log2FC[i] < 0,
        "-",
        "0"
      )
    )
  )

  direction_pattern <- paste0(
    direction_172274,
    direction_179277,
    direction_231409
  )

  meta_results$direction_pattern[i] <- direction_pattern

  available_signs <- sign(
    gene_data$log2FoldChange
  )

  available_signs <- available_signs[
    available_signs != 0
  ]

  if (
    length(available_signs) >= 2
  ) {

    meta_results$direction_consistent[i] <- (
      length(
        unique(
          available_signs
        )
      ) == 1
    )

  } else {

    meta_results$direction_consistent[i] <- NA
  }


  # ---------------------------------------------------------------------------
  # Random-effects model, REML
  # ---------------------------------------------------------------------------

  fit_re <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )


  # ---------------------------------------------------------------------------
  # Fixed-effect sensitivity model
  # ---------------------------------------------------------------------------

  fit_fe <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "FE",
      test = "z"
    ),
    error = function(e) NULL
  )


  # ---------------------------------------------------------------------------
  # Save random-effects model outputs
  # ---------------------------------------------------------------------------

  if (!is.null(fit_re)) {

    meta_results$RE_log2FC[i] <- as.numeric(
      fit_re$b[1]
    )

    meta_results$RE_SE[i] <- as.numeric(
      fit_re$se
    )

    meta_results$RE_CI_lower[i] <- as.numeric(
      fit_re$ci.lb
    )

    meta_results$RE_CI_upper[i] <- as.numeric(
      fit_re$ci.ub
    )

    meta_results$RE_z[i] <- as.numeric(
      fit_re$zval
    )

    meta_results$RE_pvalue[i] <- as.numeric(
      fit_re$pval
    )

    meta_results$tau2[i] <- as.numeric(
      fit_re$tau2
    )

    meta_results$I2[i] <- as.numeric(
      fit_re$I2
    )

    meta_results$Q[i] <- as.numeric(
      fit_re$QE
    )

    meta_results$Q_df[i] <- as.numeric(
      fit_re$k - fit_re$p
    )

    meta_results$Q_pvalue[i] <- as.numeric(
      fit_re$QEp
    )
  }


  # ---------------------------------------------------------------------------
  # Save fixed-effect model outputs
  # ---------------------------------------------------------------------------

  if (!is.null(fit_fe)) {

    meta_results$FE_log2FC[i] <- as.numeric(
      fit_fe$b[1]
    )

    meta_results$FE_SE[i] <- as.numeric(
      fit_fe$se
    )

    meta_results$FE_CI_lower[i] <- as.numeric(
      fit_fe$ci.lb
    )

    meta_results$FE_CI_upper[i] <- as.numeric(
      fit_fe$ci.ub
    )

    meta_results$FE_z[i] <- as.numeric(
      fit_fe$zval
    )

    meta_results$FE_pvalue[i] <- as.numeric(
      fit_fe$pval
    )
  }


  if (
    i %% 1000 == 0 ||
      i == length(genes_for_meta)
  ) {
    message(
      "Meta-analyzed ",
      i,
      " / ",
      length(genes_for_meta),
      " genes."
    )
  }
}


# =============================================================================
# 06. Remove failed meta-analysis fits
# =============================================================================

failed_re <- meta_results %>%
  filter(
    is.na(RE_pvalue)
  )

write.csv(
  failed_re,
  file.path(
    diagnostics_dir,
    "failed_random_effects_models.csv"
  ),
  row.names = FALSE
)

meta_results_valid <- meta_results %>%
  filter(
    !is.na(RE_pvalue),
    is.finite(RE_pvalue)
  )

if (nrow(meta_results_valid) == 0) {
  stop(
    "No valid random-effects meta-analysis results were produced."
  )
}


# =============================================================================
# 07. Multiple-testing correction
# =============================================================================

# Primary FDR is calculated once across the complete universe of genes that
# produced a valid meta-analysis result with k >= 2.

meta_results_valid$RE_FDR <- p.adjust(
  meta_results_valid$RE_pvalue,
  method = "BH"
)

meta_results_valid$FE_FDR <- p.adjust(
  meta_results_valid$FE_pvalue,
  method = "BH"
)


# =============================================================================
# 08. Derived interpretation fields
# =============================================================================

meta_results_valid <- meta_results_valid %>%
  mutate(
    RE_direction = case_when(
      RE_log2FC > 0 ~ "Pediatric",
      RE_log2FC < 0 ~ "Adult",
      TRUE ~ "Zero"
    ),

    RE_significant_FDR05 = (
      RE_FDR < ALPHA
    ),

    FE_significant_FDR05 = (
      FE_FDR < ALPHA
    ),

    Heterogeneity_Q_significant = (
      Q_pvalue < ALPHA
    ),

    I2_category = case_when(
      is.na(I2) ~ NA_character_,
      I2 < 25 ~ "Low",
      I2 < 50 ~ "Low_to_moderate",
      I2 < 75 ~ "Moderate_to_high",
      TRUE ~ "High"
    ),

    RE_FE_same_direction = (
      sign(RE_log2FC) ==
        sign(FE_log2FC)
    )
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )


# =============================================================================
# 09. Export complete primary meta-analysis table
# =============================================================================

write.csv(
  meta_results_valid,
  file.path(
    primary_dir,
    "meta_analysis_REM_all_genes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 10. Export genes represented in all three datasets
# =============================================================================

meta_results_k3 <- meta_results_valid %>%
  filter(
    k == 3
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  meta_results_k3,
  file.path(
    primary_dir,
    "meta_analysis_REM_k3_genes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. Export genes represented in exactly two datasets
# =============================================================================

meta_results_k2 <- meta_results_valid %>%
  filter(
    k == 2
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  meta_results_k2,
  file.path(
    primary_dir,
    "meta_analysis_REM_k2_genes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 12. Export significant random-effects results
# =============================================================================

meta_significant <- meta_results_valid %>%
  filter(
    RE_FDR < ALPHA
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  meta_significant,
  file.path(
    primary_dir,
    "meta_analysis_significant_FDR05.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 13. Export significant results with direction consistency
# =============================================================================

meta_significant_consistent <- meta_results_valid %>%
  filter(
    RE_FDR < ALPHA,
    direction_consistent %in% TRUE
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  meta_significant_consistent,
  file.path(
    primary_dir,
    "meta_analysis_significant_FDR05_direction_consistent.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 14. Export significant k=3 results
# =============================================================================

meta_significant_k3 <- meta_results_valid %>%
  filter(
    RE_FDR < ALPHA,
    k == 3
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  meta_significant_k3,
  file.path(
    primary_dir,
    "meta_analysis_significant_FDR05_k3.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. Export heterogeneity-focused tables
# =============================================================================

meta_high_heterogeneity <- meta_results_valid %>%
  filter(
    !is.na(I2),
    I2 >= 75
  ) %>%
  arrange(
    desc(I2),
    RE_FDR
  )

write.csv(
  meta_high_heterogeneity,
  file.path(
    diagnostics_dir,
    "meta_analysis_high_heterogeneity_I2_ge_75.csv"
  ),
  row.names = FALSE
)

meta_significant_heterogeneous <- meta_results_valid %>%
  filter(
    RE_FDR < ALPHA,
    (
      I2 >= 75 |
        Q_pvalue < ALPHA
    )
  ) %>%
  arrange(
    RE_FDR,
    desc(I2)
  )

write.csv(
  meta_significant_heterogeneous,
  file.path(
    diagnostics_dir,
    "significant_genes_with_heterogeneity.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 16. Random-effects vs fixed-effect comparison summary
# =============================================================================

re_fe_comparison_summary <- data.frame(
  Metric = c(
    "Genes analyzed",
    "RE significant FDR < 0.05",
    "FE significant FDR < 0.05",
    "Significant in both RE and FE",
    "RE-significant with same RE/FE direction",
    "RE-significant and direction-consistent across available datasets",
    "RE-significant k=3",
    "RE-significant k=2"
  ),
  N = c(
    nrow(meta_results_valid),
    sum(
      meta_results_valid$RE_significant_FDR05,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$FE_significant_FDR05,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_significant_FDR05 &
        meta_results_valid$FE_significant_FDR05,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_significant_FDR05 &
        meta_results_valid$RE_FE_same_direction,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_significant_FDR05 &
        meta_results_valid$direction_consistent %in% TRUE,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_significant_FDR05 &
        meta_results_valid$k == 3,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_significant_FDR05 &
        meta_results_valid$k == 2,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  re_fe_comparison_summary,
  file.path(
    diagnostics_dir,
    "RE_vs_FE_summary.csv"
  ),
  row.names = FALSE
)

print(re_fe_comparison_summary)


# =============================================================================
# 17. Heterogeneity summary
# =============================================================================

heterogeneity_summary <- data.frame(
  Metric = c(
    "Genes with I2 < 25",
    "Genes with 25 <= I2 < 50",
    "Genes with 50 <= I2 < 75",
    "Genes with I2 >= 75",
    "Genes with significant Cochran Q",
    "RE-significant genes with I2 >= 75",
    "RE-significant genes with significant Cochran Q"
  ),
  N = c(
    sum(
      meta_results_valid$I2 < 25,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$I2 >= 25 &
        meta_results_valid$I2 < 50,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$I2 >= 50 &
        meta_results_valid$I2 < 75,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$I2 >= 75,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$Q_pvalue < ALPHA,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_FDR < ALPHA &
        meta_results_valid$I2 >= 75,
      na.rm = TRUE
    ),
    sum(
      meta_results_valid$RE_FDR < ALPHA &
        meta_results_valid$Q_pvalue < ALPHA,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  heterogeneity_summary,
  file.path(
    diagnostics_dir,
    "heterogeneity_summary.csv"
  ),
  row.names = FALSE
)

print(heterogeneity_summary)


# =============================================================================
# 18. Direction-pattern summary
# =============================================================================

direction_pattern_summary <- meta_results_valid %>%
  count(
    k,
    direction_pattern,
    direction_consistent,
    name = "N_genes"
  ) %>%
  arrange(
    desc(k),
    desc(N_genes)
  )

write.csv(
  direction_pattern_summary,
  file.path(
    diagnostics_dir,
    "direction_pattern_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 19. Top results table
# =============================================================================

top_results <- meta_results_valid %>%
  dplyr::select(
    SYMBOL,
    k,

    GSE172274_log2FC,
    GSE179277_log2FC,
    GSE231409_log2FC,

    RE_log2FC,
    RE_SE,
    RE_CI_lower,
    RE_CI_upper,
    RE_pvalue,
    RE_FDR,

    tau2,
    I2,
    Q_pvalue,

    FE_log2FC,
    FE_pvalue,
    FE_FDR,

    direction_pattern,
    direction_consistent
  ) %>%
  slice_head(
    n = min(
      100,
      nrow(meta_results_valid)
    )
  )

write.csv(
  top_results,
  file.path(
    primary_dir,
    "meta_analysis_top100_by_RE_FDR.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 20. Primary meta-analysis summary
# =============================================================================

n_valid <- nrow(
  meta_results_valid
)

n_k3 <- sum(
  meta_results_valid$k == 3,
  na.rm = TRUE
)

n_k2 <- sum(
  meta_results_valid$k == 2,
  na.rm = TRUE
)

n_significant <- sum(
  meta_results_valid$RE_FDR < ALPHA,
  na.rm = TRUE
)

n_significant_consistent <- sum(
  meta_results_valid$RE_FDR < ALPHA &
    meta_results_valid$direction_consistent %in% TRUE,
  na.rm = TRUE
)

n_significant_k3 <- sum(
  meta_results_valid$RE_FDR < ALPHA &
    meta_results_valid$k == 3,
  na.rm = TRUE
)

n_significant_k2 <- sum(
  meta_results_valid$RE_FDR < ALPHA &
    meta_results_valid$k == 2,
  na.rm = TRUE
)

n_significant_high_i2 <- sum(
  meta_results_valid$RE_FDR < ALPHA &
    meta_results_valid$I2 >= 75,
  na.rm = TRUE
)

meta_summary <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  paste(
    "Project root:",
    project_root
  ),
  "",
  "Primary meta-analysis:",
  "  Effect measure: unshrunken DESeq2 log2FoldChange",
  "  Standard error: DESeq2 lfcSE",
  "  Direction: Pediatric vs Adult",
  "  Positive pooled effect: Pediatric-enriched",
  "  Model: random effects",
  "  Between-study variance estimator: REML",
  "  Test statistic: z",
  "  Multiple-testing correction: Benjamini-Hochberg",
  "  Primary FDR universe: all genes with valid k >= 2 meta-analysis",
  "",
  "Sensitivity model:",
  "  Fixed-effect inverse-variance model",
  "",
  paste(
    "Genes successfully meta-analyzed:",
    n_valid
  ),
  paste(
    "Genes with k=3:",
    n_k3
  ),
  paste(
    "Genes with k=2:",
    n_k2
  ),
  paste(
    "RE FDR < 0.05:",
    n_significant
  ),
  paste(
    "RE FDR < 0.05 and direction-consistent:",
    n_significant_consistent
  ),
  paste(
    "RE FDR < 0.05 and k=3:",
    n_significant_k3
  ),
  paste(
    "RE FDR < 0.05 and k=2:",
    n_significant_k2
  ),
  paste(
    "RE FDR < 0.05 and I2 >= 75:",
    n_significant_high_i2
  ),
  "",
  "Important interpretation:",
  "  Direction consistency and heterogeneity are descriptive robustness",
  "  characteristics and are not used as pre-analysis significance filters.",
  "",
  "Primary output:",
  "  meta_analysis_REM_all_genes.csv"
)

writeLines(
  meta_summary,
  file.path(
    primary_dir,
    "meta_analysis_summary.txt"
  )
)

cat(
  paste(
    meta_summary,
    collapse = "\n"
  ),
  "\n"
)


# =============================================================================
# 21. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_meta_analysis.txt"
  )
)

writeLines(
  c(
    paste(
      "R version:",
      R.version.string
    ),
    paste(
      "metafor version:",
      as.character(
        packageVersion(
          "metafor"
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
  ),
  file.path(
    diagnostics_dir,
    "package_versions_meta_analysis.txt"
  )
)


# =============================================================================
# 22. Final message
# =============================================================================

message("")
message("============================================================")
message("Random-effects meta-analysis completed successfully.")
message("============================================================")

message(
  "Genes successfully meta-analyzed: ",
  n_valid
)

message(
  "k = 3 genes: ",
  n_k3
)

message(
  "k = 2 genes: ",
  n_k2
)

message(
  "RE FDR < 0.05: ",
  n_significant
)

message(
  "RE FDR < 0.05 + direction consistent: ",
  n_significant_consistent
)

message(
  "RE FDR < 0.05 + k = 3: ",
  n_significant_k3
)

message(
  "RE FDR < 0.05 + I2 >= 75: ",
  n_significant_high_i2
)

message(
  "Results saved to: ",
  primary_dir
)

message("============================================================")
