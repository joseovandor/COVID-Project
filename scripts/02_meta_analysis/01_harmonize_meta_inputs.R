# =============================================================================
# 01_harmonize_meta_inputs.R
#
# Harmonization of cohort-level differential-expression inputs for meta-analysis
#
# Datasets
#   - GSE172274
#   - GSE179277
#   - GSE231409
#
# Purpose
#   This script imports the cohort-level DESeq2 outputs prepared specifically
#   for meta-analysis, validates their structure, confirms effect-direction
#   conventions, audits gene identifiers, and creates harmonized long- and
#   wide-format tables for downstream random-effects meta-analysis.
#
# Important conventions
#   - All cohort-level effects must represent Pediatric vs Adult.
#   - Positive log2FoldChange = higher expression in Pediatric samples.
#   - Inputs must contain UNSHRUNKEN DESeq2 log2FoldChange and lfcSE.
#   - No p-value, FDR, direction, or effect-size threshold is applied here.
#   - Gene harmonization is performed by SYMBOL because the source datasets use
#     different primary identifiers (ENTREZID or ENSEMBL).
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project. Paths are resolved
#   with here::here(); no machine-specific working directory is required.
#
# The script is intentionally linear and does not define custom functions.
# =============================================================================

SCRIPT_BUILD <- "META_HARMONIZATION_2026-08-17_v1"
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

dataset_ids <- c(
  "GSE172274",
  "GSE179277",
  "GSE231409"
)

expected_columns <- c(
  "Dataset",
  "Original_ID",
  "SYMBOL",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "MetaEligibleAnnotation"
)

project_root <- here::here()

input_files <- c(
  GSE172274 = here::here(
    "results",
    "deseq2",
    "GSE172274",
    "primary",
    "GSE172274_meta_input.csv"
  ),
  GSE179277 = here::here(
    "results",
    "deseq2",
    "GSE179277",
    "primary",
    "GSE179277_meta_input.csv"
  ),
  GSE231409 = here::here(
    "results",
    "deseq2",
    "GSE231409",
    "primary",
    "GSE231409_meta_input.csv"
  )
)

results_root <- here::here(
  "results",
  "meta_analysis"
)

harmonization_dir <- file.path(
  results_root,
  "harmonization"
)

diagnostics_dir <- file.path(
  harmonization_dir,
  "diagnostics"
)

dir.create(
  results_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  harmonization_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", project_root)
message("Harmonization directory: ", harmonization_dir)


# =============================================================================
# 02. Validate input files
# =============================================================================

missing_input_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "The following meta-analysis input files were not found:\n",
    paste(
      paste0(
        names(missing_input_files),
        ": ",
        missing_input_files
      ),
      collapse = "\n"
    )
  )
}

message("All three cohort-level meta-analysis inputs were found.")


# =============================================================================
# 03. Read GSE172274 meta-analysis input
# =============================================================================

meta_GSE172274 <- read.csv(
  input_files["GSE172274"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns_GSE172274 <- setdiff(
  expected_columns,
  colnames(meta_GSE172274)
)

if (length(missing_columns_GSE172274) > 0) {
  stop(
    "GSE172274 meta input is missing required columns: ",
    paste(missing_columns_GSE172274, collapse = ", ")
  )
}

if (
  any(
    meta_GSE172274$Dataset != "GSE172274",
    na.rm = TRUE
  )
) {
  stop(
    "Unexpected Dataset values were found in GSE172274_meta_input.csv."
  )
}

meta_GSE172274$Dataset <- "GSE172274"


# =============================================================================
# 04. Read GSE179277 meta-analysis input
# =============================================================================

meta_GSE179277 <- read.csv(
  input_files["GSE179277"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns_GSE179277 <- setdiff(
  expected_columns,
  colnames(meta_GSE179277)
)

if (length(missing_columns_GSE179277) > 0) {
  stop(
    "GSE179277 meta input is missing required columns: ",
    paste(missing_columns_GSE179277, collapse = ", ")
  )
}

if (
  any(
    meta_GSE179277$Dataset != "GSE179277",
    na.rm = TRUE
  )
) {
  stop(
    "Unexpected Dataset values were found in GSE179277_meta_input.csv."
  )
}

meta_GSE179277$Dataset <- "GSE179277"


# =============================================================================
# 05. Read GSE231409 meta-analysis input
# =============================================================================

meta_GSE231409 <- read.csv(
  input_files["GSE231409"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns_GSE231409 <- setdiff(
  expected_columns,
  colnames(meta_GSE231409)
)

if (length(missing_columns_GSE231409) > 0) {
  stop(
    "GSE231409 meta input is missing required columns: ",
    paste(missing_columns_GSE231409, collapse = ", ")
  )
}

if (
  any(
    meta_GSE231409$Dataset != "GSE231409",
    na.rm = TRUE
  )
) {
  stop(
    "Unexpected Dataset values were found in GSE231409_meta_input.csv."
  )
}

meta_GSE231409$Dataset <- "GSE231409"


# =============================================================================
# 06. Standardize column types
# =============================================================================

meta_GSE172274 <- meta_GSE172274 %>%
  mutate(
    Dataset = as.character(Dataset),
    Original_ID = as.character(Original_ID),
    SYMBOL = as.character(SYMBOL),
    baseMean = as.numeric(baseMean),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE),
    stat = as.numeric(stat),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj),
    MetaEligibleAnnotation = as.logical(MetaEligibleAnnotation)
  )

meta_GSE179277 <- meta_GSE179277 %>%
  mutate(
    Dataset = as.character(Dataset),
    Original_ID = as.character(Original_ID),
    SYMBOL = as.character(SYMBOL),
    baseMean = as.numeric(baseMean),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE),
    stat = as.numeric(stat),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj),
    MetaEligibleAnnotation = as.logical(MetaEligibleAnnotation)
  )

meta_GSE231409 <- meta_GSE231409 %>%
  mutate(
    Dataset = as.character(Dataset),
    Original_ID = as.character(Original_ID),
    SYMBOL = as.character(SYMBOL),
    baseMean = as.numeric(baseMean),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE),
    stat = as.numeric(stat),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj),
    MetaEligibleAnnotation = as.logical(MetaEligibleAnnotation)
  )


# =============================================================================
# 07. Cohort-level validation
# =============================================================================

# No missing SYMBOLs should remain in the meta_input files.
if (any(is.na(meta_GSE172274$SYMBOL) | meta_GSE172274$SYMBOL == "")) {
  stop("GSE172274 meta input contains missing SYMBOL values.")
}

if (any(is.na(meta_GSE179277$SYMBOL) | meta_GSE179277$SYMBOL == "")) {
  stop("GSE179277 meta input contains missing SYMBOL values.")
}

if (any(is.na(meta_GSE231409$SYMBOL) | meta_GSE231409$SYMBOL == "")) {
  stop("GSE231409 meta input contains missing SYMBOL values.")
}

# SYMBOLs must be unique within each dataset.
if (anyDuplicated(meta_GSE172274$SYMBOL) > 0) {
  stop("Duplicated SYMBOL values were found in GSE172274 meta input.")
}

if (anyDuplicated(meta_GSE179277$SYMBOL) > 0) {
  stop("Duplicated SYMBOL values were found in GSE179277 meta input.")
}

if (anyDuplicated(meta_GSE231409$SYMBOL) > 0) {
  stop("Duplicated SYMBOL values were found in GSE231409 meta input.")
}

# MetaEligibleAnnotation should remain TRUE for every exported row.
if (!all(meta_GSE172274$MetaEligibleAnnotation, na.rm = TRUE)) {
  stop("GSE172274 contains rows not marked as MetaEligibleAnnotation.")
}

if (!all(meta_GSE179277$MetaEligibleAnnotation, na.rm = TRUE)) {
  stop("GSE179277 contains rows not marked as MetaEligibleAnnotation.")
}

if (!all(meta_GSE231409$MetaEligibleAnnotation, na.rm = TRUE)) {
  stop("GSE231409 contains rows not marked as MetaEligibleAnnotation.")
}

# The effect and standard error required for inverse-variance meta-analysis
# must be finite and the SE must be strictly positive.
invalid_effect_GSE172274 <- (
  is.na(meta_GSE172274$log2FoldChange) |
    !is.finite(meta_GSE172274$log2FoldChange) |
    is.na(meta_GSE172274$lfcSE) |
    !is.finite(meta_GSE172274$lfcSE) |
    meta_GSE172274$lfcSE <= 0
)

invalid_effect_GSE179277 <- (
  is.na(meta_GSE179277$log2FoldChange) |
    !is.finite(meta_GSE179277$log2FoldChange) |
    is.na(meta_GSE179277$lfcSE) |
    !is.finite(meta_GSE179277$lfcSE) |
    meta_GSE179277$lfcSE <= 0
)

invalid_effect_GSE231409 <- (
  is.na(meta_GSE231409$log2FoldChange) |
    !is.finite(meta_GSE231409$log2FoldChange) |
    is.na(meta_GSE231409$lfcSE) |
    !is.finite(meta_GSE231409$lfcSE) |
    meta_GSE231409$lfcSE <= 0
)

if (any(invalid_effect_GSE172274)) {
  stop(
    "Invalid log2FoldChange and/or lfcSE values were found in GSE172274."
  )
}

if (any(invalid_effect_GSE179277)) {
  stop(
    "Invalid log2FoldChange and/or lfcSE values were found in GSE179277."
  )
}

if (any(invalid_effect_GSE231409)) {
  stop(
    "Invalid log2FoldChange and/or lfcSE values were found in GSE231409."
  )
}


# =============================================================================
# 08. Export per-cohort validation summary
# =============================================================================

cohort_validation_summary <- data.frame(
  Dataset = dataset_ids,
  Rows = c(
    nrow(meta_GSE172274),
    nrow(meta_GSE179277),
    nrow(meta_GSE231409)
  ),
  Unique_SYMBOLs = c(
    length(unique(meta_GSE172274$SYMBOL)),
    length(unique(meta_GSE179277$SYMBOL)),
    length(unique(meta_GSE231409$SYMBOL))
  ),
  Missing_SYMBOL = c(
    sum(is.na(meta_GSE172274$SYMBOL) | meta_GSE172274$SYMBOL == ""),
    sum(is.na(meta_GSE179277$SYMBOL) | meta_GSE179277$SYMBOL == ""),
    sum(is.na(meta_GSE231409$SYMBOL) | meta_GSE231409$SYMBOL == "")
  ),
  Duplicated_SYMBOL = c(
    sum(duplicated(meta_GSE172274$SYMBOL)),
    sum(duplicated(meta_GSE179277$SYMBOL)),
    sum(duplicated(meta_GSE231409$SYMBOL))
  ),
  Invalid_effect_or_SE = c(
    sum(invalid_effect_GSE172274),
    sum(invalid_effect_GSE179277),
    sum(invalid_effect_GSE231409)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cohort_validation_summary,
  file.path(
    diagnostics_dir,
    "cohort_validation_summary.csv"
  ),
  row.names = FALSE
)

print(cohort_validation_summary)


# =============================================================================
# 09. Combine datasets in long format
# =============================================================================

meta_long <- bind_rows(
  meta_GSE172274,
  meta_GSE179277,
  meta_GSE231409
) %>%
  arrange(
    SYMBOL,
    Dataset
  )

write.csv(
  meta_long,
  file.path(
    harmonization_dir,
    "meta_input_combined_long.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 10. Gene coverage across datasets
# =============================================================================

gene_coverage <- meta_long %>%
  distinct(
    SYMBOL,
    Dataset
  ) %>%
  count(
    SYMBOL,
    name = "k_available"
  ) %>%
  arrange(
    desc(k_available),
    SYMBOL
  )

coverage_distribution <- gene_coverage %>%
  count(
    k_available,
    name = "N_genes"
  ) %>%
  arrange(
    k_available
  )

write.csv(
  gene_coverage,
  file.path(
    harmonization_dir,
    "gene_coverage_by_dataset.csv"
  ),
  row.names = FALSE
)

write.csv(
  coverage_distribution,
  file.path(
    diagnostics_dir,
    "gene_coverage_distribution.csv"
  ),
  row.names = FALSE
)

print(coverage_distribution)


# =============================================================================
# 11. Identify genes present in all three datasets
# =============================================================================

genes_k3 <- gene_coverage %>%
  filter(
    k_available == 3
  ) %>%
  pull(
    SYMBOL
  )

genes_present_in_3 <- meta_long %>%
  filter(
    SYMBOL %in% genes_k3
  ) %>%
  arrange(
    SYMBOL,
    Dataset
  )

write.csv(
  genes_present_in_3,
  file.path(
    harmonization_dir,
    "genes_present_in_3_datasets.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 12. Identify genes present in exactly two datasets
# =============================================================================

genes_k2 <- gene_coverage %>%
  filter(
    k_available == 2
  ) %>%
  pull(
    SYMBOL
  )

genes_present_in_2 <- meta_long %>%
  filter(
    SYMBOL %in% genes_k2
  ) %>%
  arrange(
    SYMBOL,
    Dataset
  )

write.csv(
  genes_present_in_2,
  file.path(
    harmonization_dir,
    "genes_present_in_2_datasets.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 13. Identify genes unique to one dataset
# =============================================================================

genes_k1 <- gene_coverage %>%
  filter(
    k_available == 1
  ) %>%
  pull(
    SYMBOL
  )

genes_present_in_1 <- meta_long %>%
  filter(
    SYMBOL %in% genes_k1
  ) %>%
  arrange(
    SYMBOL,
    Dataset
  )

write.csv(
  genes_present_in_1,
  file.path(
    harmonization_dir,
    "genes_present_in_1_dataset.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 14. Create harmonized wide table
# =============================================================================

meta_wide <- meta_long %>%
  dplyr::select(
    SYMBOL,
    Dataset,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = c(
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    ),
    names_glue = "{Dataset}_{.value}"
  ) %>%
  left_join(
    gene_coverage,
    by = "SYMBOL"
  ) %>%
  arrange(
    desc(k_available),
    SYMBOL
  )

write.csv(
  meta_wide,
  file.path(
    harmonization_dir,
    "meta_input_harmonized_wide.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. Direction-of-effect audit
# =============================================================================

direction_audit <- meta_wide %>%
  mutate(
    GSE172274_direction = case_when(
      is.na(GSE172274_log2FoldChange) ~ NA_character_,
      GSE172274_log2FoldChange > 0 ~ "Pediatric",
      GSE172274_log2FoldChange < 0 ~ "Adult",
      TRUE ~ "Zero"
    ),
    GSE179277_direction = case_when(
      is.na(GSE179277_log2FoldChange) ~ NA_character_,
      GSE179277_log2FoldChange > 0 ~ "Pediatric",
      GSE179277_log2FoldChange < 0 ~ "Adult",
      TRUE ~ "Zero"
    ),
    GSE231409_direction = case_when(
      is.na(GSE231409_log2FoldChange) ~ NA_character_,
      GSE231409_log2FoldChange > 0 ~ "Pediatric",
      GSE231409_log2FoldChange < 0 ~ "Adult",
      TRUE ~ "Zero"
    )
  )

direction_audit$direction_pattern <- apply(
  direction_audit[
    ,
    c(
      "GSE172274_direction",
      "GSE179277_direction",
      "GSE231409_direction"
    )
  ],
  1,
  function(x) {
    paste(
      ifelse(
        is.na(x),
        "NA",
        ifelse(
          x == "Pediatric",
          "+",
          ifelse(
            x == "Adult",
            "-",
            "0"
          )
        )
      ),
      collapse = ""
    )
  }
)

write.csv(
  direction_audit,
  file.path(
    harmonization_dir,
    "direction_of_effect_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 16. Shared-gene direction consistency summary
# =============================================================================

direction_summary_k3 <- direction_audit %>%
  filter(
    k_available == 3
  ) %>%
  mutate(
    DirectionConsistent = (
      direction_pattern %in% c(
        "+++",
        "---"
      )
    )
  ) %>%
  count(
    direction_pattern,
    DirectionConsistent,
    name = "N_genes"
  ) %>%
  arrange(
    desc(N_genes)
  )

write.csv(
  direction_summary_k3,
  file.path(
    diagnostics_dir,
    "direction_consistency_summary_k3.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 17. Dataset pair overlap summary
# =============================================================================

symbols_172274 <- unique(
  meta_GSE172274$SYMBOL
)

symbols_179277 <- unique(
  meta_GSE179277$SYMBOL
)

symbols_231409 <- unique(
  meta_GSE231409$SYMBOL
)

pair_overlap_summary <- data.frame(
  Comparison = c(
    "GSE172274_vs_GSE179277",
    "GSE172274_vs_GSE231409",
    "GSE179277_vs_GSE231409",
    "All_3_datasets"
  ),
  N_shared_SYMBOLs = c(
    length(
      intersect(
        symbols_172274,
        symbols_179277
      )
    ),
    length(
      intersect(
        symbols_172274,
        symbols_231409
      )
    ),
    length(
      intersect(
        symbols_179277,
        symbols_231409
      )
    ),
    length(
      Reduce(
        intersect,
        list(
          symbols_172274,
          symbols_179277,
          symbols_231409
        )
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pair_overlap_summary,
  file.path(
    diagnostics_dir,
    "dataset_pair_overlap_summary.csv"
  ),
  row.names = FALSE
)

print(pair_overlap_summary)


# =============================================================================
# 18. Final harmonization summary
# =============================================================================

n_total_unique_symbols <- length(
  unique(
    meta_long$SYMBOL
  )
)

n_k3 <- sum(
  gene_coverage$k_available == 3
)

n_k2 <- sum(
  gene_coverage$k_available == 2
)

n_k1 <- sum(
  gene_coverage$k_available == 1
)

n_k3_consistent <- direction_audit %>%
  filter(
    k_available == 3,
    direction_pattern %in% c(
      "+++",
      "---"
    )
  ) %>%
  nrow()

harmonization_summary <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  paste(
    "Project root:",
    project_root
  ),
  "",
  "Meta-analysis convention:",
  "  Contrast in all datasets: Pediatric vs Adult",
  "  Positive log2FoldChange: higher expression in Pediatric samples",
  "  Input effects: unshrunken DESeq2 log2FoldChange",
  "  Input uncertainty: DESeq2 lfcSE",
  "  Harmonization key: SYMBOL",
  "  Statistical preselection before harmonization: none",
  "",
  "Cohort-level input sizes:",
  paste(
    "  GSE172274:",
    nrow(meta_GSE172274),
    "genes"
  ),
  paste(
    "  GSE179277:",
    nrow(meta_GSE179277),
    "genes"
  ),
  paste(
    "  GSE231409:",
    nrow(meta_GSE231409),
    "genes"
  ),
  "",
  paste(
    "Unique SYMBOLs across all datasets:",
    n_total_unique_symbols
  ),
  paste(
    "Genes available in all 3 datasets (k=3):",
    n_k3
  ),
  paste(
    "Genes available in exactly 2 datasets (k=2):",
    n_k2
  ),
  paste(
    "Genes available in exactly 1 dataset (k=1):",
    n_k1
  ),
  paste(
    "k=3 genes with fully consistent direction (+++ or ---):",
    n_k3_consistent
  ),
  "",
  "Outputs:",
  "  meta_input_combined_long.csv",
  "  meta_input_harmonized_wide.csv",
  "  gene_coverage_by_dataset.csv",
  "  genes_present_in_3_datasets.csv",
  "  genes_present_in_2_datasets.csv",
  "  genes_present_in_1_dataset.csv",
  "  direction_of_effect_audit.csv"
)

writeLines(
  harmonization_summary,
  file.path(
    harmonization_dir,
    "harmonization_summary.txt"
  )
)

cat(
  paste(
    harmonization_summary,
    collapse = "\n"
  ),
  "\n"
)


# =============================================================================
# 19. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_harmonization.txt"
  )
)

writeLines(
  c(
    paste(
      "R version:",
      R.version.string
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
    "package_versions_harmonization.txt"
  )
)


# =============================================================================
# 20. Final message
# =============================================================================

message("")
message("============================================================")
message("Meta-analysis input harmonization completed successfully.")
message("============================================================")
message(
  "Unique SYMBOLs across all datasets: ",
  n_total_unique_symbols
)
message(
  "Genes with k = 3: ",
  n_k3
)
message(
  "Genes with k = 2: ",
  n_k2
)
message(
  "Genes with k = 1: ",
  n_k1
)
message(
  "k = 3 genes with consistent direction: ",
  n_k3_consistent
)
message(
  "Results saved to: ",
  harmonization_dir
)
message("============================================================")
