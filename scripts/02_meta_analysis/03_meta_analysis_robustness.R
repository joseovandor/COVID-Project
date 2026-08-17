# =============================================================================
# 03_meta_analysis_robustness.R
#
# Robustness analysis for the random-effects meta-analysis
#
# Input
#   results/meta_analysis/primary/meta_analysis_REM_all_genes.csv
#   results/meta_analysis/harmonization/meta_input_combined_long.csv
#
# Purpose
#   This script evaluates whether statistically significant pooled effects are
#   supported consistently across studies or are driven by one cohort.
#
# Robustness components
#   1. Direction-of-effect consistency across available datasets
#   2. Leave-one-dataset-out sensitivity for genes represented in all 3 studies
#   3. Stability of pooled effect direction and nominal significance
#   4. Heterogeneity-aware descriptive classification
#
# Important
#   - Primary statistical significance remains RE_FDR < 0.05 from script 02.
#   - This script does NOT recalculate the primary FDR.
#   - Robustness labels are descriptive summaries, not replacement significance
#     criteria.
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project. Paths are resolved
#   with here::here(); no machine-specific working directory is required.
#
# The script is intentionally linear and does not define custom functions.
# =============================================================================

SCRIPT_BUILD <- "META_ROBUSTNESS_2026-08-17_v1"
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
HIGH_I2 <- 75

project_root <- here::here()

meta_results_file <- here::here(
  "results",
  "meta_analysis",
  "primary",
  "meta_analysis_REM_all_genes.csv"
)

meta_long_file <- here::here(
  "results",
  "meta_analysis",
  "harmonization",
  "meta_input_combined_long.csv"
)

results_root <- here::here(
  "results",
  "meta_analysis"
)

sensitivity_dir <- file.path(
  results_root,
  "sensitivity"
)

diagnostics_dir <- file.path(
  sensitivity_dir,
  "diagnostics"
)

dir.create(
  sensitivity_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", project_root)
message("Robustness directory: ", sensitivity_dir)


# =============================================================================
# 02. Validate inputs
# =============================================================================

if (!file.exists(meta_results_file)) {
  stop(
    "Primary meta-analysis results not found:\n",
    meta_results_file,
    "\nRun 02_random_effects_meta_analysis.R first."
  )
}

if (!file.exists(meta_long_file)) {
  stop(
    "Harmonized long-format meta input not found:\n",
    meta_long_file,
    "\nRun 01_harmonize_meta_inputs.R first."
  )
}

meta_results <- read.csv(
  meta_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_long <- read.csv(
  meta_long_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta_columns <- c(
  "SYMBOL",
  "k",
  "RE_log2FC",
  "RE_pvalue",
  "RE_FDR",
  "FE_log2FC",
  "FE_pvalue",
  "FE_FDR",
  "tau2",
  "I2",
  "Q_pvalue",
  "direction_pattern",
  "direction_consistent"
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  colnames(meta_results)
)

if (length(missing_meta_columns) > 0) {
  stop(
    "Primary meta-analysis table is missing required columns: ",
    paste(missing_meta_columns, collapse = ", ")
  )
}

required_long_columns <- c(
  "Dataset",
  "SYMBOL",
  "log2FoldChange",
  "lfcSE"
)

missing_long_columns <- setdiff(
  required_long_columns,
  colnames(meta_long)
)

if (length(missing_long_columns) > 0) {
  stop(
    "Harmonized long-format input is missing required columns: ",
    paste(missing_long_columns, collapse = ", ")
  )
}

meta_results <- meta_results %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    k = as.integer(k),
    RE_log2FC = as.numeric(RE_log2FC),
    RE_pvalue = as.numeric(RE_pvalue),
    RE_FDR = as.numeric(RE_FDR),
    FE_log2FC = as.numeric(FE_log2FC),
    FE_pvalue = as.numeric(FE_pvalue),
    FE_FDR = as.numeric(FE_FDR),
    tau2 = as.numeric(tau2),
    I2 = as.numeric(I2),
    Q_pvalue = as.numeric(Q_pvalue),
    direction_pattern = as.character(direction_pattern),
    direction_consistent = as.logical(direction_consistent)
  )

meta_long <- meta_long %>%
  mutate(
    Dataset = as.character(Dataset),
    SYMBOL = as.character(SYMBOL),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE)
  )


# =============================================================================
# 03. Significant genes from the primary random-effects model
# =============================================================================

significant_genes <- meta_results %>%
  filter(
    RE_FDR < ALPHA
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  significant_genes,
  file.path(
    sensitivity_dir,
    "significant_genes_primary_RE.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 04. Direction-consistency summary
# =============================================================================

direction_consistency_summary <- meta_results %>%
  mutate(
    PrimarySignificant = RE_FDR < ALPHA
  ) %>%
  count(
    k,
    PrimarySignificant,
    direction_pattern,
    direction_consistent,
    name = "N_genes"
  ) %>%
  arrange(
    desc(PrimarySignificant),
    desc(k),
    desc(N_genes)
  )

write.csv(
  direction_consistency_summary,
  file.path(
    diagnostics_dir,
    "direction_consistency_summary.csv"
  ),
  row.names = FALSE
)

significant_direction_summary <- significant_genes %>%
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
  significant_direction_summary,
  file.path(
    diagnostics_dir,
    "significant_direction_pattern_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 05. Prepare k=3 genes for leave-one-dataset-out analysis
# =============================================================================

genes_k3 <- meta_results %>%
  filter(
    k == 3
  ) %>%
  pull(
    SYMBOL
  )

meta_long_k3 <- meta_long %>%
  filter(
    SYMBOL %in% genes_k3
  )

if (length(genes_k3) == 0) {
  warning(
    "No k=3 genes are available. Leave-one-dataset-out analysis will be empty."
  )
}


# =============================================================================
# 06. Preallocate leave-one-out result table
# =============================================================================

leave_one_out <- data.frame(
  SYMBOL = character(),
  Omitted_dataset = character(),
  Included_datasets = character(),
  k = integer(),

  Pooled_log2FC = numeric(),
  Pooled_SE = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  z = numeric(),
  pvalue = numeric(),

  tau2 = numeric(),
  I2 = numeric(),
  Q = numeric(),
  Q_pvalue = numeric(),

  Same_direction_as_full_RE = logical(),
  Nominal_p_lt_0_05 = logical(),

  stringsAsFactors = FALSE
)


# =============================================================================
# 07. Leave one dataset out: omit GSE172274
# =============================================================================

for (current_gene in genes_k3) {

  gene_data <- meta_long_k3 %>%
    filter(
      SYMBOL == current_gene,
      Dataset != "GSE172274"
    ) %>%
    arrange(
      Dataset
    )

  fit <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )

  if (!is.null(fit)) {

    full_effect <- meta_results$RE_log2FC[
      meta_results$SYMBOL == current_gene
    ][1]

    leave_one_out <- bind_rows(
      leave_one_out,
      data.frame(
        SYMBOL = current_gene,
        Omitted_dataset = "GSE172274",
        Included_datasets = paste(
          sort(gene_data$Dataset),
          collapse = " + "
        ),
        k = nrow(gene_data),

        Pooled_log2FC = as.numeric(fit$b[1]),
        Pooled_SE = as.numeric(fit$se),
        CI_lower = as.numeric(fit$ci.lb),
        CI_upper = as.numeric(fit$ci.ub),
        z = as.numeric(fit$zval),
        pvalue = as.numeric(fit$pval),

        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2),
        Q = as.numeric(fit$QE),
        Q_pvalue = as.numeric(fit$QEp),

        Same_direction_as_full_RE = (
          sign(as.numeric(fit$b[1])) ==
            sign(full_effect)
        ),

        Nominal_p_lt_0_05 = (
          as.numeric(fit$pval) < ALPHA
        ),

        stringsAsFactors = FALSE
      )
    )
  }
}


# =============================================================================
# 08. Leave one dataset out: omit GSE179277
# =============================================================================

for (current_gene in genes_k3) {

  gene_data <- meta_long_k3 %>%
    filter(
      SYMBOL == current_gene,
      Dataset != "GSE179277"
    ) %>%
    arrange(
      Dataset
    )

  fit <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )

  if (!is.null(fit)) {

    full_effect <- meta_results$RE_log2FC[
      meta_results$SYMBOL == current_gene
    ][1]

    leave_one_out <- bind_rows(
      leave_one_out,
      data.frame(
        SYMBOL = current_gene,
        Omitted_dataset = "GSE179277",
        Included_datasets = paste(
          sort(gene_data$Dataset),
          collapse = " + "
        ),
        k = nrow(gene_data),

        Pooled_log2FC = as.numeric(fit$b[1]),
        Pooled_SE = as.numeric(fit$se),
        CI_lower = as.numeric(fit$ci.lb),
        CI_upper = as.numeric(fit$ci.ub),
        z = as.numeric(fit$zval),
        pvalue = as.numeric(fit$pval),

        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2),
        Q = as.numeric(fit$QE),
        Q_pvalue = as.numeric(fit$QEp),

        Same_direction_as_full_RE = (
          sign(as.numeric(fit$b[1])) ==
            sign(full_effect)
        ),

        Nominal_p_lt_0_05 = (
          as.numeric(fit$pval) < ALPHA
        ),

        stringsAsFactors = FALSE
      )
    )
  }
}


# =============================================================================
# 09. Leave one dataset out: omit GSE231409
# =============================================================================

for (current_gene in genes_k3) {

  gene_data <- meta_long_k3 %>%
    filter(
      SYMBOL == current_gene,
      Dataset != "GSE231409"
    ) %>%
    arrange(
      Dataset
    )

  fit <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )

  if (!is.null(fit)) {

    full_effect <- meta_results$RE_log2FC[
      meta_results$SYMBOL == current_gene
    ][1]

    leave_one_out <- bind_rows(
      leave_one_out,
      data.frame(
        SYMBOL = current_gene,
        Omitted_dataset = "GSE231409",
        Included_datasets = paste(
          sort(gene_data$Dataset),
          collapse = " + "
        ),
        k = nrow(gene_data),

        Pooled_log2FC = as.numeric(fit$b[1]),
        Pooled_SE = as.numeric(fit$se),
        CI_lower = as.numeric(fit$ci.lb),
        CI_upper = as.numeric(fit$ci.ub),
        z = as.numeric(fit$zval),
        pvalue = as.numeric(fit$pval),

        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2),
        Q = as.numeric(fit$QE),
        Q_pvalue = as.numeric(fit$QEp),

        Same_direction_as_full_RE = (
          sign(as.numeric(fit$b[1])) ==
            sign(full_effect)
        ),

        Nominal_p_lt_0_05 = (
          as.numeric(fit$pval) < ALPHA
        ),

        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(
  leave_one_out,
  file.path(
    sensitivity_dir,
    "leave_one_out_results.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 10. Summarize leave-one-out behavior by gene
# =============================================================================

if (nrow(leave_one_out) > 0) {

  leave_one_out_summary <- leave_one_out %>%
    group_by(
      SYMBOL
    ) %>%
    summarise(
      N_leave_one_out_models = n(),
      N_same_direction_as_full_RE = sum(
        Same_direction_as_full_RE,
        na.rm = TRUE
      ),
      N_nominal_p_lt_0_05 = sum(
        Nominal_p_lt_0_05,
        na.rm = TRUE
      ),
      Min_abs_pooled_log2FC = min(
        abs(Pooled_log2FC),
        na.rm = TRUE
      ),
      Max_abs_pooled_log2FC = max(
        abs(Pooled_log2FC),
        na.rm = TRUE
      ),
      Max_I2_leave_one_out = max(
        I2,
        na.rm = TRUE
      ),
      All_leave_one_out_same_direction = all(
        Same_direction_as_full_RE
      ),
      All_leave_one_out_nominal_significant = all(
        Nominal_p_lt_0_05
      ),
      .groups = "drop"
    )

} else {

  leave_one_out_summary <- data.frame(
    SYMBOL = character(),
    N_leave_one_out_models = integer(),
    N_same_direction_as_full_RE = integer(),
    N_nominal_p_lt_0_05 = integer(),
    Min_abs_pooled_log2FC = numeric(),
    Max_abs_pooled_log2FC = numeric(),
    Max_I2_leave_one_out = numeric(),
    All_leave_one_out_same_direction = logical(),
    All_leave_one_out_nominal_significant = logical(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  leave_one_out_summary,
  file.path(
    sensitivity_dir,
    "leave_one_out_summary_by_gene.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. Merge robustness information with primary results
# =============================================================================

robustness_table <- meta_results %>%
  left_join(
    leave_one_out_summary,
    by = "SYMBOL"
  )


# =============================================================================
# 12. Define descriptive robustness flags
# =============================================================================

robustness_table <- robustness_table %>%
  mutate(

    PrimarySignificant = (
      RE_FDR < ALPHA
    ),

    DirectionConsistent = (
      direction_consistent %in% TRUE
    ),

    HighHeterogeneity = (
      !is.na(I2) &
        I2 >= HIGH_I2
    ),

    SignificantQ = (
      !is.na(Q_pvalue) &
        Q_pvalue < ALPHA
    ),

    RE_FE_same_direction = (
      sign(RE_log2FC) ==
        sign(FE_log2FC)
    ),

    LeaveOneOutDirectionStable = case_when(
      k == 3 &
        !is.na(All_leave_one_out_same_direction) ~
        All_leave_one_out_same_direction,
      TRUE ~ NA
    ),

    LeaveOneOutNominallyStable = case_when(
      k == 3 &
        !is.na(All_leave_one_out_nominal_significant) ~
        All_leave_one_out_nominal_significant,
      TRUE ~ NA
    )
  )


# =============================================================================
# 13. Descriptive robustness classification
# =============================================================================

robustness_table <- robustness_table %>%
  mutate(

    RobustnessClass = case_when(

      !PrimarySignificant ~
        "Not_significant",

      PrimarySignificant &
        k == 3 &
        DirectionConsistent &
        !HighHeterogeneity &
        RE_FE_same_direction &
        LeaveOneOutDirectionStable %in% TRUE &
        LeaveOneOutNominallyStable %in% TRUE ~
        "Robust",

      PrimarySignificant &
        (
          HighHeterogeneity |
            !DirectionConsistent
        ) ~
        "Heterogeneous",

      PrimarySignificant ~
        "Supported",

      TRUE ~
        "Unclassified"
    )
  )

write.csv(
  robustness_table,
  file.path(
    sensitivity_dir,
    "robustness_classification.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 14. Export significant genes by robustness class
# =============================================================================

significant_robust <- robustness_table %>%
  filter(
    PrimarySignificant,
    RobustnessClass == "Robust"
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

significant_supported <- robustness_table %>%
  filter(
    PrimarySignificant,
    RobustnessClass == "Supported"
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

significant_heterogeneous <- robustness_table %>%
  filter(
    PrimarySignificant,
    RobustnessClass == "Heterogeneous"
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  significant_robust,
  file.path(
    sensitivity_dir,
    "significant_genes_robust.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_supported,
  file.path(
    sensitivity_dir,
    "significant_genes_supported.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_heterogeneous,
  file.path(
    sensitivity_dir,
    "significant_genes_heterogeneous.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. Dataset influence summary
# =============================================================================

if (nrow(leave_one_out) > 0) {

  dataset_influence_summary <- leave_one_out %>%
    group_by(
      Omitted_dataset
    ) %>%
    summarise(
      N_models = n(),
      N_same_direction_as_full_RE = sum(
        Same_direction_as_full_RE,
        na.rm = TRUE
      ),
      Proportion_same_direction = mean(
        Same_direction_as_full_RE,
        na.rm = TRUE
      ),
      N_nominal_p_lt_0_05 = sum(
        Nominal_p_lt_0_05,
        na.rm = TRUE
      ),
      Proportion_nominal_p_lt_0_05 = mean(
        Nominal_p_lt_0_05,
        na.rm = TRUE
      ),
      Median_abs_pooled_log2FC = median(
        abs(Pooled_log2FC),
        na.rm = TRUE
      ),
      Median_I2 = median(
        I2,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

} else {

  dataset_influence_summary <- data.frame(
    Omitted_dataset = character(),
    N_models = integer(),
    N_same_direction_as_full_RE = integer(),
    Proportion_same_direction = numeric(),
    N_nominal_p_lt_0_05 = integer(),
    Proportion_nominal_p_lt_0_05 = numeric(),
    Median_abs_pooled_log2FC = numeric(),
    Median_I2 = numeric(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  dataset_influence_summary,
  file.path(
    diagnostics_dir,
    "dataset_influence_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 16. Robustness-class summary
# =============================================================================

robustness_class_summary <- robustness_table %>%
  count(
    RobustnessClass,
    name = "N_genes"
  ) %>%
  arrange(
    desc(N_genes)
  )

write.csv(
  robustness_class_summary,
  file.path(
    diagnostics_dir,
    "robustness_class_summary.csv"
  ),
  row.names = FALSE
)

print(robustness_class_summary)


# =============================================================================
# 17. Significant-only robustness summary
# =============================================================================

significant_robustness_summary <- robustness_table %>%
  filter(
    PrimarySignificant
  ) %>%
  count(
    RobustnessClass,
    k,
    name = "N_genes"
  ) %>%
  arrange(
    RobustnessClass,
    desc(k)
  )

write.csv(
  significant_robustness_summary,
  file.path(
    diagnostics_dir,
    "significant_robustness_summary.csv"
  ),
  row.names = FALSE
)

print(significant_robustness_summary)


# =============================================================================
# 18. Primary-significant k=3 leave-one-out table
# =============================================================================

significant_k3_symbols <- meta_results %>%
  filter(
    RE_FDR < ALPHA,
    k == 3
  ) %>%
  pull(
    SYMBOL
  )

significant_k3_leave_one_out <- leave_one_out %>%
  filter(
    SYMBOL %in% significant_k3_symbols
  ) %>%
  arrange(
    SYMBOL,
    Omitted_dataset
  )

write.csv(
  significant_k3_leave_one_out,
  file.path(
    sensitivity_dir,
    "significant_k3_leave_one_out_results.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 19. Robustness summary text
# =============================================================================

n_primary_significant <- sum(
  robustness_table$PrimarySignificant,
  na.rm = TRUE
)

n_robust <- sum(
  robustness_table$RobustnessClass == "Robust",
  na.rm = TRUE
)

n_supported <- sum(
  robustness_table$RobustnessClass == "Supported",
  na.rm = TRUE
)

n_heterogeneous <- sum(
  robustness_table$RobustnessClass == "Heterogeneous",
  na.rm = TRUE
)

n_significant_k3 <- sum(
  robustness_table$PrimarySignificant &
    robustness_table$k == 3,
  na.rm = TRUE
)

n_significant_direction_consistent <- sum(
  robustness_table$PrimarySignificant &
    robustness_table$DirectionConsistent,
  na.rm = TRUE
)

n_significant_high_i2 <- sum(
  robustness_table$PrimarySignificant &
    robustness_table$HighHeterogeneity,
  na.rm = TRUE
)

n_significant_loo_direction_stable <- sum(
  robustness_table$PrimarySignificant &
    robustness_table$k == 3 &
    robustness_table$LeaveOneOutDirectionStable %in% TRUE,
  na.rm = TRUE
)

n_significant_loo_nominal_stable <- sum(
  robustness_table$PrimarySignificant &
    robustness_table$k == 3 &
    robustness_table$LeaveOneOutNominallyStable %in% TRUE,
  na.rm = TRUE
)

robustness_summary <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  paste(
    "Project root:",
    project_root
  ),
  "",
  "Primary statistical significance:",
  "  RE_FDR < 0.05 from 02_random_effects_meta_analysis.R",
  "  Primary FDR is not recalculated in this script.",
  "",
  "Robustness dimensions evaluated:",
  "  Direction consistency across available studies",
  "  Heterogeneity (I2 and Cochran Q)",
  "  Random-effects vs fixed-effect direction",
  "  Leave-one-dataset-out stability for k=3 genes",
  "",
  paste(
    "Primary significant genes:",
    n_primary_significant
  ),
  paste(
    "Primary significant k=3 genes:",
    n_significant_k3
  ),
  paste(
    "Primary significant + direction-consistent:",
    n_significant_direction_consistent
  ),
  paste(
    "Primary significant + I2 >= 75:",
    n_significant_high_i2
  ),
  paste(
    "Primary significant k=3 + stable direction in all leave-one-out models:",
    n_significant_loo_direction_stable
  ),
  paste(
    "Primary significant k=3 + nominal p < 0.05 in all leave-one-out models:",
    n_significant_loo_nominal_stable
  ),
  "",
  "Descriptive robustness classes:",
  paste(
    "  Robust:",
    n_robust
  ),
  paste(
    "  Supported:",
    n_supported
  ),
  paste(
    "  Heterogeneous:",
    n_heterogeneous
  ),
  "",
  "Important interpretation:",
  "  RobustnessClass is a descriptive synthesis of sensitivity evidence.",
  "  It does not replace the primary RE_FDR significance criterion."
)

writeLines(
  robustness_summary,
  file.path(
    sensitivity_dir,
    "robustness_summary.txt"
  )
)

cat(
  paste(
    robustness_summary,
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
    "sessionInfo_robustness.txt"
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
    "package_versions_robustness.txt"
  )
)


# =============================================================================
# 21. Final message
# =============================================================================

message("")
message("============================================================")
message("Meta-analysis robustness analysis completed successfully.")
message("============================================================")

message(
  "Primary significant genes: ",
  n_primary_significant
)

message(
  "Robust: ",
  n_robust
)

message(
  "Supported: ",
  n_supported
)

message(
  "Heterogeneous: ",
  n_heterogeneous
)

message(
  "Results saved to: ",
  sensitivity_dir
)

message("============================================================")
