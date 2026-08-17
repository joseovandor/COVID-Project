# =============================================================================
# 04_meta_analysis_figures.R
#
# Publication-ready figures and reporting tables for the transcriptomic
# random-effects meta-analysis.
#
# Inputs
#   results/meta_analysis/primary/meta_analysis_REM_all_genes.csv
#   results/meta_analysis/sensitivity/robustness_classification.csv
#   results/meta_analysis/sensitivity/leave_one_out_results.csv
#
# Figures
#   1. RE versus FE pooled effects
#   2. Pooled random-effects estimate versus I2
#   3. Direction-consistency summary among significant genes
#   4. Forest plots for primary significant genes
#   5. Leave-one-dataset-out plots for significant k=3 genes
#
# Tables
#   - Primary significant genes
#   - Robust genes
#   - Complete meta-analysis table
#
# Important
#   This script performs visualization and reporting only. It does not redefine
#   statistical significance, recalculate FDR, or select genes using new
#   inferential criteria.
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project.
# =============================================================================

SCRIPT_BUILD <- "META_FIGURES_2026-08-17_v1"
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

project_root <- here::here()

primary_file <- here::here(
  "results",
  "meta_analysis",
  "primary",
  "meta_analysis_REM_all_genes.csv"
)

robustness_file <- here::here(
  "results",
  "meta_analysis",
  "sensitivity",
  "robustness_classification.csv"
)

leave_one_out_file <- here::here(
  "results",
  "meta_analysis",
  "sensitivity",
  "leave_one_out_results.csv"
)

results_root <- here::here(
  "results",
  "meta_analysis"
)

figures_dir <- file.path(
  results_root,
  "figures"
)

forest_dir <- file.path(
  figures_dir,
  "forest"
)

leave_one_out_figures_dir <- file.path(
  figures_dir,
  "leave_one_out"
)

tables_dir <- file.path(
  results_root,
  "tables"
)

diagnostics_dir <- file.path(
  figures_dir,
  "diagnostics"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  forest_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  leave_one_out_figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", project_root)
message("Figures directory: ", figures_dir)
message("Tables directory: ", tables_dir)


# =============================================================================
# 02. Validate inputs
# =============================================================================

if (!file.exists(primary_file)) {
  stop(
    "Primary meta-analysis file not found:\n",
    primary_file,
    "\nRun 02_random_effects_meta_analysis.R first."
  )
}

if (!file.exists(robustness_file)) {
  stop(
    "Robustness classification file not found:\n",
    robustness_file,
    "\nRun 03_meta_analysis_robustness.R first."
  )
}

if (!file.exists(leave_one_out_file)) {
  stop(
    "Leave-one-out results file not found:\n",
    leave_one_out_file,
    "\nRun 03_meta_analysis_robustness.R first."
  )
}


# =============================================================================
# 03. Read results
# =============================================================================

meta_results <- read.csv(
  primary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

robustness_results <- read.csv(
  robustness_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

leave_one_out <- read.csv(
  leave_one_out_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta_columns <- c(
  "SYMBOL",
  "k",
  "GSE172274_log2FC",
  "GSE172274_SE",
  "GSE179277_log2FC",
  "GSE179277_SE",
  "GSE231409_log2FC",
  "GSE231409_SE",
  "RE_log2FC",
  "RE_SE",
  "RE_CI_lower",
  "RE_CI_upper",
  "RE_pvalue",
  "RE_FDR",
  "FE_log2FC",
  "FE_SE",
  "FE_CI_lower",
  "FE_CI_upper",
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

required_robustness_columns <- c(
  "SYMBOL",
  "PrimarySignificant",
  "RobustnessClass"
)

missing_robustness_columns <- setdiff(
  required_robustness_columns,
  colnames(robustness_results)
)

if (length(missing_robustness_columns) > 0) {
  stop(
    "Robustness table is missing required columns: ",
    paste(missing_robustness_columns, collapse = ", ")
  )
}

required_loo_columns <- c(
  "SYMBOL",
  "Omitted_dataset",
  "Included_datasets",
  "Pooled_log2FC",
  "CI_lower",
  "CI_upper",
  "pvalue"
)

missing_loo_columns <- setdiff(
  required_loo_columns,
  colnames(leave_one_out)
)

if (length(missing_loo_columns) > 0) {
  stop(
    "Leave-one-out table is missing required columns: ",
    paste(missing_loo_columns, collapse = ", ")
  )
}


# =============================================================================
# 04. Standardize types
# =============================================================================

meta_results <- meta_results %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    k = as.integer(k),
    RE_log2FC = as.numeric(RE_log2FC),
    RE_SE = as.numeric(RE_SE),
    RE_CI_lower = as.numeric(RE_CI_lower),
    RE_CI_upper = as.numeric(RE_CI_upper),
    RE_pvalue = as.numeric(RE_pvalue),
    RE_FDR = as.numeric(RE_FDR),
    FE_log2FC = as.numeric(FE_log2FC),
    FE_SE = as.numeric(FE_SE),
    FE_CI_lower = as.numeric(FE_CI_lower),
    FE_CI_upper = as.numeric(FE_CI_upper),
    FE_pvalue = as.numeric(FE_pvalue),
    FE_FDR = as.numeric(FE_FDR),
    tau2 = as.numeric(tau2),
    I2 = as.numeric(I2),
    Q_pvalue = as.numeric(Q_pvalue),
    direction_pattern = as.character(direction_pattern),
    direction_consistent = as.logical(direction_consistent)
  )

robustness_results <- robustness_results %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    PrimarySignificant = as.logical(PrimarySignificant),
    RobustnessClass = as.character(RobustnessClass)
  )

leave_one_out <- leave_one_out %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    Omitted_dataset = as.character(Omitted_dataset),
    Included_datasets = as.character(Included_datasets),
    Pooled_log2FC = as.numeric(Pooled_log2FC),
    CI_lower = as.numeric(CI_lower),
    CI_upper = as.numeric(CI_upper),
    pvalue = as.numeric(pvalue)
  )


# =============================================================================
# 05. Merge primary and robustness results
# =============================================================================

plot_results <- meta_results %>%
  left_join(
    robustness_results %>%
      dplyr::select(
        SYMBOL,
        PrimarySignificant,
        RobustnessClass
      ),
    by = "SYMBOL"
  )

plot_results$PrimarySignificant[
  is.na(plot_results$PrimarySignificant)
] <- plot_results$RE_FDR[
  is.na(plot_results$PrimarySignificant)
] < ALPHA

plot_results$RobustnessClass[
  is.na(plot_results$RobustnessClass)
] <- "Not_classified"


# =============================================================================
# 06. Export manuscript-oriented tables
# =============================================================================

table_significant <- plot_results %>%
  filter(
    RE_FDR < ALPHA
  ) %>%
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
    direction_pattern,
    direction_consistent,
    FE_log2FC,
    FE_pvalue,
    FE_FDR,
    RobustnessClass
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  table_significant,
  file.path(
    tables_dir,
    "Table_meta_significant_genes.csv"
  ),
  row.names = FALSE
)

table_robust <- plot_results %>%
  filter(
    RE_FDR < ALPHA,
    RobustnessClass == "Robust"
  ) %>%
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
    direction_pattern,
    direction_consistent,
    FE_log2FC,
    FE_pvalue,
    FE_FDR,
    RobustnessClass
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  )

write.csv(
  table_robust,
  file.path(
    tables_dir,
    "Table_meta_robust_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_results,
  file.path(
    tables_dir,
    "Supplementary_meta_all_genes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 07. Common publication theme
# =============================================================================

publication_theme <- theme_classic(
  base_size = 12
) +
  theme(
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    )
  )


# =============================================================================
# 08. Figure 1: Random-effects versus fixed-effect pooled estimates
# =============================================================================

re_fe_plot_data <- plot_results %>%
  filter(
    is.finite(RE_log2FC),
    is.finite(FE_log2FC)
  ) %>%
  mutate(
    Significance = ifelse(
      RE_FDR < ALPHA,
      "RE FDR < 0.05",
      "RE FDR >= 0.05"
    )
  )

p_re_fe <- ggplot(
  re_fe_plot_data,
  aes(
    x = FE_log2FC,
    y = RE_log2FC,
    shape = Significance
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    alpha = 0.65,
    size = 2
  ) +
  labs(
    title = "Random-effects versus fixed-effect estimates",
    x = "Fixed-effect pooled log2 fold change",
    y = "Random-effects pooled log2 fold change",
    shape = NULL
  ) +
  publication_theme

ggsave(
  filename = file.path(
    figures_dir,
    "meta_RE_vs_FE.tiff"
  ),
  plot = p_re_fe,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    figures_dir,
    "meta_RE_vs_FE.pdf"
  ),
  plot = p_re_fe,
  width = 7,
  height = 6,
  units = "in"
)


# =============================================================================
# 09. Figure 2: Random-effects estimate versus heterogeneity
# =============================================================================

heterogeneity_plot_data <- plot_results %>%
  filter(
    is.finite(RE_log2FC),
    is.finite(I2)
  ) %>%
  mutate(
    Significance = ifelse(
      RE_FDR < ALPHA,
      "RE FDR < 0.05",
      "RE FDR >= 0.05"
    )
  )

p_heterogeneity <- ggplot(
  heterogeneity_plot_data,
  aes(
    x = RE_log2FC,
    y = I2,
    shape = Significance
  )
) +
  geom_hline(
    yintercept = 75,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    alpha = 0.65,
    size = 2
  ) +
  labs(
    title = "Pooled effect and between-study heterogeneity",
    x = "Random-effects pooled log2 fold change",
    y = expression(I^2~"(%)"),
    shape = NULL
  ) +
  publication_theme

ggsave(
  filename = file.path(
    figures_dir,
    "meta_effect_vs_I2.tiff"
  ),
  plot = p_heterogeneity,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    figures_dir,
    "meta_effect_vs_I2.pdf"
  ),
  plot = p_heterogeneity,
  width = 7,
  height = 6,
  units = "in"
)


# =============================================================================
# 10. Figure 3: Direction patterns among significant genes
# =============================================================================

direction_plot_data <- plot_results %>%
  filter(
    RE_FDR < ALPHA
  ) %>%
  count(
    direction_pattern,
    name = "N_genes"
  ) %>%
  arrange(
    N_genes
  )

if (nrow(direction_plot_data) > 0) {

  direction_plot_data$direction_pattern <- factor(
    direction_plot_data$direction_pattern,
    levels = direction_plot_data$direction_pattern
  )

  p_direction <- ggplot(
    direction_plot_data,
    aes(
      x = direction_pattern,
      y = N_genes
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    coord_flip() +
    labs(
      title = "Direction patterns among significant meta-analysis genes",
      subtitle = "Order: GSE172274, GSE179277, GSE231409",
      x = "Direction pattern",
      y = "Number of genes"
    ) +
    publication_theme

  ggsave(
    filename = file.path(
      figures_dir,
      "meta_direction_consistency.tiff"
    ),
    plot = p_direction,
    width = 7,
    height = 5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      figures_dir,
      "meta_direction_consistency.pdf"
    ),
    plot = p_direction,
    width = 7,
    height = 5,
    units = "in"
  )
}


# =============================================================================
# 11. Prepare genes for forest plots
# =============================================================================

forest_genes <- plot_results %>%
  filter(
    RE_FDR < ALPHA
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  ) %>%
  pull(
    SYMBOL
  )

message(
  "Significant genes selected for forest plots: ",
  length(forest_genes)
)


# =============================================================================
# 12. Forest plots for primary significant genes
# =============================================================================

for (current_gene in forest_genes) {

  current <- plot_results %>%
    filter(
      SYMBOL == current_gene
    ) %>%
    slice(
      1
    )

  forest_data <- data.frame(
    Study = c(
      "GSE172274",
      "GSE179277",
      "GSE231409",
      "Random effects"
    ),
    Estimate = c(
      current$GSE172274_log2FC,
      current$GSE179277_log2FC,
      current$GSE231409_log2FC,
      current$RE_log2FC
    ),
    SE = c(
      current$GSE172274_SE,
      current$GSE179277_SE,
      current$GSE231409_SE,
      current$RE_SE
    ),
    stringsAsFactors = FALSE
  )

  forest_data <- forest_data %>%
    mutate(
      CI_lower = Estimate - 1.96 * SE,
      CI_upper = Estimate + 1.96 * SE
    )

  forest_data$CI_lower[
    forest_data$Study == "Random effects"
  ] <- current$RE_CI_lower

  forest_data$CI_upper[
    forest_data$Study == "Random effects"
  ] <- current$RE_CI_upper

  forest_data <- forest_data %>%
    filter(
      is.finite(Estimate),
      is.finite(CI_lower),
      is.finite(CI_upper)
    )

  forest_data$Study <- factor(
    forest_data$Study,
    levels = rev(
      c(
        "GSE172274",
        "GSE179277",
        "GSE231409",
        "Random effects"
      )
    )
  )

  p_forest <- ggplot(
    forest_data,
    aes(
      x = Estimate,
      y = Study
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      height = 0.15,
      linewidth = 0.7
    ) +
    geom_point(
      size = 3
    ) +
    labs(
      title = current_gene,
      subtitle = paste0(
        "RE FDR = ",
        format(
          current$RE_FDR,
          digits = 3,
          scientific = TRUE
        ),
        " | I2 = ",
        round(
          current$I2,
          1
        ),
        "%"
      ),
      x = "log2 fold change (Pediatric vs Adult)",
      y = NULL
    ) +
    publication_theme

  safe_gene <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    current_gene
  )

  ggsave(
    filename = file.path(
      forest_dir,
      paste0(
        "forest_",
        safe_gene,
        ".tiff"
      )
    ),
    plot = p_forest,
    width = 7,
    height = 4.5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      forest_dir,
      paste0(
        "forest_",
        safe_gene,
        ".pdf"
      )
    ),
    plot = p_forest,
    width = 7,
    height = 4.5,
    units = "in"
  )
}


# =============================================================================
# 13. Prepare significant k=3 genes for leave-one-out plots
# =============================================================================

loo_genes <- plot_results %>%
  filter(
    RE_FDR < ALPHA,
    k == 3
  ) %>%
  arrange(
    RE_FDR,
    RE_pvalue
  ) %>%
  pull(
    SYMBOL
  )

message(
  "Significant k=3 genes selected for leave-one-out plots: ",
  length(loo_genes)
)


# =============================================================================
# 14. Leave-one-dataset-out plots
# =============================================================================

for (current_gene in loo_genes) {

  current_full <- plot_results %>%
    filter(
      SYMBOL == current_gene
    ) %>%
    slice(
      1
    )

  current_loo <- leave_one_out %>%
    filter(
      SYMBOL == current_gene
    )

  if (nrow(current_loo) == 0) {
    next
  }

  loo_plot_data <- current_loo %>%
    transmute(
      Analysis = paste0(
        "Omit ",
        Omitted_dataset
      ),
      Estimate = Pooled_log2FC,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )

  full_row <- data.frame(
    Analysis = "Full REM",
    Estimate = current_full$RE_log2FC,
    CI_lower = current_full$RE_CI_lower,
    CI_upper = current_full$RE_CI_upper,
    stringsAsFactors = FALSE
  )

  loo_plot_data <- bind_rows(
    full_row,
    loo_plot_data
  )

  loo_plot_data$Analysis <- factor(
    loo_plot_data$Analysis,
    levels = rev(
      loo_plot_data$Analysis
    )
  )

  p_loo <- ggplot(
    loo_plot_data,
    aes(
      x = Estimate,
      y = Analysis
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      height = 0.15,
      linewidth = 0.7
    ) +
    geom_point(
      size = 3
    ) +
    labs(
      title = paste0(
        current_gene,
        ": leave-one-dataset-out"
      ),
      subtitle = paste0(
        "Full-model RE FDR = ",
        format(
          current_full$RE_FDR,
          digits = 3,
          scientific = TRUE
        )
      ),
      x = "Pooled log2 fold change (Pediatric vs Adult)",
      y = NULL
    ) +
    publication_theme

  safe_gene <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    current_gene
  )

  ggsave(
    filename = file.path(
      leave_one_out_figures_dir,
      paste0(
        "leave_one_out_",
        safe_gene,
        ".tiff"
      )
    ),
    plot = p_loo,
    width = 7,
    height = 4.5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      leave_one_out_figures_dir,
      paste0(
        "leave_one_out_",
        safe_gene,
        ".pdf"
      )
    ),
    plot = p_loo,
    width = 7,
    height = 4.5,
    units = "in"
  )
}


# =============================================================================
# 15. Export figure-generation summary
# =============================================================================

n_all <- nrow(
  plot_results
)

n_significant <- sum(
  plot_results$RE_FDR < ALPHA,
  na.rm = TRUE
)

n_significant_k3 <- sum(
  plot_results$RE_FDR < ALPHA &
    plot_results$k == 3,
  na.rm = TRUE
)

n_robust <- sum(
  plot_results$RE_FDR < ALPHA &
    plot_results$RobustnessClass == "Robust",
  na.rm = TRUE
)

n_supported <- sum(
  plot_results$RE_FDR < ALPHA &
    plot_results$RobustnessClass == "Supported",
  na.rm = TRUE
)

n_heterogeneous <- sum(
  plot_results$RE_FDR < ALPHA &
    plot_results$RobustnessClass == "Heterogeneous",
  na.rm = TRUE
)

figure_summary <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  paste(
    "Project root:",
    project_root
  ),
  "",
  "Reporting convention:",
  "  Primary significance: RE_FDR < 0.05",
  "  Effect direction: Pediatric vs Adult",
  "  Positive effect: Pediatric-enriched",
  "  No statistical significance criterion was recalculated in this script.",
  "",
  paste(
    "Genes in complete meta-analysis table:",
    n_all
  ),
  paste(
    "Primary significant genes:",
    n_significant
  ),
  paste(
    "Primary significant k=3 genes:",
    n_significant_k3
  ),
  paste(
    "Robust significant genes:",
    n_robust
  ),
  paste(
    "Supported significant genes:",
    n_supported
  ),
  paste(
    "Heterogeneous significant genes:",
    n_heterogeneous
  ),
  "",
  paste(
    "Forest plots generated:",
    length(forest_genes)
  ),
  paste(
    "Leave-one-out plot candidates:",
    length(loo_genes)
  ),
  "",
  "Main figures:",
  "  meta_RE_vs_FE.tiff",
  "  meta_effect_vs_I2.tiff",
  "  meta_direction_consistency.tiff",
  "",
  "Per-gene figures:",
  "  figures/forest/",
  "  figures/leave_one_out/",
  "",
  "Reporting tables:",
  "  Table_meta_significant_genes.csv",
  "  Table_meta_robust_genes.csv",
  "  Supplementary_meta_all_genes.csv"
)

writeLines(
  figure_summary,
  file.path(
    figures_dir,
    "figure_summary.txt"
  )
)

cat(
  paste(
    figure_summary,
    collapse = "\n"
  ),
  "\n"
)


# =============================================================================
# 16. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_figures.txt"
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
    "package_versions_figures.txt"
  )
)


# =============================================================================
# 17. Final message
# =============================================================================

message("")
message("============================================================")
message("Meta-analysis figure generation completed successfully.")
message("============================================================")

message(
  "Primary significant genes: ",
  n_significant
)

message(
  "Forest plots generated: ",
  length(forest_genes)
)

message(
  "Significant k=3 leave-one-out plot candidates: ",
  length(loo_genes)
)

message(
  "Figures saved to: ",
  figures_dir
)

message(
  "Tables saved to: ",
  tables_dir
)

message("============================================================")
