# Clinical covariates and age-window sensitivity analyses
#
# Cohort models are refitted on the same filtered counts as the primary
# analysis (normalized counts x size factors) with the covariates available in
# each cohort, and every scenario is meta-analysed with the primary model
# (REML, z test, BH-FDR across the eligible genes of that scenario).
#
# Cohort models
#   GSE179277  A primary refit       ~ Sex + Group
#              B viral load          ~ Sex + log10(SARS-CoV-2 RPM + 1) + Group
#              C exclude 16-20 y     ~ Sex + Group, 16 <= age < 21 removed
#              D log2(age)           ~ Sex + log2(age + 1)
#   GSE231409  A primary refit       ~ Batch + Sex + Group
#              B0 VL subset, unadj.  ~ Batch + Sex + Group  (samples with VL)
#              B viral load          ~ Batch + Sex + log10(copies + 1) + Group
#              B2 VL + symptoms      ~ Batch + Sex + log10(copies + 1) + Symptomatic + Group
#              C exclude 16-20 y     ~ Batch + Sex + Group
#              D log2(age)           ~ Batch + Sex + log2(age + 1)
#   GSE172274  A primary refit       ~ Sex + Group
#              B severity            ~ Sex + Severe(score >= 3) + Group
#              B2 severe excluded    ~ Sex + Group, adults with score >= 3 removed
#              C exclude 16-20 y     ~ Sex + Group
#              D log2(age)           ~ Sex + log2(age + 1)
#
# Meta-analysis scenarios
#   S0_primary_refit           A + A + A  (should reproduce the primary analysis)
#   S1_viral_load_severity     B(179277) + B(231409) + B(172274)
#   S1b_VL_subset_unadjusted   A(179277) + B0(231409) + A(172274)
#   S2_extended_clinical       B(179277) + B2(231409) + B2(172274)
#   S3_exclude_age_16_20       C + C + C
#   S4_log2_age_continuous     D + D + D  (log2FC per doubling of age; negative
#                              values correspond to the pediatric direction)
#
# Viral load may lie between age and the host response, so these models are
# sensitivity analyses rather than replacements of the primary model.
#
# Output: results/sensitivity/clinical_covariates/

SCRIPT_BUILD <- "CLINICAL_COVARIATES_2026-10-26"
message("Running script build: ", SCRIPT_BUILD)

# DESeq2 (and its Bioconductor dependencies) is attached FIRST so that dplyr
# verbs (count, select, rename, first, ...) are not masked.
suppressPackageStartupMessages({
  library(DESeq2)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(metafor)
})

# Force dplyr verbs even if Bioconductor packages (AnnotationDbi, DESeq2,
# S4Vectors, ...) are attached in the session and mask them.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
count <- dplyr::count
first <- dplyr::first
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
distinct <- dplyr::distinct

ALPHA <- 0.05

out_dir <- here::here("results", "sensitivity", "clinical_covariates")
fig_dir <- file.path(out_dir, "figures")
model_dir <- file.path(out_dir, "cohort_models")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

primary <- read.csv(
  here::here("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"),
  stringsAsFactors = FALSE
)
signature <- primary %>% filter(RE_FDR < ALPHA)
cat("Signature genes:", nrow(signature), "\n")


# 01. Load filtered raw counts, metadata and annotation per cohort ----

load_cohort <- function(gse) {
  base <- here::here("results", "deseq2", gse)

  norm <- read.csv(file.path(base, "primary", paste0(gse, "_normalized_counts.csv")),
                   check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(norm[[1]])
  norm_mat <- as.matrix(norm[, -1, drop = FALSE])
  rownames(norm_mat) <- ids

  sf <- read.csv(file.path(base, "primary", paste0(gse, "_size_factors.csv")), stringsAsFactors = FALSE)
  sf_vec <- setNames(sf$SizeFactor, sf$Sample)
  if (!all(colnames(norm_mat) %in% names(sf_vec))) stop(gse, ": size factors do not match count columns")

  counts_num <- sweep(norm_mat, 2, sf_vec[colnames(norm_mat)], "*")
  max_dev <- max(abs(counts_num - round(counts_num)))
  counts <- round(counts_num)
  storage.mode(counts) <- "integer"

  meta <- read.csv(file.path(base, "metadata", paste0(gse, "_sample_metadata.csv")),
                   stringsAsFactors = FALSE, check.names = FALSE)
  key <- if (all(colnames(counts) %in% meta$Sample)) "Sample" else "GEO_ID"
  meta <- meta[match(colnames(counts), meta[[key]]), , drop = FALSE]
  if (any(is.na(meta$GEO_ID))) stop(gse, ": metadata could not be matched to count columns")

  raw <- read.csv(file.path(base, "metadata", paste0(gse, "_metadata_raw.csv")),
                  stringsAsFactors = FALSE, check.names = FALSE)
  raw <- raw[match(meta$GEO_ID, raw$geo_accession), , drop = FALSE]

  annot <- read.csv(file.path(base, "primary", paste0(gse, "_meta_input.csv")),
                    stringsAsFactors = FALSE, colClasses = c(Original_ID = "character"))
  annot <- annot %>%
    filter(MetaEligibleAnnotation %in% c(TRUE, "TRUE")) %>%
    select(Original_ID, SYMBOL)

  primary_res <- read.csv(file.path(base, "primary", paste0(gse, "_meta_input.csv")),
                          stringsAsFactors = FALSE, colClasses = c(Original_ID = "character"))

  meta$Group <- factor(meta$Group, levels = c("Adult", "Pediatric"))
  meta$Sex <- factor(meta$Sex)
  rownames(meta) <- colnames(counts)

  message(gse, ": ", nrow(counts), " genes x ", ncol(counts),
          " samples; max deviation of reconstructed counts from integers = ", signif(max_dev, 3))

  list(gse = gse, counts = counts, meta = meta, raw = raw, annot = annot,
       primary_res = primary_res, max_dev = max_dev)
}

cohorts <- lapply(c("GSE172274", "GSE179277", "GSE231409"), load_cohort)
names(cohorts) <- c("GSE172274", "GSE179277", "GSE231409")

# --- Clinical covariates -------------------------------------------------------

as_num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "NA", "N/A", "na"), NA, x)))

cohorts$GSE179277$meta$ViralLoad <- as_num(cohorts$GSE179277$raw[["sars-cov2_rpm:ch1"]])
cohorts$GSE179277$meta$ViralLoad_unit <- "SARS-CoV-2 reads per million"

cohorts$GSE231409$meta$ViralLoad <- as_num(cohorts$GSE231409$raw[["vl copies:ch1"]])
cohorts$GSE231409$meta$ViralLoad_unit <- "SARS-CoV-2 copies"
cohorts$GSE231409$meta$Symptomatic <- factor(
  ifelse(cohorts$GSE231409$raw[["covid cat:ch1"]] == "POS_SX", "Yes",
         ifelse(cohorts$GSE231409$raw[["covid cat:ch1"]] == "POS_ASX", "No", NA)),
  levels = c("No", "Yes")
)
cohorts$GSE231409$meta$VaccineDoses <- as_num(cohorts$GSE231409$raw[["vaccine doses:ch1"]])
cohorts$GSE231409$meta$DaysSymptoms <- as_num(cohorts$GSE231409$raw[["timing sx:ch1"]])
cohorts$GSE231409$meta$DaysDiagnosis <- as_num(cohorts$GSE231409$raw[["timing dx:ch1"]])
cohorts$GSE231409$meta$Hospitalized <- cohorts$GSE231409$raw[["hospitalized:ch1"]]
cohorts$GSE231409$meta$Batch <- factor(cohorts$GSE231409$meta$Batch)

cohorts$GSE172274$meta$RespiratoryScore <- as_num(cohorts$GSE172274$meta$RespiratoryScore)
cohorts$GSE172274$meta$Severe <- factor(
  ifelse(cohorts$GSE172274$meta$RespiratoryScore >= 3, "Yes", "No"),
  levels = c("No", "Yes")
)

for (g in c("GSE179277", "GSE231409")) {
  vl <- cohorts[[g]]$meta$ViralLoad
  cohorts[[g]]$meta$log10VL <- log10(vl + 1)
  cohorts[[g]]$meta$log10VL_z <- as.numeric(scale(log10(vl + 1)))
}

for (g in names(cohorts)) {
  cohorts[[g]]$meta$log2Age <- log2(cohorts[[g]]$meta$Age + 1)
  cohorts[[g]]$meta$log2Age_z <- as.numeric(scale(cohorts[[g]]$meta$log2Age))
}


# 02. Covariate availability and distribution by age group ----

availability <- bind_rows(
  data.frame(Cohort = "GSE179277", Variable = c("Age", "Sex", "Sequencing platform", "Viral load (RPM)", "Respiratory severity score",
                                                  "Symptom status", "Days from symptom onset", "Vaccination", "Hospitalization", "Batch"),
             N_available = c(sum(!is.na(cohorts$GSE179277$meta$Age)), sum(!is.na(cohorts$GSE179277$meta$Sex)), nrow(cohorts$GSE179277$meta),
                             sum(!is.na(cohorts$GSE179277$meta$ViralLoad)), 0, 0, 0, 0, 0, 0)),
  data.frame(Cohort = "GSE231409", Variable = c("Age", "Sex", "Sequencing platform", "Viral load (copies)", "Respiratory severity score",
                                                  "Symptom status", "Days from symptom onset", "Vaccination", "Hospitalization", "Batch"),
             N_available = c(sum(!is.na(cohorts$GSE231409$meta$Age)), sum(!is.na(cohorts$GSE231409$meta$Sex)), nrow(cohorts$GSE231409$meta),
                             sum(!is.na(cohorts$GSE231409$meta$ViralLoad)), 0, sum(!is.na(cohorts$GSE231409$meta$Symptomatic)),
                             sum(!is.na(cohorts$GSE231409$meta$DaysSymptoms)), sum(!is.na(cohorts$GSE231409$meta$VaccineDoses)),
                             sum(cohorts$GSE231409$meta$Hospitalized %in% c("Yes", "No")), sum(!is.na(cohorts$GSE231409$meta$Batch)))),
  data.frame(Cohort = "GSE172274", Variable = c("Age", "Sex", "Sequencing platform", "Viral load", "Respiratory severity score",
                                                  "Symptom status", "Days from symptom onset", "Vaccination", "Hospitalization", "Batch"),
             N_available = c(sum(!is.na(cohorts$GSE172274$meta$Age)), sum(!is.na(cohorts$GSE172274$meta$Sex)), nrow(cohorts$GSE172274$meta),
                             0, sum(!is.na(cohorts$GSE172274$meta$RespiratoryScore)), 0, 0, 0, 0, 0))
) %>%
  group_by(Cohort) %>%
  mutate(
    N_samples = nrow(cohorts[[unique(Cohort)]]$meta),
    Availability = ifelse(N_available == 0, "Not reported",
                          ifelse(N_available == N_samples, "Complete", paste0("Partial (", N_available, "/", N_samples, ")"))),
    Used_in_models = case_when(
      Variable %in% c("Age", "Sex") ~ "Primary model",
      Variable == "Batch" & Cohort == "GSE231409" ~ "Primary model",
      grepl("Viral load", Variable) & N_available > 0 ~ "Sensitivity model",
      Variable == "Respiratory severity score" & N_available > 0 ~ "Sensitivity model",
      Variable == "Symptom status" & N_available > 0 ~ "Sensitivity model",
      Variable == "Sequencing platform" ~ "Constant (Illumina NovaSeq 6000 in all cohorts)",
      N_available > 0 ~ "Descriptive only",
      TRUE ~ "Not available"
    )
  ) %>%
  ungroup()

write.csv(availability, file.path(out_dir, "covariate_availability_by_cohort.csv"), row.names = FALSE)

clinical_by_group <- bind_rows(lapply(c("GSE179277", "GSE231409"), function(g) {
  m <- cohorts[[g]]$meta
  ad <- m$ViralLoad[m$Group == "Adult"]; pe <- m$ViralLoad[m$Group == "Pediatric"]
  data.frame(
    Cohort = g, Variable = unique(m$ViralLoad_unit),
    Adult_n = sum(!is.na(ad)), Adult_median = median(ad, na.rm = TRUE),
    Adult_Q1 = quantile(ad, 0.25, na.rm = TRUE), Adult_Q3 = quantile(ad, 0.75, na.rm = TRUE),
    Pediatric_n = sum(!is.na(pe)), Pediatric_median = median(pe, na.rm = TRUE),
    Pediatric_Q1 = quantile(pe, 0.25, na.rm = TRUE), Pediatric_Q3 = quantile(pe, 0.75, na.rm = TRUE),
    Wilcoxon_p = wilcox.test(ad, pe, exact = FALSE)$p.value
  )
}))

sev <- cohorts$GSE172274$meta
severity_table <- as.data.frame.matrix(table(sev$Group, sev$RespiratoryScore))
severity_table$Group <- rownames(severity_table)
sev_p <- fisher.test(table(sev$Group, sev$Severe))$p.value

sym <- cohorts$GSE231409$meta
symptom_table <- as.data.frame.matrix(table(sym$Group, sym$Symptomatic, useNA = "ifany"))
symptom_table$Group <- rownames(symptom_table)
sym_p <- fisher.test(table(sym$Group, sym$Symptomatic))$p.value

write.csv(clinical_by_group, file.path(out_dir, "viral_load_by_age_group.csv"), row.names = FALSE)
write.csv(severity_table, file.path(out_dir, "GSE172274_respiratory_score_by_group.csv"), row.names = FALSE)
write.csv(symptom_table, file.path(out_dir, "GSE231409_symptom_status_by_group.csv"), row.names = FALSE)
writeLines(c(
  paste("GSE172274 severe (score >= 3) vs Group, Fisher p =", signif(sev_p, 3)),
  paste("GSE231409 symptomatic vs Group, Fisher p =", signif(sym_p, 3))
), file.path(out_dir, "clinical_association_tests.txt"))
print(clinical_by_group)


# 03. Cohort-level DESeq2 models ----

fit_model <- function(cohort, label, design, keep = NULL, coef = "Group_Pediatric_vs_Adult") {
  counts <- cohort$counts
  meta <- cohort$meta
  if (!is.null(keep)) {
    keep[is.na(keep)] <- FALSE
    counts <- counts[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  vars <- all.vars(design)
  ok <- stats::complete.cases(meta[, vars, drop = FALSE])
  counts <- counts[, ok, drop = FALSE]
  meta <- droplevels(meta[ok, , drop = FALSE])

  n_adult <- sum(meta$Group == "Adult"); n_ped <- sum(meta$Group == "Pediatric")
  feasible <- n_adult >= 2 && n_ped >= 2
  for (v in vars) if (is.factor(meta[[v]]) && nlevels(meta[[v]]) < 2) feasible <- FALSE
  mm <- model.matrix(design, meta)
  if (qr(mm)$rank < ncol(mm)) feasible <- FALSE

  info <- data.frame(Cohort = cohort$gse, Model = label, Design = paste(deparse(design), collapse = ""),
                     Coefficient = coef, N = ncol(counts), N_Adult = n_adult, N_Pediatric = n_ped,
                     Feasible = feasible)
  if (!feasible) {
    message("  ", cohort$gse, " ", label, ": not estimable, skipped")
    return(list(info = info, res = NULL))
  }

  message("  ", cohort$gse, " ", label, " (n = ", ncol(counts), ")")
  dds <- DESeqDataSetFromMatrix(counts, meta, design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = coef)
  out <- data.frame(
    Original_ID = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    pvalue = res$pvalue,
    padj = res$padj,
    stringsAsFactors = FALSE
  ) %>%
    left_join(cohort$annot, by = "Original_ID") %>%
    mutate(Cohort = cohort$gse, Model = label)
  write.csv(out, file.path(model_dir, paste0(cohort$gse, "_", label, ".csv")), row.names = FALSE)
  list(info = info, res = out)
}

age_window <- function(m) !(m$Age >= 16 & m$Age < 21)

message("Fitting cohort models...")
models <- list(
  GSE179277 = list(
    A  = fit_model(cohorts$GSE179277, "A_primary_refit", ~ Sex + Group),
    B  = fit_model(cohorts$GSE179277, "B_viral_load", ~ Sex + log10VL_z + Group),
    C  = fit_model(cohorts$GSE179277, "C_exclude_age_16_20", ~ Sex + Group, keep = age_window(cohorts$GSE179277$meta)),
    D  = fit_model(cohorts$GSE179277, "D_log2_age", ~ Sex + log2Age, coef = "log2Age")
  ),
  GSE231409 = list(
    A  = fit_model(cohorts$GSE231409, "A_primary_refit", ~ Batch + Sex + Group),
    B0 = fit_model(cohorts$GSE231409, "B0_VL_subset_unadjusted", ~ Batch + Sex + Group,
                   keep = !is.na(cohorts$GSE231409$meta$ViralLoad)),
    B  = fit_model(cohorts$GSE231409, "B_viral_load", ~ Batch + Sex + log10VL_z + Group),
    B2 = fit_model(cohorts$GSE231409, "B2_viral_load_symptoms", ~ Batch + Sex + log10VL_z + Symptomatic + Group),
    C  = fit_model(cohorts$GSE231409, "C_exclude_age_16_20", ~ Batch + Sex + Group, keep = age_window(cohorts$GSE231409$meta)),
    D  = fit_model(cohorts$GSE231409, "D_log2_age", ~ Batch + Sex + log2Age, coef = "log2Age")
  ),
  GSE172274 = list(
    A  = fit_model(cohorts$GSE172274, "A_primary_refit", ~ Sex + Group),
    B  = fit_model(cohorts$GSE172274, "B_severity", ~ Sex + Severe + Group),
    B2 = fit_model(cohorts$GSE172274, "B2_severe_adults_excluded", ~ Sex + Group,
                   keep = !(cohorts$GSE172274$meta$Group == "Adult" & cohorts$GSE172274$meta$Severe == "Yes")),
    C  = fit_model(cohorts$GSE172274, "C_exclude_age_16_20", ~ Sex + Group, keep = age_window(cohorts$GSE172274$meta)),
    D  = fit_model(cohorts$GSE172274, "D_log2_age", ~ Sex + log2Age, coef = "log2Age")
  )
)

model_info <- bind_rows(lapply(models, function(x) bind_rows(lapply(x, `[[`, "info"))))
write.csv(model_info, file.path(out_dir, "cohort_model_specifications.csv"), row.names = FALSE)
print(model_info)

# Sanity check: primary refit reproduces the primary cohort estimates.
refit_check <- bind_rows(lapply(names(cohorts), function(g) {
  a <- models[[g]]$A$res
  p <- cohorts[[g]]$primary_res %>% select(Original_ID, lfc_primary = log2FoldChange)
  j <- inner_join(a, p, by = "Original_ID")
  data.frame(Cohort = g, n = nrow(j), max_abs_diff = max(abs(j$log2FoldChange - j$lfc_primary), na.rm = TRUE))
}))
print(refit_check)
write.csv(refit_check, file.path(out_dir, "primary_refit_reproduction_check.csv"), row.names = FALSE)


# 04. Meta-analysis per scenario ----

meta_scenario <- function(parts, label) {
  long <- bind_rows(lapply(parts, function(p) if (is.null(p)) NULL else p$res)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "", is.finite(log2FoldChange), is.finite(lfcSE), lfcSE > 0) %>%
    filter(SYMBOL %in% primary$SYMBOL)
  # Duplicated symbols within a cohort were excluded upstream; keep a guard.
  long <- long %>% group_by(SYMBOL, Cohort) %>% filter(n() == 1) %>% ungroup()
  split_g <- split(long, long$SYMBOL)
  split_g <- split_g[vapply(split_g, nrow, integer(1)) >= 2]
  message("  ", label, ": ", length(split_g), " genes with k >= 2")
  res <- bind_rows(lapply(names(split_g), function(g) {
    d <- split_g[[g]]
    fit <- tryCatch(
      suppressWarnings(metafor::rma(yi = d$log2FoldChange, sei = d$lfcSE, method = "REML", test = "z",
                                    control = list(stepadj = 0.5, maxiter = 1000))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    data.frame(SYMBOL = g, k = nrow(d), pooled_log2FC = as.numeric(fit$b), pooled_SE = fit$se,
               ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, pvalue = fit$pval, I2 = fit$I2)
  }))
  res$FDR <- p.adjust(res$pvalue, method = "BH")
  res$Scenario <- label
  res
}

message("Meta-analysis per scenario...")
scenarios <- list(
  S0_primary_refit = list(models$GSE179277$A, models$GSE231409$A, models$GSE172274$A),
  S1_viral_load_severity = list(models$GSE179277$B, models$GSE231409$B, models$GSE172274$B),
  S1b_VL_subset_unadjusted = list(models$GSE179277$A, models$GSE231409$B0, models$GSE172274$A),
  S2_extended_clinical = list(models$GSE179277$B, models$GSE231409$B2, models$GSE172274$B2),
  S3_exclude_age_16_20 = list(models$GSE179277$C, models$GSE231409$C, models$GSE172274$C),
  S4_log2_age_continuous = list(models$GSE179277$D, models$GSE231409$D, models$GSE172274$D)
)

scenario_results <- bind_rows(lapply(names(scenarios), function(s) meta_scenario(scenarios[[s]], s)))
write.csv(scenario_results, file.path(out_dir, "scenario_meta_analysis_all_genes.csv"), row.names = FALSE)


# 05. Signature genes across scenarios ----

sig_ref <- signature %>% select(SYMBOL, k_primary = k, RE_log2FC, RE_FDR)

sig_scen <- scenario_results %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  left_join(sig_ref, by = "SYMBOL") %>%
  mutate(
    expected_sign = ifelse(Scenario == "S4_log2_age_continuous", -sign(RE_log2FC), sign(RE_log2FC)),
    same_direction = sign(pooled_log2FC) == expected_sign,
    significant_FDR05 = FDR < ALPHA,
    nominal_p05_same_direction = same_direction & pvalue < 0.05
  )

write.csv(sig_scen, file.path(out_dir, "signature_genes_by_scenario_long.csv"), row.names = FALSE)

sig_wide <- sig_scen %>%
  select(SYMBOL, k_primary, RE_log2FC, RE_FDR, Scenario, pooled_log2FC, FDR, same_direction) %>%
  pivot_wider(names_from = Scenario, values_from = c(pooled_log2FC, FDR, same_direction),
              names_glue = "{Scenario}_{.value}") %>%
  arrange(RE_FDR)

write.csv(sig_wide, file.path(out_dir, "signature_genes_by_scenario_wide.csv"), row.names = FALSE)

scenario_summary <- sig_scen %>%
  group_by(Scenario) %>%
  summarise(
    n_signature_tested = n(),
    n_same_direction = sum(same_direction, na.rm = TRUE),
    n_same_direction_nominal_p05 = sum(nominal_p05_same_direction, na.rm = TRUE),
    n_same_direction_FDR05 = sum(same_direction & significant_FDR05, na.rm = TRUE),
    pearson_r_vs_primary = ifelse(first(Scenario) == "S4_log2_age_continuous",
                                  cor(-pooled_log2FC, RE_log2FC), cor(pooled_log2FC, RE_log2FC)),
    median_ratio_vs_primary = median(ifelse(first(Scenario) == "S4_log2_age_continuous", NA, pooled_log2FC / RE_log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    scenario_results %>% group_by(Scenario) %>%
      summarise(n_genes_meta = n(), n_genome_FDR05 = sum(FDR < ALPHA), .groups = "drop"),
    by = "Scenario"
  )

write.csv(scenario_summary, file.path(out_dir, "scenario_summary.csv"), row.names = FALSE)
print(as.data.frame(scenario_summary))

# Per-cohort comparison of the signature genes: primary refit vs adjusted
cohort_sig <- bind_rows(lapply(models, function(x) bind_rows(lapply(x, `[[`, "res")))) %>%
  filter(SYMBOL %in% signature$SYMBOL)
write.csv(cohort_sig, file.path(out_dir, "signature_genes_cohort_models_long.csv"), row.names = FALSE)


# 06. Figures ----

vl_df <- bind_rows(lapply(c("GSE179277", "GSE231409"), function(g) {
  m <- cohorts[[g]]$meta
  data.frame(Cohort = paste0(g, "\n(", unique(m$ViralLoad_unit), ")"), Group = m$Group, log10VL = m$log10VL)
})) %>% filter(!is.na(log10VL))

p_vl <- ggplot(vl_df, aes(x = Group, y = log10VL, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6) +
  geom_point(position = position_jitter(width = 0.15, seed = 1), size = 1.2, alpha = 0.7) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c(Adult = "#2C6FB7", Pediatric = "#E8B500")) +
  labs(x = NULL, y = "log10(viral load + 1)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_viral_load_by_age_group.pdf"), p_vl, width = 7, height = 4)
ggsave(file.path(fig_dir, "FigS_viral_load_by_age_group.tiff"), p_vl, width = 7, height = 4, dpi = 600, compression = "lzw")

scen_labels <- c(
  S1_viral_load_severity = "Adjusted for viral load\n(GSE179277, GSE231409) / severity (GSE172274)",
  S1b_VL_subset_unadjusted = "GSE231409 restricted to samples\nwith viral load, unadjusted",
  S2_extended_clinical = "Extended: + symptoms (GSE231409),\nsevere adults excluded (GSE172274)",
  S3_exclude_age_16_20 = "Participants aged 16-20 years\nexcluded",
  S4_log2_age_continuous = "Continuous log2(age)\n(sign reversed: pediatric direction)"
)

scatter_df <- sig_scen %>%
  filter(Scenario != "S0_primary_refit") %>%
  mutate(
    y = ifelse(Scenario == "S4_log2_age_continuous", -pooled_log2FC, pooled_log2FC),
    status = case_when(
      same_direction & significant_FDR05 ~ "Same direction, FDR < 0.05",
      same_direction & pvalue < 0.05 ~ "Same direction, p < 0.05",
      same_direction ~ "Same direction, n.s.",
      TRUE ~ "Opposite direction"
    ),
    Scenario = factor(Scenario, levels = names(scen_labels), labels = scen_labels)
  )

r_lab <- scatter_df %>% group_by(Scenario) %>%
  summarise(lab = paste0("r = ", sprintf("%.2f", cor(y, RE_log2FC)), "\n",
                         sum(same_direction), "/", n(), " same direction"), .groups = "drop")

p_sc <- ggplot(scatter_df, aes(x = RE_log2FC, y = y)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = status), size = 1.8) +
  geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ Scenario, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("Same direction, FDR < 0.05" = "#2C6FB7", "Same direction, p < 0.05" = "#7FA7D9",
                                 "Same direction, n.s." = "#BDBDBD", "Opposite direction" = "#D55E00")) +
  labs(x = "Primary pooled log2FC (Pediatric vs Adult)", y = "Sensitivity pooled estimate", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_signature_clinical_sensitivity_scatter.pdf"), p_sc, width = 11, height = 7.5)
ggsave(file.path(fig_dir, "FigS_signature_clinical_sensitivity_scatter.tiff"), p_sc, width = 11, height = 7.5, dpi = 600, compression = "lzw")

order_genes <- signature %>% arrange(RE_log2FC) %>% pull(SYMBOL)
forest_df <- bind_rows(
  signature %>% transmute(SYMBOL, Scenario = "Primary", est = RE_log2FC, lb = RE_CI_lower, ub = RE_CI_upper),
  sig_scen %>% filter(Scenario %in% c("S1_viral_load_severity", "S3_exclude_age_16_20")) %>%
    transmute(SYMBOL, Scenario, est = pooled_log2FC, lb = ci_lb, ub = ci_ub)
) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = order_genes),
         Scenario = factor(Scenario, levels = c("Primary", "S1_viral_load_severity", "S3_exclude_age_16_20"),
                           labels = c("Primary", "Viral load / severity adjusted", "Age 16-20 excluded")))

p_forest <- ggplot(forest_df, aes(y = SYMBOL, x = est, xmin = lb, xmax = ub, colour = Scenario)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.15, linewidth = 0.4) +
  scale_colour_manual(values = c("Primary" = "#222222", "Viral load / severity adjusted" = "#D55E00", "Age 16-20 excluded" = "#0072B2")) +
  labs(x = "Pooled log2FC (Pediatric vs Adult), 95% CI", y = NULL, colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "FigS_signature_clinical_sensitivity_forest.pdf"), p_forest, width = 7, height = 9.5)
ggsave(file.path(fig_dir, "FigS_signature_clinical_sensitivity_forest.tiff"), p_forest, width = 7, height = 9.5, dpi = 600, compression = "lzw")


# 07. Summary ----

writeLines(c(
  paste("Script build:", SCRIPT_BUILD),
  "Primary-refit reproduction (max |log2FC difference| per cohort):",
  paste0("  ", refit_check$Cohort, ": ", signif(refit_check$max_abs_diff, 3)),
  "",
  "Viral load by age group:",
  paste0("  ", clinical_by_group$Cohort, ": Adult median ", signif(clinical_by_group$Adult_median, 4),
         " vs Pediatric median ", signif(clinical_by_group$Pediatric_median, 4),
         " (Wilcoxon p = ", signif(clinical_by_group$Wilcoxon_p, 3), ")"),
  paste("GSE172274 severe vs group Fisher p =", signif(sev_p, 3)),
  "",
  "Signature genes across scenarios:",
  paste0("  ", scenario_summary$Scenario, ": same direction ", scenario_summary$n_same_direction, "/",
         scenario_summary$n_signature_tested, "; same direction & FDR < 0.05 ",
         scenario_summary$n_same_direction_FDR05, "; r = ", sprintf("%.3f", scenario_summary$pearson_r_vs_primary))
), file.path(out_dir, "run_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done: ", out_dir)
