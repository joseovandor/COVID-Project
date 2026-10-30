# Supplementary Tables S1-S22
#
# One Excel workbook (README sheet plus one sheet per table) and one CSV per
# table. Tables are numbered in order of first citation in the manuscript.
# Tables whose inputs are missing are skipped and reported at the end.
#
# Output: results/manuscript/supplementary_tables/

SCRIPT_BUILD <- "SUPPLEMENTARY_TABLES_2026-10-28"
message("Running script build: ", SCRIPT_BUILD)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count <- dplyr::count; first <- dplyr::first; mutate <- dplyr::mutate
summarise <- dplyr::summarise; arrange <- dplyr::arrange; distinct <- dplyr::distinct

ALPHA <- 0.05
COHORTS <- c("GSE172274", "GSE179277", "GSE231409")

out_dir <- here::here("results", "manuscript", "supplementary_tables")
csv_dir <- file.path(out_dir, "csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

p <- function(...) here::here(...)
rd <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
opt <- function(path) if (file.exists(path)) rd(path) else NULL

tables <- list()
add_table <- function(id, sheet, title, description, data, source) {
  tables[[length(tables) + 1]] <<- list(id = id, sheet = sheet, title = title, description = description,
                                        data = data, source = source)
}


# Samples, cohort summary and covariates ----

raw_get <- function(raw, col) if (col %in% colnames(raw)) raw[[col]] else rep(NA, nrow(raw))
as_num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "NA", "N/A"), NA, x)))

samples <- bind_rows(lapply(COHORTS, function(g) {
  m <- rd(p("results", "deseq2", g, "metadata", paste0(g, "_sample_metadata.csv")))
  raw <- rd(p("results", "deseq2", g, "metadata", paste0(g, "_metadata_raw.csv")))
  raw <- raw[match(m$GEO_ID, raw$geo_accession), ]
  data.frame(
    Cohort = g,
    GEO_ID = m$GEO_ID,
    Sample = if ("Sample" %in% colnames(m)) m$Sample else NA,
    Age_years = m$Age,
    Age_group = m$Group,
    Original_study_label = if ("OriginalLabel" %in% colnames(m)) m$OriginalLabel else NA,
    Sex = m$Sex,
    Specimen = coalesce(as.character(raw_get(raw, "tissue:ch1")), as.character(raw_get(raw, "sample type:ch1")),
                        as.character(raw_get(raw, "source_name_ch1"))),
    Sequencing_instrument = as.character(raw_get(raw, "instrument_model")),
    Batch = if ("Batch" %in% colnames(m)) as.character(m$Batch) else NA,
    Viral_load_RPM = as_num(raw_get(raw, "sars-cov2_rpm:ch1")),
    Viral_load_copies = as_num(raw_get(raw, "vl copies:ch1")),
    Respiratory_score = as_num(raw_get(raw, "respiratory score:ch1")),
    Symptom_status = as.character(raw_get(raw, "covid cat:ch1")),
    Hospitalized = as.character(raw_get(raw, "hospitalized:ch1")),
    Vaccine_doses = as_num(raw_get(raw, "vaccine doses:ch1")),
    Days_from_symptom_onset = as_num(raw_get(raw, "timing sx:ch1")),
    Days_from_diagnosis = as_num(raw_get(raw, "timing dx:ch1")),
    stringsAsFactors = FALSE
  )
}))

add_table("S1", "S1_samples", "Sample-level metadata of the 225 SARS-CoV-2-positive bulk RNA-seq samples",
          "One row per SARS-CoV-2-positive sample (N = 225). Age group derived from recorded age (pediatric < 18 years; adult >= 18 years). The original study label is the label in the GEO sample title; GSM5413083 (19 years) carried a pediatric label in GEO and is analysed as adult. Clinical variables are reported as deposited in GEO; NA = not available.", samples, "results/deseq2/*/metadata/")

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("Not reported")
  sprintf("%s (%s-%s); n = %d", signif(median(x), 3), signif(quantile(x, .25), 3), signif(quantile(x, .75), 3), length(x))
}
summary_tab <- samples %>%
  group_by(Cohort, Age_group) %>%
  summarise(
    N = n(),
    Age_median_IQR = med_iqr(Age_years),
    Age_range = sprintf("%s-%s", min(Age_years), max(Age_years)),
    Female_n_pct = sprintf("%d (%.0f%%)", sum(Sex == "Female"), 100 * mean(Sex == "Female")),
    Aged_16_to_20_n = sum(Age_years >= 16 & Age_years < 21),
    Viral_load_RPM_median_IQR = med_iqr(Viral_load_RPM),
    Viral_load_copies_median_IQR = med_iqr(Viral_load_copies),
    Respiratory_score_ge3_n = ifelse(all(is.na(Respiratory_score)), "Not reported", as.character(sum(Respiratory_score >= 3, na.rm = TRUE))),
    Symptomatic_n = ifelse(all(is.na(Symptom_status)), "Not reported", as.character(sum(Symptom_status == "POS_SX", na.rm = TRUE))),
    Sequencing_instrument = paste(unique(na.omit(Sequencing_instrument)), collapse = "; "),
    .groups = "drop"
  )
add_table("S2", "S2_cohort_summary", "Demographic and clinical characteristics by cohort and age group",
          "Values are median (IQR) or n (%). Not reported = variable not available in GEO for that cohort.", summary_tab, "derived from Supplementary Table S1")

cov <- opt(p("results", "sensitivity", "clinical_covariates", "covariate_availability_by_cohort.csv"))
vl <- opt(p("results", "sensitivity", "clinical_covariates", "viral_load_by_age_group.csv"))
if (!is.null(cov)) add_table("S3", "S3_covariates", "Availability of clinical and technical covariates and their use in the models",
                             "Availability of each clinical and technical covariate and its use in the statistical models (primary model, sensitivity model, or descriptive only).", cov, "results/sensitivity/clinical_covariates/covariate_availability_by_cohort.csv")


# DESeq2 models ----

spec <- bind_rows(lapply(COHORTS, function(g) rd(p("results", "deseq2", g, paste0(g, "_primary_model_specification.csv")))))
lowc <- bind_rows(lapply(COHORTS, function(g) {
  d <- opt(p("results", "deseq2", g, "diagnostics", paste0(g, "_low_count_filtering_summary.csv")))
  if (!is.null(d)) d$Dataset <- g
  d
}))
add_table("S4", "S4_DESeq2_models", "Exact cohort-level DESeq2 models",
          "Exact cohort-level DESeq2 models. The same DESeq2 results object was used for the meta-analysis (unshrunken log2FC and SE) and for pre-ranked GSEA (Wald statistic).", spec, "results/deseq2/*/*_primary_model_specification.csv")


# Meta-analysis ----

primary <- rd(p("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"))
robust <- rd(p("results", "meta_analysis", "sensitivity", "robustness_classification.csv"))
inf <- opt(p("results", "meta_analysis", "inference_methods", "signature_genes_by_inference_method.csv"))

sig_tab <- primary %>% filter(RE_FDR < ALPHA) %>%
  left_join(robust %>% select(SYMBOL, N_leave_one_out_models, All_leave_one_out_same_direction,
                              All_leave_one_out_nominal_significant, RobustnessClass), by = "SYMBOL")
if (!is.null(inf)) sig_tab <- sig_tab %>% left_join(inf %>% select(SYMBOL, EvidenceTier), by = "SYMBOL")
sig_tab <- sig_tab %>% arrange(RE_FDR)

add_table("S8", "S8_signature_genes", "Age-associated signature genes (random-effects FDR < 0.05)",
          "Signature genes (random-effects FDR < 0.05 across 18,209 eligible genes). Robustness class and evidence tier are descriptive and do not replace the primary FDR criterion (see Supplementary Table S11 for the evidence tiers).", sig_tab, "results/meta_analysis/primary + sensitivity")
add_table("S7", "S7_meta_all_genes", "Random-effects meta-analysis of all eligible genes",
          "Random-effects meta-analysis of all genes present in at least two cohorts; Benjamini-Hochberg correction across the complete eligible gene universe without preselection.", primary, "results/meta_analysis/primary/meta_analysis_REM_all_genes.csv")

loo <- rd(p("results", "meta_analysis", "sensitivity", "leave_one_out_results.csv")) %>%
  filter(SYMBOL %in% sig_tab$SYMBOL)
add_table("S9", "S9_leave_one_out", "Leave-one-cohort-out sensitivity analysis of the signature genes",
          "Two-cohort random-effects estimates after omitting each cohort, for the 33 signature genes present in three cohorts. This is a stability analysis, not an independent validation, because the omitted cohort contributed to gene selection.", loo, "results/meta_analysis/sensitivity/leave_one_out_results.csv")


# Inference methods, clinical covariates, continuous age and infection status ----

if (!is.null(inf)) {
  gw <- rd(p("results", "meta_analysis", "inference_methods", "meta_inference_genome_wide_summary.csv"))
  add_table("S11", "S11_inference_methods", "Sensitivity of the signature to the meta-analytic inference method",
            "Each signature gene under alternative inference methods. Evidence tiers: Tier 1 = k = 3, HKSJ p < .05, and Bayesian HN(0.5) 95% CrI excluding 0; Tier 2 = Bayesian CrI excluding 0 without meeting Tier 1; Tier 3 = exploratory.", inf, "results/meta_analysis/inference_methods/signature_genes_by_inference_method.csv")
  add_table("S12", "S12_inference_genome_wide", "Genome-wide number of significant genes by inference method",
            "Minimum p = smallest p value observed among the 18,209 genes; BH threshold (rank 1) = p value required for the top-ranked gene to reach FDR < .05. With k <= 3, HKSJ (t distribution with k - 1 degrees of freedom) cannot reach genome-wide FDR control.", gw, "results/meta_analysis/inference_methods/meta_inference_genome_wide_summary.csv")
  grid <- rd(p("results", "meta_analysis", "inference_methods", "robust_threshold_sensitivity.csv"))
  add_table("S10", "S10_robust_thresholds", "Alternative thresholds for the descriptive 'Robust' class",
            "Number and identity of Robust genes under every combination of maximum I2 and maximum LOO p value, with REML z inference or additionally requiring HKSJ p < .05. Primary definition: I2 < 75% and LOO p < .05 (26 genes).", grid, "results/meta_analysis/inference_methods/robust_threshold_sensitivity.csv")
}

clin_wide <- opt(p("results", "sensitivity", "clinical_covariates", "signature_genes_by_scenario_wide.csv"))
clin_sum <- opt(p("results", "sensitivity", "clinical_covariates", "scenario_summary.csv"))
clin_models <- opt(p("results", "sensitivity", "clinical_covariates", "cohort_model_specifications.csv"))
if (!is.null(clin_wide)) {
  add_table("S13", "S13_clinical_sensitivity", "Signature genes after adjustment for clinical covariates and alternative age modelling",
            "Pooled log2FC, FDR, and direction of the 44 signature genes in each sensitivity scenario: S0 = primary refit; S1 = adjustment for viral load (GSE179277, GSE231409) and respiratory severity (GSE172274); S1b = GSE231409 restricted to samples with viral load, unadjusted; S2 = extended adjustment (symptom status in GSE231409; adults with severity score >= 3 excluded in GSE172274); S3 = participants aged 16-20 years excluded; S4 = continuous log2(age + 1). VL = viral load.", clin_wide, "results/sensitivity/clinical_covariates/")
  add_table("S13b", "S13b_scenario_summary", "Summary of the clinical sensitivity scenarios", "Summary of each sensitivity scenario (S0-S4 as in Supplementary Table S13). Median ratio = median of the ratio of scenario to primary pooled log2FC.", clin_sum,
            "results/sensitivity/clinical_covariates/scenario_summary.csv")
  add_table("S13c", "S13c_scenario_models", "Cohort-level models of the clinical sensitivity scenarios", "Cohort-level DESeq2 model of each sensitivity scenario. log10VL_z = standardized log10(viral load + 1).",
            clin_models, "results/sensitivity/clinical_covariates/cohort_model_specifications.csv")
  if (!is.null(vl)) add_table("S13d", "S13d_viral_load_by_group", "Viral load by age group", "Viral load by age group; two-sided Wilcoxon rank-sum test. Q1/Q3 = first/third quartile.", vl,
                              "results/sensitivity/clinical_covariates/viral_load_by_age_group.csv")
}

ac <- opt(p("results", "sensitivity", "continuous_age_sex_interaction", "signature_genes_continuous_age_interaction_by_cohort.csv"))
if (!is.null(ac)) {
  add_table("S14", "S14_age_continuous_sex", "Continuous age and sex-by-group interaction for the signature genes",
            "Per gene and cohort: pediatric-vs-adult effect (group); effect per 1 SD of age (ageSD; concordant = opposite sign to the group effect); log2(age + 1) effect; effect after excluding ages 16-20 years; and Sex x Group interaction with genome-wide adjusted p and FDR within the signature genes.", ac, "results/sensitivity/continuous_age_sex_interaction/")
  add_table("S14b", "S14b_summary", "Summary by cohort of the continuous-age and interaction analyses", "Cohort-level summary of the continuous-age and interaction analyses, for the signature genes and genome-wide.",
            rd(p("results", "sensitivity", "continuous_age_sex_interaction", "summary_by_cohort.csv")), "")
  add_table("S14c", "S14c_sex_by_group", "Sex x age-group cell sizes", "Number of samples in each sex-by-age-group cell.",
            rd(p("results", "sensitivity", "continuous_age_sex_interaction", "sex_by_group_cell_sizes.csv")), "")
  add_table("S14d", "S14d_age_structure", "Age structure by cohort and age group", "Age (years) by cohort and age group. Q1/Q3 = first/third quartile.",
            rd(p("results", "sensitivity", "continuous_age_sex_interaction", "age_summary_by_cohort_group.csv")), "")
}

ai <- opt(p("results", "sensitivity", "age_by_infection", "signature_genes_age_by_infection.csv"))
if (!is.null(ai)) {
  add_table("S15", "S15_age_by_infection", "Age effect of the signature genes with and without SARS-CoV-2 (GSE179277)",
            "DESeq2 model ~ Sex + Group + Infection + Group:Infection in GSE179277, including SARS-CoV-2-negative individuals with no respiratory virus detected (not healthy volunteers). Age effects are pediatric vs adult within negative (neg) or positive (pos) individuals; interaction = difference between the two age effects.", ai, "results/sensitivity/age_by_infection/")
  add_table("S15b", "S15b_infection_samples", "Sample sizes of the age x infection analysis", "Sample sizes of the age x infection analysis in GSE179277.",
            rd(p("results", "sensitivity", "age_by_infection", "sample_sizes_infection_by_group_by_sex.csv")), "")
  add_table("S15c", "S15c_infection_genome", "Genome-wide summary of the age x infection analysis", "Genome-wide summary of the age x infection analysis in GSE179277.",
            rd(p("results", "sensitivity", "age_by_infection", "genome_wide_summary.csv")), "")
}


# GSEA ----

gsea <- rd(p("results", "gsea", "GSEA_all_cohorts_long.csv"))
add_table("S16", "S16_GSEA_by_cohort", "Pre-ranked GSEA results by cohort",
          "fgsea (multilevel) results on the complete DESeq2 Wald-statistic ranking of each cohort (Hallmark, Reactome, KEGG, and GO Biological Process). ES = enrichment score; log2err = expected error of the p value estimate.", gsea, "results/gsea/GSEA_all_cohorts_long.csv")

repro <- bind_rows(
  opt(p("results", "gsea", "cross_cohort", "tables", "GSEA_reproducible_significant_3of3.csv")),
  opt(p("results", "gsea", "cross_cohort", "tables", "GSEA_reproducible_supported_2of3.csv"))
) %>% distinct(Collection, pathway, .keep_all = TRUE)
pc <- opt(p("results", "gsea", "themes", "gsea_partial_conjunction_all_pathways.csv"))
if (!is.null(pc)) {
  repro <- repro %>% left_join(pc %>% select(Collection, pathway, p_PC_2, FDR_PC_2, p_PC_all, FDR_PC_all),
                               by = c("Collection", "pathway"))
}
add_table("S17", "S17_GSEA_reproducible", "Pathways enriched with a consistent direction in >= 2 cohorts",
          "Pathways enriched with the same NES direction and cohort-level FDR < .05 in at least two cohorts (descriptive criterion; it does not control the FDR of the selected set). Partial-conjunction p values for enrichment in >= 2 cohorts (PC 2) and in all cohorts (PC all), with BH correction across all pathways tested in >= 2 cohorts, provide the formal replicability analysis.", repro, "results/gsea/cross_cohort + results/gsea/themes")

themes <- opt(p("results", "gsea", "themes", "gsea_themes.csv"))
if (!is.null(themes)) {
  p2t <- rd(p("results", "gsea", "themes", "gsea_pathway_to_theme.csv"))
  lab_file <- p("config", "gsea_theme_labels.csv")
  if (file.exists(lab_file)) {
    labs <- rd(lab_file)
    themes <- themes %>% left_join(labs, by = c("representative_pathway" = "Representative_pathway")) %>%
      relocate(Label, .after = Direction)
  }
  add_table("S18", "S18_GSEA_themes", "Non-redundant GSEA themes",
            "Non-redundant themes obtained by average-linkage clustering of the overlap coefficient of pooled leading-edge genes (similarity >= 0.5). Core leading-edge genes are present in at least 50% of the theme's pathways (first 25 shown).", themes, "results/gsea/themes/gsea_themes.csv")
  add_table("S18b", "S18b_pathway_to_theme", "Assignment of reproducible pathways to themes", "Assignment of each consistently enriched pathway to its theme.", p2t,
            "results/gsea/themes/gsea_pathway_to_theme.csv")
}


# Single-cell ----

scp <- function(f) p("results", "scRNAseq_COVID_nasal_age", "pseudobulk_signature", "tables", f)
sc_status <- opt(scp("signature_gene_status_singlecell.csv"))
if (!is.null(sc_status)) {
  add_table("S20", "S20_sc_gene_status", "Status of every signature gene in the single-cell dataset",
            "Status of every signature gene in the single-cell dataset (presence in the expression matrix; number of cell types tested in the primary pseudobulk and secondary analyses).", sc_status, "signature_gene_status_singlecell.csv")
  add_table("S21", "S21_sc_pseudobulk_edgeR", "Donor-level pseudobulk analysis of the signature genes (primary)",
            "Primary analysis: raw UMI counts summed per donor and cell type (>= 20 cells per donor-cell-type combination; >= 3 donors per age group); edgeR quasi-likelihood F test, ~ Sex + Group. FDR within cell type (genome-wide and signature genes) and FDR_global across all signature gene-by-cell-type tests. Ped = pediatric donors.", rd(scp("PRIMARY_pseudobulk_edgeR_signature_genes.csv")), "PRIMARY_pseudobulk_edgeR_signature_genes.csv")
  add_table("S21b", "S21b_sc_models", "Pseudobulk models and donors per cell type", "Pseudobulk model of each cell type, number of donors, genes retained by filterByExpr, and signature genes tested or filtered.",
            rd(scp("PRIMARY_pseudobulk_edgeR_models_by_celltype.csv")), "PRIMARY_pseudobulk_edgeR_models_by_celltype.csv")
  add_table("S22", "S22_sc_mean_lognorm", "Donor-level mean log-normalized expression (secondary analysis)",
            "Secondary analysis: mean of the processed log-normalized expression per donor and cell type (not an aggregation of raw counts); two-sided Wilcoxon rank-sum test; Cliff's delta positive = higher in pediatric donors. Ped = pediatric donors.", rd(scp("SECONDARY_mean_lognorm_wilcoxon_signature_genes.csv")), "SECONDARY_mean_lognorm_wilcoxon_signature_genes.csv")
  ab <- opt(scp("celltype_abundance_edgeR.csv"))
  if (!is.null(ab)) add_table("S19", "S19_sc_abundance", "Cell-type abundance by age group",
                              "Epithelial cell counts per donor and cell type tested with edgeR quasi-likelihood models (library size = total epithelial cells per donor). Prop = proportion; Ped = pediatric donors.", ab, "celltype_abundance_edgeR.csv")
  donors <- opt(p("results", "scRNAseq_COVID_nasal_age", "pseudobulk_signature", "tables", "donors_group_sex.csv"))
  if (!is.null(donors)) add_table("S19b", "S19b_sc_donors", "Single-cell donors included", "Donors (9 adult, 18 pediatric) contributing at least one donor-cell-type combination with >= 20 cells to the pseudobulk analysis.", donors, "donors_group_sex.csv")
}


# Low-count filter sensitivity and software ----

sens <- bind_rows(lapply(COHORTS, function(g) opt(p("results", "deseq2", g, "sensitivity", paste0(g, "_sensitivity_summary.csv")))))
add_table("S6", "S6_filter_sensitivity", "Alternative low-count filtering and cohort-level sensitivity models",
          "Genes tested and significant (FDR < .05) under the primary prevalence filter and an independent CPM >= 1 filter, Pearson correlation of log2FC between filters, and number of significant genes in the continuous-age and Sex x Group interaction models (all Ensembl genes passing the primary filter).", sens, "results/deseq2/*/sensitivity/*_sensitivity_summary.csv")

si_files <- list.files(p("results"), pattern = "sessionInfo.*\\.txt$|_sessionInfo\\.txt$", recursive = TRUE, full.names = TRUE)
software <- bind_rows(lapply(si_files, function(f) {
  x <- readLines(f, warn = FALSE)
  rv <- x[grepl("^R version", x)][1]
  pk <- unique(unlist(regmatches(x, gregexpr("[A-Za-z0-9.]+_[0-9][0-9.\\-]*", x))))
  data.frame(Stage = sub(paste0(p("results"), "/"), "", f), R_version = rv, Packages = paste(sort(pk), collapse = "; "))
}))
add_table("S5", "S5_software", "Software environment of every analysis stage",
          "R version and attached/loaded packages parsed from the sessionInfo() output of each analysis stage.", software, "results/**/sessionInfo*.txt")


# Write workbook, CSVs and index ----

# Order tables by number (S1, S2, ..., S13b, ...), i.e. by first citation in the manuscript
tab_key <- vapply(tables, function(t) {
  n <- as.integer(sub("^S([0-9]+).*$", "\\1", t$id))
  sprintf("%03d%s", n, sub("^S[0-9]+", "", t$id))
}, character(1))
tables <- tables[order(tab_key)]

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#DDE7F3", border = "bottom", wrapText = TRUE)
ttl <- createStyle(textDecoration = "bold", fontSize = 12)

index <- bind_rows(lapply(tables, function(t) data.frame(
  Table = paste0("Table ", t$id), Sheet = t$sheet, Title = t$title, Description = t$description,
  Rows = nrow(t$data), Columns = ncol(t$data), Source = t$source
)))

addWorksheet(wb, "README")
writeData(wb, "README", "Supplementary Tables - Age-associated upper-airway transcriptomic signatures in SARS-CoV-2-positive individuals", startRow = 1)
addStyle(wb, "README", ttl, rows = 1, cols = 1)
writeData(wb, "README", index, startRow = 3, headerStyle = hdr)
setColWidths(wb, "README", cols = 1:7, widths = c(10, 24, 60, 90, 8, 8, 50))

for (t in tables) {
  d <- t$data
  chr <- vapply(d, is.character, logical(1))
  d[chr] <- lapply(d[chr], function(x) ifelse(nchar(x) > 32000, paste0(substr(x, 1, 31990), "..."), x))
  sh <- substr(t$sheet, 1, 31)
  addWorksheet(wb, sh)
  writeData(wb, sh, paste0("Table ", t$id, ". ", t$title), startRow = 1)
  addStyle(wb, sh, ttl, rows = 1, cols = 1)
  if (nzchar(t$description)) writeData(wb, sh, t$description, startRow = 2)
  writeData(wb, sh, d, startRow = 4, headerStyle = hdr)
  freezePane(wb, sh, firstActiveRow = 5)
  setColWidths(wb, sh, cols = seq_len(ncol(d)), widths = pmin(40, pmax(10, nchar(colnames(d)) + 2)))
  write.csv(t$data, file.path(csv_dir, paste0("Table_", t$id, "_", sub("^S[0-9]+[a-z]?_", "", t$sheet), ".csv")), row.names = FALSE)
}

saveWorkbook(wb, file.path(out_dir, "Supplementary_Tables.xlsx"), overwrite = TRUE)
write.csv(index, file.path(out_dir, "supplementary_table_index.csv"), row.names = FALSE)

expected <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22")
missing <- setdiff(expected, vapply(tables, `[[`, character(1), "id"))
if (length(missing) > 0) message("Tables not generated (source files missing): ", paste(missing, collapse = ", "))
print(index[, c("Table", "Rows", "Columns")])
message("Done: ", out_dir)
