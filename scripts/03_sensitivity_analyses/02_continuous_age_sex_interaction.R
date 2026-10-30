# Continuous age and Sex x Group interaction for the signature genes
#
# Uses the sensitivity models fitted by the DESeq2 scripts (same genes and
# samples as the primary model):
#   continuous age : ~ [Batch +] Sex + AgeScaled          (log2FC per 1 SD of age)
#   interaction    : ~ [Batch +] Sex + Group + Sex:Group  (SexMale.GroupPediatric)
# and, when available, the log2(age + 1) and 16-20-year exclusion models from
# 01_clinical_covariates.R. A concordant continuous-age effect has the opposite
# sign of the pediatric-vs-adult log2FC.
#
# Output: results/sensitivity/continuous_age_sex_interaction/

SCRIPT_BUILD <- "CONTINUOUS_AGE_SEX_INTERACTION_2026-10-26"
message("Running script build: ", SCRIPT_BUILD)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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
cohort_ids <- c("GSE172274", "GSE179277", "GSE231409")
id_column <- c(GSE172274 = "ENTREZID", GSE179277 = "ENSEMBL", GSE231409 = "ENSEMBL")

out_dir <- here::here("results", "sensitivity", "continuous_age_sex_interaction")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

primary <- read.csv(here::here("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"),
                    stringsAsFactors = FALSE)
signature <- primary %>% filter(RE_FDR < ALPHA)


# 01. Per-cohort table for the signature genes ----

read_model <- function(gse, file, prefix) {
  d <- read.csv(here::here("results", "deseq2", gse, "sensitivity", file),
                stringsAsFactors = FALSE, colClasses = setNames("character", id_column[[gse]]))
  d <- d[, c(id_column[[gse]], "log2FoldChange", "lfcSE", "pvalue", "padj")]
  colnames(d) <- c("Original_ID", paste0(prefix, c("_log2FC", "_SE", "_p", "_padj")))
  d
}

per_cohort <- bind_rows(lapply(cohort_ids, function(gse) {
  prim <- read.csv(here::here("results", "deseq2", gse, "primary", paste0(gse, "_meta_input.csv")),
                   stringsAsFactors = FALSE, colClasses = c(Original_ID = "character")) %>%
    filter(MetaEligibleAnnotation %in% c(TRUE, "TRUE")) %>%
    transmute(Original_ID, SYMBOL, group_log2FC = log2FoldChange, group_SE = lfcSE,
              group_p = pvalue, group_padj = padj)

  age <- read_model(gse, paste0(gse, "_DESeq2_continuous_age.csv"), "ageSD")
  int <- read_model(gse, paste0(gse, "_DESeq2_group_sex_interaction.csv"), "interaction")

  out <- prim %>%
    left_join(age, by = "Original_ID") %>%
    left_join(int, by = "Original_ID") %>%
    mutate(Cohort = gse)

  # optional models from script 02
  m02 <- here::here("results", "sensitivity", "clinical_covariates", "cohort_models")
  f_log <- file.path(m02, paste0(gse, "_D_log2_age.csv"))
  f_exc <- file.path(m02, paste0(gse, "_C_exclude_age_16_20.csv"))
  if (file.exists(f_log)) {
    out <- out %>% left_join(
      read.csv(f_log, stringsAsFactors = FALSE, colClasses = c(Original_ID = "character")) %>%
        transmute(Original_ID, log2age_log2FC = log2FoldChange, log2age_p = pvalue, log2age_padj = padj),
      by = "Original_ID")
  }
  if (file.exists(f_exc)) {
    out <- out %>% left_join(
      read.csv(f_exc, stringsAsFactors = FALSE, colClasses = c(Original_ID = "character")) %>%
        transmute(Original_ID, excl16_20_log2FC = log2FoldChange, excl16_20_p = pvalue, excl16_20_padj = padj),
      by = "Original_ID")
  }
  out
}))

genome_wide_counts <- per_cohort %>%
  group_by(Cohort) %>%
  summarise(
    genome_genes = n(),
    genome_group_FDR05 = sum(group_padj < ALPHA, na.rm = TRUE),
    genome_continuous_age_FDR05 = sum(ageSD_padj < ALPHA, na.rm = TRUE),
    genome_interaction_FDR05 = sum(interaction_padj < ALPHA, na.rm = TRUE),
    genome_interaction_nominal_p05 = sum(interaction_p < 0.05, na.rm = TRUE),
    genome_interaction_tested = sum(!is.na(interaction_p)),
    .groups = "drop"
  ) %>%
  mutate(genome_interaction_nominal_fraction = genome_interaction_nominal_p05 / genome_interaction_tested)

sig_cohort <- per_cohort %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  group_by(Cohort) %>%
  mutate(interaction_FDR_within_signature = p.adjust(interaction_p, method = "BH")) %>%
  ungroup() %>%
  left_join(signature %>% select(SYMBOL, k, RE_log2FC, RE_FDR), by = "SYMBOL") %>%
  mutate(
    group_direction_matches_pooled = sign(group_log2FC) == sign(RE_log2FC),
    ageSD_concordant = sign(ageSD_log2FC) == -sign(RE_log2FC),
    ageSD_concordant_p05 = ageSD_concordant & ageSD_p < 0.05
  )

if ("log2age_log2FC" %in% colnames(sig_cohort)) {
  sig_cohort <- sig_cohort %>%
    mutate(log2age_concordant = sign(log2age_log2FC) == -sign(RE_log2FC))
}
if ("excl16_20_log2FC" %in% colnames(sig_cohort)) {
  sig_cohort <- sig_cohort %>%
    mutate(excl16_20_same_direction = sign(excl16_20_log2FC) == sign(RE_log2FC),
           excl16_20_ratio = excl16_20_log2FC / group_log2FC)
}

sig_cohort <- sig_cohort %>% arrange(RE_FDR, Cohort)
write.csv(sig_cohort, file.path(out_dir, "signature_genes_continuous_age_interaction_by_cohort.csv"), row.names = FALSE)

summary_cohort <- sig_cohort %>%
  group_by(Cohort) %>%
  summarise(
    signature_genes_in_cohort = n(),
    group_same_direction_as_pooled = sum(group_direction_matches_pooled, na.rm = TRUE),
    continuous_age_concordant = sum(ageSD_concordant, na.rm = TRUE),
    continuous_age_concordant_p05 = sum(ageSD_concordant_p05, na.rm = TRUE),
    continuous_age_concordant_FDR05 = sum(ageSD_concordant & ageSD_padj < ALPHA, na.rm = TRUE),
    log2age_concordant = if ("log2age_concordant" %in% colnames(sig_cohort)) sum(log2age_concordant, na.rm = TRUE) else NA,
    excl16_20_same_direction = if ("excl16_20_same_direction" %in% colnames(sig_cohort)) sum(excl16_20_same_direction, na.rm = TRUE) else NA,
    interaction_tested = sum(!is.na(interaction_p)),
    interaction_nominal_p05 = sum(interaction_p < 0.05, na.rm = TRUE),
    interaction_genomewide_FDR05 = sum(interaction_padj < ALPHA, na.rm = TRUE),
    interaction_within_signature_FDR05 = sum(interaction_FDR_within_signature < ALPHA, na.rm = TRUE),
    interaction_genes_nominal = paste(sort(SYMBOL[interaction_p < 0.05 & !is.na(interaction_p)]), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(genome_wide_counts, by = "Cohort")

write.csv(summary_cohort, file.path(out_dir, "summary_by_cohort.csv"), row.names = FALSE)
print(as.data.frame(summary_cohort))

# Pooled concordance of the continuous-age effect (correlation across signature genes)
cor_tab <- sig_cohort %>%
  group_by(Cohort) %>%
  summarise(
    pearson_r_group_vs_minus_ageSD = cor(group_log2FC, -ageSD_log2FC, use = "complete.obs"),
    spearman_rho_group_vs_minus_ageSD = cor(group_log2FC, -ageSD_log2FC, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )
write.csv(cor_tab, file.path(out_dir, "group_vs_continuous_age_correlation.csv"), row.names = FALSE)


# 02. Sample structure: sex x group cells and age distribution ----

meta_all <- bind_rows(lapply(cohort_ids, function(gse) {
  read.csv(here::here("results", "deseq2", gse, "metadata", paste0(gse, "_sample_metadata.csv")),
           stringsAsFactors = FALSE) %>%
    transmute(Cohort = gse, GEO_ID, Age = as.numeric(Age), Sex, Group)
}))

sex_group <- meta_all %>% count(Cohort, Group, Sex) %>%
  pivot_wider(names_from = Sex, values_from = n, values_fill = 0)
write.csv(sex_group, file.path(out_dir, "sex_by_group_cell_sizes.csv"), row.names = FALSE)

age_bands <- meta_all %>%
  mutate(Band = cut(Age, c(-Inf, 5, 12, 16, 18, 20, 30, 50, 65, Inf), right = FALSE,
                    labels = c("<5", "5-11", "12-15", "16-17", "18-19", "20-29", "30-49", "50-64", ">=65"))) %>%
  count(Cohort, Group, Band) %>%
  pivot_wider(names_from = Band, values_from = n, values_fill = 0)
write.csv(age_bands, file.path(out_dir, "age_bands_by_cohort_group.csv"), row.names = FALSE)

age_summary <- meta_all %>% group_by(Cohort, Group) %>%
  summarise(n = n(), median_age = median(Age), q1 = quantile(Age, 0.25), q3 = quantile(Age, 0.75),
            min_age = min(Age), max_age = max(Age),
            n_age_16_to_20 = sum(Age >= 16 & Age < 21), .groups = "drop")
write.csv(age_summary, file.path(out_dir, "age_summary_by_cohort_group.csv"), row.names = FALSE)
print(as.data.frame(age_summary))


# 03. Figures ----

p_age <- ggplot(meta_all, aes(x = Age, y = Cohort, colour = Group)) +
  geom_vline(xintercept = 18, linetype = 2, colour = "grey40") +
  annotate("rect", xmin = 16, xmax = 21, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "black") +
  geom_point(position = position_jitter(height = 0.18, width = 0, seed = 1), size = 1.6, alpha = 0.75) +
  scale_colour_manual(values = c(Adult = "#2C6FB7", Pediatric = "#E8B500")) +
  labs(x = "Age (years)", y = NULL, colour = NULL,
       caption = "Dashed line: 18-year cut-off; shaded band: 16-20 years (excluded in sensitivity analysis)") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

ggsave(file.path(fig_dir, "FigS_age_distribution_by_cohort.pdf"), p_age, width = 8, height = 3.6)
ggsave(file.path(fig_dir, "FigS_age_distribution_by_cohort.tiff"), p_age, width = 8, height = 3.6, dpi = 600, compression = "lzw")

order_genes <- signature %>% arrange(RE_log2FC) %>% pull(SYMBOL)

tile_df <- bind_rows(
  sig_cohort %>% transmute(SYMBOL, Cohort, Model = "Group\n(Pediatric vs Adult)", z = group_log2FC / group_SE, p = group_p),
  sig_cohort %>% transmute(SYMBOL, Cohort, Model = "Continuous age\n(sign reversed)", z = -ageSD_log2FC / ageSD_SE, p = ageSD_p),
  sig_cohort %>% transmute(SYMBOL, Cohort, Model = "Sex x Group\ninteraction", z = interaction_log2FC / interaction_SE, p = interaction_p)
) %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = order_genes),
    Model = factor(Model, levels = c("Group\n(Pediatric vs Adult)", "Continuous age\n(sign reversed)", "Sex x Group\ninteraction")),
    z_plot = pmax(pmin(z, 5), -5),
    star = ifelse(!is.na(p) & p < 0.05, "*", "")
  )

p_tile <- ggplot(tile_df, aes(x = Cohort, y = SYMBOL, fill = z_plot)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_gradient2(low = "#2C6FB7", mid = "white", high = "#E8B500", midpoint = 0, limits = c(-5, 5),
                       name = "Wald z\n(+ = pediatric\ndirection)") +
  labs(x = NULL, y = NULL, caption = "* nominal p < 0.05; z truncated at +/-5") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_signature_group_continuous_interaction_heatmap.pdf"), p_tile, width = 8.5, height = 9.5)
ggsave(file.path(fig_dir, "FigS_signature_group_continuous_interaction_heatmap.tiff"), p_tile, width = 8.5, height = 9.5, dpi = 600, compression = "lzw")

p_sc <- ggplot(sig_cohort, aes(x = group_log2FC, y = -ageSD_log2FC)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = ageSD_p < 0.05), size = 1.8) +
  geom_text(data = cor_tab, aes(x = -Inf, y = Inf, label = paste0("r = ", sprintf("%.2f", pearson_r_group_vs_minus_ageSD))),
            hjust = -0.1, vjust = 1.4, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ Cohort, scales = "free") +
  scale_colour_manual(values = c(`TRUE` = "#2C6FB7", `FALSE` = "#BDBDBD"), labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = "Continuous age") +
  labs(x = "Pediatric vs Adult log2FC", y = "-(log2FC per 1 SD of age)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_signature_group_vs_continuous_age.pdf"), p_sc, width = 9, height = 3.8)
ggsave(file.path(fig_dir, "FigS_signature_group_vs_continuous_age.tiff"), p_sc, width = 9, height = 3.8, dpi = 600, compression = "lzw")

writeLines(c(
  paste("Script build:", SCRIPT_BUILD),
  paste0(summary_cohort$Cohort, ": continuous-age concordant ", summary_cohort$continuous_age_concordant, "/",
         summary_cohort$signature_genes_in_cohort, " (p<0.05: ", summary_cohort$continuous_age_concordant_p05,
         "); interaction nominal p<0.05: ", summary_cohort$interaction_nominal_p05, "/", summary_cohort$interaction_tested,
         "; interaction FDR<0.05 (within signature): ", summary_cohort$interaction_within_signature_FDR05,
         "; genome-wide interaction FDR<0.05 genes: ", summary_cohort$genome_interaction_FDR05)
), file.path(out_dir, "run_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done: ", out_dir)
