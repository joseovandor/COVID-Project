# Pediatric-vs-adult differences with and without SARS-CoV-2 (GSE179277)
#
# GSE179277 also contains upper-airway samples in which no respiratory virus
# was detected ("no_virus"); these are analysed together with the
# SARS-CoV-2-positive samples of the same cohort:
#
#   ~ Sex + Group + Infection + Group:Infection
#   Group     : Adult (reference) / Pediatric (age < 18 y, from recorded age)
#   Infection : Negative (no_virus, reference) / SARS-CoV-2 positive
#
# Effects extracted (log2FC):
#   age_neg      pediatric vs adult in SARS-CoV-2-negative samples
#   age_pos      pediatric vs adult in SARS-CoV-2-positive samples
#   interaction  difference between the two age effects
#
# The negative samples come from patients tested for SARS-CoV-2, not from
# healthy volunteers. Samples with other respiratory viruses are excluded.
#
# Output: results/sensitivity/age_by_infection/

SCRIPT_BUILD <- "AGE_BY_INFECTION_2026-10-26"
message("Running script build: ", SCRIPT_BUILD)

# DESeq2 (and its Bioconductor dependencies) is attached FIRST so that dplyr
# verbs (count, select, rename, first, ...) are not masked.
suppressPackageStartupMessages({
  library(DESeq2)
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
gse <- "GSE179277"

out_dir <- here::here("results", "sensitivity", "age_by_infection")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

primary <- read.csv(here::here("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"),
                    stringsAsFactors = FALSE)
signature <- primary %>% filter(RE_FDR < ALPHA)


# 01. Metadata ----

raw_meta <- read.csv(here::here("results", "deseq2", gse, "metadata", paste0(gse, "_metadata_raw.csv")),
                     stringsAsFactors = FALSE, check.names = FALSE)

meta <- data.frame(
  GEO_ID = raw_meta$geo_accession,
  SampleTitle = raw_meta$title,
  Sample = sub("_(Peds|Adult)$", "", raw_meta$title, ignore.case = TRUE),
  Age = suppressWarnings(as.numeric(raw_meta[["age:ch1"]])),
  SexRaw = raw_meta[["gender:ch1"]],
  InfectionRaw = raw_meta[["infection:ch1"]],
  ViralLoadRPM = suppressWarnings(as.numeric(raw_meta[["sars-cov2_rpm:ch1"]])),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Sex = case_when(toupper(SexRaw) %in% c("M", "MALE") ~ "Male",
                    toupper(SexRaw) %in% c("F", "FEMALE") ~ "Female", TRUE ~ NA_character_),
    Infection = case_when(grepl("SC2", InfectionRaw, ignore.case = TRUE) ~ "Positive",
                          grepl("no_virus", InfectionRaw, ignore.case = TRUE) ~ "Negative",
                          TRUE ~ NA_character_),
    Group = ifelse(Age < 18, "Pediatric", "Adult")
  ) %>%
  filter(!is.na(Infection), !is.na(Age), !is.na(Sex))

meta$Group <- factor(meta$Group, levels = c("Adult", "Pediatric"))
meta$Infection <- factor(meta$Infection, levels = c("Negative", "Positive"))
meta$Sex <- factor(meta$Sex)


# 02. Raw counts ----

counts_file <- list.files(here::here("data", "raw", gse), pattern = "counts.*\\.csv(\\.gz)?$", full.names = TRUE)[1]
if (is.na(counts_file)) stop("Raw count file for ", gse, " not found under data/raw/", gse)

raw_counts <- read.csv(counts_file, stringsAsFactors = FALSE, check.names = FALSE)
name_col <- if ("gene_name" %in% colnames(raw_counts)) "gene_name" else NA

# The gene-identifier column is unnamed in the GEO file: always use position 1.
raw_counts <- raw_counts[!duplicated(raw_counts[[1]]), ]
gene_ids <- as.character(raw_counts[[1]])
gene_names <- if (!is.na(name_col)) raw_counts[[name_col]] else NA

missing_samples <- setdiff(meta$Sample, colnames(raw_counts))
if (length(missing_samples) > 0) {
  message("Samples without counts (excluded): ", length(missing_samples))
  meta <- meta[!meta$Sample %in% missing_samples, ]
}

counts <- as.matrix(raw_counts[, meta$Sample])
rownames(counts) <- gene_ids
storage.mode(counts) <- "integer"
rownames(meta) <- meta$Sample

cell_sizes <- meta %>% count(Infection, Group, Sex) %>%
  pivot_wider(names_from = Sex, values_from = n, values_fill = 0) %>%
  mutate(Total = Female + Male)
write.csv(cell_sizes, file.path(out_dir, "sample_sizes_infection_by_group_by_sex.csv"), row.names = FALSE)
print(cell_sizes)

# Prevalence filter, same rule as the primary analysis: raw count >= 10 in at
# least as many samples as the smallest Group x Infection cell.
min_cell <- min(table(meta$Group, meta$Infection))
keep <- rowSums(counts >= 10) >= min_cell
counts <- counts[keep, ]
message("Genes retained: ", nrow(counts), " (count >= 10 in >= ", min_cell, " samples)")
if (nrow(counts) < 1000) stop("Too few genes retained; check the raw count import.")


# 03. DESeq2 interaction model ----

dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ Sex + Group + Infection + Group:Infection)
dds <- DESeq(dds, quiet = TRUE)
print(resultsNames(dds))

int_name <- "GroupPediatric.InfectionPositive"

res_to_df <- function(res, label) {
  data.frame(Original_ID = rownames(res), baseMean = res$baseMean,
             log2FC = res$log2FoldChange, SE = res$lfcSE, pvalue = res$pvalue, padj = res$padj,
             effect = label, stringsAsFactors = FALSE)
}

effects <- bind_rows(
  res_to_df(results(dds, name = "Group_Pediatric_vs_Adult"), "age_neg"),
  res_to_df(results(dds, contrast = list(c("Group_Pediatric_vs_Adult", int_name))), "age_pos"),
  res_to_df(results(dds, name = int_name), "interaction"),
  res_to_df(results(dds, name = "Infection_Positive_vs_Negative"), "infection_in_adults"),
  res_to_df(results(dds, contrast = list(c("Infection_Positive_vs_Negative", int_name))), "infection_in_pediatric")
)

# Gene symbols: harmonized primary annotation when available, otherwise GEO gene_name
annot <- read.csv(here::here("results", "deseq2", gse, "primary", paste0(gse, "_meta_input.csv")),
                  stringsAsFactors = FALSE, colClasses = c(Original_ID = "character")) %>%
  filter(MetaEligibleAnnotation %in% c(TRUE, "TRUE")) %>%
  select(Original_ID, SYMBOL, primary_log2FC = log2FoldChange)

name_map <- data.frame(Original_ID = gene_ids, gene_name = gene_names, stringsAsFactors = FALSE)

effects <- effects %>%
  left_join(annot, by = "Original_ID") %>%
  left_join(name_map, by = "Original_ID") %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), gene_name, SYMBOL))

write.csv(effects, file.path(out_dir, "age_by_infection_all_genes_long.csv"), row.names = FALSE)

wide <- effects %>%
  select(Original_ID, SYMBOL, primary_log2FC, baseMean, effect, log2FC, SE, pvalue, padj) %>%
  pivot_wider(names_from = effect, values_from = c(log2FC, SE, pvalue, padj), names_glue = "{effect}_{.value}")


# 04. Genome-wide comparison of age effects in negative vs positive individuals ----

gw <- wide %>% filter(is.finite(age_neg_log2FC), is.finite(age_pos_log2FC))

genome_summary <- data.frame(
  genes_tested = nrow(gw),
  age_neg_FDR05 = sum(gw$age_neg_padj < ALPHA, na.rm = TRUE),
  age_pos_FDR05 = sum(gw$age_pos_padj < ALPHA, na.rm = TRUE),
  interaction_FDR05 = sum(gw$interaction_padj < ALPHA, na.rm = TRUE),
  pearson_age_neg_vs_age_pos = cor(gw$age_neg_log2FC, gw$age_pos_log2FC),
  pearson_age_pos_vs_primary_cohort_model = cor(gw$age_pos_log2FC, gw$primary_log2FC, use = "complete.obs")
)
write.csv(genome_summary, file.path(out_dir, "genome_wide_summary.csv"), row.names = FALSE)
print(genome_summary)


# 05. Signature genes ----

sig <- wide %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  left_join(signature %>% select(SYMBOL, k, RE_log2FC, RE_FDR), by = "SYMBOL") %>%
  mutate(
    interaction_FDR_within_signature = p.adjust(interaction_pvalue, method = "BH"),
    age_neg_FDR_within_signature = p.adjust(age_neg_pvalue, method = "BH"),
    age_neg_same_direction = sign(age_neg_log2FC) == sign(RE_log2FC),
    age_pos_same_direction = sign(age_pos_log2FC) == sign(RE_log2FC),
    ratio_neg_to_pos = age_neg_log2FC / age_pos_log2FC,
    Pattern = case_when(
      interaction_FDR_within_signature < ALPHA ~ "Age effect differs by infection status (interaction FDR < 0.05)",
      age_neg_same_direction & age_neg_pvalue < 0.05 ~ "Age difference also present without SARS-CoV-2 (p < 0.05)",
      age_neg_same_direction ~ "Same direction without SARS-CoV-2, not significant",
      TRUE ~ "Not observed without SARS-CoV-2 (opposite direction)"
    )
  ) %>%
  arrange(Pattern, RE_FDR)

write.csv(sig, file.path(out_dir, "signature_genes_age_by_infection.csv"), row.names = FALSE)

sig_summary <- sig %>% count(Pattern, name = "n_genes") %>%
  mutate(genes = vapply(Pattern, function(p) paste(sort(sig$SYMBOL[sig$Pattern == p]), collapse = ";"), character(1)))
write.csv(sig_summary, file.path(out_dir, "signature_pattern_summary.csv"), row.names = FALSE)
print(as.data.frame(sig_summary[, 1:2]))

sig_cor <- cor(sig$age_neg_log2FC, sig$age_pos_log2FC, use = "complete.obs")
missing_sig <- setdiff(signature$SYMBOL, sig$SYMBOL)


# 06. Figures ----

p_sig <- ggplot(sig, aes(x = age_pos_log2FC, y = age_neg_log2FC)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = age_neg_log2FC - 1.96 * age_neg_SE, ymax = age_neg_log2FC + 1.96 * age_neg_SE),
                colour = "grey75", width = 0) +
  geom_point(aes(colour = Pattern), size = 2.2) +
  geom_text(aes(label = SYMBOL), size = 2.4, vjust = -0.8, check_overlap = TRUE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.3,
           label = paste0("r = ", sprintf("%.2f", sig_cor), " (", nrow(sig), " signature genes)")) +
  scale_colour_manual(values = c(
    "Age difference also present without SARS-CoV-2 (p < 0.05)" = "#2C6FB7",
    "Same direction without SARS-CoV-2, not significant" = "#9EBCDA",
    "Not observed without SARS-CoV-2 (opposite direction)" = "#BDBDBD",
    "Age effect differs by infection status (interaction FDR < 0.05)" = "#D55E00")) +
  labs(x = "Pediatric vs Adult log2FC, SARS-CoV-2-positive",
       y = "Pediatric vs Adult log2FC, SARS-CoV-2-negative\n(95% CI)", colour = NULL,
       title = "GSE179277: age effect with and without SARS-CoV-2") +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "FigS_age_effect_negative_vs_positive_signature.pdf"), p_sig, width = 7, height = 7.5)
ggsave(file.path(fig_dir, "FigS_age_effect_negative_vs_positive_signature.tiff"), p_sig, width = 7, height = 7.5, dpi = 600, compression = "lzw")

p_gw <- ggplot(gw, aes(x = age_pos_log2FC, y = age_neg_log2FC)) +
  geom_bin2d(bins = 120) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey30") +
  scale_fill_gradient(low = "#DCE6F2", high = "#08306B", trans = "log10") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
  annotate("text", x = -4, y = 4, hjust = 0, vjust = 1, size = 3.3,
           label = paste0("r = ", sprintf("%.2f", genome_summary$pearson_age_neg_vs_age_pos), " (", nrow(gw), " genes)")) +
  labs(x = "Pediatric vs Adult log2FC, SARS-CoV-2-positive", y = "Pediatric vs Adult log2FC, SARS-CoV-2-negative", fill = "Genes") +
  theme_bw(base_size = 10)

ggsave(file.path(fig_dir, "FigS_age_effect_negative_vs_positive_genome_wide.pdf"), p_gw, width = 6, height = 5)
ggsave(file.path(fig_dir, "FigS_age_effect_negative_vs_positive_genome_wide.tiff"), p_gw, width = 6, height = 5, dpi = 600, compression = "lzw")

# Expression of NTS, SCRN1 and the signature genes with the strongest age effect
norm <- counts(dds, normalized = TRUE)

# Sample-level expression of every signature gene (used by scripts/06_figures_tables)
sig_expr_long <- bind_rows(lapply(sig$Original_ID, function(id) {
  data.frame(SYMBOL = sig$SYMBOL[sig$Original_ID == id], Sample = colnames(norm),
             log2_norm_count = log2(norm[id, ] + 1),
             Group = as.character(meta[colnames(norm), "Group"]),
             Infection = as.character(meta[colnames(norm), "Infection"]),
             Sex = as.character(meta[colnames(norm), "Sex"]),
             Age = meta[colnames(norm), "Age"])
}))
write.csv(sig_expr_long, file.path(out_dir, "signature_genes_expression_by_sample_long.csv"), row.names = FALSE)
show_genes <- unique(c("NTS", "SCRN1", sig %>% arrange(age_pos_pvalue) %>% pull(SYMBOL) %>% head(6)))
show_ids <- sig$Original_ID[match(show_genes, sig$SYMBOL)]
show_ids <- show_ids[!is.na(show_ids)]

expr_df <- bind_rows(lapply(show_ids, function(id) {
  data.frame(SYMBOL = sig$SYMBOL[sig$Original_ID == id], expr = log2(norm[id, ] + 1),
             Group = meta[colnames(norm), "Group"], Infection = meta[colnames(norm), "Infection"])
})) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = show_genes),
         Infection = factor(Infection, levels = c("Negative", "Positive"),
                            labels = c("SARS-CoV-2 negative", "SARS-CoV-2 positive")),
         Group = factor(Group, levels = c("Pediatric", "Adult")))

p_box <- ggplot(expr_df, aes(x = Infection, y = expr, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(aes(colour = Group), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1),
             size = 0.7, alpha = 0.6) +
  facet_wrap(~ SYMBOL, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Pediatric = "#E8B500", Adult = "#2C6FB7")) +
  scale_colour_manual(values = c(Pediatric = "#B08800", Adult = "#1D4E85")) +
  labs(x = NULL, y = "log2(normalized count + 1)", fill = NULL, colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1),
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_selected_genes_by_age_and_infection.pdf"), p_box, width = 10, height = 5.5)
ggsave(file.path(fig_dir, "FigS_selected_genes_by_age_and_infection.tiff"), p_box, width = 10, height = 5.5, dpi = 600, compression = "lzw")


# 07. Summary ----

writeLines(c(
  paste("Script build:", SCRIPT_BUILD),
  paste("Samples:", nrow(meta), "| genes tested:", nrow(counts)),
  capture.output(print(as.data.frame(cell_sizes))),
  "",
  paste("Genome-wide: age effect FDR<0.05 in negatives =", genome_summary$age_neg_FDR05,
        "; in positives =", genome_summary$age_pos_FDR05, "; age x infection FDR<0.05 =", genome_summary$interaction_FDR05),
  paste("Genome-wide correlation of age effects (negative vs positive): r =", sprintf("%.3f", genome_summary$pearson_age_neg_vs_age_pos)),
  "",
  paste("Signature genes evaluated:", nrow(sig), "; not in this model:", paste(missing_sig, collapse = ", ")),
  paste("Signature correlation of age effects (negative vs positive): r =", sprintf("%.3f", sig_cor)),
  paste0("  ", sig_summary$Pattern, ": ", sig_summary$n_genes)
), file.path(out_dir, "run_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done: ", out_dir)
