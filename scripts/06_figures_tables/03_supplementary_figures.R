# Supplementary Figures S1-S17
#
# Copies the figures written by the analysis scripts under standard names
# (FigS01_..., numbered in order of first citation in the manuscript), draws
# the single-cell secondary heatmap and writes an index with titles and legends.
# Figures whose source files do not exist yet are listed as missing.
#
# Output: results/manuscript/supplementary_figures/

SCRIPT_BUILD <- "SUPPLEMENTARY_FIGURES_2026-10-28"
message("Running script build: ", SCRIPT_BUILD)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise

COL_PED <- "#E8B500"
COL_ADULT <- "#2C6FB7"

out_dir <- here::here("results", "manuscript", "supplementary_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

R <- function(...) here::here("results", ...)
scp <- function(f) R("scRNAseq_COVID_nasal_age", "pseudobulk_signature", "tables", f)
new_dir <- R("scRNAseq_COVID_nasal_age", "pseudobulk_signature", "figures")
dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)


# 01. Additional single-cell supplementary figures (from saved tables) ----

if (file.exists(scp("SECONDARY_mean_lognorm_wilcoxon_signature_genes.csv"))) {
  sec <- read.csv(scp("SECONDARY_mean_lognorm_wilcoxon_signature_genes.csv"), stringsAsFactors = FALSE)
  prim <- read.csv(R("meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"), stringsAsFactors = FALSE)
  ord <- prim %>% filter(RE_FDR < 0.05) %>% arrange(RE_log2FC) %>% pull(SYMBOL)
  sec <- sec %>% mutate(SYMBOL = factor(SYMBOL, levels = ord[ord %in% SYMBOL]),
                        star = ifelse(FDR_global < 0.05, "*", ""))
  p <- ggplot(sec, aes(cell_type, SYMBOL, fill = cliffs_delta)) +
    geom_tile(colour = "white") + geom_text(aes(label = star), size = 3, vjust = 0.75) +
    scale_fill_gradient2(low = COL_ADULT, mid = "white", high = COL_PED, limits = c(-1, 1), name = "Cliff's delta\n(+ = pediatric)") +
    labs(x = NULL, y = NULL, caption = "Donor-level mean log-normalized expression; Wilcoxon rank-sum; * FDR (global) < 0.05") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "italic"), panel.grid = element_blank())
  ggsave(file.path(new_dir, "sc_secondary_mean_lognorm_heatmap.pdf"), p, width = 7, height = 8.5)
  ggsave(file.path(new_dir, "sc_secondary_mean_lognorm_heatmap.tiff"), p, width = 7, height = 8.5, dpi = 600, compression = "lzw")
}

if (file.exists(scp("celltype_proportions_by_donor.csv"))) {
  prop <- read.csv(scp("celltype_proportions_by_donor.csv"), stringsAsFactors = FALSE)
  ab <- read.csv(scp("celltype_abundance_edgeR.csv"), stringsAsFactors = FALSE)
  prop <- prop %>% filter(cell_type %in% ab$cell_type) %>%
    left_join(ab %>% select(cell_type, FDR), by = "cell_type") %>%
    mutate(Group = factor(ifelse(Group == "Ped", "Pediatric", "Adult"), levels = c("Pediatric", "Adult")),
           panel = sprintf("%s\nFDR = %.2g", cell_type, FDR))
  p <- ggplot(prop, aes(Group, 100 * proportion, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), size = 0.8) +
    facet_wrap(~ panel, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = c(Pediatric = COL_PED, Adult = COL_ADULT), guide = "none") +
    labs(x = NULL, y = "% of donor's cells", caption = "edgeR quasi-likelihood on cell counts (offset: total cells per donor)") +
    theme_bw(base_size = 8) + theme(strip.text = element_text(size = 6.5))
  ggsave(file.path(new_dir, "sc_celltype_abundance.pdf"), p, width = 10, height = 6)
  ggsave(file.path(new_dir, "sc_celltype_abundance.tiff"), p, width = 10, height = 6, dpi = 600, compression = "lzw")
}


# 02. Supplementary figure list ----

FS <- function(...) here::here("results", "manuscript", "figures_for_supplement", ...)

sf <- list(
  # Figure S1
  list("Age_distribution", R("sensitivity", "continuous_age_sex_interaction", "figures", "FigS_age_distribution_by_cohort"),
       "Age Distribution of the Three Bulk RNA-seq Cohorts",
       "Each point represents one SARS-CoV-2-positive participant (n = 225), plotted by age in years (x-axis) and cohort (y-axis; GSE172274, n = 18; GSE179277, n = 83; GSE231409, n = 124); points are jittered vertically to reduce overlap. Point colour indicates the age group assigned from recorded age (yellow, pediatric, < 18 years; blue, adult, >= 18 years). The vertical dashed line marks the 18-year cut-off, and the grey shaded band marks the 16-20-year window excluded in the corresponding sensitivity analysis. In GSE231409, 9 of the 20 adults were aged 18.3-19.5 years. Source data: Supplementary Tables S1 and S14d."),
  # Figure S2
  list("Robust_thresholds", R("meta_analysis", "inference_methods", "figures", "FigS_robust_threshold_grid"),
       "Number of Robust Genes Under Alternative Classification Thresholds",
       "Heatmaps of the number of genes classified as Robust (number printed in each cell; darker blue = more genes) as a function of the maximum I2 allowed (y-axis: 75%, 50%, 25%) and of the maximum leave-one-cohort-out p value required in every leave-one-out model (x-axis: 0.10, 0.05, 0.01). Right panel: primary definition (REML with z test, FDR < 0.05); left panel: additionally requiring a nominal Hartung-Knapp-Sidik-Jonkman p < 0.05. The primary threshold combination (I2 < 75%, leave-one-out p < 0.05) yields 26 Robust genes. Source data: Supplementary Table S10."),
  # Figure S3
  list("Signature_retention", FS("SuppFig_signature_retention"),
       "Retention of the 44-Gene Signature Across Meta-Analytic Inference Methods and Clinical Sensitivity Analyses",
       "Horizontal bars show the number of the 44 signature genes (x-axis) retained under each alternative analysis (y-axis). The grey bar indicates genes with the same direction of effect as in the primary analysis, and the overlaid blue bar indicates genes that were also statistically supported; the number of supported genes is printed at the end of each blue bar. Upper panel (Inference): REML with Wald z test (primary), DerSimonian-Laird and Paule-Mandel estimators with z tests (support = FDR < 0.05 across 18,209 genes), REML with the Hartung-Knapp-Sidik-Jonkman adjustment (support = nominal p < 0.05), and Bayesian random-effects models with half-normal(0.5), half-Cauchy(0.5), and half-normal(1.0) priors on tau (support = 95% credible interval excluding zero). Lower panel (Clinical/age): adjustment for viral load and respiratory severity; further adjustment for symptom status with exclusion of adults with severity score >= 3; exclusion of participants aged 16-20 years; and continuous modelling of age as log2(age + 1) (support = FDR < 0.05). The dotted vertical line marks 44 genes. Source data: Supplementary Tables S11, S13, and S13b."),
  # Figure S4
  list("Inference_methods", R("meta_analysis", "inference_methods", "figures", "FigS_meta_inference_signature_tiles"),
       "Support for Each Signature Gene Across Meta-Analytic Inference Methods",
       "Tile plot of the 44 signature genes (rows) by inference method (columns). Blue tiles indicate that the gene was supported by the method and light grey tiles that it was not. Support was defined as FDR < 0.05 across all 18,209 eligible genes for the z-based methods (REML, DerSimonian-Laird [DL], Paule-Mandel [PM]) and the fixed-effect model; nominal p < 0.05 for the Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment (t distribution with k - 1 degrees of freedom; standard and ad hoc variants, and PM with HKSJ); and a 95% credible interval excluding zero for the Bayesian random-effects models with half-normal (HN) or half-Cauchy (HC) priors on tau with the indicated scale. Rows are grouped by the number of contributing cohorts (k = 2, upper block; k = 3, lower block) and ordered by pooled effect. Source data: Supplementary Table S11."),
  # Figure S5
  list("Inference_intervals", R("meta_analysis", "inference_methods", "figures", "FigS_meta_inference_intervals"),
       "Pooled Estimates and 95% Intervals of the Signature Genes Under Alternative Inference Methods",
       "For each signature gene (y-axis), the point shows the pooled log2 fold change (pediatric vs adult; x-axis) and the horizontal line its 95% interval, coloured by method: REML with Wald z test (primary; black), REML with the Hartung-Knapp-Sidik-Jonkman adjustment (orange), its ad hoc variant (yellow), and the Bayesian random-effects model with a half-normal(0.5) prior on tau (blue; posterior median and 95% credible interval). The vertical dashed line marks zero (no age difference); intervals that do not cross it indicate support. Genes are grouped by the number of contributing cohorts (k = 2, upper panel; k = 3, lower panel). Intervals were truncated at +/-6 for display; the very wide HKSJ intervals of k = 2 genes reflect the t distribution with one degree of freedom. Source data: Supplementary Table S11."),
  # Figure S6
  list("Viral_load", R("sensitivity", "clinical_covariates", "figures", "FigS_viral_load_by_age_group"),
       "SARS-CoV-2 Viral Load by Age Group",
       "Viral load expressed as log10(viral load + 1) in GSE179277 (left; SARS-CoV-2 reads per million; 46 adults and 37 pediatric samples) and GSE231409 (right; viral copies; available for 12 of 20 adults and 80 of 104 pediatric samples). GSE172274 did not report viral load. Boxes show the median (central line) and interquartile range; whiskers extend to the most extreme values within 1.5 times the interquartile range; box fill indicates the age group (blue, adult; yellow, pediatric), and each point is one sample. Differences between age groups were not significant (two-sided Wilcoxon rank-sum test, p = 0.32 and p = 0.70; Supplementary Table S13d)."),
  # Figure S7
  list("Clinical_sensitivity_scatter", R("sensitivity", "clinical_covariates", "figures", "FigS_signature_clinical_sensitivity_scatter"),
       "Pooled Effects of the Signature Genes Under Clinical Adjustment and Alternative Age Modelling",
       "Each panel compares the primary pooled log2 fold change (x-axis) of the 44 signature genes with the pooled estimate obtained in one sensitivity analysis (y-axis): adjustment for viral load (GSE179277, GSE231409) and respiratory severity (GSE172274); GSE231409 restricted to samples with viral load, without adjustment; extended adjustment (symptom status in GSE231409 and exclusion of adults with severity score >= 3 in GSE172274); exclusion of participants aged 16-20 years; and continuous log2(age + 1) (sign reversed so that positive values indicate the pediatric direction; different scale). Each point is one gene, coloured by the result of the sensitivity analysis: dark blue, same direction and FDR < 0.05; light blue, same direction and nominal p < 0.05; grey, same direction and not significant; orange, opposite direction. The dashed line is the identity line, and grey lines mark zero. The Pearson correlation (r) with the primary estimates and the number of genes with the same direction are shown in each panel. Source data: Supplementary Tables S13 and S13b."),
  # Figure S8
  list("Clinical_sensitivity_forest", R("sensitivity", "clinical_covariates", "figures", "FigS_signature_clinical_sensitivity_forest"),
       "Forest Plot of the Signature Genes in the Primary, Clinically Adjusted, and Age-Window-Excluded Analyses",
       "For each signature gene (y-axis), points and horizontal lines show the random-effects pooled log2 fold change (pediatric vs adult; x-axis) and its 95% confidence interval in the primary analysis (black), after adjustment for viral load and respiratory severity (orange), and after exclusion of participants aged 16-20 years (blue). Genes are ordered by the primary pooled effect; the vertical dashed line marks zero. Positive values indicate higher expression in pediatric samples. Source data: Supplementary Table S13."),
  # Figure S9
  list("Continuous_age_interaction", R("sensitivity", "continuous_age_sex_interaction", "figures", "FigS_signature_group_continuous_interaction_heatmap"),
       "Pediatric-Versus-Adult, Continuous-Age, and Sex-by-Group Effects of the Signature Genes in Each Cohort",
       "Heatmap of Wald z statistics for each signature gene (rows) in each cohort (columns) for three DESeq2 models: left, pediatric vs adult group effect; centre, continuous age per standard deviation (sign reversed so that positive values indicate higher expression at younger ages); right, Sex x Group interaction. Yellow indicates positive z values (pediatric direction) and blue negative z values (adult direction); colours are truncated at +/-5. Asterisks mark nominal p < 0.05. Empty (white) cells indicate that the gene was not tested in that cohort. Source data: Supplementary Tables S14 and S14b."),
  # Figure S10
  list("Group_vs_continuous_age", R("sensitivity", "continuous_age_sex_interaction", "figures", "FigS_signature_group_vs_continuous_age"),
       "Concordance Between the Dichotomous and Continuous Age Effects of the Signature Genes",
       "For each cohort (panels), each point is one signature gene, with the pediatric-vs-adult log2 fold change on the x-axis and the negative of the log2 fold change per 1 standard deviation of age on the y-axis (so that both axes are positive for genes with higher expression in younger individuals). Point colour indicates whether the continuous-age effect was nominally significant (blue, p < 0.05; grey, not significant). Grey lines mark zero, and the Pearson correlation (r) between the two effects is shown in each panel (r = 0.93, 0.99, and 0.94 in GSE172274, GSE179277, and GSE231409). Source data: Supplementary Table S14."),
  # Figure S11
  list("Age_by_infection", FS("SuppFig_age_by_infection"),
       "Pediatric-Versus-Adult Differences of the Signature Genes in SARS-CoV-2-Negative and SARS-CoV-2-Positive Individuals of GSE179277",
       "(A) For each of the 40 evaluable signature genes, the pediatric-vs-adult log2 fold change among SARS-CoV-2-positive individuals (x-axis; 46 adults and 37 pediatric) is plotted against the corresponding log2 fold change among SARS-CoV-2-negative individuals with no respiratory virus detected (y-axis; 81 adults and 34 pediatric), estimated with the DESeq2 model ~ Sex + Group + Infection + Group:Infection. Vertical grey lines are 95% confidence intervals of the effect in negative individuals. Points are coloured by pattern: dark blue, same direction and p < 0.05 also without SARS-CoV-2; light blue, same direction but not significant; grey, opposite direction; orange, infection-dependent (age x infection interaction FDR < 0.05 within the signature). Shaded quadrants indicate concordant pediatric-enriched (upper right, yellow) and adult-enriched (lower left, blue) effects; the dashed line is the identity line, and r is the Pearson correlation. Selected genes are labelled. (B) Number of signature genes in each pattern. (C) Expression (log2[normalized count + 1]) of six selected genes by infection status (x-axis) and age group (yellow, pediatric; blue, adult); boxes show the median and interquartile range, whiskers extend to 1.5 times the interquartile range, and points are individual samples. Source data: Supplementary Tables S15 and S15b."),
  # Figure S12
  list("Age_by_infection_genome", R("sensitivity", "age_by_infection", "figures", "FigS_age_effect_negative_vs_positive_genome_wide"),
       "Genome-Wide Pediatric-Versus-Adult Differences With and Without SARS-CoV-2 in GSE179277",
       "Two-dimensional histogram of the pediatric-vs-adult log2 fold change in SARS-CoV-2-positive individuals (x-axis) and in SARS-CoV-2-negative individuals (y-axis) for all 15,979 genes tested. Each tile is a bin, coloured by the number of genes it contains (log scale; darker = more genes). The dashed line is the identity line, solid grey lines mark zero, and r is the Pearson correlation across all genes. Source data: Supplementary Table S15c."),
  # Figure S13
  list("GSEA_overlap_pediatric", R("gsea", "themes", "figures", "FigS_gsea_leading_edge_overlap_Pediatric"),
       "Leading-Edge Overlap Among the Pediatric-Enriched Consistently Enriched Pathways",
       "Symmetric heatmap of the overlap coefficient (|A n B| / min[|A|, |B|]) between the leading-edge genes of each pair of the 46 pediatric-enriched pathways, pooled across cohorts. Colour intensity (white to gold, 0 to 1) indicates the overlap; the diagonal equals 1. Pathways (rows and columns, same order) are ordered by average-linkage hierarchical clustering of 1 - overlap; blocks of high overlap correspond to the non-redundant themes defined at overlap >= 0.5 (Supplementary Tables S18 and S18b). The largest block corresponds to the replication-dependent histone theme."),
  # Figure S14
  list("GSEA_overlap_adult", R("gsea", "themes", "figures", "FigS_gsea_leading_edge_overlap_Adult"),
       "Leading-Edge Overlap Among the Adult-Enriched Consistently Enriched Pathways",
       "Symmetric heatmap of the overlap coefficient between the pooled leading-edge genes of each pair of the 17 adult-enriched pathways (white to dark blue, 0 to 1). Pathways are ordered by average-linkage hierarchical clustering as in Supplementary Figure S13; themes were defined at overlap >= 0.5 (Supplementary Tables S18 and S18b)."),
  # Figure S15
  list("GSEA_curves", R("gsea", "figures_final", "supplementary_curves"),
       "Enrichment Curves of Selected Pathways in Each Cohort",
       "Pre-ranked GSEA enrichment plots in GSE172274, GSE179277, and GSE231409 (left to right) for (a) antimicrobial humoral immune response mediated by antimicrobial peptide (GO:BP), (b) bile acid metabolic process (GO:BP), (c) humoral immune response (GO:BP), (d) long-chain fatty acid metabolic process (GO:BP), and (e) xenobiotic metabolic process (GO:BP). In each plot, the x-axis is the rank of the gene in the cohort's DESeq2 Wald-statistic ranking (left, most pediatric-enriched; right, most adult-enriched); the green line is the running enrichment score; black vertical ticks mark the positions of the gene-set members; and the red dashed lines mark the maximum and minimum enrichment scores. The normalized enrichment score (NES) and Benjamini-Hochberg-adjusted p value (FDR) are given below each plot; positive NES indicates pediatric enrichment. Cohort labels are coloured by the direction of the pathway (orange, pediatric-enriched; blue, adult-enriched). Source data: Supplementary Table S16."),
  # Figure S16
  list("Single_cell_donor_level", FS("SuppFig_sc_donor_level_NTS_SCRN1"),
       "Donor-Level Pseudobulk Expression of NTS and SCRN1 in Secretory Cells",
       "Each point is one donor: pseudobulk expression (log2 counts per million of raw UMI counts summed over the donor's Secretory cells) of NTS (left) and SCRN1 (right) in pediatric (yellow; n = 16) and adult (blue; n = 9) donors. Boxes show the median and interquartile range, and whiskers extend to 1.5 times the interquartile range. The log2 fold change (pediatric vs adult) and FDR_global from the edgeR quasi-likelihood model (~ Sex + Group; Benjamini-Hochberg correction across all signature gene-by-cell-type tests) and the number of donors are shown above each panel. Source data: Supplementary Table S21."),
  # Figure S17
  list("Single_cell_secondary_heatmap", file.path(new_dir, "sc_secondary_mean_lognorm_heatmap"),
       "Secondary Single-Cell Analysis: Donor-Level Mean Log-Normalized Expression",
       "Heatmap of Cliff's delta for each signature gene (rows) in each cell type (columns), comparing the donor-level mean of the processed log-normalized expression between pediatric and adult donors (yellow, positive delta = higher in pediatric donors; blue, negative delta = higher in adult donors; range -1 to 1). Blank cells indicate gene-cell-type combinations that were not tested (fewer than 3 donors per age group with >= 20 cells). Asterisks mark FDR_global < 0.05 (two-sided Wilcoxon rank-sum test with Benjamini-Hochberg correction across all tests). This analysis averages log-normalized values and is not an aggregation of raw counts. Source data: Supplementary Table S22.")
)

index <- list()
for (i in seq_along(sf)) {
  e <- sf[[i]]
  id <- sprintf("FigS%02d", i)
  src <- e[[2]]
  copied <- character(0)
  if (dir.exists(src)) {
    files <- list.files(src, pattern = "\\.(pdf|tiff)$", full.names = TRUE)
    stems <- unique(tools::file_path_sans_ext(basename(files)))
    for (j in seq_along(stems)) for (ext in c("pdf", "tiff")) {
      f <- file.path(src, paste0(stems[j], ".", ext))
      if (file.exists(f)) {
        dest <- file.path(out_dir, paste0(id, letters[j], "_", e[[1]], ".", ext))
        file.copy(f, dest, overwrite = TRUE); copied <- c(copied, basename(dest))
      }
    }
  } else {
    for (ext in c("pdf", "tiff")) {
      f <- paste0(src, ".", ext)
      if (file.exists(f)) {
        dest <- file.path(out_dir, paste0(id, "_", e[[1]], ".", ext))
        file.copy(f, dest, overwrite = TRUE); copied <- c(copied, basename(dest))
      }
    }
  }
  index[[i]] <- data.frame(Figure = sub("FigS0?", "Figure S", id), File_id = id, Title = e[[3]], Legend = e[[4]],
                           Status = ifelse(length(copied) > 0, "ok", "missing"),
                           Files = paste(copied, collapse = "; "))
}
index <- bind_rows(index)
write.csv(index, file.path(out_dir, "supplementary_figure_index.csv"), row.names = FALSE)
print(index[, c("Figure", "Title", "Status")])
if (any(index$Status == "missing")) message("Missing figures (run the corresponding analysis first): ",
                                             paste(index$Figure[index$Status == "missing"], collapse = ", "))
message("Done: ", out_dir)
