# Figures 1-4 and the data behind every panel
#
#   Figure 1  PCA of each cohort; forest plot of the signature genes
#   Figure 2  leave-one-cohort-out analysis (k = 3 genes)
#   Figure 3  cross-cohort GSEA, one row per non-redundant theme
#   Figure 4  single cell: UMAP, epithelial composition and expression of the
#             signature genes by cell type and age group
#
# Also writes three figures used in the supplement (signature retention across
# analyses, age effect with and without SARS-CoV-2, donor-level NTS and SCRN1).
#
# Output: results/manuscript/main_figures/, figures_for_supplement/, source_data/

SCRIPT_BUILD <- "MAIN_FIGURES_2026-10-28"
message("Running script build: ", SCRIPT_BUILD)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count <- dplyr::count; first <- dplyr::first; mutate <- dplyr::mutate
summarise <- dplyr::summarise; arrange <- dplyr::arrange; distinct <- dplyr::distinct

ALPHA <- 0.05
COHORT_ORDER <- c("GSE172274", "GSE179277", "GSE231409")

# Palette of the first-round figures
COHORT_COLOR <- c(GSE172274 = "#f0027f", GSE179277 = "#bf5b17", GSE231409 = "#7fc97f")
COHORT_SHAPE <- c(GSE172274 = 16, GSE179277 = 17, GSE231409 = 15)
DIV_PED <- "#FCCE24"
DIV_ADULT <- "#3E76BC"
PCA_PED <- "#E9A23B"
PCA_ADULT <- "#3B6EA8"
SC_PED <- "#E9C46A"
SC_ADULT <- "#457B9D"
CONCORDANT_COLOR <- "#7fc97f"
DISCORDANT_COLOR <- "#D95F0E"
HEAT_COLS <- c("#3B6FB6", "#6F91CF", "#B8C9E6", "#F7F7F7", "#F6E8B1", "#F5D66E", "#FFD43B")
TAB10 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
           "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5")
EPITHELIAL <- c("Basal", "Basal cycling", "Ciliated", "Deuterosomal", "Multiciliated", "Club", "Duct",
                "Goblet", "Goblet inflammatory", "Hillock", "Hillock cycling", "Rare epi",
                "Secretory", "Squamous", "Transit epi")

fig_dir <- here::here("results", "manuscript", "main_figures")
sup_dir <- here::here("results", "manuscript", "figures_for_supplement")
src_dir <- here::here("results", "manuscript", "source_data")
for (d in c(fig_dir, sup_dir, src_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, w, h, dir = fig_dir) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  ggsave(file.path(dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", device = "tiff")
  message("  saved ", name)
}
src <- function(d, name) write.csv(d, file.path(src_dir, paste0(name, ".csv")), row.names = FALSE)
P <- function(...) here::here("results", ...)
rd <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
fmt_fdr <- function(x) ifelse(is.na(x), "", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
stars <- function(fdr) dplyr::case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**", fdr < 0.05 ~ "*", TRUE ~ "")
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 24))

primary <- rd(P("meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"))
signature <- primary %>% filter(RE_FDR < ALPHA) %>%
  mutate(direction = factor(ifelse(RE_log2FC > 0, "Pediatric-enriched", "Adult-enriched"),
                            levels = c("Pediatric-enriched", "Adult-enriched")))


# FIGURE 1 ----
message("Figure 1")

pca_plots <- lapply(seq_along(COHORT_ORDER), function(i) {
  g <- COHORT_ORDER[i]
  vst <- rd(P("deseq2", g, "primary", paste0(g, "_VST_matrix.csv")))
  m <- as.matrix(vst[, -1]); rownames(m) <- vst[[1]]
  rv <- apply(m, 1, var)
  top <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pc <- prcomp(t(m[top, ]))
  pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2))
  meta <- rd(P("deseq2", g, "metadata", paste0(g, "_sample_metadata.csv")))
  key <- if (all(colnames(m) %in% meta$Sample)) "Sample" else "GEO_ID"
  meta <- meta[match(colnames(m), meta[[key]]), ]
  d <- data.frame(Cohort = g, Sample = colnames(m), PC1 = pc$x[, 1], PC2 = pc$x[, 2], Group = meta$Group)
  src(d, paste0("Figure1_", LETTERS[i], "_PCA_", g))
  ggplot(d, aes(PC1, PC2, colour = Group)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_point(size = 2.6, alpha = 0.95) +
    scale_colour_manual(values = c(Adult = PCA_ADULT, Pediatric = PCA_PED), breaks = c("Adult", "Pediatric"), name = NULL) +
    labs(title = g, x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance")) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", size = 16), legend.text = element_text(size = 15))
})

df <- signature %>%
  group_by(direction) %>%
  arrange(abs(RE_log2FC), RE_FDR, .by_group = TRUE) %>%
  mutate(y = row_number()) %>%
  ungroup() %>%
  mutate(effect_ci = sprintf("%.2f [%.2f, %.2f]", RE_log2FC, RE_CI_lower, RE_CI_upper), fdr_text = fmt_fdr(RE_FDR))

cohort_long <- df %>%
  select(SYMBOL, direction, y, all_of(paste0(COHORT_ORDER, "_log2FC")), all_of(paste0(COHORT_ORDER, "_SE"))) %>%
  pivot_longer(-c(SYMBOL, direction, y), names_to = c("cohort", ".value"), names_pattern = "(GSE\\d+)_(log2FC|SE)") %>%
  filter(!is.na(log2FC), !is.na(SE)) %>%
  mutate(cohort = factor(cohort, levels = COHORT_ORDER),
         lo = log2FC - 1.96 * SE, hi = log2FC + 1.96 * SE,
         y_dodge = y + c(GSE172274 = 0.20, GSE179277 = 0, GSE231409 = -0.20)[as.character(cohort)])

src(df %>% select(SYMBOL, direction, k, RE_log2FC, RE_CI_lower, RE_CI_upper, RE_FDR,
                  all_of(paste0(COHORT_ORDER, "_log2FC")), all_of(paste0(COHORT_ORDER, "_SE"))), "Figure1_D_forest")

forest_abs <- min(ceiling(max(abs(c(cohort_long$lo, cohort_long$hi, df$RE_CI_lower, df$RE_CI_upper)), na.rm = TRUE)), 8)
cohort_long <- cohort_long %>% mutate(lo = pmax(lo, -forest_abs), hi = pmin(hi, forest_abs))
GENE_X <- -forest_abs - 0.85
X_EFF <- forest_abs + 0.95
X_FDR <- forest_abs + 6.10
X_MIN <- -forest_abs - 4.80
X_MAX <- forest_abs + 8.20
header_df <- df %>% group_by(direction) %>% summarise(y_header = max(y) + 1.2, .groups = "drop")
row_bg <- df %>% filter(y %% 2 == 0)

p_forest <- ggplot() +
  geom_rect(data = row_bg, aes(xmin = -Inf, xmax = Inf, ymin = y - 0.47, ymax = y + 0.47), fill = "#F7F7F7", colour = NA) +
  geom_vline(xintercept = 0, colour = "#8A8A8A", linewidth = 0.5, linetype = "22") +
  geom_text(data = df, aes(x = GENE_X, y = y, label = SYMBOL), hjust = 1, fontface = "italic", size = 5.6, colour = "#111111") +
  geom_errorbarh(data = cohort_long, aes(xmin = lo, xmax = hi, y = y_dodge, colour = cohort), height = 0.07, linewidth = 0.6, alpha = 0.68) +
  geom_point(data = cohort_long, aes(x = log2FC, y = y_dodge, colour = cohort, shape = cohort), size = 2.8, alpha = 0.96) +
  geom_errorbarh(data = df, aes(xmin = RE_CI_lower, xmax = RE_CI_upper, y = y), height = 0.16, linewidth = 1.0, colour = "#151515") +
  geom_point(data = df, aes(x = RE_log2FC, y = y, fill = direction), shape = 23, size = 5.0, stroke = 0.9, colour = "#111111") +
  geom_vline(xintercept = forest_abs + 0.5, colour = "#D0D0D0", linewidth = 0.48) +
  geom_text(data = df, aes(x = X_EFF, y = y, label = effect_ci), hjust = 0, size = 5.1, colour = "#222222") +
  geom_text(data = df, aes(x = X_FDR, y = y, label = fdr_text), hjust = 0, size = 5.1, colour = "#222222") +
  geom_text(data = header_df, aes(x = X_EFF, y = y_header, label = "RE log2FC [95% CI]"), hjust = 0, size = 5.5, fontface = "bold") +
  geom_text(data = header_df, aes(x = X_FDR, y = y_header, label = "FDR"), hjust = 0, size = 5.5, fontface = "bold") +
  facet_wrap(~ direction, ncol = 1, scales = "free_y",
             labeller = as_labeller(c("Pediatric-enriched" = "Pediatric-enriched genes", "Adult-enriched" = "Adult-enriched genes"))) +
  scale_colour_manual(name = "Cohort", values = COHORT_COLOR, breaks = COHORT_ORDER) +
  scale_shape_manual(name = "Cohort", values = COHORT_SHAPE, breaks = COHORT_ORDER) +
  scale_fill_manual(values = c("Pediatric-enriched" = DIV_PED, "Adult-enriched" = DIV_ADULT), guide = "none") +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.65, 1.75))) +
  scale_x_continuous(limits = c(X_MIN, X_MAX), breaks = pretty(c(-forest_abs, forest_abs), n = 7), expand = expansion(mult = 0)) +
  labs(x = expression(log[2] * " fold change (Pediatric vs Adult)"), y = NULL) +
  theme_classic(base_size = 16) +
  theme(strip.background = element_rect(fill = "#EFEFEF", colour = NA),
        strip.text = element_text(size = 18, face = "bold", margin = margin(t = 8, b = 8)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_text(size = 14), legend.position = "right",
        legend.title = element_text(face = "bold", size = 16), legend.text = element_text(size = 15),
        panel.spacing.y = unit(0.8, "cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5)))

pca_row <- (pca_plots[[1]] | pca_plots[[2]] | pca_plots[[3]]) + plot_layout(guides = "collect") & theme(legend.position = "top")
fig1 <- wrap_elements(full = pca_row) / p_forest + plot_layout(heights = c(1, 5.2)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig1, "Figure1_cohorts_and_signature", 15, 24)


# FIGURE 2 ----
message("Figure 2")

lin_ccc <- function(x, y) 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
loo <- rd(P("meta_analysis", "sensitivity", "leave_one_out_results.csv"))
scatter_df <- bind_rows(lapply(COHORT_ORDER, function(omit) {
  loo %>% filter(Omitted_dataset == omit, SYMBOL %in% signature$SYMBOL[signature$k == 3]) %>%
    left_join(signature %>% select(SYMBOL, held_out = all_of(paste0(omit, "_log2FC"))), by = "SYMBOL") %>%
    filter(!is.na(Pooled_log2FC), !is.na(held_out)) %>%
    transmute(SYMBOL, panel = omit, discovery = Pooled_log2FC, held_out, concordant = sign(Pooled_log2FC) == sign(held_out))
})) %>% mutate(panel = factor(panel, levels = COHORT_ORDER))
src(scatter_df, "Figure2_ABC_leave_one_cohort_out")

lim <- ceiling(max(abs(c(scatter_df$discovery, scatter_df$held_out))) * 1.10)
panel_stats <- scatter_df %>% group_by(panel) %>%
  summarise(r = cor(discovery, held_out), ccc = lin_ccc(discovery, held_out), conc = 100 * mean(concordant), .groups = "drop") %>%
  mutate(label = sprintf("Pearson r = %.2f\nLin CCC = %.2f\nDirection concordance = %.0f%%", r, ccc, conc))
discordant <- scatter_df %>% filter(!concordant)

p_loo_panel <- function(g, tag) {
  d <- scatter_df %>% filter(panel == g)
  s <- panel_stats %>% filter(panel == g)
  ggplot(d, aes(discovery, held_out)) +
    annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim, fill = DIV_PED, alpha = 0.10) +
    annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0, fill = DIV_ADULT, alpha = 0.10) +
    annotate("text", x = lim * 0.35, y = lim * 0.77, label = "Pediatric-enriched", size = 5.6, fontface = "bold", colour = "#8A6D00") +
    annotate("text", x = -lim * 0.50, y = -lim * 0.77, label = "Adult-enriched", size = 5.6, fontface = "bold", colour = "#244B7A") +
    geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.75) +
    geom_vline(xintercept = 0, colour = "#B8B8B8", linewidth = 0.75) +
    geom_abline(slope = 1, intercept = 0, colour = "#686868", linewidth = 1, linetype = "22") +
    geom_point(aes(fill = concordant, shape = concordant), size = 5.6, stroke = 0.75, colour = "white", alpha = 0.95) +
    geom_text(data = discordant %>% filter(panel == g), aes(label = SYMBOL), nudge_x = 0.22, nudge_y = -0.18, hjust = 0,
              size = 5, fontface = "italic", colour = "#8C3B00") +
    annotate("text", x = -lim * 0.92, y = lim * 0.91, hjust = 0, vjust = 1, size = 5.4, lineheight = 1.12, label = s$label) +
    scale_fill_manual(name = NULL, breaks = c("TRUE", "FALSE"), values = c(`TRUE` = CONCORDANT_COLOR, `FALSE` = DISCORDANT_COLOR),
                      labels = c(`TRUE` = "Concordant direction", `FALSE` = "Discordant direction"), drop = FALSE) +
    scale_shape_manual(name = NULL, breaks = c("TRUE", "FALSE"), values = c(`TRUE` = 21, `FALSE` = 24),
                       labels = c(`TRUE` = "Concordant direction", `FALSE` = "Discordant direction"), drop = FALSE) +
    facet_wrap(~ panel, labeller = as_labeller(setNames(paste("Held out:", COHORT_ORDER), COHORT_ORDER))) +
    scale_x_continuous(limits = c(-lim, lim), breaks = pretty(c(-lim, lim), n = 5), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(-lim, lim), breaks = pretty(c(-lim, lim), n = 5), expand = expansion(mult = 0)) +
    labs(x = if (g == COHORT_ORDER[2]) expression("LOO pooled " * log[2] * "FC (remaining two cohorts)") else NULL,
         y = if (g == COHORT_ORDER[1]) expression("Observed " * log[2] * "FC in held-out cohort") else NULL) +
    theme_classic(base_size = 18) +
    theme(aspect.ratio = 1, strip.background = element_rect(fill = "#ECECEC", colour = NA),
          strip.text = element_text(face = "bold", size = 20), legend.position = "bottom", legend.text = element_text(size = 16.5)) +
    guides(fill = guide_legend(override.aes = list(size = 6.2)), shape = guide_legend(override.aes = list(size = 6.2)))
}
scatter_df$concordant <- factor(scatter_df$concordant, levels = c(TRUE, FALSE))
discordant$concordant <- factor(discordant$concordant, levels = c(TRUE, FALSE))
fig2 <- (p_loo_panel(COHORT_ORDER[1]) | p_loo_panel(COHORT_ORDER[2]) | p_loo_panel(COHORT_ORDER[3])) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & tag_theme & theme(legend.position = "bottom")
save_fig(fig2, "Figure2_leave_one_cohort_out", 16.5, 7.4)


# FIGURE 3 ----
message("Figure 3")

themes_file <- P("gsea", "themes", "gsea_themes.csv")
if (file.exists(themes_file)) {
  themes <- rd(themes_file)
  gsea <- rd(P("gsea", "GSEA_all_cohorts_long.csv"))
  conc <- rd(P("gsea", "cross_cohort", "tables", "GSEA_pathway_concordance_complete.csv"))
  nice <- function(x) {
    x <- tolower(gsub("_", " ", sub("^(GOBP|REACTOME|KEGG|HALLMARK)_", "", x)))
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  themes$Label <- nice(themes$representative_pathway)
  lab_file <- here::here("config", "gsea_theme_labels.csv")
  if (file.exists(lab_file)) {
    lo <- rd(lab_file)
    i <- match(lo$Representative_pathway, themes$representative_pathway)
    themes$Label[i[!is.na(i)]] <- lo$Label[!is.na(i)]
  }
  themes <- themes %>%
    left_join(conc %>% select(Collection, pathway, N_LE_shared_all3, Shared_LE_all3),
              by = c("representative_collection" = "Collection", "representative_pathway" = "pathway")) %>%
    mutate(Label = vapply(Label, function(x) paste(strwrap(x, 40), collapse = "\n"), character(1)),
           N_LE_shared_all3 = ifelse(is.na(N_LE_shared_all3), 0, N_LE_shared_all3))

  tiles <- gsea %>%
    inner_join(themes %>% select(Theme, Direction, Label, n_pathways, representative_pathway, representative_collection),
               by = c("pathway" = "representative_pathway", "Collection" = "representative_collection")) %>%
    mutate(star = stars(padj), Dataset = factor(Dataset, levels = COHORT_ORDER))
  src(tiles %>% select(Theme, Direction, Label, pathway, n_pathways, Dataset, NES, padj), "Figure3_AC_GSEA_NES")
  src(themes %>% select(Theme, Direction, Label, representative_pathway, n_pathways, N_LE_shared_all3, Shared_LE_all3),
      "Figure3_BD_shared_leading_edge")

  row_order <- function(dir) themes %>% filter(Direction == dir) %>% arrange(n_pathways, abs(mean_abs_NES)) %>% pull(Label)
  heat <- function(dir) {
    d <- tiles %>% filter(Direction == dir) %>% mutate(Label = factor(Label, levels = row_order(dir)))
    ggplot(d, aes(Dataset, Label, fill = NES)) +
      geom_tile(colour = "white", linewidth = 1.4) +
      geom_text(aes(label = star), size = 6.5, vjust = 0.75) +
      scale_fill_gradient2(low = DIV_ADULT, mid = "#F7F7F7", high = DIV_PED, midpoint = 0, limits = c(-3, 3),
                           breaks = -3:3, oob = scales::squish, name = "NES") +
      scale_x_discrete(position = "top") +
      labs(x = NULL, y = NULL, title = paste0(dir, "-enriched pathways")) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(size = 16, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"), plot.title = element_text(hjust = 0.5, size = 19),
            legend.key.width = unit(2, "cm"))
  }
  lolli <- function(dir, col) {
    d <- themes %>% filter(Direction == dir) %>% mutate(Label = factor(Label, levels = row_order(dir)))
    ggplot(d, aes(N_LE_shared_all3, Label)) +
      geom_segment(aes(x = 0, xend = N_LE_shared_all3, yend = Label), colour = col, linewidth = 1.3) +
      geom_point(colour = col, size = 5) +
      geom_text(aes(label = N_LE_shared_all3), hjust = -0.7, size = 5) +
      scale_x_continuous(limits = c(0, max(10, max(d$N_LE_shared_all3) + 2)), expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Leading-edge genes shared across all three cohorts", y = NULL, title = paste0(dir, "-enriched pathways")) +
      theme_minimal(base_size = 16) +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 14, colour = "black"), plot.title = element_text(hjust = 0.5, size = 19))
  }
  fig3 <- (heat("Pediatric") + lolli("Pediatric", "#F2C200") +
             heat("Adult") + lolli("Adult", DIV_ADULT)) +
    plot_layout(ncol = 2, widths = c(1.25, 1), guides = "collect",
                heights = c(sum(themes$Direction == "Pediatric"), sum(themes$Direction == "Adult"))) +
    plot_annotation(tag_levels = "A",
                    caption = "Each row is the representative pathway of a non-redundant theme (pathways grouped by leading-edge overlap; number of grouped pathways in Supplementary Table S18).\n* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001. NES, normalized enrichment score.") &
    tag_theme & theme(legend.position = "bottom", plot.caption = element_text(size = 12.5, hjust = 0))
  save_fig(fig3, "Figure3_cross_cohort_GSEA", 20, 4 + 0.62 * nrow(themes))
} else message("  GSEA themes not found; run scripts/04_gsea/04_themes_partial_conjunction.R")


# FIGURE 4 (single cell) ----
message("Figure 4")

scp <- function(f) P("scRNAseq_COVID_nasal_age", "pseudobulk_signature", "tables", f)
comp <- function(f) P("scRNAseq_COVID_nasal_age", "cellular_composition", f)

if (file.exists(scp("PRIMARY_pseudobulk_edgeR_signature_genes.csv")) &&
    file.exists(scp("SECONDARY_mean_lognorm_signature_donor_celltype_long.csv"))) {

  pb <- rd(scp("PRIMARY_pseudobulk_edgeR_signature_genes.csv"))
  ml <- rd(scp("SECONDARY_mean_lognorm_signature_donor_celltype_long.csv"))
  ml$Group <- ifelse(ml$Group %in% c("Ped", "Pediatric"), "Pediatric", "Adult")
  cts <- intersect(EPITHELIAL, unique(pb$cell_type))
  ct_cols <- setNames(TAB10[seq_along(cts)], cts)

  # A. UMAP
  p4a <- NULL
  if (file.exists(comp("Yoshida_epithelial_UMAP_coordinates.csv.gz"))) {
    um <- read.csv(gzfile(comp("Yoshida_epithelial_UMAP_coordinates.csv.gz")), stringsAsFactors = FALSE)
    um <- um[um$cell_type %in% cts, ]
    um$cell_type <- factor(um$cell_type, levels = cts)
    set.seed(1)
    um <- um[sample(nrow(um)), ]
    p4a <- ggplot(um, aes(UMAP1, UMAP2, colour = cell_type)) +
      geom_point(size = 0.12, alpha = 0.7, stroke = 0) +
      scale_colour_manual(values = ct_cols, name = "Cell type") +
      guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
      labs(title = "Yoshida et al. Airway COVID-19") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(size = 13, hjust = 0.3), legend.title = element_text(size = 12),
            legend.text = element_text(size = 11))
  }

  # B. Composition (first-round style)
  p4b <- NULL
  if (file.exists(comp("Yoshida_epithelial_proportions_by_donor.csv"))) {
    pr <- rd(comp("Yoshida_epithelial_proportions_by_donor.csv"))
    gcol <- intersect(c("Group_plot", "Group"), colnames(pr))[1]
    pr$Group <- factor(ifelse(pr[[gcol]] %in% c("Ped", "Pediatric"), "Pediatric", "Adult"), levels = c("Pediatric", "Adult"))
    pr <- pr[pr$cell_type %in% cts, ]
    pr$cell_type <- factor(pr$cell_type, levels = cts)
    sm <- pr %>% group_by(cell_type, Group) %>%
      summarise(mean = mean(proportion, na.rm = TRUE), se = sd(proportion, na.rm = TRUE) / sqrt(sum(!is.na(proportion))),
                .groups = "drop")
    ab <- if (file.exists(scp("celltype_abundance_edgeR.csv"))) rd(scp("celltype_abundance_edgeR.csv")) else NULL
    star_df <- if (!is.null(ab)) ab %>% filter(cell_type %in% cts, FDR < ALPHA) %>%
      transmute(cell_type = factor(cell_type, levels = cts), star = stars(FDR)) %>%
      left_join(pr %>% group_by(cell_type) %>% summarise(ymax = max(proportion, na.rm = TRUE), .groups = "drop"), by = "cell_type") else NULL
    src(sm, "Figure4_B_composition")
    p4b <- ggplot(sm, aes(cell_type, mean, fill = Group)) +
      geom_col(position = position_dodge(0.8), width = 0.75) +
      geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se), position = position_dodge(0.8), width = 0.25, linewidth = 0.4) +
      geom_point(data = pr, aes(cell_type, proportion, colour = Group),
                 position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1), size = 1.4,
                 alpha = 0.75, inherit.aes = FALSE, show.legend = FALSE) +
      { if (!is.null(star_df) && nrow(star_df) > 0) geom_text(data = star_df, aes(cell_type, ymax + 0.03, label = star), inherit.aes = FALSE, size = 5) } +
      scale_fill_manual(values = c(Pediatric = SC_PED, Adult = SC_ADULT), name = NULL) +
      scale_colour_manual(values = c(Pediatric = SC_PED, Adult = SC_ADULT)) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
      labs(x = NULL, y = "Cell proportion") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5), legend.position = "right")
  }

  # C. Heatmap of relative expression (first-round layout)
  ped_genes <- signature %>% filter(RE_log2FC > 0) %>% arrange(desc(RE_log2FC)) %>% pull(SYMBOL)
  ad_genes <- signature %>% filter(RE_log2FC < 0) %>% arrange(RE_log2FC) %>% pull(SYMBOL)
  hm <- ml %>% filter(cell_type %in% cts) %>%
    group_by(SYMBOL, cell_type, Group) %>% summarise(mean_expr = mean(mean_lognorm, na.rm = TRUE), .groups = "drop") %>%
    group_by(SYMBOL) %>% mutate(z = (mean_expr - mean(mean_expr)) / sd(mean_expr)) %>% ungroup() %>%
    mutate(z = pmax(pmin(ifelse(is.finite(z), z, 0), 3), -3))
  st <- pb %>% filter(cell_type %in% cts, FDR_signature_global < ALPHA) %>%
    transmute(SYMBOL, cell_type, Group = ifelse(logFC > 0, "Pediatric", "Adult"), star = stars(FDR_signature_global))
  hm <- hm %>% left_join(st, by = c("SYMBOL", "cell_type", "Group")) %>%
    mutate(star = ifelse(is.na(star), "", star),
           block = factor(ifelse(SYMBOL %in% ped_genes, "Pediatric-enriched", "Adult-enriched"),
                          levels = c("Pediatric-enriched", "Adult-enriched")),
           SYMBOL = factor(SYMBOL, levels = rev(c(ped_genes, ad_genes))),
           cell_type = factor(cell_type, levels = cts),
           Group = factor(Group, levels = c("Pediatric", "Adult")))
  src(hm, "Figure4_C_relative_expression")

  bar <- expand.grid(cell_type = factor(cts, levels = cts), Group = factor(c("Pediatric", "Adult"), levels = c("Pediatric", "Adult")))
  p4c_bar <- ggplot(bar, aes(Group, 1, fill = Group)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(. ~ cell_type) +
    scale_fill_manual(values = c(Pediatric = SC_PED, Adult = SC_ADULT), guide = "none") +
    theme_void() +
    theme(strip.text = element_text(size = 8, margin = margin(b = 2, t = 2)),
          strip.background = element_rect(fill = "grey95", colour = "grey60"), panel.spacing.x = unit(0.12, "cm"))

  p4c <- ggplot(hm, aes(Group, SYMBOL, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = star), size = 4, vjust = 0.75) +
    facet_grid(block ~ cell_type, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = HEAT_COLS, limits = c(-3, 3), breaks = -3:3, name = "Relative expression") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5),
          axis.text.y = element_text(face = "italic", size = 9.5, colour = "black"),
          strip.text.y.left = element_text(angle = 90, size = 11), strip.placement = "outside",
          strip.background.y = element_rect(fill = "#EAF0F7", colour = NA),
          panel.spacing.x = unit(0.12, "cm"), panel.spacing.y = unit(0.25, "cm"),
          legend.position = "bottom", legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.35, "cm"))

  p4c_full <- p4c_bar / p4c + plot_layout(heights = c(0.035, 1))
  fig4 <- if (!is.null(p4a) && !is.null(p4b)) {
    ((p4a | p4b) + plot_layout(widths = c(1, 1.15))) / wrap_elements(full = p4c_full) + plot_layout(heights = c(1, 2.6))
  } else wrap_elements(full = p4c_full)
  fig4 <- fig4 + plot_annotation(tag_levels = "A",
                                 caption = "C: mean donor-level log-normalized expression scaled per gene. * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001 (raw-count pseudobulk, edgeR; BH across all signature-gene x cell-type tests),\nshown in the age group with higher expression.") &
    theme(plot.tag = element_text(face = "bold", size = 20), plot.caption = element_text(size = 9.5, hjust = 0))
  save_fig(fig4, "Figure4_single_cell", 11, 15)
} else message("  single-cell results not found; run scripts/05_single_cell (01-03)")


# FIGURES FOR THE SUPPLEMENT ----
message("Supplementary composite figures")

# Signature retention across inference methods and clinical models
retention <- list()
inf_file <- P("meta_analysis", "inference_methods", "signature_genes_by_inference_method.csv")
if (file.exists(inf_file)) {
  inf <- rd(inf_file)
  same <- function(col) sum(sign(inf[[col]]) == sign(inf$RE_log2FC), na.rm = TRUE)
  retention$a <- data.frame(
    Block = "Inference",
    Criterion = c("REML + Wald z (primary; FDR < 0.05)", "DerSimonian-Laird + z (FDR < 0.05)", "Paule-Mandel + z (FDR < 0.05)",
                  "REML + Hartung-Knapp (p < 0.05)", "Bayesian, half-normal(0.5) prior (95% CrI)",
                  "Bayesian, half-Cauchy(0.5) prior (95% CrI)", "Bayesian, half-normal(1.0) prior (95% CrI)"),
    Same = c(same("REML_z_est"), same("DL_z_est"), same("PM_z_est"), same("REML_HKSJ_est"),
             same("Bayes_HN0.5_post_median"), same("Bayes_HC0.5_post_median"), same("Bayes_HN1.0_post_median")),
    Supported = c(sum(inf$REML_z_FDR < ALPHA), sum(inf$DL_z_FDR < ALPHA), sum(inf$PM_z_FDR < ALPHA), sum(inf$REML_HKSJ_pval < 0.05),
                  sum(inf$Bayes_HN0.5_cri_excludes_0), sum(inf$Bayes_HC0.5_cri_excludes_0), sum(inf$Bayes_HN1.0_cri_excludes_0)))
}
sc_file <- P("sensitivity", "clinical_covariates", "scenario_summary.csv")
if (file.exists(sc_file)) {
  sc <- rd(sc_file)
  lab <- c(S1_viral_load_severity = "Adjusted for viral load / severity",
           S2_extended_clinical = "Further adjusted for symptoms; severe adults excluded",
           S3_exclude_age_16_20 = "Participants aged 16-20 years excluded",
           S4_log2_age_continuous = "Age modelled continuously (log2 age)")
  sc <- sc[sc$Scenario %in% names(lab), ]
  retention$b <- data.frame(Block = "Clinical / age", Criterion = lab[sc$Scenario],
                            Same = sc$n_same_direction, Supported = sc$n_same_direction_FDR05)
}
if (length(retention) > 0) {
  ret <- bind_rows(retention) %>%
    mutate(Criterion = factor(Criterion, levels = rev(Criterion)), Block = factor(Block, levels = c("Inference", "Clinical / age")))
  src(ret, "SuppFig_signature_retention")
  n_sig <- nrow(signature)
  p_ret <- ggplot(ret, aes(y = Criterion)) +
    geom_col(aes(x = Same, fill = "Same direction as primary"), width = 0.72) +
    geom_col(aes(x = Supported, fill = "Statistically supported"), width = 0.72) +
    geom_text(aes(x = Supported, label = Supported), hjust = 1.25, colour = "white", size = 4.5, fontface = "bold") +
    geom_vline(xintercept = n_sig, linetype = "22", colour = "#686868") +
    facet_grid(Block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c("Same direction as primary" = "#D9D9D9", "Statistically supported" = DIV_ADULT), name = NULL) +
    scale_x_continuous(limits = c(0, n_sig), breaks = seq(0, n_sig, by = 10), expand = expansion(mult = c(0, 0.01))) +
    labs(x = paste0("Signature genes (of ", n_sig, ")"), y = NULL) +
    theme_classic(base_size = 13) +
    theme(strip.placement = "outside", strip.background = element_rect(fill = "#ECECEC", colour = NA),
          strip.text.y.left = element_text(angle = 90, face = "bold"), legend.position = "bottom")
  save_fig(p_ret, "SuppFig_signature_retention", 11, 6, sup_dir)
}

# Age effect with and without SARS-CoV-2
inf4_file <- P("sensitivity", "age_by_infection", "signature_genes_age_by_infection.csv")
if (file.exists(inf4_file)) {
  inf4 <- rd(inf4_file)
  pshort <- c("Age difference also present without SARS-CoV-2 (p < 0.05)" = "Present without SARS-CoV-2",
              "Same direction without SARS-CoV-2, not significant" = "Same direction, not significant",
              "Not observed without SARS-CoV-2 (opposite direction)" = "Opposite direction",
              "Age effect differs by infection status (interaction FDR < 0.05)" = "Infection-dependent")
  pcol <- c("Present without SARS-CoV-2" = DIV_ADULT, "Same direction, not significant" = "#9EBCDA",
            "Opposite direction" = "#BDBDBD", "Infection-dependent" = DISCORDANT_COLOR)
  inf4$Pattern_short <- factor(pshort[inf4$Pattern], levels = names(pcol))
  r5 <- cor(inf4$age_neg_log2FC, inf4$age_pos_log2FC)
  src(inf4, "SuppFig_age_by_infection")
  lab_genes <- inf4 %>% filter(SYMBOL %in% c("NTS", "SCRN1") | interaction_FDR_within_signature < ALPHA |
                                 (age_neg_pvalue < 0.001 & age_neg_same_direction))
  pa <- ggplot(inf4, aes(age_pos_log2FC, age_neg_log2FC)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = DIV_PED, alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = DIV_ADULT, alpha = 0.10) +
    geom_hline(yintercept = 0, colour = "#B8B8B8") + geom_vline(xintercept = 0, colour = "#B8B8B8") +
    geom_abline(slope = 1, intercept = 0, colour = "#686868", linetype = "22") +
    geom_errorbar(aes(ymin = age_neg_log2FC - 1.96 * age_neg_SE, ymax = age_neg_log2FC + 1.96 * age_neg_SE), colour = "grey75", width = 0) +
    geom_point(aes(fill = Pattern_short), shape = 21, size = 3.8, colour = "white", stroke = 0.5) +
    geom_text(data = lab_genes, aes(label = SYMBOL), size = 3.4, vjust = -0.9, fontface = "italic", check_overlap = TRUE) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 4.5, label = sprintf("r = %.2f", r5)) +
    scale_fill_manual(values = pcol, name = NULL, drop = FALSE) +
    labs(x = expression("Pediatric vs Adult " * log[2] * "FC, SARS-CoV-2-positive"),
         y = expression("Pediatric vs Adult " * log[2] * "FC (95% CI), SARS-CoV-2-negative")) +
    theme_classic(base_size = 13) + theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, override.aes = list(size = 4.5)))
  pb5 <- ggplot(inf4 %>% count(Pattern_short), aes(n, Pattern_short, fill = Pattern_short)) +
    geom_col(width = 0.7) + geom_text(aes(label = n), hjust = -0.4, size = 4.5) +
    scale_fill_manual(values = pcol, guide = "none") + scale_y_discrete(limits = rev(names(pcol))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Signature genes", y = NULL) + theme_classic(base_size = 13)
  pc5 <- NULL
  ex_file <- P("sensitivity", "age_by_infection", "signature_genes_expression_by_sample_long.csv")
  if (file.exists(ex_file)) {
    ex <- rd(ex_file)
    pick1 <- inf4 %>% filter(Pattern_short == "Present without SARS-CoV-2") %>% arrange(age_neg_pvalue) %>% head(2) %>% pull(SYMBOL)
    pick2 <- inf4 %>% filter(Pattern_short == "Infection-dependent") %>% arrange(interaction_pvalue) %>% head(2) %>% pull(SYMBOL)
    show <- unique(c("NTS", "SCRN1", pick1, pick2))
    ex <- ex %>% filter(SYMBOL %in% show) %>%
      mutate(SYMBOL = factor(SYMBOL, levels = show), Group = factor(Group, levels = c("Pediatric", "Adult")),
             Infection = factor(Infection, levels = c("Negative", "Positive"), labels = c("SARS-CoV-2\nnegative", "SARS-CoV-2\npositive")))
    pc5 <- ggplot(ex, aes(Infection, log2_norm_count, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7, position = position_dodge(0.8)) +
      geom_point(aes(colour = Group), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1), size = 0.7, alpha = 0.6) +
      facet_wrap(~ SYMBOL, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = c(Pediatric = DIV_PED, Adult = DIV_ADULT), name = NULL) +
      scale_colour_manual(values = c(Pediatric = "#9C7A00", Adult = "#1D4E85"), guide = "none") +
      labs(x = NULL, y = expression(log[2] * "(normalized count + 1)")) +
      theme_classic(base_size = 12) + theme(legend.position = "top", strip.text = element_text(face = "bold.italic"))
  }
  f5 <- (pa | pb5) + plot_layout(widths = c(1.2, 1))
  if (!is.null(pc5)) f5 <- f5 / pc5 + plot_layout(heights = c(1.3, 1))
  f5 <- f5 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 18))
  save_fig(f5, "SuppFig_age_by_infection", 14, 11, sup_dir)
}

# Donor-level NTS and SCRN1 in secretory cells
if (file.exists(scp("PRIMARY_pseudobulk_log2CPM_signature_donor_celltype_long.csv")) &&
    file.exists(scp("PRIMARY_pseudobulk_edgeR_signature_genes.csv"))) {
  pb <- rd(scp("PRIMARY_pseudobulk_edgeR_signature_genes.csv"))
  dl <- rd(scp("PRIMARY_pseudobulk_log2CPM_signature_donor_celltype_long.csv")) %>%
    filter(SYMBOL %in% c("NTS", "SCRN1"), cell_type == "Secretory") %>%
    mutate(Group = factor(ifelse(Group %in% c("Ped", "Pediatric"), "Pediatric", "Adult"), levels = c("Pediatric", "Adult")))
  ann <- pb %>% filter(SYMBOL %in% c("NTS", "SCRN1"), cell_type == "Secretory") %>%
    transmute(SYMBOL, lab = sprintf("log2FC = %.2f\nFDR = %.3f\nn = %d / %d donors", logFC, FDR_signature_global, n_Ped, n_Adult))
  src(dl, "SuppFig_sc_donor_level_NTS_SCRN1")
  pdl <- ggplot(dl, aes(Group, log2CPM, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.6) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), size = 2.2, shape = 21, colour = "black") +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = lab), vjust = 1.15, size = 4, inherit.aes = FALSE) +
    facet_wrap(~ SYMBOL, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.38))) +
    scale_fill_manual(values = c(Pediatric = SC_PED, Adult = SC_ADULT), guide = "none") +
    labs(x = NULL, y = "Secretory-cell pseudobulk log2 CPM\n(one point per donor)") +
    theme_classic(base_size = 13) +
    theme(strip.text = element_text(face = "bold.italic", size = 14), strip.background = element_rect(fill = "#ECECEC", colour = NA))
  save_fig(pdl, "SuppFig_sc_donor_level_NTS_SCRN1", 7, 5.5, sup_dir)
}

writeLines(capture.output(sessionInfo()), file.path(fig_dir, "sessionInfo_main_figures.txt"))
message("Done: ", fig_dir)
