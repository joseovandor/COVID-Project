# Non-redundant GSEA themes and partial-conjunction test
#
# A. Pathways enriched in the same direction with FDR < 0.05 in >= 2 cohorts
#    (from 02_cross_cohort_concordance.R) are grouped by the overlap
#    coefficient |A n B| / min(|A|, |B|) of their leading-edge genes, pooled
#    across the cohorts where they were significant. Average-linkage
#    clustering, separately for pediatric and adult pathways, cut at 0.5.
#    Theme counts with other overlap and Jaccard cut-offs are also reported.
#    Core genes are present in >= 50% of the pathways of a theme.
# B. Partial conjunction (Benjamini & Heller, 2008) for enrichment in at least
#    u cohorts: p_PC(u/n) = (n - u + 1) * p_(u) for each direction, times 2.
#    A cohort enriched in the opposite direction contributes p = 1.
#    BH-FDR across all pathways tested in >= 2 cohorts.
#
# Output: results/gsea/themes/

SCRIPT_BUILD <- "GSEA_THEMES_2026-10-26"
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
SIM_CUT <- 0.5
CORE_FRACTION <- 0.5

out_dir <- here::here("results", "gsea", "themes")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gsea <- read.csv(here::here("results", "gsea", "GSEA_all_cohorts_long.csv"), stringsAsFactors = FALSE)

rep_files <- here::here("results", "gsea", "cross_cohort", "tables",
                        c("GSEA_reproducible_significant_3of3.csv", "GSEA_reproducible_supported_2of3.csv"))
reproducible <- bind_rows(lapply(rep_files, read.csv, stringsAsFactors = FALSE)) %>%
  distinct(Collection, pathway, .keep_all = TRUE) %>%
  mutate(
    Direction = ifelse(N_FDR05_Pediatric >= N_FDR05_Adult, "Pediatric", "Adult"),
    key = paste(Collection, pathway, sep = "::")
  )

primary <- read.csv(here::here("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"),
                    stringsAsFactors = FALSE)
signature_genes <- primary$SYMBOL[primary$RE_FDR < ALPHA]

cat("Reproducible pathways:", nrow(reproducible), "\n")
print(table(reproducible$Direction, reproducible$N_FDR05))


# Part A. Leading-edge-based themes ----

gsea$key <- paste(gsea$Collection, gsea$pathway, sep = "::")

le_sets <- lapply(seq_len(nrow(reproducible)), function(i) {
  rows <- gsea %>%
    filter(key == reproducible$key[i], padj < ALPHA, EnrichmentDirection == reproducible$Direction[i])
  unique(unlist(strsplit(rows$leadingEdge, ";", fixed = TRUE)))
})
names(le_sets) <- reproducible$key
reproducible$N_LE_pooled <- vapply(le_sets, length, integer(1))

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
overlap_coef <- function(a, b) length(intersect(a, b)) / min(length(a), length(b))

pathway_themes <- list()
theme_table <- list()
jac_long <- list()

for (dir_i in c("Pediatric", "Adult")) {
  keys <- reproducible$key[reproducible$Direction == dir_i]
  if (length(keys) == 0) next
  J <- matrix(1, length(keys), length(keys), dimnames = list(keys, keys))
  O <- J
  if (length(keys) > 1) {
    for (a in seq_along(keys)) for (b in seq_along(keys)) {
      J[a, b] <- jaccard(le_sets[[keys[a]]], le_sets[[keys[b]]])
      O[a, b] <- overlap_coef(le_sets[[keys[a]]], le_sets[[keys[b]]])
    }
    hc <- hclust(as.dist(1 - O), method = "average")
    cl <- cutree(hc, h = 1 - SIM_CUT)
    ord <- hc$labels[hc$order]
    hc_j <- hclust(as.dist(1 - J), method = "average")
    n_clusters_sensitivity <- c(
      vapply(c(0.3, 0.5, 0.7), function(t) length(unique(cutree(hc, h = 1 - t))), numeric(1)),
      vapply(c(0.10, 0.25, 0.50), function(t) length(unique(cutree(hc_j, h = 1 - t))), numeric(1))
    )
  } else {
    cl <- setNames(1, keys); ord <- keys; n_clusters_sensitivity <- rep(1, 6)
  }
  message(dir_i, ": ", length(keys), " pathways -> ", length(unique(cl)), " themes (overlap >= ", SIM_CUT,
          "); overlap 0.3/0.5/0.7: ", paste(n_clusters_sensitivity[1:3], collapse = "/"),
          "; Jaccard 0.10/0.25/0.50: ", paste(n_clusters_sensitivity[4:6], collapse = "/"))

  # renumber clusters in dendrogram order
  cl_ord <- unique(cl[ord])
  cl <- setNames(match(cl, cl_ord), names(cl))
  theme_id <- setNames(paste0(substr(dir_i, 1, 3), "_T", sprintf("%02d", cl)), names(cl))

  pathway_themes[[dir_i]] <- reproducible %>%
    filter(key %in% keys) %>%
    mutate(Theme = theme_id[key], dendrogram_order = match(key, ord))

  jac_long[[dir_i]] <- as.data.frame(as.table(O)) %>%
    rename(key_a = Var1, key_b = Var2, Similarity = Freq) %>%
    mutate(Direction = dir_i,
           key_a = factor(key_a, levels = ord), key_b = factor(key_b, levels = ord),
           theme_a = theme_id[as.character(key_a)], theme_b = theme_id[as.character(key_b)])

  theme_table[[dir_i]] <- data.frame(
    Direction = dir_i,
    n_pathways = length(keys),
    themes_overlap_0.3 = n_clusters_sensitivity[1],
    themes_overlap_0.5_primary = n_clusters_sensitivity[2],
    themes_overlap_0.7 = n_clusters_sensitivity[3],
    themes_jaccard_0.10 = n_clusters_sensitivity[4],
    themes_jaccard_0.25 = n_clusters_sensitivity[5],
    themes_jaccard_0.50 = n_clusters_sensitivity[6]
  )
}

pathway_themes <- bind_rows(pathway_themes)

themes <- pathway_themes %>%
  group_by(Theme, Direction) %>%
  arrange(desc(N_FDR05), desc(abs(Mean_NES)), .by_group = TRUE) %>%
  summarise(
    n_pathways = n(),
    representative_pathway = first(pathway),
    representative_collection = first(Collection),
    max_cohorts_FDR05 = max(N_FDR05),
    mean_abs_NES = mean(abs(Mean_NES)),
    pathways = paste(pathway, collapse = "; "),
    keys = list(key),
    .groups = "drop"
  )

themes$core_leading_edge <- vapply(themes$keys, function(k) {
  g <- unlist(le_sets[k])
  freq <- sort(table(g), decreasing = TRUE)
  core <- names(freq)[freq >= max(1, ceiling(CORE_FRACTION * length(k)))]
  paste(head(core, 25), collapse = ";")
}, character(1))
themes$n_core_genes <- vapply(strsplit(themes$core_leading_edge, ";"), function(x) sum(nzchar(x)), integer(1))
themes$signature_genes_in_leading_edge <- vapply(themes$keys, function(k) {
  paste(sort(intersect(unique(unlist(le_sets[k])), signature_genes)), collapse = ";")
}, character(1))
themes$keys <- NULL
themes <- themes %>% arrange(Direction, desc(n_pathways))

write.csv(themes, file.path(out_dir, "gsea_themes.csv"), row.names = FALSE)
write.csv(pathway_themes %>% select(Theme, Direction, Collection, pathway, N_cohorts_tested, N_FDR05, Mean_NES, Min_NES, Max_NES, N_LE_pooled) %>%
            arrange(Direction, Theme, desc(N_FDR05)),
          file.path(out_dir, "gsea_pathway_to_theme.csv"), row.names = FALSE)
write.csv(bind_rows(theme_table), file.path(out_dir, "gsea_theme_number_by_similarity_threshold.csv"), row.names = FALSE)

print(themes %>% select(Theme, Direction, n_pathways, representative_pathway, n_core_genes) %>% as.data.frame())


# Part B. Partial-conjunction replicability test ----

pc <- gsea %>%
  filter(!is.na(pval), !is.na(NES)) %>%
  group_by(Collection, pathway) %>%
  filter(n() >= 2) %>%
  summarise(
    n_cohorts = n(),
    p_ped = list(sort(ifelse(NES > 0, pval, 1))),
    p_adult = list(sort(ifelse(NES < 0, pval, 1))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    pc2_ped = min(1, (n_cohorts - 2 + 1) * p_ped[2]),
    pc2_adult = min(1, (n_cohorts - 2 + 1) * p_adult[2]),
    p_PC_2 = min(1, 2 * min(pc2_ped, pc2_adult)),
    dir_PC_2 = ifelse(pc2_ped <= pc2_adult, "Pediatric", "Adult"),
    p_PC_all = ifelse(n_cohorts == 3, min(1, 2 * min(p_ped[3], p_adult[3])), NA_real_),
    dir_PC_all = ifelse(n_cohorts == 3, ifelse(p_ped[3] <= p_adult[3], "Pediatric", "Adult"), NA_character_)
  ) %>%
  ungroup() %>%
  select(-p_ped, -p_adult) %>%
  mutate(
    FDR_PC_2 = p.adjust(p_PC_2, method = "BH"),
    FDR_PC_all = p.adjust(p_PC_all, method = "BH")
  ) %>%
  mutate(key = paste(Collection, pathway, sep = "::")) %>%
  left_join(reproducible %>% select(key, Direction_descriptive = Direction, N_FDR05), by = "key") %>%
  mutate(descriptive_reproducible = !is.na(Direction_descriptive)) %>%
  arrange(FDR_PC_2)

write.csv(pc, file.path(out_dir, "gsea_partial_conjunction_all_pathways.csv"), row.names = FALSE)

pc_summary <- data.frame(
  pathways_tested_in_2plus_cohorts = nrow(pc),
  descriptive_rule_2of3 = sum(pc$descriptive_reproducible),
  PC_2_FDR05 = sum(pc$FDR_PC_2 < ALPHA),
  PC_2_FDR05_and_descriptive = sum(pc$FDR_PC_2 < ALPHA & pc$descriptive_reproducible),
  PC_all3_FDR05 = sum(pc$FDR_PC_all < ALPHA, na.rm = TRUE)
)
write.csv(pc_summary, file.path(out_dir, "gsea_partial_conjunction_summary.csv"), row.names = FALSE)
print(pc_summary)

themes_pc <- pathway_themes %>%
  left_join(pc %>% select(key, FDR_PC_2, FDR_PC_all), by = "key") %>%
  group_by(Theme) %>%
  summarise(n_pathways_PC_2_FDR05 = sum(FDR_PC_2 < ALPHA, na.rm = TRUE),
            min_FDR_PC_2 = min(FDR_PC_2, na.rm = TRUE), .groups = "drop")
themes <- themes %>% left_join(themes_pc, by = "Theme")
write.csv(themes, file.path(out_dir, "gsea_themes.csv"), row.names = FALSE)


# Figures ----

short_name <- function(x) {
  x <- sub("^(GOBP|REACTOME|KEGG|HALLMARK)_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  ifelse(nchar(x) > 55, paste0(substr(x, 1, 52), "..."), x)
}

jac_df <- bind_rows(jac_long) %>%
  mutate(label_a = factor(short_name(sub(".*::", "", as.character(key_a))),
                          levels = unique(short_name(sub(".*::", "", levels(key_a))))),
         label_b = factor(short_name(sub(".*::", "", as.character(key_b))),
                          levels = unique(short_name(sub(".*::", "", levels(key_b))))))

for (dir_i in unique(jac_df$Direction)) {
  d <- jac_df %>% filter(Direction == dir_i)
  n_p <- length(unique(d$key_a))
  p_j <- ggplot(d, aes(x = label_a, y = label_b, fill = Similarity)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = ifelse(dir_i == "Pediatric", "#B08800", "#1D4E85"), limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "Overlap\ncoefficient", title = paste0(dir_i, "-enriched reproducible pathways: leading-edge overlap")) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank())
  size_in <- max(5, 0.16 * n_p + 3)
  ggsave(file.path(fig_dir, paste0("FigS_gsea_leading_edge_overlap_", dir_i, ".pdf")), p_j, width = size_in + 1.5, height = size_in)
  ggsave(file.path(fig_dir, paste0("FigS_gsea_leading_edge_overlap_", dir_i, ".tiff")), p_j, width = size_in + 1.5, height = size_in,
         dpi = 400, compression = "lzw")
}

rep_rows <- themes %>% select(Theme, Direction, n_pathways, representative_pathway, representative_collection)
dot_df <- gsea %>%
  inner_join(rep_rows, by = c("pathway" = "representative_pathway", "Collection" = "representative_collection")) %>%
  mutate(
    theme_label = paste0(short_name(pathway), "  [", n_pathways, " pathway", ifelse(n_pathways > 1, "s", ""), "]"),
    neglog10padj = pmin(-log10(padj), 10)
  )

dot_df$theme_label <- factor(dot_df$theme_label,
                             levels = unique(dot_df$theme_label[order(dot_df$Direction, dot_df$n_pathways)]))

p_dot <- ggplot(dot_df, aes(x = Dataset, y = theme_label)) +
  geom_point(aes(size = neglog10padj, colour = NES)) +
  geom_point(data = dot_df %>% filter(padj < ALPHA), aes(size = neglog10padj), shape = 21, colour = "black", stroke = 0.5) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient2(low = "#2C6FB7", mid = "grey90", high = "#E8B500", midpoint = 0) +
  scale_size_continuous(range = c(1, 6), name = "-log10 FDR") +
  labs(x = NULL, y = NULL, caption = "Representative pathway per theme; black outline: cohort FDR < 0.05") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_gsea_themes_representative.pdf"), p_dot, width = 8.5, height = max(4, 0.32 * nrow(themes) + 2))
ggsave(file.path(fig_dir, "FigS_gsea_themes_representative.tiff"), p_dot, width = 8.5, height = max(4, 0.32 * nrow(themes) + 2),
       dpi = 600, compression = "lzw")

writeLines(c(
  paste("Script build:", SCRIPT_BUILD),
  paste("Reproducible pathways (descriptive rule):", nrow(reproducible)),
  paste0("  ", names(table(reproducible$Direction)), ": ", table(reproducible$Direction)),
  paste("Themes (overlap coefficient >=", SIM_CUT, "):", nrow(themes)),
  paste0("  ", themes$Theme, " (", themes$n_pathways, " pathways): ", themes$representative_pathway),
  "",
  "Partial conjunction (>= 2 cohorts, BH-FDR):",
  capture.output(print(pc_summary))
), file.path(out_dir, "run_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done: ", out_dir)
