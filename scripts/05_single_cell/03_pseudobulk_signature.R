# Donor-level analysis of the signature genes in the single-cell data
#
#   Primary   : raw UMI counts summed per donor x cell type (pseudobulk),
#               edgeR quasi-likelihood, ~ Sex + Group (~ Group if sex is not
#               estimable)
#   Secondary : mean log-normalized expression per donor x cell type, two-sided
#               Wilcoxon test and Cliff's delta
# Every signature gene gets a status (absent from the matrix / not expressed /
# filtered / tested). Epithelial cell-type abundance is tested with edgeR QL
# (offset: log total epithelial cells per donor).
# Positive logFC or Cliff's delta = higher in pediatric donors.
#
# Run after 01_import_yoshida.R in the same R session; it uses sce,
# sce_age_filt, file_h5ad, genes_primary_sc, signature_table, min_cells,
# min_donors and results_dir.
# Output: results/scRNAseq_COVID_nasal_age/pseudobulk_signature/

library(SingleCellExperiment)
library(SummarizedExperiment)
library(DelayedArray)
library(DelayedMatrixStats)
library(scuttle)
library(edgeR)
library(Matrix)
library(rhdf5)
library(dplyr)
library(tidyr)
library(ggplot2)

required_objects <- c(
  "sce",
  "sce_age",
  "sce_age_filt",
  "file_h5ad",
  "genes_primary_sc",
  "signature_table",
  "min_cells",
  "min_donors",
  "results_dir"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "Run 01_import_yoshida.R first in this session. Missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

out_dir <- file.path(results_dir, "pseudobulk_signature")
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)


# 1. CELL / DONOR ANNOTATION ----

meta_cells <- as.data.frame(colData(sce_age_filt))
meta_cells$donor_id <- as.character(meta_cells$donor_id)
meta_cells$cell_type <- as.character(meta_cells$Cell_type_annotation_level2)
meta_cells$Group <- as.character(meta_cells$Group)

sex_column <- intersect(c("sex", "Sex", "gender", "Gender"), colnames(meta_cells))[1]

if (is.na(sex_column)) {
  meta_cells$Sex <- NA_character_
  message("No sex column found in the single-cell metadata.")
} else {
  meta_cells$Sex <- tolower(as.character(meta_cells[[sex_column]]))
  meta_cells$Sex <- ifelse(
    meta_cells$Sex %in% c("male", "m"), "Male",
    ifelse(meta_cells$Sex %in% c("female", "f"), "Female", NA)
  )
}

meta_cells$aggregate_id <- paste(meta_cells$donor_id, meta_cells$cell_type, sep = "___")

donor_info <- meta_cells %>%
  dplyr::distinct(donor_id, Group, Sex)

if (anyDuplicated(donor_info$donor_id) > 0) {
  stop("A donor has more than one Group or Sex value. Inspect the metadata.")
}

write.csv(
  donor_info,
  file.path(out_dir, "tables", "donors_group_sex.csv"),
  row.names = FALSE
)

aggregate_info <- meta_cells %>%
  dplyr::count(aggregate_id, donor_id, cell_type, Group, Sex, name = "n_cells")

# Every aggregate already has >= min_cells cells (filter applied in script 00).
stopifnot(all(aggregate_info$n_cells >= min_cells))

donors_per_celltype <- aggregate_info %>%
  dplyr::count(cell_type, Group, name = "n_donors") %>%
  tidyr::pivot_wider(names_from = Group, values_from = n_donors, values_fill = 0)

if (!"Adult" %in% colnames(donors_per_celltype)) donors_per_celltype$Adult <- 0
if (!"Ped" %in% colnames(donors_per_celltype)) donors_per_celltype$Ped <- 0

donors_per_celltype <- donors_per_celltype %>%
  dplyr::mutate(
    testable = Adult >= min_donors & Ped >= min_donors
  )

write.csv(
  donors_per_celltype,
  file.path(out_dir, "tables", "donors_per_celltype_after_cell_filter.csv"),
  row.names = FALSE
)

testable_celltypes <- donors_per_celltype$cell_type[donors_per_celltype$testable]

cat("\nTestable cell types:", length(testable_celltypes), "\n")
print(donors_per_celltype, n = Inf)


# 2-3. RAW COUNTS AND PSEUDOBULK (SUM PER DONOR x CELL TYPE)
#
# CELLxGENE H5AD files keep the raw UMI counts in raw.X when X
# holds normalized values (as in this dataset: X is log-normalized,
# raw.X contains integer counts with the same dimensions).
# zellkonverter's R reader does not import raw.X, so raw.X is read
# directly from the HDF5 file (CSR: cells x genes) in blocks of
# cells and summed per donor x cell type on the fly, without
# loading the full matrix into memory.

raw_source <- NA_character_
pb_counts <- NULL

if ("counts" %in% assayNames(sce_age_filt)) {

  raw_source <- "assay 'counts'"
  pb_summary <- scuttle::summarizeAssayByGroup(assay(sce_age_filt, "counts"),
                                               ids = meta_cells$aggregate_id, statistics = "sum")
  pb_counts <- as.matrix(assay(pb_summary, "sum"))
  rownames(pb_counts) <- rownames(sce_age_filt)

} else {

  h5_has <- function() {
    ls <- rhdf5::h5ls(file_h5ad, recursive = TRUE)
    paste0(ifelse(ls$group == "/", "", ls$group), "/", ls$name)
  }
  h5_objects <- h5_has()
  if (!"/raw/X/data" %in% h5_objects) {
    stop("No raw counts found: the H5AD has no 'counts' assay and no raw/X. ",
         "Report only the secondary analysis (mean log-normalized expression).")
  }

  raw_source <- "raw.X (read directly from the H5AD with rhdf5)"
  message("Reading raw.X from ", file_h5ad, " in blocks...")

  # The name of the index column is stored in the "_index" attribute of each
  # dataframe group (e.g. obs/index in this file, raw/var/_index).
  h5_index <- function(group) {
    col <- rhdf5::h5readAttributes(file_h5ad, group)[["_index"]]
    if (is.null(col)) col <- "_index"
    as.character(rhdf5::h5read(file_h5ad, paste0(group, "/", col)))
  }
  obs_index <- h5_index("obs")
  raw_gene_index <- h5_index("raw/var")
  indptr <- as.numeric(rhdf5::h5read(file_h5ad, "raw/X/indptr"))
  n_genes <- length(raw_gene_index)

  cell_rows <- match(colnames(sce_age_filt), obs_index)
  if (anyNA(cell_rows)) stop(sum(is.na(cell_rows)), " cells of sce_age_filt not found in obs/_index")

  agg_levels <- unique(meta_cells$aggregate_id)
  cell_group <- match(meta_cells$aggregate_id, agg_levels)
  ord <- order(cell_rows)
  cell_rows <- cell_rows[ord]
  cell_group <- cell_group[ord]

  pb_counts <- matrix(0, n_genes, length(agg_levels), dimnames = list(raw_gene_index, agg_levels))

  block_rows <- 10000
  windows <- split(seq_along(cell_rows), floor((cell_rows - 1) / block_rows))
  probe_ok <- TRUE

  for (w in seq_along(windows)) {
    sel <- windows[[w]]
    r0 <- min(cell_rows[sel]); r1 <- max(cell_rows[sel])
    start <- indptr[r0] + 1                 # indptr is 0-based
    n_val <- indptr[r1 + 1] - indptr[r0]
    if (n_val <= 0) next
    vals <- as.numeric(rhdf5::h5read(file_h5ad, "raw/X/data", start = start, count = n_val))
    idx <- as.integer(rhdf5::h5read(file_h5ad, "raw/X/indices", start = start, count = n_val))
    if (w == 1) probe_ok <- all(vals[seq_len(min(5000, length(vals)))] ==
                                  round(vals[seq_len(min(5000, length(vals)))]))
    # transpose of the CSR block = CSC genes x cells
    block <- Matrix::sparseMatrix(
      i = idx + 1L,
      p = as.integer(indptr[r0:(r1 + 1)] - indptr[r0]),
      x = vals,
      dims = c(n_genes, r1 - r0 + 1)
    )
    block <- block[, cell_rows[sel] - r0 + 1, drop = FALSE]
    indicator <- Matrix::sparseMatrix(i = seq_along(sel), j = cell_group[sel], x = 1,
                                      dims = c(length(sel), length(agg_levels)))
    pb_counts <- pb_counts + as.matrix(block %*% indicator)
    if (w %% 5 == 0 || w == length(windows)) message("  block ", w, "/", length(windows))
  }

  if (!probe_ok) stop("raw.X does not contain integer counts.")
}

if (!all(pb_counts == round(pb_counts))) stop("Pseudobulk counts are not integer-valued.")
storage.mode(pb_counts) <- "integer"

cat("\nRaw counts source:", raw_source, "\n")
cat("Pseudobulk matrix:", nrow(pb_counts), "genes x", ncol(pb_counts), "donor x cell-type aggregates\n")

# Gene symbols for the raw-count genes.
raw_gene_ids <- sub("\\..*$", "", rownames(pb_counts))
sce_gene_ids <- sub("\\..*$", "", rownames(sce))
raw_symbols <- rowData(sce)$SYMBOL[match(raw_gene_ids, sce_gene_ids)]
if ("/raw/var/feature_name" %in% (if (exists("h5_objects")) h5_objects else character(0))) {
  raw_feature <- tryCatch(as.character(rhdf5::h5read(file_h5ad, "raw/var/feature_name")), error = function(e) NULL)
  # categorical columns are stored as codes + categories
  if (is.null(raw_feature) || length(raw_feature) != nrow(pb_counts)) {
    codes <- tryCatch(rhdf5::h5read(file_h5ad, "raw/var/feature_name/codes"), error = function(e) NULL)
    cats <- tryCatch(rhdf5::h5read(file_h5ad, "raw/var/feature_name/categories"), error = function(e) NULL)
    raw_feature <- if (!is.null(codes) && !is.null(cats)) as.character(cats[codes + 1]) else NULL
  }
  if (!is.null(raw_feature)) raw_symbols <- ifelse(is.na(raw_symbols) | raw_symbols == "", raw_feature, raw_symbols)
}
cat("Signature genes found among raw-count genes:", sum(genes_primary_sc %in% raw_symbols), "/",
    length(genes_primary_sc), "\n")

pb_meta <- aggregate_info[match(colnames(pb_counts), aggregate_info$aggregate_id), ]
stopifnot(identical(pb_meta$aggregate_id, colnames(pb_counts)))

saveRDS(
  list(counts = pb_counts, meta = pb_meta, symbols = raw_symbols, source = raw_source),
  file.path(out_dir, "tables", "pseudobulk_raw_counts_donor_celltype.rds")
)


# 4. PRIMARY ANALYSIS: edgeR QUASI-LIKELIHOOD PER CELL TYPE ----

signature_ids <- which(raw_symbols %in% genes_primary_sc)

edger_results_list <- list()
edger_design_list <- list()

for (ct in testable_celltypes) {

  cols <- which(pb_meta$cell_type == ct)
  meta_ct <- pb_meta[cols, ]
  meta_ct$Group <- factor(meta_ct$Group, levels = c("Adult", "Ped"))

  y <- edgeR::DGEList(
    counts = pb_counts[, cols, drop = FALSE],
    group = meta_ct$Group
  )

  sex_ok <- all(!is.na(meta_ct$Sex)) && length(unique(meta_ct$Sex)) == 2

  if (sex_ok) {
    design_ct <- model.matrix(~ Sex + Group, data = meta_ct)
    design_label <- "~ Sex + Group"
    if (qr(design_ct)$rank < ncol(design_ct)) {
      design_ct <- model.matrix(~ Group, data = meta_ct)
      design_label <- "~ Group"
    }
  } else {
    design_ct <- model.matrix(~ Group, data = meta_ct)
    design_label <- "~ Group"
  }

  keep_genes <- edgeR::filterByExpr(y, design = design_ct)
  y_f <- y[keep_genes, seq_len(ncol(y)), keep.lib.sizes = FALSE]
  y_f <- edgeR::calcNormFactors(y_f)
  y_f <- edgeR::estimateDisp(y_f, design_ct)
  fit <- edgeR::glmQLFit(y_f, design_ct, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "GroupPed")

  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tt$gene_id <- rownames(tt)
  tt$SYMBOL <- raw_symbols[match(tt$gene_id, rownames(pb_counts))]
  tt$FDR_genomewide_celltype <- tt$FDR

  sig_rows <- tt[tt$SYMBOL %in% genes_primary_sc, , drop = FALSE]

  if (nrow(sig_rows) > 0) {
    sig_rows$cell_type <- ct
    sig_rows$n_Adult <- sum(meta_ct$Group == "Adult")
    sig_rows$n_Ped <- sum(meta_ct$Group == "Ped")
    sig_rows$design <- design_label
    edger_results_list[[ct]] <- sig_rows
  }

  # Signature genes present in the matrix but removed by filterByExpr
  filtered_symbols <- setdiff(
    raw_symbols[signature_ids],
    sig_rows$SYMBOL
  )

  edger_design_list[[ct]] <- data.frame(
    cell_type = ct,
    design = design_label,
    n_Adult = sum(meta_ct$Group == "Adult"),
    n_Ped = sum(meta_ct$Group == "Ped"),
    genes_after_filterByExpr = nrow(y_f),
    signature_genes_tested = nrow(sig_rows),
    signature_genes_filtered = paste(sort(filtered_symbols), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

edger_results <- dplyr::bind_rows(edger_results_list)
edger_design <- dplyr::bind_rows(edger_design_list)

edger_results <- edger_results %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(FDR_signature_celltype = p.adjust(PValue, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_signature_global = p.adjust(PValue, method = "BH"),
    sc_direction = ifelse(logFC > 0, "Pediatric", "Adult")
  ) %>%
  dplyr::left_join(
    signature_table %>% dplyr::select(SYMBOL, k, RE_log2FC, RE_FDR, RE_direction),
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(direction_concordant_with_bulk = sc_direction == RE_direction) %>%
  dplyr::select(
    SYMBOL, k, cell_type, n_Adult, n_Ped, design,
    logFC, logCPM, F, PValue,
    FDR_genomewide_celltype, FDR_signature_celltype, FDR_signature_global,
    sc_direction, RE_log2FC, RE_FDR, RE_direction, direction_concordant_with_bulk
  ) %>%
  dplyr::arrange(FDR_signature_global, PValue)

write.csv(
  edger_results,
  file.path(out_dir, "tables", "PRIMARY_pseudobulk_edgeR_signature_genes.csv"),
  row.names = FALSE
)

write.csv(
  edger_design,
  file.path(out_dir, "tables", "PRIMARY_pseudobulk_edgeR_models_by_celltype.csv"),
  row.names = FALSE
)

cat("\nPRIMARY (raw-count pseudobulk, edgeR QL), top results:\n")
print(head(as.data.frame(edger_results), 20))


# 5. SECONDARY ANALYSIS: MEAN LOG-NORMALIZED EXPRESSION ----
#    per donor x cell type, Wilcoxon + Cliff's delta

x_ids <- which(rowData(sce_age_filt)$SYMBOL %in% genes_primary_sc)

mean_summary <- scuttle::summarizeAssayByGroup(
  assay(sce_age_filt, "X")[x_ids, , drop = FALSE],
  ids = meta_cells$aggregate_id,
  statistics = "mean"
)

mean_lognorm <- as.matrix(assay(mean_summary, "mean"))
rownames(mean_lognorm) <- rowData(sce_age_filt)$SYMBOL[x_ids]

mean_long <- as.data.frame(mean_lognorm) %>%
  tibble::rownames_to_column("SYMBOL") %>%
  tidyr::pivot_longer(-SYMBOL, names_to = "aggregate_id", values_to = "mean_lognorm") %>%
  dplyr::left_join(aggregate_info, by = "aggregate_id")

wilcox_results_list <- list()

for (ct in testable_celltypes) {
  for (gene in unique(mean_long$SYMBOL)) {

    tmp <- mean_long[mean_long$cell_type == ct & mean_long$SYMBOL == gene, ]
    ped <- tmp$mean_lognorm[tmp$Group == "Ped"]
    adult <- tmp$mean_lognorm[tmp$Group == "Adult"]

    if (length(ped) < min_donors || length(adult) < min_donors) next
    if (length(unique(c(ped, adult))) < 2) next  # constant values: not testable

    test <- wilcox.test(ped, adult, exact = FALSE)

    # Cliff's delta: P(Ped > Adult) - P(Ped < Adult); positive = Pediatric higher
    cliffs_delta <- mean(sign(outer(ped, adult, "-")))

    wilcox_results_list[[length(wilcox_results_list) + 1]] <- data.frame(
      SYMBOL = gene,
      cell_type = ct,
      n_Adult = length(adult),
      n_Ped = length(ped),
      median_Adult = median(adult),
      median_Ped = median(ped),
      cliffs_delta = cliffs_delta,
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  }
}

wilcox_results <- dplyr::bind_rows(wilcox_results_list) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(FDR_celltype = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    sc_direction = ifelse(cliffs_delta > 0, "Pediatric", ifelse(cliffs_delta < 0, "Adult", "None"))
  ) %>%
  dplyr::left_join(
    signature_table %>% dplyr::select(SYMBOL, k, RE_log2FC, RE_FDR, RE_direction),
    by = "SYMBOL"
  ) %>%
  dplyr::mutate(direction_concordant_with_bulk = sc_direction == RE_direction) %>%
  dplyr::arrange(FDR_global, p_value)

write.csv(
  wilcox_results,
  file.path(out_dir, "tables", "SECONDARY_mean_lognorm_wilcoxon_signature_genes.csv"),
  row.names = FALSE
)

cat("\nSECONDARY (mean log-normalized, Wilcoxon), top results:\n")
print(head(as.data.frame(wilcox_results), 20))


# 6. STATUS OF EVERY SIGNATURE GENE ----

in_x_matrix <- genes_primary_sc %in% rowData(sce)$SYMBOL
in_raw_matrix <- genes_primary_sc %in% raw_symbols

expressed_cells <- setNames(rep(0, length(genes_primary_sc)), genes_primary_sc)
if (length(x_ids) > 0) {
  n_expr <- DelayedMatrixStats::rowSums2(assay(sce_age_filt, "X")[x_ids, , drop = FALSE] > 0)
  n_expr <- tapply(n_expr, rowData(sce_age_filt)$SYMBOL[x_ids], sum)
  expressed_cells[names(n_expr)] <- n_expr
}

gene_status <- signature_table %>%
  dplyr::select(SYMBOL, k, RE_log2FC, RE_FDR, RE_direction) %>%
  dplyr::mutate(
    in_singlecell_matrix = SYMBOL %in% rowData(sce)$SYMBOL | SYMBOL %in% raw_symbols,
    n_cells_expressing_filtered = as.numeric(expressed_cells[SYMBOL]),
    n_celltypes_tested_primary = vapply(
      SYMBOL, function(g) sum(edger_results$SYMBOL == g), integer(1)
    ),
    n_celltypes_tested_secondary = vapply(
      SYMBOL, function(g) sum(wilcox_results$SYMBOL == g), integer(1)
    ),
    status = dplyr::case_when(
      !in_singlecell_matrix ~ "Absent from single-cell expression matrix",
      n_cells_expressing_filtered == 0 ~ "Present in matrix, not detected in any retained cell",
      n_celltypes_tested_primary == 0 & n_celltypes_tested_secondary == 0 ~
        "Present, excluded by analysis filters (filterByExpr / constant values)",
      TRUE ~ "Tested"
    )
  ) %>%
  dplyr::arrange(status, k, SYMBOL)

write.csv(
  gene_status,
  file.path(out_dir, "tables", "signature_gene_status_singlecell.csv"),
  row.names = FALSE
)

cat("\nSignature-gene status:\n")
print(table(gene_status$status))
print(gene_status %>% dplyr::filter(status != "Tested"), n = Inf)


# 6b. DONOR-LEVEL EXPRESSION TABLES FOR THE FIGURE SCRIPTS ----

pb_logcpm <- edgeR::cpm(edgeR::DGEList(pb_counts), log = TRUE, prior.count = 1)

sig_rows_pb <- which(raw_symbols %in% genes_primary_sc)
pb_sig_long <- as.data.frame(pb_logcpm[sig_rows_pb, , drop = FALSE])
pb_sig_long$SYMBOL <- raw_symbols[sig_rows_pb]
pb_sig_long <- pb_sig_long %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  tidyr::pivot_longer(-SYMBOL, names_to = "aggregate_id", values_to = "log2CPM") %>%
  dplyr::left_join(aggregate_info, by = "aggregate_id")

write.csv(
  pb_sig_long,
  file.path(out_dir, "tables", "PRIMARY_pseudobulk_log2CPM_signature_donor_celltype_long.csv"),
  row.names = FALSE
)

write.csv(
  mean_long,
  file.path(out_dir, "tables", "SECONDARY_mean_lognorm_signature_donor_celltype_long.csv"),
  row.names = FALSE
)


# 6c. CELL-TYPE ABUNDANCE (edgeR quasi-likelihood) ----
# Counts of cells per donor x epithelial cell type among all retained
# COVID+ nasal-cavity cells (sce_age, before the >= 20-cell filter).

abund_meta <- as.data.frame(colData(sce_age))
abund_meta$donor_id <- as.character(abund_meta$donor_id)
abund_meta$cell_type <- as.character(abund_meta$Cell_type_annotation_level2)

# Epithelial compartment only (same denominator as the composition panel of
# 03_UMAP.R): proportions are relative to each donor's epithelial cells.
epithelial_celltypes <- c("Basal", "Basal cycling", "Ciliated", "Deuterosomal", "Multiciliated", "Club",
                          "Duct", "Goblet", "Goblet inflammatory", "Hillock", "Hillock cycling",
                          "Rare epi", "Secretory", "Squamous", "Transit epi")
abund_meta <- abund_meta[abund_meta$cell_type %in% epithelial_celltypes, ]

abund_counts <- as.matrix(table(abund_meta$cell_type, abund_meta$donor_id))

donor_tab <- abund_meta %>%
  dplyr::distinct(donor_id, Group) %>%
  dplyr::mutate(Group = factor(as.character(Group), levels = c("Adult", "Ped")))
donor_tab <- donor_tab[match(colnames(abund_counts), donor_tab$donor_id), ]

y_ab <- edgeR::DGEList(counts = abund_counts, lib.size = colSums(abund_counts))
design_ab <- model.matrix(~ Group, data = donor_tab)
keep_ab <- edgeR::filterByExpr(y_ab, design = design_ab)
y_ab <- y_ab[keep_ab, seq_len(ncol(y_ab)), keep.lib.sizes = TRUE]
y_ab <- edgeR::estimateDisp(y_ab, design_ab, trend.method = "none")
fit_ab <- edgeR::glmQLFit(y_ab, design_ab, robust = TRUE, abundance.trend = FALSE)
res_ab <- edgeR::topTags(edgeR::glmQLFTest(fit_ab, coef = "GroupPed"), n = Inf)$table
res_ab$cell_type <- rownames(res_ab)

prop_tab <- as.data.frame(t(t(abund_counts) / colSums(abund_counts)))
colnames(prop_tab) <- c("cell_type", "donor_id", "proportion")
prop_tab <- prop_tab %>%
  dplyr::left_join(donor_tab %>% dplyr::mutate(donor_id = as.character(donor_id)), by = "donor_id")

abundance_results <- res_ab %>%
  dplyr::select(cell_type, logFC, logCPM, F, PValue, FDR) %>%
  dplyr::left_join(
    prop_tab %>%
      dplyr::group_by(cell_type) %>%
      dplyr::summarise(
        median_prop_Adult = median(proportion[Group == "Adult"]),
        median_prop_Ped = median(proportion[Group == "Ped"]),
        .groups = "drop"
      ),
    by = "cell_type"
  ) %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "Pediatric", "Adult")) %>%
  dplyr::arrange(FDR)

write.csv(abundance_results, file.path(out_dir, "tables", "celltype_abundance_edgeR.csv"), row.names = FALSE)
write.csv(prop_tab, file.path(out_dir, "tables", "celltype_proportions_by_donor.csv"), row.names = FALSE)

cat("\nCell-type abundance (edgeR QL), FDR < 0.05:", sum(abundance_results$FDR < 0.05), "\n")
print(head(abundance_results, 10))


# 7. DONOR-LEVEL PLOTS ----
# Genes plotted: NTS and SCRN1, plus any gene significant at
# FDR_signature_global < 0.05 (primary) or FDR_global < 0.05
# (secondary). One point = one donor.

plot_pairs <- dplyr::bind_rows(
  edger_results %>%
    dplyr::filter(FDR_signature_global < 0.05) %>%
    dplyr::select(SYMBOL, cell_type),
  wilcox_results %>%
    dplyr::filter(FDR_global < 0.05) %>%
    dplyr::select(SYMBOL, cell_type),
  data.frame(SYMBOL = c("NTS", "SCRN1"), cell_type = "Secretory")
) %>%
  dplyr::distinct() %>%
  dplyr::filter(cell_type %in% testable_celltypes)

# log2-CPM from the raw-count pseudobulk (pb_logcpm, section 6b)

donor_plot_data <- dplyr::bind_rows(
  lapply(seq_len(nrow(plot_pairs)), function(i) {
    g <- plot_pairs$SYMBOL[i]
    ct <- plot_pairs$cell_type[i]
    row_id <- which(raw_symbols == g)[1]
    cols <- which(pb_meta$cell_type == ct)
    if (is.na(row_id) || length(cols) == 0) return(NULL)
    data.frame(
      SYMBOL = g,
      cell_type = ct,
      donor_id = pb_meta$donor_id[cols],
      Group = pb_meta$Group[cols],
      value = as.numeric(pb_logcpm[row_id, cols]),
      measure = "Pseudobulk log2-CPM (summed raw counts)",
      stringsAsFactors = FALSE
    )
  }),
  mean_long %>%
    dplyr::inner_join(plot_pairs, by = c("SYMBOL", "cell_type")) %>%
    dplyr::transmute(
      SYMBOL, cell_type, donor_id, Group,
      value = mean_lognorm,
      measure = "Mean log-normalized expression"
    )
)

donor_plot_data$Group <- factor(
  ifelse(donor_plot_data$Group == "Ped", "Pediatric", "Adult"),
  levels = c("Pediatric", "Adult")
)
donor_plot_data$panel <- paste0(donor_plot_data$SYMBOL, " (", donor_plot_data$cell_type, ")")

write.csv(
  donor_plot_data,
  file.path(out_dir, "tables", "donor_level_plot_data.csv"),
  row.names = FALSE
)

p_donor <- ggplot(donor_plot_data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 2, shape = 21) +
  facet_grid(measure ~ panel, scales = "free_y") +
  scale_fill_manual(values = c(Pediatric = "#E8B500", Adult = "#2C6FB7")) +
  labs(x = NULL, y = "Expression (one point per donor)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

ggsave(
  file.path(out_dir, "figures", "donor_level_expression_significant_genes.pdf"),
  p_donor,
  width = 3 + 2.6 * length(unique(donor_plot_data$panel)),
  height = 6.5
)

ggsave(
  file.path(out_dir, "figures", "donor_level_expression_significant_genes.tiff"),
  p_donor,
  width = 3 + 2.6 * length(unique(donor_plot_data$panel)),
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)


# 8. RUN SUMMARY ----

run_summary <- c(
  paste("Signature genes read from meta-analysis:", length(genes_primary_sc)),
  paste("  k = 3:", sum(signature_table$k == 3), " k = 2:", sum(signature_table$k == 2)),
  paste("Raw counts source:", raw_source),
  paste("Minimum cells per donor x cell type:", min_cells),
  paste("Minimum donors per group and cell type:", min_donors),
  paste("Testable cell types:", paste(testable_celltypes, collapse = ", ")),
  paste("Gene status:", paste(names(table(gene_status$status)), table(gene_status$status), sep = " = ", collapse = "; ")),
  paste("PRIMARY tests:", nrow(edger_results),
        "| significant FDR_signature_global < 0.05:", sum(edger_results$FDR_signature_global < 0.05)),
  paste("SECONDARY tests:", nrow(wilcox_results),
        "| significant FDR_global < 0.05:", sum(wilcox_results$FDR_global < 0.05)),
  paste("Cell-type abundance (edgeR QL) FDR < 0.05:", sum(abundance_results$FDR < 0.05),
        paste0("(", paste(abundance_results$cell_type[abundance_results$FDR < 0.05], collapse = ", "), ")")),
  paste("R:", R.version.string),
  paste("edgeR:", packageVersion("edgeR"), "| scuttle:", packageVersion("scuttle"))
)

writeLines(run_summary, file.path(out_dir, "run_summary.txt"))
cat("\n", paste(run_summary, collapse = "\n"), "\n")
