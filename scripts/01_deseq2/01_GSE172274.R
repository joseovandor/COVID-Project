# ============================================================================
# 01_GSE172274_deseq2.R
# Cohort-level RNA-seq analysis for GSE172274
#
# IMPORTANT:
# This version is intentionally written as a LINEAR SCRIPT.
# No custom functions are created.
#
# Run from the COVID-Project RStudio Project.
# The script uses here::here(), so no absolute Windows path is required.
# ============================================================================


# ============================================================================
# 00. Packages
# ============================================================================

required_packages <- c(
  "GEOquery",
  "DESeq2",
  "tidyverse",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "apeglm",
  "edgeR",
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
  library(GEOquery)
  library(DESeq2)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(apeglm)
  library(edgeR)
  library(here)
})


# ============================================================================
# 01. Configuration
# ============================================================================

dataset_id <- "GSE172274"

AGE_CUTOFF <- 18
MIN_COUNT <- 10
MIN_CPM <- 1
ALPHA <- 0.05

raw_counts_filename <- "GSE172274_raw_counts_GRCh38.p13_NCBI.tsv"

project_root <- here::here()

data_dir <- here::here("data", "raw", dataset_id)

results_root <- here::here("results", "deseq2", dataset_id)
metadata_dir <- file.path(results_root, "metadata")
primary_dir <- file.path(results_root, "primary")
sensitivity_dir <- file.path(results_root, "sensitivity")
diagnostics_dir <- file.path(results_root, "diagnostics")
figures_dir <- file.path(results_root, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(primary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_root)
message("Results directory: ", results_root)


# ============================================================================
# 02. Download GEO metadata
# ============================================================================

message("Downloading GEO metadata for ", dataset_id, "...")

gse <- GEOquery::getGEO(
  dataset_id,
  GSEMatrix = TRUE
)

if (length(gse) < 1) {
  stop("No ExpressionSet object was returned by GEOquery.")
}

if (length(gse) > 1) {
  warning(
    dataset_id,
    " returned ",
    length(gse),
    " ExpressionSet objects. The first object will be used."
  )
}

gse_data <- gse[[1]]
metadata_raw <- Biobase::pData(gse_data)

write.csv(
  metadata_raw,
  file.path(metadata_dir, paste0(dataset_id, "_metadata_raw.csv")),
  row.names = FALSE
)


# ============================================================================
# 03. Inspect metadata columns
# ============================================================================

metadata_columns <- colnames(metadata_raw)

writeLines(
  metadata_columns,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_metadata_columns.txt")
  )
)

# Expected column names based on the current GEO metadata.
# If GEO changes these names, the script will stop rather than silently guess.

required_metadata_columns <- c(
  "geo_accession",
  "title",
  "diagnosis:ch1",
  "age:ch1",
  "Sex:ch1",
  "tissue:ch1"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  metadata_columns
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata_columns, collapse = ", "),
    "\nInspect ",
    file.path(
      diagnostics_dir,
      paste0(dataset_id, "_metadata_columns.txt")
    )
  )
}

# Respiratory score is useful if available, but it is not required to run DESeq2.
respiratory_score_column <- grep(
  "respiratory.*score|severity.*score|severity",
  metadata_columns,
  ignore.case = TRUE,
  value = TRUE
)

if (length(respiratory_score_column) > 0) {
  respiratory_score_column <- respiratory_score_column[1]
} else {
  respiratory_score_column <- NA_character_
}


# ============================================================================
# 04. Curate metadata
# ============================================================================

sample_metadata <- data.frame(
  GEO_ID = as.character(metadata_raw$geo_accession),
  Sample = as.character(metadata_raw$title),
  Diagnosis = as.character(metadata_raw[["diagnosis:ch1"]]),
  Age = suppressWarnings(as.numeric(metadata_raw[["age:ch1"]])),
  Sex = as.character(metadata_raw[["Sex:ch1"]]),
  Tissue = as.character(metadata_raw[["tissue:ch1"]]),
  stringsAsFactors = FALSE
)

if (!is.na(respiratory_score_column)) {
  sample_metadata$RespiratoryScore <- as.character(
    metadata_raw[[respiratory_score_column]]
  )
} else {
  sample_metadata$RespiratoryScore <- NA_character_
}

sample_metadata$Sex <- trimws(sample_metadata$Sex)

sample_metadata$Sex <- ifelse(
  tolower(sample_metadata$Sex) %in% c("male", "m"),
  "Male",
  ifelse(
    tolower(sample_metadata$Sex) %in% c("female", "f"),
    "Female",
    NA
  )
)

sample_metadata$Group <- ifelse(
  !is.na(sample_metadata$Age) & sample_metadata$Age < AGE_CUTOFF,
  "Pediatric",
  ifelse(
    !is.na(sample_metadata$Age) & sample_metadata$Age >= AGE_CUTOFF,
    "Adult",
    NA
  )
)

sample_metadata$IsPooled <- grepl(
  "POOL",
  sample_metadata$Sample,
  ignore.case = TRUE
)


# ============================================================================
# 05. Document diagnosis labels
# ============================================================================

diagnosis_summary <- sample_metadata %>%
  count(Diagnosis, name = "N")

write.csv(
  diagnosis_summary,
  file.path(metadata_dir, paste0(dataset_id, "_diagnosis_summary.csv")),
  row.names = FALSE
)

print(diagnosis_summary)

# No automatic "Positive"/"Negative" recoding is performed.
# GEO diagnosis labels are preserved exactly as provided.


# ============================================================================
# 06. Exclude samples lacking age/group/sex for primary adjusted model
# ============================================================================

excluded_metadata <- sample_metadata[
  is.na(sample_metadata$Age) |
    is.na(sample_metadata$Group) |
    is.na(sample_metadata$Sex),
  ,
  drop = FALSE
]

write.csv(
  excluded_metadata,
  file.path(metadata_dir, paste0(dataset_id, "_excluded_samples_primary.csv")),
  row.names = FALSE
)

sample_metadata_primary <- sample_metadata[
  !is.na(sample_metadata$Age) &
    !is.na(sample_metadata$Group) &
    !is.na(sample_metadata$Sex),
  ,
  drop = FALSE
]

sample_metadata_primary$Group <- factor(
  sample_metadata_primary$Group,
  levels = c("Adult", "Pediatric")
)

sample_metadata_primary$Sex <- factor(
  sample_metadata_primary$Sex
)

if (nrow(sample_metadata_primary) == 0) {
  stop("No samples remained after metadata quality control.")
}


# ============================================================================
# 07. Load raw counts
# ============================================================================

counts_file <- file.path(
  data_dir,
  raw_counts_filename
)

if (!file.exists(counts_file)) {
  stop(
    "Raw-count file not found:\n",
    counts_file,
    "\n\nPlace the file here before running the script:\n",
    data_dir
  )
}

raw_counts <- readr::read_tsv(
  counts_file,
  show_col_types = FALSE
)

if (!"GeneID" %in% colnames(raw_counts)) {
  stop("The raw-count file does not contain a GeneID column.")
}

raw_counts$GeneID <- as.character(raw_counts$GeneID)

if (anyDuplicated(raw_counts$GeneID) > 0) {
  duplicated_gene_ids <- unique(
    raw_counts$GeneID[
      duplicated(raw_counts$GeneID)
    ]
  )

  writeLines(
    duplicated_gene_ids,
    file.path(
      diagnostics_dir,
      paste0(dataset_id, "_duplicated_raw_GeneID.txt")
    )
  )

  stop(
    "Duplicated GeneID values were found in the raw count matrix."
  )
}

rownames(raw_counts) <- raw_counts$GeneID
raw_counts$GeneID <- NULL

counts_matrix <- as.matrix(raw_counts)
storage.mode(counts_matrix) <- "numeric"

if (any(is.na(counts_matrix))) {
  stop("NA values were detected in the raw count matrix.")
}

if (any(counts_matrix < 0)) {
  stop("Negative values were detected in the raw count matrix.")
}

if (any(abs(counts_matrix - round(counts_matrix)) > 1e-8)) {
  stop(
    "Non-integer values were detected. ",
    "DESeq2 requires raw integer gene-level counts."
  )
}

counts_matrix <- round(counts_matrix)


# ============================================================================
# 08. Match metadata and counts
# ============================================================================

samples_in_both <- intersect(
  sample_metadata_primary$GEO_ID,
  colnames(counts_matrix)
)

missing_from_counts <- setdiff(
  sample_metadata_primary$GEO_ID,
  colnames(counts_matrix)
)

missing_from_metadata <- setdiff(
  colnames(counts_matrix),
  sample_metadata_primary$GEO_ID
)

writeLines(
  missing_from_counts,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_samples_missing_from_counts.txt")
  )
)

writeLines(
  missing_from_metadata,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_count_columns_not_used.txt")
  )
)

sample_metadata_primary <- sample_metadata_primary[
  sample_metadata_primary$GEO_ID %in% samples_in_both,
  ,
  drop = FALSE
]

sample_metadata_primary <- sample_metadata_primary[
  match(samples_in_both, sample_metadata_primary$GEO_ID),
  ,
  drop = FALSE
]

counts_matrix <- counts_matrix[
  ,
  sample_metadata_primary$GEO_ID,
  drop = FALSE
]

rownames(sample_metadata_primary) <- sample_metadata_primary$GEO_ID

stopifnot(
  identical(
    colnames(counts_matrix),
    rownames(sample_metadata_primary)
  )
)


# ============================================================================
# 09. Cohort composition
# ============================================================================

cohort_summary <- sample_metadata_primary %>%
  count(
    Group,
    Sex,
    IsPooled,
    name = "N"
  )

write.csv(
  sample_metadata_primary,
  file.path(metadata_dir, paste0(dataset_id, "_sample_metadata.csv")),
  row.names = FALSE
)

write.csv(
  cohort_summary,
  file.path(metadata_dir, paste0(dataset_id, "_cohort_summary.csv")),
  row.names = FALSE
)

print(cohort_summary)

group_sex_table <- table(
  sample_metadata_primary$Group,
  sample_metadata_primary$Sex
)

capture.output(
  group_sex_table,
  file = file.path(
    diagnostics_dir,
    paste0(dataset_id, "_group_by_sex_table.txt")
  )
)


# ============================================================================
# 10. Primary low-count filtering
# ============================================================================

group_sizes <- table(sample_metadata_primary$Group)

if (length(group_sizes) != 2) {
  stop("Both Adult and Pediatric groups must be represented.")
}

min_group_size <- min(group_sizes)

keep_primary <- rowSums(
  counts_matrix >= MIN_COUNT
) >= min_group_size

counts_primary <- counts_matrix[
  keep_primary,
  ,
  drop = FALSE
]


# ============================================================================
# 11. CPM-based sensitivity filter
# ============================================================================

cpm_matrix <- edgeR::cpm(counts_matrix)

keep_cpm <- rowSums(
  cpm_matrix >= MIN_CPM
) >= min_group_size

counts_cpm <- counts_matrix[
  keep_cpm,
  ,
  drop = FALSE
]

filtering_summary <- data.frame(
  Dataset = dataset_id,
  Total_samples = ncol(counts_matrix),
  Adult_samples = sum(sample_metadata_primary$Group == "Adult"),
  Pediatric_samples = sum(sample_metadata_primary$Group == "Pediatric"),
  Minimum_group_size = min_group_size,
  Genes_before_filtering = nrow(counts_matrix),
  Primary_MIN_COUNT = MIN_COUNT,
  Primary_genes_after_filtering = nrow(counts_primary),
  Primary_genes_removed = nrow(counts_matrix) - nrow(counts_primary),
  Sensitivity_MIN_CPM = MIN_CPM,
  CPM_genes_after_filtering = nrow(counts_cpm),
  CPM_genes_removed = nrow(counts_matrix) - nrow(counts_cpm)
)

write.csv(
  filtering_summary,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_low_count_filtering_summary.csv")
  ),
  row.names = FALSE
)

print(filtering_summary)


# ============================================================================
# 12. Determine whether sex adjustment is feasible
# ============================================================================

sex_adjustment_feasible <- (
  nlevels(droplevels(sample_metadata_primary$Sex)) >= 2
)

if (sex_adjustment_feasible) {
  primary_design <- ~ Sex + Group
} else {
  primary_design <- ~ Group

  warning(
    "Sex adjustment is not feasible. ",
    "Primary model will use ~ Group."
  )
}


# ============================================================================
# 13. Primary DESeq2 model
# ============================================================================

dds_primary <- DESeqDataSetFromMatrix(
  countData = counts_primary,
  colData = sample_metadata_primary,
  design = primary_design
)

dds_primary <- DESeq(
  dds_primary,
  test = "Wald"
)

print(resultsNames(dds_primary))

res_primary <- results(
  dds_primary,
  contrast = c(
    "Group",
    "Pediatric",
    "Adult"
  ),
  alpha = ALPHA,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

res_primary_df <- as.data.frame(res_primary)
res_primary_df$ENTREZID <- rownames(res_primary_df)
rownames(res_primary_df) <- NULL

res_primary_df <- res_primary_df %>%
  relocate(ENTREZID)


# ============================================================================
# 14. Annotation diagnostics
# ============================================================================

annotation_list <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = as.character(res_primary_df$ENTREZID),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "CharacterList"
)

annotation_diagnostics <- data.frame(
  ENTREZID = names(annotation_list),
  N_SYMBOLS = lengths(annotation_list),
  SYMBOLS = vapply(
    annotation_list,
    function(x) paste(as.character(x), collapse = ";"),
    character(1)
  ),
  stringsAsFactors = FALSE
)

annotation_diagnostics$AnnotationStatus <- ifelse(
  annotation_diagnostics$N_SYMBOLS == 0,
  "No_symbol",
  ifelse(
    annotation_diagnostics$N_SYMBOLS == 1,
    "Unique_symbol",
    "Ambiguous_symbol"
  )
)

write.csv(
  annotation_diagnostics,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_annotation_diagnostics.csv")
  ),
  row.names = FALSE
)

unique_annotation <- annotation_diagnostics[
  annotation_diagnostics$N_SYMBOLS == 1,
  c("ENTREZID", "SYMBOLS"),
  drop = FALSE
]

colnames(unique_annotation)[2] <- "SYMBOL"

res_primary_annotated <- merge(
  res_primary_df,
  unique_annotation,
  by = "ENTREZID",
  all.x = TRUE,
  sort = FALSE
)

res_primary_annotated <- res_primary_annotated[
  match(res_primary_df$ENTREZID, res_primary_annotated$ENTREZID),
  ,
  drop = FALSE
]

duplicate_symbols <- res_primary_annotated %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  count(
    SYMBOL,
    name = "N_ENTREZ_IDS"
  ) %>%
  filter(
    N_ENTREZ_IDS > 1
  )

write.csv(
  duplicate_symbols,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_duplicated_symbols.csv")
  ),
  row.names = FALSE
)

res_primary_annotated$MetaEligibleAnnotation <- (
  !is.na(res_primary_annotated$SYMBOL) &
    res_primary_annotated$SYMBOL != "" &
    !(res_primary_annotated$SYMBOL %in% duplicate_symbols$SYMBOL)
)


# ============================================================================
# 15. Export COMPLETE primary DESeq2 results
# ============================================================================

write.csv(
  res_primary_annotated,
  file.path(
    primary_dir,
    paste0(dataset_id, "_DESeq2_complete.csv")
  ),
  row.names = FALSE
)


# ============================================================================
# 16. Export meta-analysis input
# ============================================================================

# No filtering by p-value, FDR, direction, or effect magnitude.

meta_input <- res_primary_annotated %>%
  filter(
    MetaEligibleAnnotation,
    !is.na(log2FoldChange),
    !is.na(lfcSE),
    is.finite(log2FoldChange),
    is.finite(lfcSE),
    lfcSE > 0
  ) %>%
  transmute(
    Dataset = dataset_id,
    Original_ID = ENTREZID,
    SYMBOL,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    MetaEligibleAnnotation
  )

write.csv(
  meta_input,
  file.path(
    primary_dir,
    paste0(dataset_id, "_meta_input.csv")
  ),
  row.names = FALSE
)


# ============================================================================
# 17. Shrinkage for visualization only
# ============================================================================

expected_coef <- "Group_Pediatric_vs_Adult"

if (expected_coef %in% resultsNames(dds_primary)) {

  res_shrunk <- lfcShrink(
    dds_primary,
    coef = expected_coef,
    type = "apeglm"
  )

  res_shrunk_df <- as.data.frame(res_shrunk)
  res_shrunk_df$ENTREZID <- rownames(res_shrunk_df)
  rownames(res_shrunk_df) <- NULL

  res_shrunk_df <- res_shrunk_df %>%
    relocate(ENTREZID)

  res_shrunk_df <- merge(
    res_shrunk_df,
    unique_annotation,
    by = "ENTREZID",
    all.x = TRUE,
    sort = FALSE
  )

  write.csv(
    res_shrunk_df,
    file.path(
      primary_dir,
      paste0(dataset_id, "_DESeq2_shrunken.csv")
    ),
    row.names = FALSE
  )
}


# ============================================================================
# 18. Normalized counts
# ============================================================================

normalized_counts <- counts(
  dds_primary,
  normalized = TRUE
)

normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$ENTREZID <- rownames(normalized_counts_df)
rownames(normalized_counts_df) <- NULL

normalized_counts_df <- normalized_counts_df %>%
  relocate(ENTREZID)

write.csv(
  normalized_counts_df,
  file.path(
    primary_dir,
    paste0(dataset_id, "_normalized_counts.csv")
  ),
  row.names = FALSE
)


# ============================================================================
# 19. Size factors
# ============================================================================

size_factor_table <- data.frame(
  Sample = names(sizeFactors(dds_primary)),
  SizeFactor = sizeFactors(dds_primary)
)

write.csv(
  size_factor_table,
  file.path(
    primary_dir,
    paste0(dataset_id, "_size_factors.csv")
  ),
  row.names = FALSE
)


# ============================================================================
# 20. Dispersion diagnostics
# ============================================================================

tiff(
  file.path(
    figures_dir,
    paste0(dataset_id, "_dispersion_estimates.tiff")
  ),
  res = 600,
  width = 2500,
  height = 2000,
  compression = "lzw"
)

plotDispEsts(dds_primary)

dev.off()

dispersion_fit_type <- attr(
  dispersionFunction(dds_primary),
  "fitType"
)

if (is.null(dispersion_fit_type)) {
  dispersion_fit_type <- NA_character_
}


# ============================================================================
# 21. Independent filtering
# ============================================================================

res_metadata <- S4Vectors::metadata(res_primary)

filter_threshold <- res_metadata$filterThreshold

if (!is.null(filter_threshold)) {
  independent_filter_threshold <- as.numeric(filter_threshold)
} else {
  independent_filter_threshold <- NA_real_
}


# ============================================================================
# 22. Cook's distance diagnostics
# ============================================================================

cooks_matrix <- assays(dds_primary)[["cooks"]]

if (!is.null(cooks_matrix)) {

  max_cooks <- apply(
    cooks_matrix,
    1,
    max,
    na.rm = TRUE
  )

  cooks_diagnostics <- data.frame(
    ENTREZID = names(max_cooks),
    MaxCooksDistance = max_cooks
  )

  write.csv(
    cooks_diagnostics,
    file.path(
      diagnostics_dir,
      paste0(dataset_id, "_cooks_distance.csv")
    ),
    row.names = FALSE
  )
}


# ============================================================================
# 23. PCA
# ============================================================================

vsd <- vst(
  dds_primary,
  blind = TRUE
)

pca <- prcomp(
  t(assay(vsd)),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- 100 * (
  pca$sdev^2 / sum(pca$sdev^2)
)

pca_scores <- as.data.frame(pca$x)
pca_scores$GEO_ID <- rownames(pca_scores)
rownames(pca_scores) <- NULL

pca_scores <- merge(
  pca_scores,
  sample_metadata_primary,
  by = "GEO_ID",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  pca_scores,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_PCA_scores.csv")
  ),
  row.names = FALSE
)

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Group
  )
) +
  geom_point(size = 3) +
  ggtitle(dataset_id) +
  xlab(
    paste0(
      "PC1: ",
      round(pca_variance[1], 1),
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      round(pca_variance[2], 1),
      "% variance"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank()
  ) +
  labs(
    colour = "Age group"
  )

ggsave(
  filename = file.path(
    figures_dir,
    paste0(dataset_id, "_PCA.tiff")
  ),
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

if (sex_adjustment_feasible) {

  pca_pc1_model <- lm(
    PC1 ~ Age + Sex,
    data = pca_scores
  )

  pca_pc2_model <- lm(
    PC2 ~ Age + Sex,
    data = pca_scores
  )

} else {

  pca_pc1_model <- lm(
    PC1 ~ Age,
    data = pca_scores
  )

  pca_pc2_model <- lm(
    PC2 ~ Age,
    data = pca_scores
  )
}

writeLines(
  c(
    "===== PC1 =====",
    capture.output(summary(pca_pc1_model)),
    "",
    "===== PC2 =====",
    capture.output(summary(pca_pc2_model))
  ),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_PCA_age_models.txt")
  )
)


# ============================================================================
# 24. Sensitivity: continuous age
# ============================================================================

metadata_age <- sample_metadata_primary
metadata_age$AgeScaled <- as.numeric(
  scale(metadata_age$Age)
)

if (sex_adjustment_feasible) {
  age_design <- ~ Sex + AgeScaled
} else {
  age_design <- ~ AgeScaled
}

dds_age <- DESeqDataSetFromMatrix(
  countData = counts_primary,
  colData = metadata_age,
  design = age_design
)

dds_age <- DESeq(
  dds_age,
  test = "Wald"
)

age_coef <- grep(
  "AgeScaled",
  resultsNames(dds_age),
  value = TRUE
)

if (length(age_coef) == 1) {

  res_age <- results(
    dds_age,
    name = age_coef,
    alpha = ALPHA,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )

  res_age_df <- as.data.frame(res_age)
  res_age_df$ENTREZID <- rownames(res_age_df)
  rownames(res_age_df) <- NULL

  res_age_df <- merge(
    res_age_df,
    unique_annotation,
    by = "ENTREZID",
    all.x = TRUE,
    sort = FALSE
  )

  write.csv(
    res_age_df,
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_DESeq2_continuous_age.csv")
    ),
    row.names = FALSE
  )
}


# ============================================================================
# 25. Sensitivity: CPM filter
# ============================================================================

dds_cpm <- DESeqDataSetFromMatrix(
  countData = counts_cpm,
  colData = sample_metadata_primary,
  design = primary_design
)

dds_cpm <- DESeq(
  dds_cpm,
  test = "Wald"
)

res_cpm <- results(
  dds_cpm,
  contrast = c(
    "Group",
    "Pediatric",
    "Adult"
  ),
  alpha = ALPHA,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

res_cpm_df <- as.data.frame(res_cpm)
res_cpm_df$ENTREZID <- rownames(res_cpm_df)
rownames(res_cpm_df) <- NULL

res_cpm_df <- merge(
  res_cpm_df,
  unique_annotation,
  by = "ENTREZID",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  res_cpm_df,
  file.path(
    sensitivity_dir,
    paste0(dataset_id, "_DESeq2_CPM_filter.csv")
  ),
  row.names = FALSE
)


# ============================================================================
# 26. Sensitivity: exclude pooled samples
# ============================================================================

metadata_no_pool <- sample_metadata_primary[
  !sample_metadata_primary$IsPooled,
  ,
  drop = FALSE
]

if (
  nrow(metadata_no_pool) >= 4 &&
    length(unique(metadata_no_pool$Group)) == 2
) {

  rownames(metadata_no_pool) <- metadata_no_pool$GEO_ID

  counts_no_pool_unfiltered <- counts_matrix[
    ,
    metadata_no_pool$GEO_ID,
    drop = FALSE
  ]

  no_pool_group_sizes <- table(metadata_no_pool$Group)
  no_pool_min_group_size <- min(no_pool_group_sizes)

  keep_no_pool <- rowSums(
    counts_no_pool_unfiltered >= MIN_COUNT
  ) >= no_pool_min_group_size

  counts_no_pool <- counts_no_pool_unfiltered[
    keep_no_pool,
    ,
    drop = FALSE
  ]

  no_pool_sex_feasible <- (
    nlevels(droplevels(metadata_no_pool$Sex)) >= 2
  )

  if (no_pool_sex_feasible) {
    no_pool_design <- ~ Sex + Group
  } else {
    no_pool_design <- ~ Group
  }

  dds_no_pool <- DESeqDataSetFromMatrix(
    countData = counts_no_pool,
    colData = metadata_no_pool,
    design = no_pool_design
  )

  dds_no_pool <- DESeq(
    dds_no_pool,
    test = "Wald"
  )

  res_no_pool <- results(
    dds_no_pool,
    contrast = c(
      "Group",
      "Pediatric",
      "Adult"
    ),
    alpha = ALPHA,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )

  res_no_pool_df <- as.data.frame(res_no_pool)
  res_no_pool_df$ENTREZID <- rownames(res_no_pool_df)
  rownames(res_no_pool_df) <- NULL

  res_no_pool_df <- merge(
    res_no_pool_df,
    unique_annotation,
    by = "ENTREZID",
    all.x = TRUE,
    sort = FALSE
  )

  write.csv(
    res_no_pool_df,
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_DESeq2_without_pools.csv")
    ),
    row.names = FALSE
  )

  no_pool_summary <- data.frame(
    Dataset = dataset_id,
    Samples_without_pools = ncol(counts_no_pool),
    Adult_samples_without_pools =
      sum(metadata_no_pool$Group == "Adult"),
    Pediatric_samples_without_pools =
      sum(metadata_no_pool$Group == "Pediatric"),
    Genes_after_filtering_without_pools =
      nrow(counts_no_pool)
  )

  write.csv(
    no_pool_summary,
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_without_pools_summary.csv")
    ),
    row.names = FALSE
  )

} else {

  writeLines(
    "No-pool sensitivity analysis was not estimable.",
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_without_pools_NOT_ESTIMABLE.txt")
    )
  )
}


# ============================================================================
# 27. Sensitivity: Group x Sex interaction
# ============================================================================

interaction_feasible <- (
  nlevels(droplevels(sample_metadata_primary$Sex)) == 2 &&
    nlevels(droplevels(sample_metadata_primary$Group)) == 2 &&
    all(group_sex_table > 0)
)

if (interaction_feasible) {

  dds_interaction <- DESeqDataSetFromMatrix(
    countData = counts_primary,
    colData = sample_metadata_primary,
    design = ~ Sex + Group + Sex:Group
  )

  dds_interaction <- DESeq(
    dds_interaction,
    test = "Wald"
  )

  interaction_coef <- grep(
    "Group.*Sex|Sex.*Group",
    resultsNames(dds_interaction),
    value = TRUE
  )

  if (length(interaction_coef) == 1) {

    res_interaction <- results(
      dds_interaction,
      name = interaction_coef,
      alpha = ALPHA,
      independentFiltering = TRUE,
      cooksCutoff = TRUE
    )

    res_interaction_df <- as.data.frame(res_interaction)
    res_interaction_df$ENTREZID <- rownames(res_interaction_df)
    rownames(res_interaction_df) <- NULL

    res_interaction_df <- merge(
      res_interaction_df,
      unique_annotation,
      by = "ENTREZID",
      all.x = TRUE,
      sort = FALSE
    )

    write.csv(
      res_interaction_df,
      file.path(
        sensitivity_dir,
        paste0(dataset_id, "_DESeq2_group_sex_interaction.csv")
      ),
      row.names = FALSE
    )

    writeLines(
      paste(
        "Interaction coefficient:",
        interaction_coef
      ),
      file.path(
        sensitivity_dir,
        paste0(dataset_id, "_group_sex_interaction_info.txt")
      )
    )

  } else {

    writeLines(
      "Interaction coefficient could not be identified uniquely.",
      file.path(
        sensitivity_dir,
        paste0(dataset_id, "_group_sex_interaction_NOT_ESTIMABLE.txt")
      )
    )
  }

} else {

  writeLines(
    "Group x Sex interaction was not estimable because at least one Group x Sex cell was empty.",
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_group_sex_interaction_NOT_ESTIMABLE.txt")
    )
  )
}


# ============================================================================
# 28. Workflow summary
# ============================================================================

workflow_summary <- c(
  paste("Dataset:", dataset_id),
  paste("Project root:", project_root),
  paste(
    "Age-group definition:",
    paste0(
      "Pediatric < ",
      AGE_CUTOFF,
      " years; Adult >= ",
      AGE_CUTOFF,
      " years"
    )
  ),
  paste(
    "Primary contrast:",
    "Pediatric vs Adult"
  ),
  paste(
    "Positive log2FoldChange:",
    "higher expression in Pediatric samples"
  ),
  paste(
    "Primary samples:",
    nrow(sample_metadata_primary)
  ),
  paste(
    "Adult samples:",
    sum(sample_metadata_primary$Group == "Adult")
  ),
  paste(
    "Pediatric samples:",
    sum(sample_metadata_primary$Group == "Pediatric")
  ),
  paste(
    "Pooled libraries:",
    sum(sample_metadata_primary$IsPooled, na.rm = TRUE)
  ),
  paste(
    "Genes before low-count filtering:",
    nrow(counts_matrix)
  ),
  paste(
    "Primary low-count rule:",
    paste0(
      "raw count >= ",
      MIN_COUNT,
      " in at least ",
      min_group_size,
      " samples"
    )
  ),
  paste(
    "Genes after primary filtering:",
    nrow(counts_primary)
  ),
  paste(
    "CPM sensitivity rule:",
    paste0(
      "CPM >= ",
      MIN_CPM,
      " in at least ",
      min_group_size,
      " samples"
    )
  ),
  paste(
    "Genes after CPM filtering:",
    nrow(counts_cpm)
  ),
  paste(
    "Primary design formula:",
    paste(
      deparse(primary_design),
      collapse = ""
    )
  ),
  paste(
    "Sex adjustment feasible:",
    sex_adjustment_feasible
  ),
  paste(
    "DESeq2 test:",
    "Wald test"
  ),
  paste(
    "FDR correction:",
    "Benjamini-Hochberg"
  ),
  paste(
    "Independent filtering:",
    "enabled"
  ),
  paste(
    "Independent filtering threshold:",
    independent_filter_threshold
  ),
  paste(
    "Dispersion fit type:",
    dispersion_fit_type
  ),
  paste(
    "Log2FC shrinkage:",
    "apeglm for visualization only"
  ),
  paste(
    "Meta-analysis estimates:",
    "unshrunken DESeq2 log2FC and lfcSE"
  ),
  paste(
    "Meta-analysis preselection:",
    "no p-value, FDR, direction, or effect-magnitude filtering"
  ),
  paste(
    "Unique mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "Unique_symbol"
    )
  ),
  paste(
    "Ambiguous mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "Ambiguous_symbol"
    )
  ),
  paste(
    "Missing SYMBOL mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "No_symbol"
    )
  ),
  paste(
    "Genes eligible for later SYMBOL-based harmonization:",
    nrow(meta_input)
  )
)

writeLines(
  workflow_summary,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_workflow_summary.txt")
  )
)


# ============================================================================
# 29. Reproducibility
# ============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_sessionInfo.txt")
  )
)

writeLines(
  c(
    paste("R version:", R.version.string),
    paste("DESeq2 version:", as.character(packageVersion("DESeq2"))),
    paste("GEOquery version:", as.character(packageVersion("GEOquery"))),
    paste("edgeR version:", as.character(packageVersion("edgeR"))),
    paste("apeglm version:", as.character(packageVersion("apeglm"))),
    paste("AnnotationDbi version:", as.character(packageVersion("AnnotationDbi"))),
    paste("org.Hs.eg.db version:", as.character(packageVersion("org.Hs.eg.db"))),
    paste("here version:", as.character(packageVersion("here")))
  ),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_package_versions.txt")
  )
)


# ============================================================================
# 30. Final summary
# ============================================================================

message("")
message("============================================================")
message(dataset_id, " analysis completed successfully.")
message("============================================================")
message("Project root: ", project_root)
message("Primary samples: ", nrow(sample_metadata_primary))
message(
  "Adult / Pediatric: ",
  sum(sample_metadata_primary$Group == "Adult"),
  " / ",
  sum(sample_metadata_primary$Group == "Pediatric")
)
message(
  "Genes before filtering: ",
  nrow(counts_matrix)
)
message(
  "Genes after primary filtering: ",
  nrow(counts_primary)
)
message(
  "Genes eligible for later meta-analysis harmonization: ",
  nrow(meta_input)
)
message(
  "Results saved to: ",
  results_root
)
message("============================================================")
