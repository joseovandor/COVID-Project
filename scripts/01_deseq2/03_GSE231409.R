# =============================================================================
# 03_GSE231409_deseq2.R
#
# Differential expression analysis of GSE231409
#
# This script reproduces the cohort-level RNA-seq analysis used for the
# cross-cohort comparison of pediatric and adult SARS-CoV-2-positive samples.
#
# Primary analysis
#   - Samples restricted to SARS-CoV-2-positive nasal/nasopharyngeal swabs.
#   - Pediatric: age < 18 years.
#   - Adult: age >= 18 years.
#   - DESeq2 design: ~ Batch + Sex + Group, when estimable.
#   - Contrast: Pediatric vs Adult.
#   - Positive log2FoldChange: higher expression in Pediatric samples.
#   - Low-count filter: raw count >= 10 in at least the size of the smaller
#     age group.
#
# Sensitivity analyses
#   - Age modeled as a continuous variable.
#   - Alternative prevalence filter: CPM >= 1 in at least the size of the
#     smaller age group.
#   - Group-by-sex interaction, when estimable.
#
# Downstream outputs
#   - Complete cohort-level DESeq2 results.
#   - Unshrunken log2FC and SE for cross-cohort meta-analysis.
#   - Genome-wide ranked statistics for pre-ranked GSEA.
#   - Normalized-count and VST matrices.
#   - QC, PCA, annotation, batch, and sensitivity diagnostics.
#
# Repository use
#   Run from the root of the COVID-Project RStudio Project. Paths are resolved
#   with here::here(); no machine-specific working directory is required.
#
# The script is intentionally linear and does not define custom functions.
# =============================================================================

SCRIPT_BUILD <- "GSE231409_REPO_FINAL_2026-08-17_v2"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

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


# =============================================================================
# 01. Configuration
# =============================================================================

dataset_id <- "GSE231409"

AGE_CUTOFF <- 18
MIN_COUNT <- 10
MIN_CPM <- 1
ALPHA <- 0.05

raw_counts_filename <- "GSE231409_BRAVE_RNASeq_counts.csv.gz"

raw_counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/",
  "GSE231409/suppl/",
  "GSE231409%5FBRAVE%5FRNASeq%5Fcounts%2Ecsv%2Egz"
)

project_root <- here::here()

data_dir <- here::here("data", "raw", dataset_id)

results_root <- here::here("results", "deseq2", dataset_id)
metadata_dir <- file.path(results_root, "metadata")
primary_dir <- file.path(results_root, "primary")
sensitivity_dir <- file.path(results_root, "sensitivity")
diagnostics_dir <- file.path(results_root, "diagnostics")
figures_dir <- file.path(results_root, "figures")
sensitivity_figures_dir <- file.path(figures_dir, "sensitivity")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(primary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_root)
message("Results directory: ", results_root)


# =============================================================================
# 02. Download GEO metadata
# =============================================================================

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
  file.path(
    metadata_dir,
    paste0(dataset_id, "_metadata_raw.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 03. Inspect metadata columns
# =============================================================================

metadata_columns <- colnames(metadata_raw)

writeLines(
  metadata_columns,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_metadata_columns.txt")
  )
)

required_metadata_columns <- c(
  "geo_accession",
  "title",
  "covid:ch1",
  "age:ch1",
  "Sex:ch1",
  "tissue:ch1",
  "batch:ch1"
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


# =============================================================================
# 04. Curate metadata
# =============================================================================

sample_metadata <- data.frame(
  GEO_ID = as.character(metadata_raw$geo_accession),
  Sample = as.character(metadata_raw$title),
  CovidRaw = as.character(metadata_raw[["covid:ch1"]]),
  Age = suppressWarnings(as.numeric(metadata_raw[["age:ch1"]])),
  Sex = as.character(metadata_raw[["Sex:ch1"]]),
  Tissue = as.character(metadata_raw[["tissue:ch1"]]),
  Batch = as.character(metadata_raw[["batch:ch1"]]),
  stringsAsFactors = FALSE
)

sample_metadata$Sex <- trimws(sample_metadata$Sex)
sample_metadata$Tissue <- trimws(sample_metadata$Tissue)
sample_metadata$Batch <- trimws(sample_metadata$Batch)

sample_metadata$Sex <- ifelse(
  tolower(sample_metadata$Sex) %in% c("male", "m"),
  "Male",
  ifelse(
    tolower(sample_metadata$Sex) %in% c("female", "f"),
    "Female",
    NA
  )
)

sample_metadata$InfectionStatus <- ifelse(
  grepl(
    "Positive",
    sample_metadata$CovidRaw,
    ignore.case = TRUE
  ),
  "SC2_Pos",
  ifelse(
    grepl(
      "Negative",
      sample_metadata$CovidRaw,
      ignore.case = TRUE
    ),
    "SC2_Neg",
    NA
  )
)

sample_metadata$Group <- ifelse(
  !is.na(sample_metadata$Age) &
    sample_metadata$Age < AGE_CUTOFF,
  "Pediatric",
  ifelse(
    !is.na(sample_metadata$Age) &
      sample_metadata$Age >= AGE_CUTOFF,
    "Adult",
    NA
  )
)


# =============================================================================
# 05. Document cohort labels
# =============================================================================

infection_summary <- sample_metadata %>%
  count(
    InfectionStatus,
    name = "N"
  )

tissue_summary <- sample_metadata %>%
  count(
    Tissue,
    name = "N"
  )

batch_summary_all <- sample_metadata %>%
  count(
    Batch,
    name = "N"
  )

write.csv(
  infection_summary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_infection_summary.csv")
  ),
  row.names = FALSE
)

write.csv(
  tissue_summary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_tissue_summary.csv")
  ),
  row.names = FALSE
)

write.csv(
  batch_summary_all,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_batch_summary_all_samples.csv")
  ),
  row.names = FALSE
)

print(infection_summary)
print(tissue_summary)
print(batch_summary_all)


# =============================================================================
# 06. Restrict cohort using the original analysis logic
# =============================================================================

# The cohort is defined exactly as in the original GSE231409 analysis:
#   1. retain recognized SARS-CoV-2 status;
#   2. retain Nasal swab / Nasopharyngeal swab samples;
#   3. define Pediatric as age < 18 and Adult as age >= 18;
#   4. retain SARS-CoV-2-positive Pediatric and Adult samples.
#
# Sex and Batch are NOT used to remove samples at this stage. They are retained
# as model covariates and validated after the cohort has been defined.

eligible_tissues <- c(
  "Nasal swab",
  "Nasopharyngeal swab"
)

sample_metadata_original_logic <- sample_metadata[
  sample_metadata$InfectionStatus %in% c(
    "SC2_Pos",
    "SC2_Neg"
  ),
  ,
  drop = FALSE
]

sample_metadata_original_logic <- sample_metadata_original_logic[
  sample_metadata_original_logic$Tissue %in% eligible_tissues,
  ,
  drop = FALSE
]

sample_metadata_original_logic$CombinedStatus <- paste(
  sample_metadata_original_logic$InfectionStatus,
  sample_metadata_original_logic$Group,
  sep = "_"
)

sample_metadata_primary <- sample_metadata_original_logic[
  sample_metadata_original_logic$CombinedStatus %in% c(
    "SC2_Pos_Pediatric",
    "SC2_Pos_Adult"
  ),
  ,
  drop = FALSE
]

excluded_metadata <- sample_metadata[
  !(sample_metadata$Sample %in% sample_metadata_primary$Sample),
  ,
  drop = FALSE
]

write.csv(
  excluded_metadata,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_excluded_samples_primary.csv")
  ),
  row.names = FALSE
)

sample_metadata_primary$Group <- factor(
  sample_metadata_primary$Group,
  levels = c(
    "Adult",
    "Pediatric"
  )
)

sample_metadata_primary$Sex <- factor(
  sample_metadata_primary$Sex
)

sample_metadata_primary$Batch <- factor(
  sample_metadata_primary$Batch
)

if (nrow(sample_metadata_primary) == 0) {
  stop(
    "No eligible SARS-CoV-2-positive nasal/nasopharyngeal samples remained after cohort selection."
  )
}

# Sex and Batch are model covariates, not cohort-selection filters.
# Do not silently remove samples for missing covariates. Instead, document and
# stop if the selected cohort contains missing values because DESeq2 cannot fit
# a design containing NA covariates.

missing_covariates <- sample_metadata_primary[
  is.na(sample_metadata_primary$Sex) |
    is.na(sample_metadata_primary$Batch) |
    sample_metadata_primary$Batch == "",
  c(
    "GEO_ID",
    "Sample",
    "Age",
    "Sex",
    "Batch",
    "Group",
    "Tissue",
    "InfectionStatus"
  ),
  drop = FALSE
]

write.csv(
  missing_covariates,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_missing_model_covariates.csv")
  ),
  row.names = FALSE
)

if (nrow(missing_covariates) > 0) {
  stop(
    "The cohort defined by the original selection logic contains samples with ",
    "missing Sex and/or Batch values. No samples were removed automatically. ",
    "Inspect GSE231409_missing_model_covariates.csv before fitting DESeq2."
  )
}


# =============================================================================
# 07. Download raw counts if absent
# =============================================================================

counts_file <- file.path(
  data_dir,
  raw_counts_filename
)

if (!file.exists(counts_file)) {
  
  message("Raw-count file not found locally.")
  message("Downloading from NCBI GEO...")
  message("Destination: ", counts_file)
  
  download.file(
    url = raw_counts_url,
    destfile = counts_file,
    mode = "wb",
    method = "auto",
    quiet = FALSE
  )
  
  if (!file.exists(counts_file)) {
    stop(
      "The raw-count download did not create the expected file:\n",
      counts_file
    )
  }
  
  downloaded_size <- file.info(counts_file)$size
  
  if (is.na(downloaded_size) || downloaded_size == 0) {
    stop(
      "The downloaded raw-count file is empty:\n",
      counts_file
    )
  }
  
  message(
    "Download completed. File size: ",
    round(downloaded_size / 1024^2, 2),
    " MB"
  )
  
} else {
  
  message(
    "Raw-count file already exists. Download skipped: ",
    counts_file
  )
}


# =============================================================================
# 08. Load raw count matrix
# =============================================================================

raw_counts <- read.csv(
  counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

writeLines(
  colnames(raw_counts),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_raw_count_columns.txt")
  )
)

# The original analysis uses target_id as the gene identifier.
if ("target_id" %in% colnames(raw_counts)) {
  
  gene_id_column <- "target_id"
  
} else {
  
  gene_id_column <- colnames(raw_counts)[1]
  
  warning(
    "Expected identifier column 'target_id' was not found. ",
    "Using the first column as the gene identifier: ",
    gene_id_column
  )
}

raw_counts[[gene_id_column]] <- as.character(
  raw_counts[[gene_id_column]]
)

if (anyDuplicated(raw_counts[[gene_id_column]]) > 0) {
  
  duplicated_gene_ids <- unique(
    raw_counts[[gene_id_column]][
      duplicated(
        raw_counts[[gene_id_column]]
      )
    ]
  )
  
  writeLines(
    duplicated_gene_ids,
    file.path(
      diagnostics_dir,
      paste0(dataset_id, "_duplicated_raw_gene_ids.txt")
    )
  )
  
  stop(
    "Duplicated gene identifiers were found in the raw count matrix."
  )
}

rownames(raw_counts) <- raw_counts[[gene_id_column]]
raw_counts[[gene_id_column]] <- NULL

counts_matrix <- as.matrix(
  raw_counts
)

storage.mode(counts_matrix) <- "numeric"

if (any(is.na(counts_matrix))) {
  stop(
    "NA values were detected in the raw count matrix."
  )
}

if (any(counts_matrix < 0)) {
  stop(
    "Negative values were detected in the raw count matrix."
  )
}

if (
  any(
    abs(
      counts_matrix -
      round(counts_matrix)
    ) > 1e-8
  )
) {
  stop(
    "Non-integer values were detected. ",
    "DESeq2 requires raw integer gene-level counts."
  )
}

counts_matrix <- round(
  counts_matrix
)

message(
  "Raw count matrix loaded successfully: ",
  nrow(counts_matrix),
  " genes x ",
  ncol(counts_matrix),
  " samples"
)


# =============================================================================
# 09. Match metadata and counts
# =============================================================================

samples_in_both <- intersect(
  sample_metadata_primary$Sample,
  colnames(counts_matrix)
)

missing_from_counts <- setdiff(
  sample_metadata_primary$Sample,
  colnames(counts_matrix)
)

count_columns_not_used <- setdiff(
  colnames(counts_matrix),
  sample_metadata_primary$Sample
)

writeLines(
  missing_from_counts,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_samples_missing_from_counts.txt")
  )
)

writeLines(
  count_columns_not_used,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_count_columns_not_used.txt")
  )
)

if (length(samples_in_both) == 0) {
  stop(
    "No sample identifiers matched between GEO metadata and the raw count matrix."
  )
}

sample_metadata_primary <- sample_metadata_primary[
  sample_metadata_primary$Sample %in% samples_in_both,
  ,
  drop = FALSE
]

sample_metadata_primary <- sample_metadata_primary[
  match(
    samples_in_both,
    sample_metadata_primary$Sample
  ),
  ,
  drop = FALSE
]

counts_matrix <- counts_matrix[
  ,
  sample_metadata_primary$Sample,
  drop = FALSE
]

rownames(sample_metadata_primary) <- sample_metadata_primary$Sample

stopifnot(
  identical(
    colnames(counts_matrix),
    rownames(sample_metadata_primary)
  )
)


# =============================================================================
# 10. Cohort composition and batch structure
# =============================================================================

cohort_summary <- sample_metadata_primary %>%
  count(
    Group,
    Sex,
    name = "N"
  )

age_summary <- sample_metadata_primary %>%
  group_by(
    Group
  ) %>%
  summarise(
    N = n(),
    Age_min = min(Age, na.rm = TRUE),
    Age_median = median(Age, na.rm = TRUE),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_max = max(Age, na.rm = TRUE),
    .groups = "drop"
  )

batch_summary <- sample_metadata_primary %>%
  count(
    Batch,
    Group,
    Sex,
    name = "N"
  )

write.csv(
  sample_metadata_primary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_sample_metadata.csv")
  ),
  row.names = FALSE
)

write.csv(
  cohort_summary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_cohort_summary.csv")
  ),
  row.names = FALSE
)

write.csv(
  age_summary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_age_summary.csv")
  ),
  row.names = FALSE
)

write.csv(
  batch_summary,
  file.path(
    metadata_dir,
    paste0(dataset_id, "_batch_summary.csv")
  ),
  row.names = FALSE
)

print(cohort_summary)
print(age_summary)
print(batch_summary)

group_sex_table <- table(
  sample_metadata_primary$Group,
  sample_metadata_primary$Sex
)

group_batch_table <- table(
  sample_metadata_primary$Group,
  sample_metadata_primary$Batch
)

capture.output(
  group_sex_table,
  file = file.path(
    diagnostics_dir,
    paste0(dataset_id, "_group_by_sex_table.txt")
  )
)

capture.output(
  group_batch_table,
  file = file.path(
    diagnostics_dir,
    paste0(dataset_id, "_group_by_batch_table.txt")
  )
)


# =============================================================================
# 11. Primary low-count filtering
# =============================================================================

group_sizes <- table(
  sample_metadata_primary$Group
)

if (length(group_sizes) != 2) {
  stop(
    "Both Adult and Pediatric groups must be represented."
  )
}

min_group_size <- min(
  group_sizes
)

keep_primary <- rowSums(
  counts_matrix >= MIN_COUNT
) >= min_group_size

counts_primary <- counts_matrix[
  keep_primary,
  ,
  drop = FALSE
]


# =============================================================================
# 12. CPM-based sensitivity filter
# =============================================================================

cpm_matrix <- edgeR::cpm(
  counts_matrix
)

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
  Adult_samples = sum(
    sample_metadata_primary$Group == "Adult"
  ),
  Pediatric_samples = sum(
    sample_metadata_primary$Group == "Pediatric"
  ),
  Minimum_group_size = min_group_size,
  Genes_before_filtering = nrow(counts_matrix),
  Primary_MIN_COUNT = MIN_COUNT,
  Primary_genes_after_filtering = nrow(counts_primary),
  Primary_genes_removed =
    nrow(counts_matrix) - nrow(counts_primary),
  Sensitivity_MIN_CPM = MIN_CPM,
  CPM_genes_after_filtering = nrow(counts_cpm),
  CPM_genes_removed =
    nrow(counts_matrix) - nrow(counts_cpm)
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


# =============================================================================
# 13. Determine primary design
# =============================================================================

sex_adjustment_feasible <- (
  nlevels(
    droplevels(
      sample_metadata_primary$Sex
    )
  ) >= 2
)

batch_adjustment_feasible <- (
  nlevels(
    droplevels(
      sample_metadata_primary$Batch
    )
  ) >= 2
)

if (
  batch_adjustment_feasible &&
  sex_adjustment_feasible
) {
  
  primary_design <- ~ Batch + Sex + Group
  
} else if (
  batch_adjustment_feasible &&
  !sex_adjustment_feasible
) {
  
  primary_design <- ~ Batch + Group
  
  warning(
    "Sex adjustment is not feasible. Primary model will use ~ Batch + Group."
  )
  
} else if (
  !batch_adjustment_feasible &&
  sex_adjustment_feasible
) {
  
  primary_design <- ~ Sex + Group
  
  warning(
    "Batch adjustment is not feasible. Primary model will use ~ Sex + Group."
  )
  
} else {
  
  primary_design <- ~ Group
  
  warning(
    "Neither Batch nor Sex adjustment is feasible. Primary model will use ~ Group."
  )
}

primary_design_matrix <- model.matrix(
  primary_design,
  data = sample_metadata_primary
)

primary_design_full_rank <- (
  qr(primary_design_matrix)$rank ==
    ncol(primary_design_matrix)
)

writeLines(
  c(
    paste(
      "Primary design:",
      paste(
        deparse(primary_design),
        collapse = ""
      )
    ),
    paste(
      "Design matrix rank:",
      qr(primary_design_matrix)$rank
    ),
    paste(
      "Design matrix columns:",
      ncol(primary_design_matrix)
    ),
    paste(
      "Full rank:",
      primary_design_full_rank
    )
  ),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_design_rank_check.txt")
  )
)

if (!primary_design_full_rank) {
  stop(
    "Primary DESeq2 design matrix is not full rank. ",
    "Inspect Group x Batch and Group x Sex diagnostics before proceeding."
  )
}


# =============================================================================
# 14. Primary DESeq2 model
# =============================================================================

dds_primary <- DESeqDataSetFromMatrix(
  countData = counts_primary,
  colData = sample_metadata_primary,
  design = primary_design
)

dds_primary <- DESeq(
  dds_primary,
  test = "Wald"
)

print(
  resultsNames(
    dds_primary
  )
)

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

res_primary_df <- as.data.frame(
  res_primary
)

res_primary_df$ENSEMBL <- rownames(
  res_primary_df
)

rownames(
  res_primary_df
) <- NULL

res_primary_df <- res_primary_df %>%
  relocate(
    ENSEMBL
  )

res_primary_df$ENSEMBL_BASE <- sub(
  "\\..*$",
  "",
  res_primary_df$ENSEMBL
)


# =============================================================================
# 15. Annotation diagnostics
# =============================================================================

annotation_list <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = as.character(
    res_primary_df$ENSEMBL_BASE
  ),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "CharacterList"
)

annotation_diagnostics <- data.frame(
  ENSEMBL_BASE = names(annotation_list),
  N_SYMBOLS = lengths(annotation_list),
  SYMBOLS = vapply(
    annotation_list,
    paste,
    character(1),
    collapse = ";"
  ),
  stringsAsFactors = FALSE
)

annotation_diagnostics$SYMBOLS[
  annotation_diagnostics$SYMBOLS %in% c(
    "",
    "NA",
    "NaN"
  )
] <- NA_character_

annotation_diagnostics$AnnotationStatus <- ifelse(
  is.na(
    annotation_diagnostics$SYMBOLS
  ),
  "No_symbol",
  ifelse(
    annotation_diagnostics$N_SYMBOLS > 1 |
      grepl(
        ";",
        annotation_diagnostics$SYMBOLS,
        fixed = TRUE
      ),
    "Ambiguous_symbol",
    "Unique_symbol"
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
  annotation_diagnostics$AnnotationStatus ==
    "Unique_symbol",
  c(
    "ENSEMBL_BASE",
    "SYMBOLS"
  ),
  drop = FALSE
]

colnames(
  unique_annotation
)[2] <- "SYMBOL"

res_primary_annotated <- merge(
  res_primary_df,
  unique_annotation,
  by = "ENSEMBL_BASE",
  all.x = TRUE,
  sort = FALSE
)

res_primary_annotated <- res_primary_annotated[
  match(
    res_primary_df$ENSEMBL,
    res_primary_annotated$ENSEMBL
  ),
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
    name = "N_ENSEMBL_IDS"
  ) %>%
  filter(
    N_ENSEMBL_IDS > 1
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
  !is.na(
    res_primary_annotated$SYMBOL
  ) &
    res_primary_annotated$SYMBOL != "" &
    !(
      res_primary_annotated$SYMBOL %in%
        duplicate_symbols$SYMBOL
    )
)


# =============================================================================
# 16. Export complete primary DESeq2 results
# =============================================================================

write.csv(
  res_primary_annotated,
  file.path(
    primary_dir,
    paste0(dataset_id, "_DESeq2_complete.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 17. Export meta-analysis input
# =============================================================================

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
    Original_ID = ENSEMBL,
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


# =============================================================================
# 18. Shrinkage for visualization and cohort-level reporting
# =============================================================================

expected_coef <- "Group_Pediatric_vs_Adult"

if (
  expected_coef %in%
  resultsNames(dds_primary)
) {
  
  res_shrunk <- lfcShrink(
    dds_primary,
    coef = expected_coef,
    type = "apeglm"
  )
  
  res_shrunk_df <- as.data.frame(
    res_shrunk
  )
  
  res_shrunk_df$ENSEMBL <- rownames(
    res_shrunk_df
  )
  
  rownames(
    res_shrunk_df
  ) <- NULL
  
  res_shrunk_df$ENSEMBL_BASE <- sub(
    "\\..*$",
    "",
    res_shrunk_df$ENSEMBL
  )
  
  res_shrunk_df <- merge(
    res_shrunk_df,
    unique_annotation,
    by = "ENSEMBL_BASE",
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
  
} else {
  
  warning(
    "Expected coefficient not found for apeglm shrinkage: ",
    expected_coef
  )
}


# =============================================================================
# 19. Normalized counts
# =============================================================================

normalized_counts <- counts(
  dds_primary,
  normalized = TRUE
)

normalized_counts_df <- as.data.frame(
  normalized_counts
)

normalized_counts_df$ENSEMBL <- rownames(
  normalized_counts_df
)

rownames(
  normalized_counts_df
) <- NULL

normalized_counts_df <- normalized_counts_df %>%
  relocate(
    ENSEMBL
  )

write.csv(
  normalized_counts_df,
  file.path(
    primary_dir,
    paste0(dataset_id, "_normalized_counts.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 20. Export matrices and ranked statistics for pathway analysis
# =============================================================================

vsd_export <- vst(
  dds_primary,
  blind = TRUE
)

vst_matrix <- as.data.frame(
  assay(vsd_export)
)

vst_matrix$ENSEMBL <- rownames(
  vst_matrix
)

rownames(
  vst_matrix
) <- NULL

vst_matrix <- vst_matrix %>%
  relocate(
    ENSEMBL
  )

write.csv(
  vst_matrix,
  file.path(
    primary_dir,
    paste0(dataset_id, "_VST_matrix.csv")
  ),
  row.names = FALSE
)

gsea_ranked <- res_primary_annotated %>%
  filter(
    MetaEligibleAnnotation,
    !is.na(stat),
    is.finite(stat)
  ) %>%
  transmute(
    Dataset = dataset_id,
    Original_ID = ENSEMBL,
    SYMBOL,
    ranking_statistic = stat,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  arrange(
    desc(
      ranking_statistic
    )
  )

write.csv(
  gsea_ranked,
  file.path(
    primary_dir,
    paste0(dataset_id, "_GSEA_ranked_genes.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 21. Size factors
# =============================================================================

size_factor_table <- data.frame(
  Sample = names(
    sizeFactors(
      dds_primary
    )
  ),
  SizeFactor = sizeFactors(
    dds_primary
  )
)

write.csv(
  size_factor_table,
  file.path(
    primary_dir,
    paste0(dataset_id, "_size_factors.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 22. Dispersion diagnostics
# =============================================================================

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

plotDispEsts(
  dds_primary
)

dev.off()

dispersion_fit_type <- attr(
  dispersionFunction(
    dds_primary
  ),
  "fitType"
)

if (
  is.null(
    dispersion_fit_type
  )
) {
  dispersion_fit_type <- NA_character_
}


# =============================================================================
# 23. Independent filtering
# =============================================================================

res_metadata <- S4Vectors::metadata(
  res_primary
)

filter_threshold <- res_metadata$filterThreshold

if (
  !is.null(
    filter_threshold
  )
) {
  independent_filter_threshold <- as.numeric(
    filter_threshold
  )
} else {
  independent_filter_threshold <- NA_real_
}


# =============================================================================
# 24. Cook's distance diagnostics
# =============================================================================

cooks_matrix <- assays(
  dds_primary
)[["cooks"]]

if (
  !is.null(
    cooks_matrix
  )
) {
  
  max_cooks <- apply(
    cooks_matrix,
    1,
    max,
    na.rm = TRUE
  )
  
  cooks_ensembl_ids <- rownames(
    dds_primary
  )
  
  if (
    is.null(
      cooks_ensembl_ids
    ) ||
    length(
      cooks_ensembl_ids
    ) != length(
      max_cooks
    )
  ) {
    stop(
      "Could not match Cook's-distance values to DESeq2 gene identifiers."
    )
  }
  
  cooks_diagnostics <- data.frame(
    ENSEMBL = cooks_ensembl_ids,
    MaxCooksDistance = as.numeric(
      max_cooks
    ),
    stringsAsFactors = FALSE
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


# =============================================================================
# 25. PCA
# =============================================================================

rld <- vst(
  dds_primary,
  blind = TRUE
)

pcaData <- plotPCA(
  rld,
  intgroup = c(
    "Group"
  ),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(
    pcaData,
    "percentVar"
  )
)

tiff(
  file.path(
    figures_dir,
    "PCA_R_GSE231409_Hurst_Kelly_Eng.tiff"
  ),
  res = 600,
  width = 3700,
  height = 2000,
  compression = "lzw"
)

ggplot(
  pcaData,
  aes(
    PC1,
    PC2
  )
) +
  geom_point(
    aes(
      colour = Group
    ),
    size = 3
  ) +
  ggtitle(
    "GSE231409"
  ) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.4,
      face = "bold"
    ),
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(
      size = 13
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  ) +
  scale_color_manual(
    values = c(
      "#3E76BC",
      "#FCCE24"
    ),
    labels = c(
      "SARS-CoV-2 + Adult",
      "SARS-CoV-2 + Pediatric"
    )
  ) +
  labs(
    color = "Status"
  )

dev.off()


# Quantitative assessment of age, sex, and batch associations with PC1/PC2.

pca_scores <- as.data.frame(
  pcaData
)

pca_scores$Sample <- rownames(
  pcaData
)

rownames(
  pca_scores
) <- NULL

pca_scores <- merge(
  pca_scores,
  sample_metadata_primary %>%
    dplyr::select(
      Sample,
      Age,
      Sex,
      Batch
    ),
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)

if (
  !all(
    c(
      "Group",
      "Age",
      "Sex",
      "Batch"
    ) %in%
    colnames(
      pca_scores
    )
  )
) {
  stop(
    "PCA metadata merge failed: Group, Age, Sex, and Batch must be present."
  )
}

write.csv(
  pca_scores,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_PCA_scores.csv")
  ),
  row.names = FALSE
)

if (
  batch_adjustment_feasible &&
  sex_adjustment_feasible
) {
  
  pca_pc1_model <- lm(
    PC1 ~ Age + Sex + Batch,
    data = pca_scores
  )
  
  pca_pc2_model <- lm(
    PC2 ~ Age + Sex + Batch,
    data = pca_scores
  )
  
} else if (
  batch_adjustment_feasible
) {
  
  pca_pc1_model <- lm(
    PC1 ~ Age + Batch,
    data = pca_scores
  )
  
  pca_pc2_model <- lm(
    PC2 ~ Age + Batch,
    data = pca_scores
  )
  
} else if (
  sex_adjustment_feasible
) {
  
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
    capture.output(
      summary(
        pca_pc1_model
      )
    ),
    "",
    "===== PC2 =====",
    capture.output(
      summary(
        pca_pc2_model
      )
    )
  ),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_PCA_age_sex_batch_models.txt")
  )
)


# =============================================================================
# 26. Expression and PCA diagnostics
# =============================================================================

expression_prevalence <- data.frame(
  ENSEMBL = rownames(
    counts_matrix
  ),
  N_samples_count_ge_10 = rowSums(
    counts_matrix >= MIN_COUNT
  ),
  Proportion_samples_count_ge_10 = rowMeans(
    counts_matrix >= MIN_COUNT
  ),
  stringsAsFactors = FALSE
)

write.csv(
  expression_prevalence,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_expression_prevalence.csv")
  ),
  row.names = FALSE
)

tiff(
  file.path(
    sensitivity_figures_dir,
    "Expression_prevalence_GSE231409.tiff"
  ),
  res = 600,
  width = 3700,
  height = 2000,
  compression = "lzw"
)

ggplot(
  expression_prevalence,
  aes(
    x = N_samples_count_ge_10
  )
) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    closed = "left",
    fill = "#3E76BC",
    colour = "black"
  ) +
  geom_vline(
    xintercept = min_group_size,
    linetype = 2,
    linewidth = 0.8
  ) +
  ggtitle(
    "GSE231409"
  ) +
  xlab(
    "Number of samples with raw count >= 10"
  ) +
  ylab(
    "Number of genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.4,
      face = "bold"
    ),
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )

dev.off()


raw_count_long <- as.data.frame(
  counts_matrix
)

raw_count_long$ENSEMBL <- rownames(
  raw_count_long
)

raw_count_long <- raw_count_long %>%
  pivot_longer(
    cols = -ENSEMBL,
    names_to = "Sample",
    values_to = "RawCount"
  ) %>%
  left_join(
    sample_metadata_primary %>%
      dplyr::select(
        Sample,
        Group
      ),
    by = "Sample"
  )

write.csv(
  raw_count_long,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_raw_count_distribution_long.csv")
  ),
  row.names = FALSE
)

tiff(
  file.path(
    sensitivity_figures_dir,
    "Raw_count_distribution_GSE231409.tiff"
  ),
  res = 600,
  width = 3700,
  height = 2000,
  compression = "lzw"
)

ggplot(
  raw_count_long,
  aes(
    x = log10(
      RawCount + 1
    ),
    colour = Group
  )
) +
  geom_density(
    linewidth = 1,
    adjust = 1
  ) +
  ggtitle(
    "GSE231409"
  ) +
  xlab(
    "log10(raw count + 1)"
  ) +
  ylab(
    "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.4,
      face = "bold"
    ),
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(
      size = 13
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  ) +
  scale_color_manual(
    values = c(
      "#3E76BC",
      "#FCCE24"
    ),
    labels = c(
      "SARS-CoV-2 + Adult",
      "SARS-CoV-2 + Pediatric"
    )
  ) +
  labs(
    color = "Status"
  )

dev.off()


# =============================================================================
# 27. Sensitivity: continuous age
# =============================================================================

metadata_age <- sample_metadata_primary

metadata_age$AgeScaled <- as.numeric(
  scale(
    metadata_age$Age
  )
)

if (
  batch_adjustment_feasible &&
  sex_adjustment_feasible
) {
  
  age_design <- ~ Batch + Sex + AgeScaled
  
} else if (
  batch_adjustment_feasible
) {
  
  age_design <- ~ Batch + AgeScaled
  
} else if (
  sex_adjustment_feasible
) {
  
  age_design <- ~ Sex + AgeScaled
  
} else {
  
  age_design <- ~ AgeScaled
}

age_design_matrix <- model.matrix(
  age_design,
  data = metadata_age
)

age_design_full_rank <- (
  qr(age_design_matrix)$rank ==
    ncol(age_design_matrix)
)

if (age_design_full_rank) {
  
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
    resultsNames(
      dds_age
    ),
    value = TRUE
  )
  
  if (
    length(
      age_coef
    ) == 1
  ) {
    
    res_age <- results(
      dds_age,
      name = age_coef,
      alpha = ALPHA,
      independentFiltering = TRUE,
      cooksCutoff = TRUE
    )
    
    res_age_df <- as.data.frame(
      res_age
    )
    
    res_age_df$ENSEMBL <- rownames(
      res_age_df
    )
    
    rownames(
      res_age_df
    ) <- NULL
    
    res_age_df$ENSEMBL_BASE <- sub(
      "\\..*$",
      "",
      res_age_df$ENSEMBL
    )
    
    res_age_df <- merge(
      res_age_df,
      unique_annotation,
      by = "ENSEMBL_BASE",
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
  
} else {
  
  writeLines(
    "Continuous-age sensitivity was not estimable because the design matrix was not full rank.",
    file.path(
      sensitivity_dir,
      paste0(dataset_id, "_continuous_age_NOT_ESTIMABLE.txt")
    )
  )
}


# =============================================================================
# 28. Sensitivity: CPM filter
# =============================================================================

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

res_cpm_df <- as.data.frame(
  res_cpm
)

res_cpm_df$ENSEMBL <- rownames(
  res_cpm_df
)

rownames(
  res_cpm_df
) <- NULL

res_cpm_df$ENSEMBL_BASE <- sub(
  "\\..*$",
  "",
  res_cpm_df$ENSEMBL
)

res_cpm_df <- merge(
  res_cpm_df,
  unique_annotation,
  by = "ENSEMBL_BASE",
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

low_count_comparison <- res_primary_annotated %>%
  dplyr::select(
    ENSEMBL,
    SYMBOL,
    log2FC_primary = log2FoldChange
  ) %>%
  inner_join(
    res_cpm_df %>%
      dplyr::select(
        ENSEMBL,
        log2FC_CPM = log2FoldChange
      ),
    by = "ENSEMBL"
  ) %>%
  filter(
    !is.na(
      log2FC_primary
    ),
    !is.na(
      log2FC_CPM
    ),
    is.finite(
      log2FC_primary
    ),
    is.finite(
      log2FC_CPM
    )
  )

write.csv(
  low_count_comparison,
  file.path(
    sensitivity_dir,
    paste0(dataset_id, "_primary_vs_CPM_filter_comparison.csv")
  ),
  row.names = FALSE
)

tiff(
  file.path(
    sensitivity_figures_dir,
    "Primary_vs_CPM_filter_GSE231409.tiff"
  ),
  res = 600,
  width = 3700,
  height = 2000,
  compression = "lzw"
)

ggplot(
  low_count_comparison,
  aes(
    x = log2FC_primary,
    y = log2FC_CPM
  )
) +
  geom_point(
    colour = "#3E76BC",
    alpha = 0.45,
    size = 1.4
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    colour = "black"
  ) +
  ggtitle(
    "GSE231409"
  ) +
  xlab(
    "Primary filter: log2 fold change"
  ) +
  ylab(
    "CPM sensitivity filter: log2 fold change"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.4,
      face = "bold"
    ),
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )

dev.off()


# =============================================================================
# 29. Sensitivity: Group x Sex interaction
# =============================================================================

interaction_feasible <- (
  nlevels(
    droplevels(
      sample_metadata_primary$Sex
    )
  ) == 2 &&
    nlevels(
      droplevels(
        sample_metadata_primary$Group
      )
    ) == 2 &&
    all(
      group_sex_table > 0
    )
)

if (
  interaction_feasible
) {
  
  if (
    batch_adjustment_feasible
  ) {
    interaction_design <- ~ Batch + Sex + Group + Sex:Group
  } else {
    interaction_design <- ~ Sex + Group + Sex:Group
  }
  
  interaction_design_matrix <- model.matrix(
    interaction_design,
    data = sample_metadata_primary
  )
  
  interaction_design_full_rank <- (
    qr(interaction_design_matrix)$rank ==
      ncol(interaction_design_matrix)
  )
  
  if (
    interaction_design_full_rank
  ) {
    
    dds_interaction <- DESeqDataSetFromMatrix(
      countData = counts_primary,
      colData = sample_metadata_primary,
      design = interaction_design
    )
    
    dds_interaction <- DESeq(
      dds_interaction,
      test = "Wald"
    )
    
    interaction_coef <- grep(
      "Group.*Sex|Sex.*Group",
      resultsNames(
        dds_interaction
      ),
      value = TRUE
    )
    
    if (
      length(
        interaction_coef
      ) == 1
    ) {
      
      res_interaction <- results(
        dds_interaction,
        name = interaction_coef,
        alpha = ALPHA,
        independentFiltering = TRUE,
        cooksCutoff = TRUE
      )
      
      res_interaction_df <- as.data.frame(
        res_interaction
      )
      
      res_interaction_df$ENSEMBL <- rownames(
        res_interaction_df
      )
      
      rownames(
        res_interaction_df
      ) <- NULL
      
      res_interaction_df$ENSEMBL_BASE <- sub(
        "\\..*$",
        "",
        res_interaction_df$ENSEMBL
      )
      
      res_interaction_df <- merge(
        res_interaction_df,
        unique_annotation,
        by = "ENSEMBL_BASE",
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
      "Group x Sex interaction was not estimable because the design matrix was not full rank.",
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


# =============================================================================
# 30. Sensitivity summary
# =============================================================================

primary_significant_ids <- res_primary_annotated$ENSEMBL[
  !is.na(
    res_primary_annotated$padj
  ) &
    res_primary_annotated$padj < ALPHA
]

cpm_significant_ids <- res_cpm_df$ENSEMBL[
  !is.na(
    res_cpm_df$padj
  ) &
    res_cpm_df$padj < ALPHA
]

primary_cpm_correlation <- cor(
  low_count_comparison$log2FC_primary,
  low_count_comparison$log2FC_CPM,
  use = "complete.obs",
  method = "pearson"
)

continuous_age_significant <- if (
  exists(
    "res_age_df"
  )
) {
  sum(
    !is.na(
      res_age_df$padj
    ) &
      res_age_df$padj < ALPHA
  )
} else {
  NA_integer_
}

interaction_significant <- if (
  exists(
    "res_interaction_df"
  )
) {
  sum(
    !is.na(
      res_interaction_df$padj
    ) &
      res_interaction_df$padj < ALPHA
  )
} else {
  NA_integer_
}

sensitivity_summary <- data.frame(
  Dataset = dataset_id,
  Primary_genes_tested = nrow(
    res_primary_annotated
  ),
  Primary_FDR_lt_0_05 = length(
    primary_significant_ids
  ),
  CPM_genes_tested = nrow(
    res_cpm_df
  ),
  CPM_FDR_lt_0_05 = length(
    cpm_significant_ids
  ),
  Primary_CPM_shared_significant = length(
    intersect(
      primary_significant_ids,
      cpm_significant_ids
    )
  ),
  Primary_CPM_log2FC_Pearson_r =
    primary_cpm_correlation,
  Continuous_age_FDR_lt_0_05 =
    continuous_age_significant,
  Group_by_sex_interaction_FDR_lt_0_05 =
    interaction_significant,
  stringsAsFactors = FALSE
)

write.csv(
  sensitivity_summary,
  file.path(
    sensitivity_dir,
    paste0(dataset_id, "_sensitivity_summary.csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 31. Workflow summary
# =============================================================================

workflow_summary <- c(
  
  paste(
    "Dataset:",
    dataset_id
  ),
  
  paste(
    "Project root:",
    project_root
  ),
  
  paste(
    "Population:",
    "SARS-CoV-2-positive nasal/nasopharyngeal swab samples selected using the original analysis logic"
  ),
  
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
    nrow(
      sample_metadata_primary
    )
  ),
  
  paste(
    "Adult samples:",
    sum(
      sample_metadata_primary$Group == "Adult"
    )
  ),
  
  paste(
    "Pediatric samples:",
    sum(
      sample_metadata_primary$Group == "Pediatric"
    )
  ),
  
  paste(
    "Batch levels:",
    nlevels(
      droplevels(
        sample_metadata_primary$Batch
      )
    )
  ),
  
  paste(
    "Genes before low-count filtering:",
    nrow(
      counts_matrix
    )
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
    nrow(
      counts_primary
    )
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
    nrow(
      counts_cpm
    )
  ),
  
  paste(
    "Primary design formula:",
    paste(
      deparse(
        primary_design
      ),
      collapse = ""
    )
  ),
  
  paste(
    "Primary design full rank:",
    primary_design_full_rank
  ),
  
  paste(
    "Batch adjustment feasible:",
    batch_adjustment_feasible
  ),
  
  paste(
    "Sex adjustment feasible:",
    sex_adjustment_feasible
  ),
  
  paste(
    "Cohort selection:",
    "Sex and Batch were not used to exclude samples; both were retained only as model covariates"
  ),
  
  paste(
    "Batch handling:",
    "included as a DESeq2 design covariate; no pre-correction of raw counts"
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
    "apeglm for cohort-level visualization/reporting only"
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
    "Unique ENSEMBL-to-SYMBOL mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "Unique_symbol",
      na.rm = TRUE
    )
  ),
  
  paste(
    "Ambiguous ENSEMBL-to-SYMBOL mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "Ambiguous_symbol",
      na.rm = TRUE
    )
  ),
  
  paste(
    "Missing ENSEMBL-to-SYMBOL mappings:",
    sum(
      annotation_diagnostics$AnnotationStatus ==
        "No_symbol",
      na.rm = TRUE
    )
  ),
  
  paste(
    "Genes eligible for later SYMBOL-based harmonization:",
    nrow(
      meta_input
    )
  ),
  
  paste(
    "GSEA ranking:",
    "complete DESeq2 Wald statistic among genes with valid unique SYMBOL mappings"
  ),
  
  paste(
    "GSEA ranked genes exported:",
    nrow(
      gsea_ranked
    )
  ),
  
  paste(
    "VST matrix exported:",
    paste0(
      nrow(
        vst_matrix
      ),
      " genes x ",
      ncol(
        vst_matrix
      ) - 1,
      " samples"
    )
  ),
  
  paste(
    "Primary-vs-CPM log2FC Pearson correlation:",
    round(
      primary_cpm_correlation,
      6
    )
  ),
  
  paste(
    "Primary significant genes retained under CPM sensitivity:",
    paste0(
      length(
        intersect(
          primary_significant_ids,
          cpm_significant_ids
        )
      ),
      " / ",
      length(
        primary_significant_ids
      )
    )
  )
)

writeLines(
  workflow_summary,
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_workflow_summary.txt")
  )
)


# =============================================================================
# 32. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    paste0(dataset_id, "_sessionInfo.txt")
  )
)

writeLines(
  c(
    paste(
      "R version:",
      R.version.string
    ),
    paste(
      "DESeq2 version:",
      as.character(
        packageVersion(
          "DESeq2"
        )
      )
    ),
    paste(
      "GEOquery version:",
      as.character(
        packageVersion(
          "GEOquery"
        )
      )
    ),
    paste(
      "edgeR version:",
      as.character(
        packageVersion(
          "edgeR"
        )
      )
    ),
    paste(
      "apeglm version:",
      as.character(
        packageVersion(
          "apeglm"
        )
      )
    ),
    paste(
      "AnnotationDbi version:",
      as.character(
        packageVersion(
          "AnnotationDbi"
        )
      )
    ),
    paste(
      "org.Hs.eg.db version:",
      as.character(
        packageVersion(
          "org.Hs.eg.db"
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
    paste0(dataset_id, "_package_versions.txt")
  )
)


# =============================================================================
# 33. Final summary
# =============================================================================

message("")
message("============================================================")
message(
  dataset_id,
  " analysis completed successfully."
)
message("============================================================")

message(
  "Project root: ",
  project_root
)

message(
  "Primary samples: ",
  nrow(
    sample_metadata_primary
  )
)

message(
  "Adult / Pediatric: ",
  sum(
    sample_metadata_primary$Group == "Adult"
  ),
  " / ",
  sum(
    sample_metadata_primary$Group == "Pediatric"
  )
)

message(
  "Batch levels: ",
  nlevels(
    droplevels(
      sample_metadata_primary$Batch
    )
  )
)

message(
  "Genes before filtering: ",
  nrow(
    counts_matrix
  )
)

message(
  "Genes after primary filtering: ",
  nrow(
    counts_primary
  )
)

message(
  "Genes eligible for later meta-analysis harmonization: ",
  nrow(
    meta_input
  )
)

message(
  "Results saved to: ",
  results_root
)

message("============================================================")
