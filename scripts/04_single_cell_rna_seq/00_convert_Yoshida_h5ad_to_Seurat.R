# =============================================================================
# 00_convert_Yoshida_h5ad_to_Seurat.R
#
# Convert the Yoshida et al. CELLxGENE airway H5AD object to a Seurat object.
#
# This script performs ONLY:
#   1. download of the H5AD file;
#   2. inspection of AnnData / SingleCellExperiment structure;
#   3. extraction of raw counts;
#   4. conversion to a Seurat object;
#   5. preservation of CELLxGENE cell metadata;
#   6. preservation of available dimensional reductions when compatible;
#   7. export of the Seurat object as RDS.
#
# No filtering, normalization, clustering, differential expression,
# pseudobulk analysis, or biological interpretation is performed here.
#
# Dataset:
#   https://datasets.cellxgene.cziscience.com/
#   f9efb73e-f116-46b5-a775-d4233e758024.h5ad
#
# Run from the root of the COVID-Project RStudio Project.
#
# No custom functions are defined in this script.
# =============================================================================

SCRIPT_BUILD <- "YOSHIDA_H5AD_TO_SEURAT_2026-08-17_v1"
message("Running script build: ", SCRIPT_BUILD)


# =============================================================================
# 00. Packages
# =============================================================================

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "here"
)

bioc_packages <- c(
  "zellkonverter",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {

  stop(
    "Missing CRAN packages: ",
    paste(
      missing_cran,
      collapse = ", "
    ),
    "\nInstall them before running the script."
  )
}

if (length(missing_bioc) > 0) {

  stop(
    "Missing Bioconductor packages: ",
    paste(
      missing_bioc,
      collapse = ", "
    ),
    "\nInstall them with BiocManager::install()."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(here)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})


# =============================================================================
# 01. Paths
# =============================================================================

dataset_url <- paste0(
  "https://datasets.cellxgene.cziscience.com/",
  "f9efb73e-f116-46b5-a775-d4233e758024.h5ad"
)

raw_dir <- here::here(
  "data",
  "raw",
  "single_cell",
  "Yoshida_airway"
)

processed_dir <- here::here(
  "data",
  "processed",
  "single_cell",
  "Yoshida_airway"
)

diagnostics_dir <- here::here(
  "results",
  "single_cell",
  "Yoshida_airway",
  "conversion_diagnostics"
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostics_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

h5ad_file <- file.path(
  raw_dir,
  "Yoshida_airway_CELLxGENE.h5ad"
)

seurat_file <- file.path(
  processed_dir,
  "Yoshida_airway_Seurat.rds"
)


# =============================================================================
# 02. Download H5AD file
# =============================================================================

if (!file.exists(h5ad_file)) {

  message("Downloading H5AD object...")

  download.file(
    url = dataset_url,
    destfile = h5ad_file,
    mode = "wb",
    method = "libcurl"
  )
}

if (!file.exists(h5ad_file)) {

  stop(
    "H5AD download failed:\n",
    h5ad_file
  )
}

h5ad_info <- file.info(
  h5ad_file
)

write.csv(
  data.frame(
    File = h5ad_file,
    Size_bytes = h5ad_info$size,
    Size_GB = h5ad_info$size / 1024^3,
    stringsAsFactors = FALSE
  ),
  file.path(
    diagnostics_dir,
    "H5AD_file_information.csv"
  ),
  row.names = FALSE
)

message(
  "H5AD size: ",
  round(
    h5ad_info$size / 1024^3,
    2
  ),
  " GB"
)


# =============================================================================
# 03. Read H5AD as SingleCellExperiment
# =============================================================================

# zellkonverter provides the bridge from AnnData/H5AD to R.
# The Python reader is used because it follows AnnData directly.
#
# use_hdf5 = FALSE is intentional for the conversion step because the raw
# expression matrix must ultimately become a Seurat counts layer.

message("Reading H5AD...")

sce <- zellkonverter::readH5AD(
  file = h5ad_file,
  use_hdf5 = FALSE,
  reader = "python"
)

message(
  "SingleCellExperiment dimensions: ",
  nrow(sce),
  " genes x ",
  ncol(sce),
  " cells"
)


# =============================================================================
# 04. Inspect imported object structure
# =============================================================================

writeLines(
  assayNames(sce),
  file.path(
    diagnostics_dir,
    "SCE_assay_names.txt"
  )
)

writeLines(
  altExpNames(sce),
  file.path(
    diagnostics_dir,
    "SCE_altExp_names.txt"
  )
)

writeLines(
  reducedDimNames(sce),
  file.path(
    diagnostics_dir,
    "SCE_reducedDim_names.txt"
  )
)

write.csv(
  data.frame(
    Metadata_column = colnames(
      colData(sce)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    diagnostics_dir,
    "CELLxGENE_metadata_columns.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    rowData(sce)
  ),
  file.path(
    diagnostics_dir,
    "SCE_feature_metadata.csv"
  ),
  row.names = TRUE
)


# =============================================================================
# 05. Extract cell metadata
# =============================================================================

cell_metadata <- as.data.frame(
  colData(sce)
)

rownames(
  cell_metadata
) <- colnames(
  sce
)

write.csv(
  cell_metadata,
  file.path(
    diagnostics_dir,
    "CELLxGENE_cell_metadata.csv"
  ),
  row.names = TRUE
)


# =============================================================================
# 06. Locate raw counts
# =============================================================================

raw_sce <- NULL
raw_assay_name <- NULL
raw_source <- NULL

# CELLxGENE H5AD files commonly store normalized X in the main AnnData matrix
# and raw counts in AnnData.raw. zellkonverter commonly imports AnnData.raw as
# an alternative experiment called "raw".

if (
  "raw" %in%
    altExpNames(
      sce
    )
) {

  raw_sce <- altExp(
    sce,
    "raw"
  )

  if (
    "X" %in%
      assayNames(
        raw_sce
      )
  ) {

    raw_assay_name <- "X"

  } else if (
    "counts" %in%
      assayNames(
        raw_sce
      )
  ) {

    raw_assay_name <- "counts"

  } else {

    raw_assay_name <- assayNames(
      raw_sce
    )[1]
  }

  raw_source <- "altExp(raw)"

} else if (
  "counts" %in%
    assayNames(
      sce
    )
) {

  raw_sce <- sce
  raw_assay_name <- "counts"
  raw_source <- "main SCE counts assay"

} else if (
  "raw" %in%
    assayNames(
      sce
    )
) {

  raw_sce <- sce
  raw_assay_name <- "raw"
  raw_source <- "main SCE raw assay"

} else {

  stop(
    "A raw-count matrix could not be identified automatically.\n",
    "Inspect:\n",
    file.path(
      diagnostics_dir,
      "SCE_assay_names.txt"
    ),
    "\n",
    file.path(
      diagnostics_dir,
      "SCE_altExp_names.txt"
    ),
    "\nDo not create a Seurat object from normalized X until raw counts are identified."
  )
}

raw_counts <- assay(
  raw_sce,
  raw_assay_name
)

message(
  "Raw count source: ",
  raw_source
)

message(
  "Raw count assay: ",
  raw_assay_name
)


# =============================================================================
# 07. Verify raw counts
# =============================================================================

check_gene_n <- min(
  250,
  nrow(
    raw_counts
  )
)

check_cell_n <- min(
  250,
  ncol(
    raw_counts
  )
)

raw_check <- as.matrix(
  raw_counts[
    seq_len(
      check_gene_n
    ),
    seq_len(
      check_cell_n
    ),
    drop = FALSE
  ]
)

raw_check_values <- as.numeric(
  raw_check
)

raw_check_values <- raw_check_values[
  is.finite(
    raw_check_values
  )
]

integer_like_fraction <- mean(
  abs(
    raw_check_values -
      round(
        raw_check_values
      )
  ) < 1e-8
)

minimum_raw_value <- min(
  raw_check_values
)

write.csv(
  data.frame(
    Raw_source = raw_source,
    Raw_assay = raw_assay_name,
    Integer_like_fraction = integer_like_fraction,
    Minimum_value = minimum_raw_value,
    stringsAsFactors = FALSE
  ),
  file.path(
    diagnostics_dir,
    "raw_count_validation.csv"
  ),
  row.names = FALSE
)

if (
  integer_like_fraction < 0.99 ||
    minimum_raw_value < 0
) {

  stop(
    "The selected matrix does not appear to contain raw non-negative integer counts.\n",
    "Inspect raw_count_validation.csv before conversion."
  )
}


# =============================================================================
# 08. Harmonize feature names for Seurat
# =============================================================================

raw_feature_metadata <- as.data.frame(
  rowData(
    raw_sce
  )
)

raw_feature_metadata$ORIGINAL_ROWNAME <- rownames(
  raw_sce
)

symbol_candidates <- c(
  "feature_name",
  "gene_symbol",
  "symbol",
  "gene_name"
)

available_symbol_columns <- symbol_candidates[
  symbol_candidates %in%
    colnames(
      raw_feature_metadata
    )
]

if (
  length(
    available_symbol_columns
  ) > 0
) {

  raw_feature_metadata$SEURAT_FEATURE <- as.character(
    raw_feature_metadata[
      [
        available_symbol_columns[1]
      ]
    ]
  )

} else {

  raw_feature_metadata$SEURAT_FEATURE <- raw_feature_metadata$ORIGINAL_ROWNAME
}

missing_feature_name <- (
  is.na(
    raw_feature_metadata$SEURAT_FEATURE
  ) |
    raw_feature_metadata$SEURAT_FEATURE == ""
)

raw_feature_metadata$SEURAT_FEATURE[
  missing_feature_name
] <- raw_feature_metadata$ORIGINAL_ROWNAME[
  missing_feature_name
]

raw_feature_metadata$SEURAT_FEATURE <- make.unique(
  raw_feature_metadata$SEURAT_FEATURE
)

rownames(
  raw_counts
) <- raw_feature_metadata$SEURAT_FEATURE

write.csv(
  raw_feature_metadata,
  file.path(
    diagnostics_dir,
    "Seurat_feature_name_mapping.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 09. Ensure cells match between counts and metadata
# =============================================================================

if (
  !identical(
    colnames(
      raw_counts
    ),
    rownames(
      cell_metadata
    )
  )
) {

  common_cells <- intersect(
    colnames(
      raw_counts
    ),
    rownames(
      cell_metadata
    )
  )

  if (
    length(
      common_cells
    ) == 0
  ) {

    stop(
      "Raw-count cell names and metadata cell names do not overlap."
    )
  }

  raw_counts <- raw_counts[
    ,
    common_cells,
    drop = FALSE
  ]

  cell_metadata <- cell_metadata[
    common_cells,
    ,
    drop = FALSE
  ]
}


# =============================================================================
# 10. Convert raw counts to sparse matrix if needed
# =============================================================================

if (
  !inherits(
    raw_counts,
    "sparseMatrix"
  )
) {

  raw_counts <- as(
    raw_counts,
    "dgCMatrix"
  )
}


# =============================================================================
# 11. Create Seurat object
# =============================================================================

message("Creating Seurat object...")

yoshida <- CreateSeuratObject(
  counts = raw_counts,
  assay = "RNA",
  project = "Yoshida_Airway_COVID19",
  meta.data = cell_metadata,
  min.cells = 0,
  min.features = 0
)

DefaultAssay(
  yoshida
) <- "RNA"

message(
  "Seurat object dimensions: ",
  nrow(yoshida),
  " genes x ",
  ncol(yoshida),
  " cells"
)


# =============================================================================
# 12. Preserve original feature identifiers
# =============================================================================

feature_metadata_seurat <- raw_feature_metadata[
  match(
    rownames(
      yoshida
    ),
    raw_feature_metadata$SEURAT_FEATURE
  ),
  ,
  drop = FALSE
]

rownames(
  feature_metadata_seurat
) <- rownames(
  yoshida
)

yoshida[
  [
    "RNA"
  ]
]@meta.features <- feature_metadata_seurat


# =============================================================================
# 13. Preserve imported dimensional reductions when available
# =============================================================================

available_reductions <- reducedDimNames(
  sce
)

if (
  length(
    available_reductions
  ) > 0
) {

  reduction_inventory <- data.frame(
    Original_reduction = available_reductions,
    Imported_to_Seurat = FALSE,
    Seurat_name = NA_character_,
    stringsAsFactors = FALSE
  )

  for (reduction_i in available_reductions) {

    embeddings_i <- reducedDim(
      sce,
      reduction_i
    )

    if (
      nrow(
        embeddings_i
      ) != ncol(
        sce
      )
    ) {

      next
    }

    embeddings_i <- as.matrix(
      embeddings_i
    )

    rownames(
      embeddings_i
    ) <- colnames(
      sce
    )

    embeddings_i <- embeddings_i[
      colnames(
        yoshida
      ),
      ,
      drop = FALSE
    ]

    seurat_reduction_name <- tolower(
      reduction_i
    )

    seurat_reduction_name <- gsub(
      "[^A-Za-z0-9]",
      "",
      seurat_reduction_name
    )

    if (
      seurat_reduction_name == ""
    ) {

      next
    }

    reduction_key <- paste0(
      toupper(
        seurat_reduction_name
      ),
      "_"
    )

    yoshida[
      [
        seurat_reduction_name
      ]
    ] <- CreateDimReducObject(
      embeddings = embeddings_i,
      key = reduction_key,
      assay = "RNA"
    )

    reduction_inventory$Imported_to_Seurat[
      reduction_inventory$Original_reduction ==
        reduction_i
    ] <- TRUE

    reduction_inventory$Seurat_name[
      reduction_inventory$Original_reduction ==
        reduction_i
    ] <- seurat_reduction_name
  }

  write.csv(
    reduction_inventory,
    file.path(
      diagnostics_dir,
      "Seurat_reduction_import_summary.csv"
    ),
    row.names = FALSE
  )
}


# =============================================================================
# 14. Record source information inside Seurat metadata
# =============================================================================

yoshida@misc$source <- list(
  study = "Yoshida et al.",
  dataset = "CELLxGENE airway",
  h5ad_url = dataset_url,
  original_file = h5ad_file,
  conversion_script = SCRIPT_BUILD,
  raw_count_source = raw_source,
  raw_count_assay = raw_assay_name
)


# =============================================================================
# 15. Validate converted Seurat object
# =============================================================================

validation_table <- data.frame(
  Metric = c(
    "Genes",
    "Cells",
    "Metadata_columns",
    "Assays",
    "Reductions"
  ),
  Value = c(
    nrow(
      yoshida
    ),
    ncol(
      yoshida
    ),
    ncol(
      yoshida[[]]
    ),
    paste(
      Assays(
        yoshida
      ),
      collapse = ";"
    ),
    paste(
      Reductions(
        yoshida
      ),
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  validation_table,
  file.path(
    diagnostics_dir,
    "Seurat_conversion_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  yoshida[[]],
  file.path(
    diagnostics_dir,
    "Seurat_cell_metadata.csv"
  ),
  row.names = TRUE
)


# =============================================================================
# 16. Save Seurat object
# =============================================================================

message("Saving Seurat RDS...")

saveRDS(
  yoshida,
  file = seurat_file,
  compress = FALSE
)

if (!file.exists(seurat_file)) {

  stop(
    "The Seurat object was not saved successfully."
  )
}

seurat_file_info <- file.info(
  seurat_file
)

write.csv(
  data.frame(
    File = seurat_file,
    Size_bytes = seurat_file_info$size,
    Size_GB = seurat_file_info$size / 1024^3,
    stringsAsFactors = FALSE
  ),
  file.path(
    diagnostics_dir,
    "Seurat_RDS_file_information.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 17. Reproducibility
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    diagnostics_dir,
    "sessionInfo_conversion.txt"
  )
)


# =============================================================================
# 18. Final summary
# =============================================================================

summary_text <- c(
  paste(
    "Script build:",
    SCRIPT_BUILD
  ),
  "",
  "H5AD to Seurat conversion completed.",
  "",
  paste(
    "Genes:",
    nrow(
      yoshida
    )
  ),
  paste(
    "Cells:",
    ncol(
      yoshida
    )
  ),
  paste(
    "Raw-count source:",
    raw_source
  ),
  paste(
    "Raw-count assay:",
    raw_assay_name
  ),
  paste(
    "Seurat assays:",
    paste(
      Assays(
        yoshida
      ),
      collapse = ", "
    )
  ),
  paste(
    "Imported reductions:",
    paste(
      Reductions(
        yoshida
      ),
      collapse = ", "
    )
  ),
  "",
  paste(
    "Saved Seurat object:",
    seurat_file
  ),
  "",
  "No biological filtering or downstream Seurat analysis was performed."
)

writeLines(
  summary_text,
  file.path(
    diagnostics_dir,
    "conversion_summary.txt"
  )
)

cat(
  paste(
    summary_text,
    collapse = "\n"
  ),
  "\n"
)

message("")
message("============================================================")
message("Yoshida H5AD -> Seurat conversion completed.")
message("============================================================")
message("Saved object: ", seurat_file)
message("============================================================")
