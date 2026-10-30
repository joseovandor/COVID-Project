# Yoshida et al. (2022) nasal single-cell data: import and donor-level summary
#
# COVID+ donors, nasal cavity, pediatric vs adult; the donor is the
# experimental unit. The SingleCellExperiment keeps all genes (Ensembl IDs as
# rownames, symbols in rowData). The signature genes are read from the
# meta-analysis results and extracted after the >= 20 cells per donor x cell
# type filter. This script summarises the log-normalized expression; the
# raw-count pseudobulk analysis is in 03_pseudobulk_signature.R.
#
# Input : data/sc_rna_seq/Yoshida_airway.h5ad
# Output: results/scRNAseq_COVID_nasal_age/

library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)
library(DelayedArray)
library(DelayedMatrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

options(timeout = 3600)


# 1. FILE ----

url <- paste0(
  "https://datasets.cellxgene.cziscience.com/",
  "f9efb73e-f116-46b5-a775-d4233e758024.h5ad"
)

file_h5ad <- "data/sc_rna_seq/Yoshida_airway.h5ad"


dir.create(
  "data/sc_rna_seq",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. DECIDE WHETHER sce NEEDS TO BE RELOADED ----

reload_sce <- FALSE


if (exists("sce")) {
  
  cat(
    "\n========================================\n",
    "The 'sce' object already exists in the environment.\n",
    "Checking whether it can be reused...\n",
    "========================================\n"
  )
  
  
  if (!"X" %in% assayNames(sce)) {
    
    cat(
      "\nThe object does not contain assay X.\n",
      "It will be reloaded.\n"
    )
    
    reload_sce <- TRUE
    
  } else {
    
    X_current <- assay(
      sce,
      "X"
    )
    
    
    cat(
      "\nCurrent class of X:\n"
    )
    
    print(
      class(X_current)
    )
    
    
    # If it is a dgCMatrix we can check whether it is empty ----
    # without summing the whole matrix
    
    if (inherits(X_current, "dgCMatrix")) {
      
      cat(
        "\nValues stored in X:",
        length(X_current@x),
        "\n"
      )
      
      
      if (length(X_current@x) == 0) {
        
        cat(
          "\nThe existing X matrix is empty.\n",
          "It will be removed and re-imported.\n"
        )
        
        reload_sce <- TRUE
        
      } else {
        
        cat(
          "\nThe existing X matrix contains data.\n",
          "sce will be reused.\n"
        )
      }
      
    } else {
      
      # For HDF5-backed matrices we do not run sum(X)
      # because it would be unnecessarily heavy.
      
      cat(
        "\nX is not an empty dgCMatrix.\n",
        "Specific genes will be checked later.\n"
      )
    }
  }
  
} else {
  
  cat(
    "\nThe sce object does not exist in the environment.\n"
  )
  
  reload_sce <- TRUE
}


# 3. DOWNLOAD ONLY IF IT DOES NOT EXIST ----

if (reload_sce) {
  
  if (!file.exists(file_h5ad)) {
    
    cat(
      "\nH5AD file not found.\n",
      "Downloading...\n"
    )
    
    
    download.file(
      url,
      destfile = file_h5ad,
      mode = "wb",
      method = "libcurl"
    )
    
  } else {
    
    cat(
      "\nThe H5AD file already exists on disk.\n",
      "It will not be downloaded again.\n"
    )
  }
  
  
  # 4. REMOVE LARGE, BROKEN OBJECTS ----
  
  if (exists("sce")) {
    rm(sce)
  }
  
  if (exists("X_current")) {
    rm(X_current)
  }
  
  if (exists("X")) {
    rm(X)
  }
  
  invisible(
    gc()
  )
  
  
  # 5. LOAD WITHOUT MATERIALIZING THE WHOLE MATRIX IN RAM ----
  
  cat(
    "\n========================================\n",
    "Reading H5AD with HDF5 backing...\n",
    "reader = 'R'\n",
    "use_hdf5 = TRUE\n",
    "========================================\n"
  )
  
  
  sce <- zellkonverter::readH5AD(
    file_h5ad,
    reader = "R",
    use_hdf5 = TRUE
  )
}


# 6. CHECK OBJECT ----

cat(
  "\n========================================\n",
  "SCE OBJECT\n",
  "========================================\n"
)

print(sce)


cat(
  "\nAvailable assays:\n"
)

print(
  assayNames(sce)
)


cat(
  "\nDimensions:\n"
)

print(
  dim(sce)
)


if (!"X" %in% assayNames(sce)) {
  
  stop(
    "\nERROR: the object does not contain an assay named X."
  )
}


cat(
  "\nClass of X:\n"
)

print(
  class(
    assay(sce, "X")
  )
)


# 7. MAP ENSEMBL -> SYMBOL ----
# sce STILL CONTAINS THE WHOLE TRANSCRIPTOME

if (!"ENSEMBL" %in% colnames(rowData(sce))) {
  
  rowData(sce)$ENSEMBL <- rownames(
    sce
  )
}


if (!"SYMBOL" %in% colnames(rowData(sce))) {
  
  cat(
    "\nMapping ENSEMBL -> SYMBOL...\n"
  )
  
  
  ensembl_clean <- sub(
    "\\..*$",
    "",
    rowData(sce)$ENSEMBL
  )
  
  
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  
  symbols <- unname(
    symbols
  )

  # CELLxGENE objects store the curated gene symbol in var$feature_name.
  # Use it when org.Hs.eg.db returns no symbol (e.g. lncRNAs, pseudogenes,
  # olfactory receptors), so that absence from the matrix is not confused
  # with a failed identifier mapping.
  if ("feature_name" %in% colnames(rowData(sce))) {
    feature_symbols <- as.character(rowData(sce)$feature_name)
    symbols <- ifelse(
      is.na(symbols) | symbols == "",
      feature_symbols,
      symbols
    )
  }

  rowData(sce)$SYMBOL <- symbols
  
} else {
  
  cat(
    "\nThe SYMBOL column already exists.\n"
  )
}


# 8. GENES OF INTEREST ----

signature_file <- here::here(
  "results",
  "meta_analysis",
  "primary",
  "meta_analysis_significant_FDR05.csv"
)

if (!file.exists(signature_file)) {
  stop(
    "Signature file not found: ", signature_file,
    "\nRun scripts/02_meta_analysis/02_random_effects_model.R first."
  )
}

signature_table <- read.csv(
  signature_file,
  stringsAsFactors = FALSE
)

genes_primary_sc <- unique(
  as.character(signature_table$SYMBOL)
)

cat(
  "\nSignature genes read from the meta-analysis:",
  length(genes_primary_sc),
  "\n"
)


# 9. CHECK PRESENCE ----

presence_primary <- data.frame(
  SYMBOL = genes_primary_sc,
  PRESENT =
    genes_primary_sc %in%
    rowData(sce)$SYMBOL
)


cat(
  "\nGenes present:\n"
)


print(
  presence_primary,
  row.names = FALSE
)


cat(
  "\nTotal:",
  sum(
    presence_primary$PRESENT
  ),
  "de",
  length(
    genes_primary_sc
  ),
  "\n"
)


# 10. ENSEMBL INDICES OF THE SIGNATURE GENES ----

idx_primary <- which(
  !is.na(
    rowData(sce)$SYMBOL
  ) &
    rowData(sce)$SYMBOL %in%
    genes_primary_sc
)


gene_annotation_primary <- data.frame(
  
  ENSEMBL =
    rownames(sce)[
      idx_primary
    ],
  
  SYMBOL =
    rowData(sce)$SYMBOL[
      idx_primary
    ],
  
  stringsAsFactors = FALSE
)


cat(
  "\nFeatures corresponding to the genes of interest:\n"
)


print(
  gene_annotation_primary,
  row.names = FALSE
)


# 11. CHECK THAT THE MATRIX ACTUALLY CONTAINS EXPRESSION ----
# We ONLY read the signature-gene rows.
# We do not sum the ~32 thousand genes.

X_check <- assay(
  sce,
  "X"
)[
  idx_primary,
  ,
  drop = FALSE
]


cat(
  "\nChecking expression of the signature genes...\n"
)


n_cells_gene <- DelayedMatrixStats::rowSums2(
  X_check != 0
)


total_gene_expression <- DelayedMatrixStats::rowSums2(
  X_check
)


gene_expression_check <- data.frame(
  
  ENSEMBL =
    rownames(sce)[
      idx_primary
    ],
  
  SYMBOL =
    rowData(sce)$SYMBOL[
      idx_primary
    ],
  
  n_cells_expressing =
    as.numeric(
      n_cells_gene
    ),
  
  total_expression =
    as.numeric(
      total_gene_expression
    ),
  
  stringsAsFactors = FALSE
)


print(
  gene_expression_check,
  row.names = FALSE
)


if (
  all(
    gene_expression_check$n_cells_expressing == 0
  )
) {
  
  stop(
    paste0(
      "\nERROR: the signature genes still show zero expression.\n",
      "The X matrix was not imported correctly.\n",
      "No artificial results will be generated."
    )
  )
}


cat(
  "\nThe matrix contains real expression.\n",
  "Continuing with the analysis.\n"
)


# Free temporary check

rm(
  X_check
)

invisible(
  gc()
)


# 12. DIRECTORIES ----

results_dir <- "results/scRNAseq_COVID_nasal_age"


dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    results_dir,
    "tables"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    results_dir,
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    results_dir,
    "figures",
    "per_gene"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  gene_expression_check,
  file.path(
    results_dir,
    "tables",
    "primary_genes_initial_expression_check.csv"
  ),
  row.names = FALSE
)


# 13. METADATA ----

meta <- as.data.frame(
  colData(sce)
)


cat(
  "\n========================================\n",
  "METADATA\n",
  "========================================\n"
)


cat(
  "\nGroup:\n"
)

print(
  table(
    meta$Group,
    useNA = "ifany"
  )
)


cat(
  "\nCOVID_status:\n"
)

print(
  table(
    meta$COVID_status,
    useNA = "ifany"
  )
)


cat(
  "\nTissue:\n"
)

print(
  table(
    meta$tissue,
    useNA = "ifany"
  )
)


# 14. FILTER CELLS
#
# ONLY:
# - Adult / Ped
# - COVID+
# - nasal cavity
#
# sce retains ALL genes

keep <- (
  meta$Group %in% c(
    "Adult",
    "Ped"
  ) &
    meta$COVID_status == "COVID+" &
    meta$tissue == "nasal cavity" &
    !is.na(
      meta$donor_id
    ) &
    !is.na(
      meta$Cell_type_annotation_level2
    )
)


cat(
  "\nCOVID+ nasal cells Adult/Ped:",
  sum(keep),
  "\n"
)


if (sum(keep) == 0) {
  
  stop(
    "\nNo cells were found with the established filters."
  )
}


sce_age <- sce[, keep]


# 15. FACTORS ----

sce_age$Group <- factor(
  sce_age$Group,
  levels = c(
    "Ped",
    "Adult"
  )
)


sce_age$donor_id <- droplevels(
  factor(
    sce_age$donor_id
  )
)


sce_age$Cell_type_annotation_level2 <-
  droplevels(
    factor(
      sce_age$Cell_type_annotation_level2
    )
  )


# 16. CELL SUMMARY ----

cat(
  "\n========================================\n",
  "COVID+ + NASAL CAVITY\n",
  "========================================\n"
)


cat(
  "\nGenes:",
  nrow(sce_age),
  "\n"
)


cat(
  "Cells:",
  ncol(sce_age),
  "\n"
)


cat(
  "\nCells per group:\n"
)


print(
  table(
    sce_age$Group
  )
)


# 17. DONORS ----

donors <- as.data.frame(
  colData(sce_age)
) %>%
  
  dplyr::select(
    donor_id,
    Group
  ) %>%
  
  dplyr::distinct()


cat(
  "\nTotal donors:",
  nrow(donors),
  "\n"
)


cat(
  "\nDonors per group:\n"
)


print(
  table(
    donors$Group
  )
)


write.csv(
  donors,
  file.path(
    results_dir,
    "tables",
    "donors_included.csv"
  ),
  row.names = FALSE
)


# 18. CELLS PER DONOR x CELL TYPE ----

cell_numbers <- as.data.frame(
  colData(sce_age)
) %>%
  
  dplyr::count(
    donor_id,
    Group,
    Cell_type_annotation_level2,
    name = "n_cells"
  )


write.csv(
  cell_numbers,
  file.path(
    results_dir,
    "tables",
    "cell_numbers_by_donor_celltype.csv"
  ),
  row.names = FALSE
)


# 19. MINIMUM 20 CELLS PER DONOR x CELL TYPE ----

min_cells <- 20


valid_combinations <- cell_numbers %>%
  
  dplyr::filter(
    n_cells >= min_cells
  ) %>%
  
  dplyr::mutate(
    combination = paste(
      donor_id,
      Cell_type_annotation_level2,
      sep = "___"
    )
  )


cell_combination <- paste(
  sce_age$donor_id,
  sce_age$Cell_type_annotation_level2,
  sep = "___"
)


keep2 <- cell_combination %in%
  valid_combinations$combination


sce_age_filt <- sce_age[, keep2]


cat(
  "\nCells before the filter:",
  ncol(sce_age),
  "\n"
)


cat(
  "Cells after >=20:",
  ncol(sce_age_filt),
  "\n"
)


if (ncol(sce_age_filt) == 0) {
  
  stop(
    "\nNo cells remained after the >=20 filter."
  )
}


# 20. FILTERED METADATA ----

meta_filt <- as.data.frame(
  colData(sce_age_filt)
)


# 21. NOW EXTRACT ONLY THE SIGNATURE-GENE EXPRESSION ROWS
#
# IMPORTANT:
# sce_age_filt still retains all genes.
#
# Only X_primary is reduced because we only need
# to study those genes.

X_primary <- assay(
  sce_age_filt,
  "X"
)[
  idx_primary,
  ,
  drop = FALSE
]


cat(
  "\nDimensions of X_primary:\n"
)


print(
  dim(X_primary)
)


cat(
  "\nExpected: signature genes present in the matrix x filtered cells.\n"
)


# 22. DONOR x CELL TYPE ----

group_id <- paste(
  meta_filt$donor_id,
  meta_filt$Cell_type_annotation_level2,
  sep = "___"
)


groups <- unique(
  group_id
)


cat(
  "\nNumber of donor x cell type aggregates:",
  length(groups),
  "\n"
)


if (length(groups) == 0) {
  
  stop(
    "\nDonor x cell type aggregates could not be generated."
  )
}


# 23. MEAN EXPRESSION ----
# Each final column = DONOR x CELL TYPE

mean_expr_list <- lapply(
  groups,
  function(g) {
    
    idx <- which(
      group_id == g
    )
    
    
    values <- DelayedMatrixStats::rowMeans2(
      X_primary[, idx, drop = FALSE],
      na.rm = TRUE
    )
    
    
    as.numeric(
      values
    )
  }
)


mean_expr <- do.call(
  cbind,
  mean_expr_list
)


mean_expr <- matrix(
  mean_expr,
  nrow = nrow(X_primary),
  ncol = length(groups)
)


rownames(mean_expr) <- rownames(
  X_primary
)


colnames(mean_expr) <- groups


cat(
  "\nDimension of mean_expr:\n"
)

print(
  dim(mean_expr)
)


# 24. FRACTION OF CELLS WITH EXPRESSION > 0 ----

fraction_expr_list <- lapply(
  groups,
  function(g) {
    
    idx <- which(
      group_id == g
    )
    
    
    values <- DelayedMatrixStats::rowMeans2(
      X_primary[, idx, drop = FALSE] > 0,
      na.rm = TRUE
    )
    
    
    as.numeric(
      values
    )
  }
)


fraction_expr <- do.call(
  cbind,
  fraction_expr_list
)


fraction_expr <- matrix(
  fraction_expr,
  nrow = nrow(X_primary),
  ncol = length(groups)
)


rownames(fraction_expr) <- rownames(
  X_primary
)


colnames(fraction_expr) <- groups


# 25. AGGREGATE METADATA ----

aggregate_meta <- data.frame(
  aggregate_id = groups,
  stringsAsFactors = FALSE
)


split_ids <- strsplit(
  aggregate_meta$aggregate_id,
  "___",
  fixed = TRUE
)


aggregate_meta$donor_id <- sapply(
  split_ids,
  `[`,
  1
)


aggregate_meta$cell_type <- sapply(
  split_ids,
  `[`,
  2
)


donor_group <- meta_filt %>%
  
  dplyr::select(
    donor_id,
    Group
  ) %>%
  
  dplyr::distinct()


donor_group$donor_id <- as.character(
  donor_group$donor_id
)


aggregate_meta <- aggregate_meta %>%
  
  dplyr::left_join(
    donor_group,
    by = "donor_id"
  )


aggregate_meta$Group <- factor(
  aggregate_meta$Group,
  levels = c(
    "Ped",
    "Adult"
  )
)


# 26. ANNOTATION OF SELECTED GENES ----

gene_annotation <- data.frame(
  
  ENSEMBL =
    rownames(
      X_primary
    ),
  
  SYMBOL =
    rowData(sce)$SYMBOL[
      match(
        rownames(X_primary),
        rownames(sce)
      )
    ],
  
  stringsAsFactors = FALSE
)


# 27. LONG FORMAT ----

expr_primary <- data.frame(
  
  ENSEMBL = rep(
    rownames(mean_expr),
    times = ncol(mean_expr)
  ),
  
  aggregate_id = rep(
    colnames(mean_expr),
    each = nrow(mean_expr)
  ),
  
  mean_expression =
    as.vector(
      mean_expr
    ),
  
  fraction_expressing =
    as.vector(
      fraction_expr
    ),
  
  stringsAsFactors = FALSE
)


expr_primary <- expr_primary %>%
  
  dplyr::left_join(
    gene_annotation,
    by = "ENSEMBL"
  ) %>%
  
  dplyr::left_join(
    aggregate_meta,
    by = "aggregate_id"
  )


expr_primary$Group <- factor(
  expr_primary$Group,
  levels = c(
    "Ped",
    "Adult"
  )
)


# 28. DONORS PER CELL TYPE ----

donor_counts <- aggregate_meta %>%
  
  dplyr::distinct(
    donor_id,
    Group,
    cell_type
  ) %>%
  
  dplyr::count(
    cell_type,
    Group,
    name = "n_donors"
  )


cat(
  "\nDonors per cell type and group:\n"
)


print(
  tibble::as_tibble(
    donor_counts
  ),
  n = Inf
)


write.csv(
  donor_counts,
  file.path(
    results_dir,
    "tables",
    "donors_by_celltype_group.csv"
  ),
  row.names = FALSE
)


# 29. MINIMUM 3 DONORS IN BOTH GROUPS ----

min_donors <- 3


valid_celltypes <- donor_counts %>%
  
  dplyr::filter(
    n_donors >= min_donors
  ) %>%
  
  dplyr::group_by(
    cell_type
  ) %>%
  
  dplyr::summarise(
    n_groups =
      dplyr::n_distinct(Group),
    .groups = "drop"
  ) %>%
  
  dplyr::filter(
    n_groups == 2
  ) %>%
  
  dplyr::pull(
    cell_type
  )


cat(
  "\nComparable cell types:",
  length(valid_celltypes),
  "\n"
)


print(
  valid_celltypes
)


if (length(valid_celltypes) == 0) {
  
  stop(
    "\nNo cell types have >=3 donors in both groups."
  )
}


# 30. FILTER VALID CELL TYPES ----

expr_primary <- expr_primary %>%
  
  dplyr::filter(
    cell_type %in%
      valid_celltypes
  )


# 31. SUMMARY ACROSS DONORS ----

expression_summary <- expr_primary %>%
  
  dplyr::group_by(
    Group,
    cell_type,
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    mean_expression =
      mean(
        mean_expression,
        na.rm = TRUE
      ),
    
    median_expression =
      median(
        mean_expression,
        na.rm = TRUE
      ),
    
    mean_fraction_expressing =
      mean(
        fraction_expressing,
        na.rm = TRUE
      ),
    
    n_donors =
      dplyr::n_distinct(
        donor_id
      ),
    
    .groups = "drop"
  )


# 32. RANKING: ----
# WHERE IS EACH GENE MOST EXPRESSED?

expression_rank <- expression_summary %>%
  
  dplyr::group_by(
    Group,
    SYMBOL
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      mean_expression
    ),
    .by_group = TRUE
  ) %>%
  
  dplyr::mutate(
    expression_rank =
      dplyr::row_number()
  ) %>%
  
  dplyr::ungroup()


top3_celltypes <- expression_rank %>%
  
  dplyr::filter(
    expression_rank <= 3
  )


cat(
  "\n========================================\n",
  "TOP 3 CELL TYPES PER GENE\n",
  "========================================\n"
)


print(
  tibble::as_tibble(
    top3_celltypes
  ),
  n = Inf
)


# 33. ADULT - PEDIATRIC DIFFERENCE ----
# Not called logFC because X is already processed data

age_difference <- expression_summary %>%
  
  dplyr::select(
    Group,
    cell_type,
    SYMBOL,
    mean_expression
  ) %>%
  
  tidyr::pivot_wider(
    names_from = Group,
    values_from = mean_expression
  ) %>%
  
  dplyr::filter(
    !is.na(Ped),
    !is.na(Adult)
  ) %>%
  
  dplyr::mutate(
    
    Adult_minus_Ped =
      Adult - Ped,
    
    abs_difference =
      abs(
        Adult_minus_Ped
      ),
    
    direction =
      dplyr::case_when(
        
        Adult_minus_Ped > 0 ~
          "Adult",
        
        Adult_minus_Ped < 0 ~
          "Pediatric",
        
        TRUE ~
          "Similar"
      )
  )


# 34. WILCOXON ACROSS DONORS ----

test_results_list <- list()

counter <- 1


for (ct in valid_celltypes) {
  
  for (gene in genes_primary_sc) {
    
    tmp <- expr_primary %>%
      
      dplyr::filter(
        cell_type == ct,
        SYMBOL == gene
      )
    
    
    ped_values <- tmp$mean_expression[
      tmp$Group == "Ped"
    ]
    
    
    adult_values <- tmp$mean_expression[
      tmp$Group == "Adult"
    ]
    
    
    if (
      length(ped_values) >= min_donors &&
      length(adult_values) >= min_donors
    ) {
      
      test <- wilcox.test(
        adult_values,
        ped_values,
        exact = FALSE
      )
      
      
      test_results_list[[counter]] <- data.frame(
        
        SYMBOL = gene,
        
        cell_type = ct,
        
        n_Ped =
          length(
            ped_values
          ),
        
        n_Adult =
          length(
            adult_values
          ),
        
        mean_Ped =
          mean(
            ped_values,
            na.rm = TRUE
          ),
        
        mean_Adult =
          mean(
            adult_values,
            na.rm = TRUE
          ),
        
        median_Ped =
          median(
            ped_values,
            na.rm = TRUE
          ),
        
        median_Adult =
          median(
            adult_values,
            na.rm = TRUE
          ),
        
        Adult_minus_Ped =
          mean(
            adult_values,
            na.rm = TRUE
          ) -
          mean(
            ped_values,
            na.rm = TRUE
          ),
        
        p_value =
          test$p.value,
        
        stringsAsFactors = FALSE
      )
      
      
      counter <- counter + 1
    }
  }
}


test_results <- dplyr::bind_rows(
  test_results_list
)


# 35. FDR ----

if (nrow(test_results) > 0) {
  
  test_results <- test_results %>%
    
    dplyr::mutate(
      
      FDR =
        p.adjust(
          p_value,
          method = "BH"
        ),
      
      direction =
        dplyr::case_when(
          
          Adult_minus_Ped > 0 ~
            "Adult",
          
          Adult_minus_Ped < 0 ~
            "Pediatric",
          
          TRUE ~
            "Similar"
        ),
      
      abs_difference =
        abs(
          Adult_minus_Ped
        )
    ) %>%
    
    dplyr::arrange(
      FDR
    )
}


cat(
  "\n========================================\n",
  "ADULT vs PEDIATRIC RESULTS\n",
  "========================================\n"
)


print(
  tibble::as_tibble(
    test_results
  ),
  n = Inf
)


# 36. SAVE TABLES ----

write.csv(
  expression_summary,
  file.path(
    results_dir,
    "tables",
    "primary_genes_expression_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  expression_rank,
  file.path(
    results_dir,
    "tables",
    "primary_genes_celltype_expression_rank.csv"
  ),
  row.names = FALSE
)


write.csv(
  top3_celltypes,
  file.path(
    results_dir,
    "tables",
    "primary_genes_top3_celltypes.csv"
  ),
  row.names = FALSE
)


write.csv(
  age_difference,
  file.path(
    results_dir,
    "tables",
    "primary_genes_Adult_minus_Pediatric.csv"
  ),
  row.names = FALSE
)


write.csv(
  test_results,
  file.path(
    results_dir,
    "tables",
    "primary_genes_Adult_vs_Pediatric_Wilcoxon.csv"
  ),
  row.names = FALSE
)


# 37. FIGURE A ----
# Color = mean expression across donors
# Size = % of positive cells

pA <- ggplot(
  expression_summary,
  aes(
    x = cell_type,
    y = SYMBOL
  )
) +
  
  geom_point(
    aes(
      size =
        mean_fraction_expressing * 100,
      
      color =
        mean_expression
    )
  ) +
  
  facet_wrap(
    ~ Group,
    ncol = 1
  ) +
  
  labs(
    
    title =
      "Cell-type-specific expression during COVID-19",
    
    subtitle =
      "Nasal samples from pediatric and adult donors",
    
    x =
      "Cell type",
    
    y =
      "Gene",
    
    color =
      "Mean\nexpression",
    
    size =
      "% expressing"
  ) +
  
  theme_bw(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
    
    panel.grid =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      )
  )


print(pA)


# 38. FIGURE B ----

pB <- ggplot(
  age_difference,
  aes(
    x = cell_type,
    y = SYMBOL,
    fill = Adult_minus_Ped
  )
) +
  
  geom_tile() +
  
  scale_fill_gradient2(
    midpoint = 0
  ) +
  
  labs(
    
    title =
      "Adult versus pediatric expression during COVID-19",
    
    subtitle =
      "Nasal samples; positive values indicate higher expression in adults",
    
    x =
      "Cell type",
    
    y =
      "Gene",
    
    fill =
      "Adult - Ped"
  ) +
  
  theme_bw(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
    
    panel.grid =
      element_blank()
  )


print(pB)


# 39. SAVE FIGURES ----

ggsave(
  file.path(
    results_dir,
    "figures",
    "Figure_A_celltype_expression_dotplot.pdf"
  ),
  pA,
  width = 15,
  height = 10
)


ggsave(
  file.path(
    results_dir,
    "figures",
    "Figure_A_celltype_expression_dotplot.tiff"
  ),
  pA,
  width = 15,
  height = 10,
  dpi = 600,
  compression = "lzw"
)


ggsave(
  file.path(
    results_dir,
    "figures",
    "Figure_B_Adult_minus_Pediatric.pdf"
  ),
  pB,
  width = 15,
  height = 9
)


ggsave(
  file.path(
    results_dir,
    "figures",
    "Figure_B_Adult_minus_Pediatric.tiff"
  ),
  pB,
  width = 15,
  height = 9,
  dpi = 600,
  compression = "lzw"
)


# 40. FINAL SUMMARY ----

cat(
  "\n========================================\n",
  "ORTHOGONAL VALIDATION COMPLETE\n",
  "========================================\n"
)


cat(
  "Condition: COVID+\n"
)


cat(
  "Tissue: nasal cavity\n"
)


cat(
  "Comparison: Pediatric vs Adult\n"
)


cat(
  "Experimental unit: donor\n"
)


cat(
  "Genes of interest:",
  length(
    unique(
      expr_primary$SYMBOL
    )
  ),
  "\n"
)


cat(
  "Final cells:",
  ncol(
    sce_age_filt
  ),
  "\n"
)


cat(
  "Donors:",
  length(
    unique(
      sce_age_filt$donor_id
    )
  ),
  "\n"
)


cat(
  "Tipos celulares comparables:",
  length(
    valid_celltypes
  ),
  "\n"
)


cat(
  "Minimum cells per donor x cell type:",
  min_cells,
  "\n"
)


cat(
  "Minimum donors per group:",
  min_donors,
  "\n"
)


cat(
  "\nResults saved to:\n",
  results_dir,
  "\n"
)
