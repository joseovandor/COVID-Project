# ============================================================
# VALIDACIÓN ORTOGONAL SINGLE-CELL
#
# Dataset: Yoshida Airway
# Condición: COVID+
# Tejido: nasal cavity
# Comparación: Pediatric vs Adult
#
# Datos procesados / normalizados
# Unidad experimental = DONANTE
#
# IMPORTANTE:
# - sce conserva TODOS los genes
# - ENSEMBL permanece como rownames
# - SYMBOL se almacena en rowData
# - Los 33 genes se extraen SOLO después del filtro de células
# - No usamos edgeR/DESeq2 porque no tenemos raw counts
# ============================================================


# ============================================================
# 0. PAQUETES
# ============================================================

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


# ============================================================
# 1. ARCHIVO
# ============================================================

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


# ============================================================
# 2. DECIDIR SI HAY QUE RECARGAR sce
# ============================================================

reload_sce <- FALSE


if (exists("sce")) {
  
  cat(
    "\n========================================\n",
    "El objeto 'sce' ya existe en el ambiente.\n",
    "Comprobando si puede reutilizarse...\n",
    "========================================\n"
  )
  
  
  if (!"X" %in% assayNames(sce)) {
    
    cat(
      "\nEl objeto no contiene assay X.\n",
      "Se recargará.\n"
    )
    
    reload_sce <- TRUE
    
  } else {
    
    X_current <- assay(
      sce,
      "X"
    )
    
    
    cat(
      "\nClase actual de X:\n"
    )
    
    print(
      class(X_current)
    )
    
    
    # --------------------------------------------------------
    # Si es dgCMatrix podemos comprobar si está vacía
    # sin sumar toda la matriz
    # --------------------------------------------------------
    
    if (inherits(X_current, "dgCMatrix")) {
      
      cat(
        "\nValores almacenados en X:",
        length(X_current@x),
        "\n"
      )
      
      
      if (length(X_current@x) == 0) {
        
        cat(
          "\nLa matriz X existente está vacía.\n",
          "Se eliminará y volverá a importar.\n"
        )
        
        reload_sce <- TRUE
        
      } else {
        
        cat(
          "\nLa matriz X existente contiene datos.\n",
          "Se reutilizará sce.\n"
        )
      }
      
    } else {
      
      # Para matrices respaldadas en HDF5 no hacemos sum(X)
      # porque sería innecesariamente pesado.
      
      cat(
        "\nX no es una dgCMatrix vacía.\n",
        "Se comprobarán genes específicos más adelante.\n"
      )
    }
  }
  
} else {
  
  cat(
    "\nEl objeto sce no existe en el ambiente.\n"
  )
  
  reload_sce <- TRUE
}


# ============================================================
# 3. DESCARGAR SOLO SI NO EXISTE
# ============================================================

if (reload_sce) {
  
  if (!file.exists(file_h5ad)) {
    
    cat(
      "\nArchivo H5AD no encontrado.\n",
      "Descargando...\n"
    )
    
    
    download.file(
      url,
      destfile = file_h5ad,
      mode = "wb",
      method = "libcurl"
    )
    
  } else {
    
    cat(
      "\nEl archivo H5AD ya existe en disco.\n",
      "No se descargará nuevamente.\n"
    )
  }
  
  
  # ==========================================================
  # 4. ELIMINAR OBJETOS GRANDES DEFECTUOSOS
  # ==========================================================
  
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
  
  
  # ==========================================================
  # 5. CARGAR SIN MATERIALIZAR TODA LA MATRIZ EN RAM
  # ==========================================================
  
  cat(
    "\n========================================\n",
    "Leyendo H5AD con respaldo HDF5...\n",
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


# ============================================================
# 6. COMPROBAR OBJETO
# ============================================================

cat(
  "\n========================================\n",
  "OBJETO SCE\n",
  "========================================\n"
)

print(sce)


cat(
  "\nAssays disponibles:\n"
)

print(
  assayNames(sce)
)


cat(
  "\nDimensiones:\n"
)

print(
  dim(sce)
)


if (!"X" %in% assayNames(sce)) {
  
  stop(
    "\nERROR: el objeto no contiene un assay llamado X."
  )
}


cat(
  "\nClase de X:\n"
)

print(
  class(
    assay(sce, "X")
  )
)


# ============================================================
# 7. MAPEAR ENSEMBL -> SYMBOL
#
# sce SIGUE CONTENIENDO TODO EL TRANSCRIPTOMA
# ============================================================

if (!"ENSEMBL" %in% colnames(rowData(sce))) {
  
  rowData(sce)$ENSEMBL <- rownames(
    sce
  )
}


if (!"SYMBOL" %in% colnames(rowData(sce))) {
  
  cat(
    "\nMapeando ENSEMBL -> SYMBOL...\n"
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
  
  
  rowData(sce)$SYMBOL <- unname(
    symbols
  )
  
} else {
  
  cat(
    "\nLa columna SYMBOL ya existe.\n"
  )
}


# ============================================================
# 8. GENES DE INTERÉS
# ============================================================

genes_primary_sc <- c(
  "DOC2B",
  "BCAP29",
  "SLC9A3",
  "KANK1",
  "CLU",
  "RNF183",
  "BEST3",
  "ERVW-1",
  "SPRR2A",
  "CDK14",
  "NTS",
  "MMP12",
  "C3",
  "SLC22A15",
  "CXCL6",
  "C4orf19",
  "RILPL2",
  "SPECC1",
  "TMCO6",
  "SNX16",
  "SCRN1",
  "ZNF774",
  "DPYSL3",
  "NEU4",
  "STARD5",
  "PDZK1IP1",
  "KRT6A",
  "NIBAN1",
  "CDC14A",
  "PRKAG2",
  "FGD4",
  "ECHS1",
  "SPRR2D"
)


# ============================================================
# 9. COMPROBAR PRESENCIA
# ============================================================

presence_primary <- data.frame(
  SYMBOL = genes_primary_sc,
  PRESENT =
    genes_primary_sc %in%
    rowData(sce)$SYMBOL
)


cat(
  "\nGenes presentes:\n"
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


# ============================================================
# 10. ÍNDICES ENSEMBL DE LOS 33 GENES
# ============================================================

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
  "\nFeatures correspondientes a los genes de interés:\n"
)


print(
  gene_annotation_primary,
  row.names = FALSE
)


# ============================================================
# 11. COMPROBAR QUE LA MATRIZ REALMENTE TIENE EXPRESIÓN
#
# SOLO leemos 33 filas.
# No sumamos los 32 mil genes.
# ============================================================

X_check <- assay(
  sce,
  "X"
)[
  idx_primary,
  ,
  drop = FALSE
]


cat(
  "\nComprobando expresión de los 33 genes...\n"
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
      "\nERROR: los 33 genes continúan con expresión cero.\n",
      "La matriz X no fue importada correctamente.\n",
      "No se generarán resultados artificiales."
    )
  )
}


cat(
  "\nLa matriz contiene expresión real.\n",
  "Continuando con el análisis.\n"
)


# Liberar comprobación temporal

rm(
  X_check
)

invisible(
  gc()
)


# ============================================================
# 12. DIRECTORIOS
# ============================================================

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


# ============================================================
# 13. METADATA
# ============================================================

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


# ============================================================
# 14. FILTRAR CÉLULAS
#
# SOLO:
# - Adult / Ped
# - COVID+
# - nasal cavity
#
# sce conserva TODOS los genes
# ============================================================

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
  "\nCélulas COVID+ nasales Adult/Ped:",
  sum(keep),
  "\n"
)


if (sum(keep) == 0) {
  
  stop(
    "\nNo se encontraron células con los filtros establecidos."
  )
}


sce_age <- sce[, keep]


# ============================================================
# 15. FACTORES
# ============================================================

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


# ============================================================
# 16. RESUMEN DE CÉLULAS
# ============================================================

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
  "Células:",
  ncol(sce_age),
  "\n"
)


cat(
  "\nCélulas por grupo:\n"
)


print(
  table(
    sce_age$Group
  )
)


# ============================================================
# 17. DONANTES
# ============================================================

donors <- as.data.frame(
  colData(sce_age)
) %>%
  
  dplyr::select(
    donor_id,
    Group
  ) %>%
  
  dplyr::distinct()


cat(
  "\nDonantes totales:",
  nrow(donors),
  "\n"
)


cat(
  "\nDonantes por grupo:\n"
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


# ============================================================
# 18. CÉLULAS POR DONOR x CELL TYPE
# ============================================================

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


# ============================================================
# 19. MÍNIMO 20 CÉLULAS POR DONOR x CELL TYPE
# ============================================================

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
  "\nCélulas antes del filtro:",
  ncol(sce_age),
  "\n"
)


cat(
  "Células después de >=20:",
  ncol(sce_age_filt),
  "\n"
)


if (ncol(sce_age_filt) == 0) {
  
  stop(
    "\nNo quedaron células después del filtro >=20."
  )
}


# ============================================================
# 20. METADATA FILTRADA
# ============================================================

meta_filt <- as.data.frame(
  colData(sce_age_filt)
)


# ============================================================
# 21. AHORA SÍ EXTRAER SOLO LAS 33 FILAS DE EXPRESIÓN
#
# IMPORTANTE:
# sce_age_filt sigue conservando todos los genes.
#
# Solo X_primary es reducido porque únicamente necesitamos
# estudiar esos genes.
# ============================================================

X_primary <- assay(
  sce_age_filt,
  "X"
)[
  idx_primary,
  ,
  drop = FALSE
]


cat(
  "\nDimensiones de X_primary:\n"
)


print(
  dim(X_primary)
)


cat(
  "\nEsperado: ~33 genes x células filtradas.\n"
)


# ============================================================
# 22. DONOR x CELL TYPE
# ============================================================

group_id <- paste(
  meta_filt$donor_id,
  meta_filt$Cell_type_annotation_level2,
  sep = "___"
)


groups <- unique(
  group_id
)


cat(
  "\nNúmero de agregados donor x cell type:",
  length(groups),
  "\n"
)


if (length(groups) == 0) {
  
  stop(
    "\nNo se pudieron generar agregados donor x cell type."
  )
}


# ============================================================
# 23. EXPRESIÓN MEDIA
#
# Cada columna final = DONANTE x CELL TYPE
# ============================================================

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
  "\nDimensión mean_expr:\n"
)

print(
  dim(mean_expr)
)


# ============================================================
# 24. FRACCIÓN DE CÉLULAS CON EXPRESIÓN > 0
# ============================================================

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


# ============================================================
# 25. METADATA DE AGREGADOS
# ============================================================

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


# ============================================================
# 26. ANOTACIÓN DE LOS GENES SELECCIONADOS
# ============================================================

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


# ============================================================
# 27. FORMATO LARGO
# ============================================================

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


# ============================================================
# 28. DONANTES POR CELL TYPE
# ============================================================

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
  "\nDonantes por cell type y grupo:\n"
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


# ============================================================
# 29. MÍNIMO 3 DONANTES EN AMBOS GRUPOS
# ============================================================

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
  "\nTipos celulares comparables:",
  length(valid_celltypes),
  "\n"
)


print(
  valid_celltypes
)


if (length(valid_celltypes) == 0) {
  
  stop(
    "\nNo existen cell types con >=3 donantes en ambos grupos."
  )
}


# ============================================================
# 30. FILTRAR TIPOS CELULARES VÁLIDOS
# ============================================================

expr_primary <- expr_primary %>%
  
  dplyr::filter(
    cell_type %in%
      valid_celltypes
  )


# ============================================================
# 31. RESUMEN ENTRE DONANTES
# ============================================================

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


# ============================================================
# 32. RANKING:
# ¿DÓNDE SE EXPRESA MÁS CADA GEN?
# ============================================================

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
  "TOP 3 CELL TYPES POR GEN\n",
  "========================================\n"
)


print(
  tibble::as_tibble(
    top3_celltypes
  ),
  n = Inf
)


# ============================================================
# 33. DIFERENCIA ADULT - PEDIATRIC
#
# NO se llama logFC porque X está procesada
# ============================================================

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


# ============================================================
# 34. WILCOXON ENTRE DONANTES
# ============================================================

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


# ============================================================
# 35. FDR
# ============================================================

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
  "RESULTADOS ADULT vs PEDIATRIC\n",
  "========================================\n"
)


print(
  tibble::as_tibble(
    test_results
  ),
  n = Inf
)


# ============================================================
# 36. GUARDAR TABLAS
# ============================================================

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


# ============================================================
# 37. FIGURA A
#
# Color = expresión media entre donantes
# Tamaño = % de células positivas
# ============================================================

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


# ============================================================
# 38. FIGURA B
# ============================================================

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


# ============================================================
# 39. GUARDAR FIGURAS
# ============================================================

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


# ============================================================
# 40. RESUMEN FINAL
# ============================================================

cat(
  "\n========================================\n",
  "VALIDACIÓN ORTOGONAL FINALIZADA\n",
  "========================================\n"
)


cat(
  "Condición: COVID+\n"
)


cat(
  "Tejido: nasal cavity\n"
)


cat(
  "Comparación: Pediatric vs Adult\n"
)


cat(
  "Unidad experimental: donor\n"
)


cat(
  "Genes de interés:",
  length(
    unique(
      expr_primary$SYMBOL
    )
  ),
  "\n"
)


cat(
  "Células finales:",
  ncol(
    sce_age_filt
  ),
  "\n"
)


cat(
  "Donantes:",
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
  "Mínimo células donor x cell type:",
  min_cells,
  "\n"
)


cat(
  "Mínimo donantes por grupo:",
  min_donors,
  "\n"
)


cat(
  "\nResultados guardados en:\n",
  results_dir,
  "\n"
)