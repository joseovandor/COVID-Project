# ============================================================
# COMPLEXHEATMAP FINAL
#
# Epithelial cell-type localization
# Pediatric | Adult emparejados
#
# Filas:
# Pediatric-enriched
# Adult-enriched
#
# Pediatric = amarillo
# Adult     = azul
#
# Expresión relativa:
# Azul -> blanco -> amarillo
#
# Significancia:
# *** = FDR < 0.001
# **  = FDR < 0.01
# *   = FDR < 0.05
#
# Sin marcas para FDR >= 0.05
#
# TIFF:
# 600 dpi
# compression = "lzw"
# ============================================================


# ============================================================
# 0. PAQUETES
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(tibble)


# ============================================================
# 1. LIMITAR Z-SCORE
# ============================================================

expression_scaled_plot <- expression_scaled


expression_scaled_plot[
  expression_scaled_plot > 3
] <- 3


expression_scaled_plot[
  expression_scaled_plot < -3
] <- -3


# ============================================================
# 2. PEDIATRIC-ENRICHED / ADULT-ENRICHED
# ============================================================

gene_direction <- age_difference %>%
  
  dplyr::group_by(
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    mean_age_difference =
      mean(
        Adult_minus_Ped,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    enrichment = dplyr::case_when(
      
      mean_age_difference < 0 ~
        "Pediatric-enriched",
      
      mean_age_difference > 0 ~
        "Adult-enriched",
      
      TRUE ~
        "Similar"
    )
  )


# ============================================================
# 3. ANOTACIÓN DE FILAS
# ============================================================

row_annotation_df <- data.frame(
  
  SYMBOL =
    rownames(
      expression_scaled_plot
    ),
  
  stringsAsFactors = FALSE
  
) %>%
  
  dplyr::left_join(
    gene_direction,
    by = "SYMBOL"
  ) %>%
  
  dplyr::filter(
    
    enrichment %in% c(
      "Pediatric-enriched",
      "Adult-enriched"
    )
  )


# ============================================================
# 4. ORDEN DE GENES
# ============================================================

genes_pediatric <- gene_order[
  gene_order %in%
    row_annotation_df$SYMBOL[
      row_annotation_df$enrichment ==
        "Pediatric-enriched"
    ]
]


genes_adult <- gene_order[
  gene_order %in%
    row_annotation_df$SYMBOL[
      row_annotation_df$enrichment ==
        "Adult-enriched"
    ]
]


final_gene_order <- c(
  genes_pediatric,
  genes_adult
)


expression_scaled_plot <- expression_scaled_plot[
  final_gene_order,
  ,
  drop = FALSE
]


# ============================================================
# 5. ROW SPLIT
# ============================================================

row_split <- gene_direction$enrichment[
  match(
    rownames(expression_scaled_plot),
    gene_direction$SYMBOL
  )
]


row_split <- factor(
  row_split,
  levels = c(
    "Pediatric-enriched",
    "Adult-enriched"
  )
)


# ============================================================
# 6. INFORMACIÓN DE COLUMNAS
# ============================================================

column_celltype <- sub(
  "___.*$",
  "",
  colnames(expression_scaled_plot)
)


column_age <- sub(
  "^.*___",
  "",
  colnames(expression_scaled_plot)
)


column_celltype <- factor(
  column_celltype,
  levels = celltypes_present
)


column_age <- factor(
  column_age,
  levels = c(
    "Pediatric",
    "Adult"
  )
)


# ============================================================
# 7. COLORES
# ============================================================

age_colors <- c(
  
  "Pediatric" = "#E9C46A",
  
  "Adult" = "#457B9D"
)


# ============================================================
# 8. PALETA DE EXPRESIÓN
# ============================================================

col_fun <- circlize::colorRamp2(
  
  c(
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3
  ),
  
  c(
    "#3B6FB6",
    "#6F91CF",
    "#B8C9E6",
    "#F7F7F7",
    "#F6E8B1",
    "#F5D66E",
    "#FFD43B"
  )
)


# ============================================================
# 9. MATRIZ DE SIGNIFICANCIA
#
# Solo FDR < 0.05
#
# *** = FDR < 0.001
# **  = FDR < 0.01
# *   = FDR < 0.05
#
# Se muestra solamente en la columna Adult.
# ============================================================

sig_matrix_final <- matrix(
  
  "",
  
  nrow =
    nrow(
      expression_scaled_plot
    ),
  
  ncol =
    ncol(
      expression_scaled_plot
    )
)


rownames(
  sig_matrix_final
) <- rownames(
  expression_scaled_plot
)


colnames(
  sig_matrix_final
) <- colnames(
  expression_scaled_plot
)


# ============================================================
# 10. SELECCIONAR RESULTADOS SIGNIFICATIVOS
# ============================================================

sig_results <- test_results %>%
  
  dplyr::filter(
    
    cell_type %in%
      celltypes_present,
    
    !is.na(FDR),
    
    FDR < 0.05
  ) %>%
  
  dplyr::mutate(
    
    mark = dplyr::case_when(
      
      FDR < 0.001 ~ "***",
      
      FDR < 0.01 ~ "**",
      
      FDR < 0.05 ~ "*",
      
      TRUE ~ ""
    )
  )


# ============================================================
# 11. MOSTRAR RESULTADOS SIGNIFICATIVOS
#
# CORREGIDO:
# convertir a tibble antes de usar n = Inf
# ============================================================

cat(
  "\n========================================\n",
  "SIGNIFICANT RESULTS SHOWN IN HEATMAP\n",
  "========================================\n"
)


sig_results_print <- sig_results %>%
  
  dplyr::select(
    SYMBOL,
    cell_type,
    FDR,
    mark
  ) %>%
  
  tibble::as_tibble()


print(
  sig_results_print,
  n = Inf
)


# ============================================================
# 12. COLOCAR ASTERISCOS EN COLUMNA ADULT
# ============================================================

if (nrow(sig_results) > 0) {
  
  for (i in seq_len(nrow(sig_results))) {
    
    gene_i <-
      sig_results$SYMBOL[i]
    
    
    celltype_i <-
      sig_results$cell_type[i]
    
    
    mark_i <-
      sig_results$mark[i]
    
    
    adult_column <- paste0(
      celltype_i,
      "___Adult"
    )
    
    
    if (
      gene_i %in%
      rownames(
        sig_matrix_final
      ) &&
      adult_column %in%
      colnames(
        sig_matrix_final
      )
    ) {
      
      sig_matrix_final[
        gene_i,
        adult_column
      ] <- mark_i
    }
  }
}


# ============================================================
# 13. ANOTACIÓN SUPERIOR
# ============================================================

top_anno <- HeatmapAnnotation(
  
  Age =
    column_age,
  
  
  CellType =
    anno_block(
      
      gp =
        gpar(
          
          fill =
            "#F2F2F2",
          
          col =
            "black",
          
          lwd =
            0.7
        ),
      
      
      labels =
        celltypes_present,
      
      
      labels_gp =
        gpar(
          
          fontsize =
            13,
          
          fontface =
            "plain",
          
          col =
            "black"
        ),
      
      
      height =
        unit(
          12,
          "mm"
        )
    ),
  
  
  col = list(
    Age = age_colors
  ),
  
  
  show_legend = c(
    Age = FALSE
  ),
  
  
  show_annotation_name =
    FALSE,
  
  
  simple_anno_size =
    unit(
      5,
      "mm"
    ),
  
  
  gap =
    unit(
      1.5,
      "mm"
    )
)


# ============================================================
# 14. LABELS DE COLUMNAS
# ============================================================

column_labels <- as.character(
  column_age
)


# ============================================================
# 15. HEATMAP
# ============================================================

ht <- Heatmap(
  
  expression_scaled_plot,
  
  
  name =
    "Relative expression",
  
  col =
    col_fun,
  
  
  # ==========================================================
  # FILAS
  # ==========================================================
  
  cluster_rows =
    FALSE,
  
  
  show_row_dend =
    FALSE,
  
  
  show_row_names =
    TRUE,
  
  
  row_names_side =
    "left",
  
  
  row_names_gp =
    gpar(
      
      fontsize =
        15,
      
      fontface =
        "italic",
      
      col =
        "black"
    ),
  
  
  row_split =
    row_split,
  
  
  cluster_row_slices =
    FALSE,
  
  
  row_gap =
    unit(
      5,
      "mm"
    ),
  
  
  row_title_rot =
    90,
  
  
  row_title_gp =
    gpar(
      
      fontsize =
        16,
      
      fontface =
        "plain",
      
      col =
        "black"
    ),
  
  
  # ==========================================================
  # COLUMNAS
  # ==========================================================
  
  cluster_columns =
    FALSE,
  
  
  cluster_column_slices =
    FALSE,
  
  
  column_split =
    column_celltype,
  
  
  column_gap =
    unit(
      2.5,
      "mm"
    ),
  
  
  column_labels =
    column_labels,
  
  
  show_column_names =
    TRUE,
  
  
  column_names_side =
    "bottom",
  
  
  column_names_rot =
    45,
  
  
  column_names_gp =
    gpar(
      
      fontsize =
        13,
      
      fontface =
        "plain",
      
      col =
        "black"
    ),
  
  
  column_title =
    NULL,
  
  
  # ==========================================================
  # ANOTACIÓN SUPERIOR
  # ==========================================================
  
  top_annotation =
    top_anno,
  
  
  # ==========================================================
  # BORDES
  # ==========================================================
  
  rect_gp =
    gpar(
      
      col =
        "black",
      
      lwd =
        0.35
    ),
  
  
  # ==========================================================
  # ASTERISCOS DE SIGNIFICANCIA
  # ==========================================================
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    mark <- sig_matrix_final[
      i,
      j
    ]
    
    
    if (
      !is.na(mark) &&
      mark != ""
    ) {
      
      grid.text(
        
        mark,
        
        x =
          x,
        
        y =
          y,
        
        gp =
          gpar(
            
            fontsize =
              16,
            
            fontface =
              "plain",
            
            col =
              "black"
          )
      )
    }
  },
  
  
  # ==========================================================
  # LEYENDA HORIZONTAL ABAJO
  # ==========================================================
  
  heatmap_legend_param = list(
    
    title =
      "Relative expression",
    
    
    at = c(
      -3,
      -2,
      -1,
      0,
      1,
      2,
      3
    ),
    
    
    labels = c(
      "-3",
      "-2",
      "-1",
      "0",
      "1",
      "2",
      "3"
    ),
    
    
    color_bar =
      "continuous",
    
    
    direction =
      "horizontal",
    
    
    legend_width =
      unit(
        75,
        "mm"
      ),
    
    
    legend_height =
      unit(
        6,
        "mm"
      ),
    
    
    title_position =
      "leftcenter",
    
    
    title_gp =
      gpar(
        
        fontsize =
          13,
        
        fontface =
          "plain",
        
        col =
          "black"
      ),
    
    
    labels_gp =
      gpar(
        
        fontsize =
          12,
        
        fontface =
          "plain",
        
        col =
          "black"
      )
  ),
  
  
  border =
    FALSE
)


# ============================================================
# 16. MOSTRAR EN RSTUDIO
# ============================================================

draw(
  
  ht,
  
  
  heatmap_legend_side =
    "bottom",
  
  
  annotation_legend_side =
    "right",
  
  
  merge_legends =
    FALSE,
  
  
  padding =
    unit(
      c(
        8,
        8,
        20,
        8
      ),
      "mm"
    )
)


# ============================================================
# 17. GUARDAR TIFF
#
# 600 dpi
# compression LZW
# ============================================================

tiff(
  
  filename =
    file.path(
      figure_dir,
      "Epithelial_Pediatric_Adult_ComplexHeatmap_significant.tiff"
    ),
  
  
  width =
    18,
  
  
  height =
    14,
  
  
  units =
    "in",
  
  
  res =
    600,
  
  
  compression =
    "lzw"
)


draw(
  
  ht,
  
  
  heatmap_legend_side =
    "bottom",
  
  
  annotation_legend_side =
    "right",
  
  
  merge_legends =
    FALSE,
  
  
  padding =
    unit(
      c(
        8,
        8,
        20,
        8
      ),
      "mm"
    )
)


dev.off()


# ============================================================
# 18. RESUMEN
# ============================================================

cat(
  "\n========================================\n",
  "HEATMAP FINAL GENERADO\n",
  "========================================\n"
)


cat(
  "\n*** = FDR < 0.001\n"
)


cat(
  "**  = FDR < 0.01\n"
)


cat(
  "*   = FDR < 0.05\n"
)


cat(
  "\nFigura guardada en:\n",
  file.path(
    figure_dir,
    "Epithelial_Pediatric_Adult_ComplexHeatmap_significant.tiff"
  ),
  "\n"
)