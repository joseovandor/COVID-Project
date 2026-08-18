# ============================================================
# YOSHIDA AIRWAY COVID-19
#
# A. UMAP de tipos celulares epiteliales utilizados
# B. Proporciones celulares Pediatric vs Adult
#
# UMAP:
# - COVID+
# - nasal cavity
# - Pediatric + Adult juntos
# - solo tipos epiteliales válidos
# - paleta Scanpy / Matplotlib tab10
#
# PROPORCIONES:
# - calculadas por DONANTE
# - ANTES del filtro >=20 células donor x cell type
# - Wilcoxon Pediatric vs Adult
# - corrección BH/FDR
#
# Estilo:
# - similar a Scanpy / Nature
# - leyendas a la derecha
# - textos grandes
# - sin títulos redundantes
# ============================================================


# ============================================================
# 0. PAQUETES
# ============================================================

library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)


# ============================================================
# 1. TIPOS CELULARES EPITELIALES
# ============================================================

epithelial_celltypes <- c(
  "Basal",
  "Basal cycling",
  "Multiciliated",
  "Club",
  "Duct",
  "Goblet",
  "Goblet inflammatory",
  "Hillock",
  "Rare epi",
  "Secretory",
  "Squamous",
  "Transit epi"
)


# ============================================================
# 2. SOLO EPITELIALES QUE FUERON VÁLIDOS
#    EN EL ANÁLISIS PRINCIPAL
#
# valid_celltypes ya contiene los tipos celulares
# que tuvieron >=3 donantes en ambos grupos.
# ============================================================

epithelial_valid <- epithelial_celltypes[
  epithelial_celltypes %in% valid_celltypes
]


cat(
  "\n========================================\n",
  "EPITHELIAL CELL TYPES INCLUDED\n",
  "========================================\n"
)


print(
  epithelial_valid
)


cat(
  "\nNúmero de tipos celulares:",
  length(epithelial_valid),
  "\n"
)


# ============================================================
# 3. PALETA SCANPY / MATPLOTLIB TAB10
#
# Colores categóricos clásicos de Matplotlib
# ============================================================

scanpy_tab10 <- c(
  "#1F77B4",
  "#FF7F0E",
  "#2CA02C",
  "#D62728",
  "#9467BD",
  "#8C564B",
  "#E377C2",
  "#7F7F7F",
  "#BCBD22",
  "#17BECF"
)


# Si tienes más de 10 poblaciones:
# usaríamos tab20.
# Para este caso deberían ser <=10.

if (length(epithelial_valid) > length(scanpy_tab10)) {
  
  stop(
    "Hay más de 10 categorías. Usa una paleta tab20."
  )
}


celltype_colors <- scanpy_tab10[
  seq_along(
    epithelial_valid
  )
]


names(
  celltype_colors
) <- epithelial_valid


# ============================================================
# PARTE A
# UMAP
# ============================================================


# ============================================================
# 4. CÉLULAS PARA UMAP
#
# sce_age_filt contiene:
# - COVID+
# - nasal cavity
# - Pediatric / Adult
# - combinaciones donor x cell type >=20 células
# ============================================================

keep_umap <- sce_age_filt$Cell_type_annotation_level2 %in%
  epithelial_valid


sce_umap_epi <- sce_age_filt[
  ,
  keep_umap
]


cat(
  "\nCélulas incluidas en UMAP:",
  ncol(sce_umap_epi),
  "\n"
)


# ============================================================
# 5. COMPROBAR UMAP
# ============================================================

if (!"X_umap" %in% reducedDimNames(sce_umap_epi)) {
  
  stop(
    "No existe X_umap en reducedDims(sce_umap_epi)."
  )
}


# ============================================================
# 6. EXTRAER COORDENADAS
# ============================================================

umap_coordinates <- reducedDim(
  sce_umap_epi,
  "X_umap"
)


umap_df <- data.frame(
  
  UMAP1 =
    umap_coordinates[, 1],
  
  UMAP2 =
    umap_coordinates[, 2],
  
  cell_type =
    as.character(
      sce_umap_epi$Cell_type_annotation_level2
    ),
  
  Group =
    as.character(
      sce_umap_epi$Group
    ),
  
  donor_id =
    as.character(
      sce_umap_epi$donor_id
    ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 7. ORDEN DE TIPOS CELULARES
# ============================================================

umap_df$cell_type <- factor(
  
  umap_df$cell_type,
  
  levels =
    epithelial_valid
)


# ============================================================
# 8. UMAP ESTILO SCANPY / NATURE
#
# Un solo UMAP:
# Pediatric + Adult juntos
# ============================================================

p_umap <- ggplot(
  
  umap_df,
  
  aes(
    x = UMAP1,
    y = UMAP2,
    color = cell_type
  )
  
) +
  
  # ----------------------------------------------------------
# Puntos:
# relativamente pequeños para conservar estructura,
# pero visibles a 600 dpi
# ----------------------------------------------------------

geom_point(
  size = 0.85,
  alpha = 0.85,
  stroke = 0
) +
  
  # ----------------------------------------------------------
# Paleta Scanpy / Matplotlib
# ----------------------------------------------------------

scale_color_manual(
  values = celltype_colors,
  drop = FALSE
) +
  
  # ----------------------------------------------------------
# Misma escala X/Y
# ----------------------------------------------------------

coord_equal() +
  
  # ----------------------------------------------------------
# Solo título principal
# ----------------------------------------------------------

labs(
  
  title =
    "Yoshida et al. Airway COVID-19",
  
  x =
    NULL,
  
  y =
    NULL,
  
  color =
    "Cell type"
) +
  
  # ----------------------------------------------------------
# Estilo limpio
# ----------------------------------------------------------

theme_void(
  base_size = 16
) +
  
  theme(
    
    # --------------------------------------------------------
    # TÍTULO
    # --------------------------------------------------------
    
    plot.title =
      element_text(
        size = 21,
        face = "plain",
        hjust = 0.5,
        color = "black",
        margin = margin(
          b = 12
        )
      ),
    
    # --------------------------------------------------------
    # LEYENDA A LA DERECHA
    # --------------------------------------------------------
    
    legend.position =
      "right",
    
    legend.title =
      element_text(
        size = 16,
        face = "plain",
        color = "black"
      ),
    
    legend.text =
      element_text(
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    legend.key.height =
      unit(
        7,
        "mm"
      ),
    
    legend.spacing.y =
      unit(
        1,
        "mm"
      ),
    
    legend.box.margin =
      margin(
        l = 8
      ),
    
    plot.margin =
      margin(
        10,
        10,
        10,
        10
      )
  ) +
  
  # ----------------------------------------------------------
# Círculos grandes en leyenda
# ----------------------------------------------------------

guides(
  
  color =
    guide_legend(
      
      override.aes = list(
        size = 5.5,
        alpha = 1
      ),
      
      title.position =
        "top"
    )
)


print(
  p_umap
)


# ============================================================
# PARTE B
# PROPORCIONES CELULARES
# ============================================================


# ============================================================
# 9. IMPORTANTE:
#
# Para proporciones usamos sce_age
#
# NO sce_age_filt
#
# porque sce_age_filt ya eliminó combinaciones
# donor x cell type con menos de 20 células.
#
# Eso podría sesgar la composición.
#
# sce_age contiene:
# COVID+
# nasal cavity
# Pediatric / Adult
# antes de ese filtro.
# ============================================================

prop_meta <- as.data.frame(
  colData(sce_age)
) %>%
  
  dplyr::filter(
    
    Cell_type_annotation_level2 %in%
      epithelial_valid
  ) %>%
  
  dplyr::transmute(
    
    donor_id =
      as.character(
        donor_id
      ),
    
    Group =
      as.character(
        Group
      ),
    
    cell_type =
      as.character(
        Cell_type_annotation_level2
      )
  ) %>%
  
  dplyr::mutate(
    
    Group_plot =
      dplyr::case_when(
        
        Group == "Ped" ~
          "Pediatric",
        
        Group == "Adult" ~
          "Adult",
        
        TRUE ~
          Group
      )
  )


# ============================================================
# 10. CONTAR CÉLULAS POR DONANTE
# ============================================================

prop_donor <- prop_meta %>%
  
  dplyr::count(
    donor_id,
    Group_plot,
    cell_type,
    name = "n_cells"
  )


# ============================================================
# 11. AGREGAR CEROS
#
# Si un donor no presenta un tipo celular,
# debe entrar como 0 y no desaparecer.
# ============================================================

donor_groups <- prop_meta %>%
  
  dplyr::distinct(
    donor_id,
    Group_plot
  )


complete_design <- tidyr::crossing(
  
  donor_groups,
  
  cell_type =
    epithelial_valid
)


prop_donor <- complete_design %>%
  
  dplyr::left_join(
    
    prop_donor,
    
    by = c(
      "donor_id",
      "Group_plot",
      "cell_type"
    )
  ) %>%
  
  dplyr::mutate(
    
    n_cells =
      tidyr::replace_na(
        n_cells,
        0L
      )
  )


# ============================================================
# 12. PROPORCIÓN POR DONANTE
#
# Dentro del compartimento epitelial
# ============================================================

prop_donor <- prop_donor %>%
  
  dplyr::group_by(
    donor_id,
    Group_plot
  ) %>%
  
  dplyr::mutate(
    
    total_epithelial =
      sum(
        n_cells
      ),
    
    proportion =
      ifelse(
        
        total_epithelial > 0,
        
        n_cells /
          total_epithelial,
        
        NA_real_
      )
  ) %>%
  
  dplyr::ungroup()


prop_donor$Group_plot <- factor(
  
  prop_donor$Group_plot,
  
  levels = c(
    "Pediatric",
    "Adult"
  )
)


prop_donor$cell_type <- factor(
  
  prop_donor$cell_type,
  
  levels =
    epithelial_valid
)


# ============================================================
# 13. WILCOXON POR TIPO CELULAR
#
# Unidad experimental = DONANTE
# ============================================================

prop_tests_list <- list()

counter <- 1


for (ct in epithelial_valid) {
  
  
  tmp <- prop_donor %>%
    
    dplyr::filter(
      cell_type == ct
    )
  
  
  ped_values <- tmp$proportion[
    tmp$Group_plot ==
      "Pediatric"
  ]
  
  
  adult_values <- tmp$proportion[
    tmp$Group_plot ==
      "Adult"
  ]
  
  
  if (
    sum(!is.na(ped_values)) >= 3 &&
    sum(!is.na(adult_values)) >= 3
  ) {
    
    
    test <- wilcox.test(
      
      adult_values,
      
      ped_values,
      
      exact = FALSE
    )
    
    
    prop_tests_list[[counter]] <- data.frame(
      
      cell_type =
        ct,
      
      n_Pediatric =
        sum(
          !is.na(
            ped_values
          )
        ),
      
      n_Adult =
        sum(
          !is.na(
            adult_values
          )
        ),
      
      mean_Pediatric =
        mean(
          ped_values,
          na.rm = TRUE
        ),
      
      mean_Adult =
        mean(
          adult_values,
          na.rm = TRUE
        ),
      
      median_Pediatric =
        median(
          ped_values,
          na.rm = TRUE
        ),
      
      median_Adult =
        median(
          adult_values,
          na.rm = TRUE
        ),
      
      Adult_minus_Pediatric =
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
      
      stringsAsFactors =
        FALSE
    )
    
    
    counter <- counter + 1
  }
}


prop_tests <- dplyr::bind_rows(
  prop_tests_list
)


# ============================================================
# 14. FDR
# ============================================================

prop_tests <- prop_tests %>%
  
  dplyr::mutate(
    
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      ),
    
    significance =
      dplyr::case_when(
        
        FDR < 0.001 ~
          "***",
        
        FDR < 0.01 ~
          "**",
        
        FDR < 0.05 ~
          "*",
        
        TRUE ~
          "ns"
      )
  ) %>%
  
  dplyr::arrange(
    FDR
  )


cat(
  "\n========================================\n",
  "CELLULAR COMPOSITION TESTS\n",
  "========================================\n"
)


print(
  tibble::as_tibble(
    prop_tests
  ),
  n = Inf
)


# ============================================================
# 15. RESUMEN PARA BARRAS
# ============================================================

prop_summary <- prop_donor %>%
  
  dplyr::group_by(
    Group_plot,
    cell_type
  ) %>%
  
  dplyr::summarise(
    
    mean_proportion =
      mean(
        proportion,
        na.rm = TRUE
      ),
    
    se =
      sd(
        proportion,
        na.rm = TRUE
      ) /
      sqrt(
        sum(
          !is.na(
            proportion
          )
        )
      ),
    
    .groups =
      "drop"
  )


# ============================================================
# 16. COLORES PEDIATRIC / ADULT
#
# Mantener los colores utilizados en el paper
# ============================================================

group_colors <- c(
  
  "Pediatric" =
    "#E9C46A",
  
  "Adult" =
    "#457B9D"
)


# ============================================================
# 17. POSICIONES PARA SIGNIFICANCIA
#
# Solo se mostrará si FDR < 0.05
# ============================================================

annotation_prop <- prop_donor %>%
  
  dplyr::group_by(
    cell_type
  ) %>%
  
  dplyr::summarise(
    
    y_max =
      max(
        proportion,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::left_join(
    prop_tests,
    by = "cell_type"
  ) %>%
  
  dplyr::mutate(
    
    y_annotation =
      y_max +
      0.05
  )


# ============================================================
# 18. GRÁFICA DE PROPORCIONES
#
# Barras = media entre donantes
# Puntos = cada donante
#
# SIN título
# SIN subtítulo
#
# Leyenda a la derecha
# ============================================================

p_prop <- ggplot(
  
  prop_summary,
  
  aes(
    x = cell_type,
    y = mean_proportion,
    fill = Group_plot
  )
  
) +
  
  # ----------------------------------------------------------
# BARRAS
# ----------------------------------------------------------

geom_col(
  
  position =
    position_dodge(
      width = 0.72
    ),
  
  width =
    0.62,
  
  alpha =
    0.82
) +
  
  # ----------------------------------------------------------
# ERROR ESTÁNDAR
# ----------------------------------------------------------

geom_errorbar(
  
  aes(
    
    ymin =
      pmax(
        mean_proportion - se,
        0
      ),
    
    ymax =
      mean_proportion + se
  ),
  
  position =
    position_dodge(
      width = 0.72
    ),
  
  width =
    0.16,
  
  linewidth =
    0.65,
  
  color =
    "black"
) +
  
  # ----------------------------------------------------------
# DONANTES INDIVIDUALES
# ----------------------------------------------------------

geom_point(
  
  data =
    prop_donor,
  
  aes(
    x = cell_type,
    y = proportion,
    color = Group_plot
  ),
  
  inherit.aes =
    FALSE,
  
  position =
    position_jitterdodge(
      jitter.width = 0.055,
      dodge.width = 0.72
    ),
  
  size =
    3,
  
  alpha =
    0.78
) +
  
  # ----------------------------------------------------------
# SIGNIFICANCIA
#
# Solo si FDR < 0.05
# ----------------------------------------------------------

geom_text(
  
  data =
    annotation_prop %>%
    
    dplyr::filter(
      !is.na(FDR),
      FDR < 0.05
    ),
  
  aes(
    x = cell_type,
    y = y_annotation,
    label = significance
  ),
  
  inherit.aes =
    FALSE,
  
  size =
    6
) +
  
  # ----------------------------------------------------------
# COLORES
# ----------------------------------------------------------

scale_fill_manual(
  values = group_colors
) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  # ----------------------------------------------------------
# EJE Y
# ----------------------------------------------------------

scale_y_continuous(
  
  labels =
    scales::percent_format(
      accuracy = 1
    ),
  
  expand =
    expansion(
      mult = c(
        0,
        0.10
      )
    )
) +
  
  # ----------------------------------------------------------
# SIN TÍTULO / SUBTÍTULO
# ----------------------------------------------------------

labs(
  
  x =
    NULL,
  
  y =
    "Cell proportion",
  
  fill =
    NULL,
  
  color =
    NULL
) +
  
  # ----------------------------------------------------------
# ESTILO
# ----------------------------------------------------------

theme_classic(
  base_size = 16
) +
  
  theme(
    
    # --------------------------------------------------------
    # EJE X
    # --------------------------------------------------------
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    # --------------------------------------------------------
    # EJE Y
    # --------------------------------------------------------
    
    axis.text.y =
      element_text(
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    axis.title.y =
      element_text(
        size = 16,
        face = "plain",
        color = "black",
        margin = margin(
          r = 8
        )
      ),
    
    axis.line =
      element_line(
        linewidth = 0.7,
        color = "black"
      ),
    
    axis.ticks =
      element_line(
        linewidth = 0.6,
        color = "black"
      ),
    
    # --------------------------------------------------------
    # LEYENDA A LA DERECHA
    # --------------------------------------------------------
    
    legend.position =
      "right",
    
    legend.text =
      element_text(
        size = 14,
        color = "black"
      ),
    
    legend.key.height =
      unit(
        7,
        "mm"
      ),
    
    legend.spacing.y =
      unit(
        1,
        "mm"
      ),
    
    plot.margin =
      margin(
        10,
        10,
        10,
        10
      )
  ) +
  
  guides(
    
    fill =
      guide_legend(
        order = 1,
        override.aes = list(
          alpha = 1
        )
      ),
    
    color =
      guide_legend(
        order = 1,
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
  )


print(
  p_prop
)


# ============================================================
# 19. DIRECTORIO
# ============================================================

composition_dir <- file.path(
  results_dir,
  "cellular_composition"
)


dir.create(
  composition_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 20. GUARDAR TABLA ESTADÍSTICA
# ============================================================

write.csv(
  
  prop_tests,
  
  file.path(
    composition_dir,
    "Yoshida_epithelial_cell_composition_Wilcoxon.csv"
  ),
  
  row.names =
    FALSE
)


# ============================================================
# 21. GUARDAR UMAP
# ============================================================

ggsave(
  
  filename =
    file.path(
      composition_dir,
      "Yoshida_nasal_epithelial_UMAP_publication.tiff"
    ),
  
  plot =
    p_umap,
  
  width =
    11,
  
  height =
    8,
  
  dpi =
    600,
  
  compression =
    "lzw"
)


# ============================================================
# 22. GUARDAR PROPORCIONES
# ============================================================

ggsave(
  
  filename =
    file.path(
      composition_dir,
      "Yoshida_epithelial_cell_composition_publication.tiff"
    ),
  
  plot =
    p_prop,
  
  width =
    11,
  
  height =
    7,
  
  dpi =
    600,
  
  compression =
    "lzw"
)


# ============================================================
# 23. FIGURA COMBINADA
#
# UMAP | PROPORCIONES
# ============================================================

final_composition <- p_umap |
  p_prop


print(
  final_composition
)


# ============================================================
# 24. GUARDAR FIGURA COMBINADA
# ============================================================

ggsave(
  
  filename =
    file.path(
      composition_dir,
      "Yoshida_UMAP_and_epithelial_composition_publication.tiff"
    ),
  
  plot =
    final_composition,
  
  width =
    19,
  
  height =
    8,
  
  dpi =
    600,
  
  compression =
    "lzw"
)


# ============================================================
# 25. RESUMEN
# ============================================================

cat(
  "\n========================================\n",
  "FIGURAS FINALIZADAS\n",
  "========================================\n"
)


cat(
  "\nUMAP:\n",
  file.path(
    composition_dir,
    "Yoshida_nasal_epithelial_UMAP_publication.tiff"
  ),
  "\n"
)


cat(
  "\nProporciones:\n",
  file.path(
    composition_dir,
    "Yoshida_epithelial_cell_composition_publication.tiff"
  ),
  "\n"
)


cat(
  "\nFigura combinada:\n",
  file.path(
    composition_dir,
    "Yoshida_UMAP_and_epithelial_composition_publication.tiff"
  ),
  "\n"
)