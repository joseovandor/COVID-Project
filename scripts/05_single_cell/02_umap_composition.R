# UMAP of the nasal epithelial cell types and cell-type proportions by age group
#
# Proportions are computed per donor before the >= 20-cell filter and compared
# with the Wilcoxon test (BH correction).
#
# Run after 01_import_yoshida.R in the same R session.
# Output: results/scRNAseq_COVID_nasal_age/cellular_composition/

library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)


# 1. EPITHELIAL CELL TYPES ----

# Level-2 epithelial labels used in the Yoshida et al. object. "Ciliated" and
# "Deuterosomal" were missing from the previous list (which used the label
# "Multiciliated"), so ciliated cells had been silently excluded.
epithelial_celltypes <- c(
  "Basal",
  "Basal cycling",
  "Ciliated",
  "Deuterosomal",
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


# 2. ONLY EPITHELIAL TYPES THAT WERE VALID ----
#    IN THE PRIMARY ANALYSIS
# valid_celltypes already contains the cell types
# that had >=3 donors in both groups.

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
  "\nNumber of cell types:",
  length(epithelial_valid),
  "\n"
)


# 3. SCANPY / MATPLOTLIB TAB10 PALETTE ----
# Classic Matplotlib categorical colors

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


# If you have more than 10 populations:
# we would use tab20.
# For this case there should be <=10.

if (length(epithelial_valid) > length(scanpy_tab10)) {
  # Matplotlib tab20 when there are more than 10 cell types
  scanpy_tab10 <- c(
    "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
    "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
    "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
    "#17BECF", "#9EDAE5"
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


# PART A ----
# UMAP


# 4. CELLS FOR UMAP
#
# sce_age_filt contains:
# - COVID+
# - nasal cavity
# - Pediatric / Adult
# - donor x cell type combinations with >=20 cells

keep_umap <- sce_age_filt$Cell_type_annotation_level2 %in%
  epithelial_valid


sce_umap_epi <- sce_age_filt[
  ,
  keep_umap
]


cat(
  "\nCells included in UMAP:",
  ncol(sce_umap_epi),
  "\n"
)


# 5. CHECK UMAP ----

if (!"X_umap" %in% reducedDimNames(sce_umap_epi)) {
  
  stop(
    "X_umap does not exist in reducedDims(sce_umap_epi)."
  )
}


# 6. EXTRACT COORDINATES ----

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


# 7. CELL TYPE ORDER ----

umap_df$cell_type <- factor(
  
  umap_df$cell_type,
  
  levels =
    epithelial_valid
)


# 8. UMAP IN SCANPY / NATURE STYLE ----
# A single UMAP:
# Pediatric + Adult together

p_umap <- ggplot(
  
  umap_df,
  
  aes(
    x = UMAP1,
    y = UMAP2,
    color = cell_type
  )
  
) +
  
  # Points: ----
# relatively small to preserve structure,
# but visible at 600 dpi

geom_point(
  size = 0.85,
  alpha = 0.85,
  stroke = 0
) +
  
  # Scanpy / Matplotlib palette ----

scale_color_manual(
  values = celltype_colors,
  drop = FALSE
) +
  
  # Same X/Y scale ----

coord_equal() +
  
  # Main title only ----

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
  
  # Clean style ----

theme_void(
  base_size = 16
) +
  
  theme(
    
    # TITLE ----
    
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
    
    # LEGEND ON THE RIGHT ----
    
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
  
  # Large circles in legend ----

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


# PART B ----
# CELL PROPORTIONS


# 9. IMPORTANT:
#
# For proportions we use sce_age
#
# NOT sce_age_filt
#
# because sce_age_filt already removed
# donor x cell type combinations with fewer than 20 cells.
#
# That could bias the composition.
#
# sce_age contains:
# COVID+
# nasal cavity
# Pediatric / Adult
# before that filter.

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


# 10. COUNT CELLS PER DONOR ----

prop_donor <- prop_meta %>%
  
  dplyr::count(
    donor_id,
    Group_plot,
    cell_type,
    name = "n_cells"
  )


# 11. ADD ZEROS ----
# If a donor does not have a given cell type,
# it must be entered as 0 rather than dropped.

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


# 12. PROPORTION PER DONOR ----
# Within the epithelial compartment

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


# 13. WILCOXON PER CELL TYPE ----
# Experimental unit = DONOR

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


# 14. FDR ----

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


# 15. SUMMARY FOR BAR CHART ----

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


# 16. PEDIATRIC / ADULT COLORS ----
# Keep the colors used in the paper

group_colors <- c(
  
  "Pediatric" =
    "#E9C46A",
  
  "Adult" =
    "#457B9D"
)


# 17. POSITIONS FOR SIGNIFICANCE ----
# Only shown if FDR < 0.05

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


# 18. PROPORTIONS PLOT
#
# Bars = mean across donors
# Points = each donor
#
# NO title
# NO subtitle
#
# Legend on the right

p_prop <- ggplot(
  
  prop_summary,
  
  aes(
    x = cell_type,
    y = mean_proportion,
    fill = Group_plot
  )
  
) +
  
  # BARS ----

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
  
  # STANDARD ERROR ----

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
  
  # INDIVIDUAL DONORS ----

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
  
  # SIGNIFICANCE ----
# Only if FDR < 0.05

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
  
  # COLORS ----

scale_fill_manual(
  values = group_colors
) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  # Y AXIS ----

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
  
  # NO TITLE / SUBTITLE ----

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
  
  # STYLE ----

theme_classic(
  base_size = 16
) +
  
  theme(
    
    # X AXIS ----
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    # Y AXIS ----
    
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
    
    # LEGEND ON THE RIGHT ----
    
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


# 19. DIRECTORY ----

composition_dir <- file.path(
  results_dir,
  "cellular_composition"
)


dir.create(
  composition_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 20. SAVE STATISTICS TABLE ----

write.csv(
  
  prop_tests,
  
  file.path(
    composition_dir,
    "Yoshida_epithelial_cell_composition_Wilcoxon.csv"
  ),
  
  row.names =
    FALSE
)


# Data behind the UMAP and composition panels (used by scripts/06_figures_tables)
write.csv(
  umap_df,
  gzfile(file.path(composition_dir, "Yoshida_epithelial_UMAP_coordinates.csv.gz")),
  row.names = FALSE
)

write.csv(
  prop_donor,
  file.path(composition_dir, "Yoshida_epithelial_proportions_by_donor.csv"),
  row.names = FALSE
)


# 21. SAVE UMAP ----

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


# 22. SAVE PROPORTIONS ----

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


# 23. COMBINED FIGURE ----
# UMAP | PROPORTIONS

final_composition <- p_umap |
  p_prop


print(
  final_composition
)


# 24. SAVE COMBINED FIGURE ----

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


# 25. SUMMARY ----

cat(
  "\n========================================\n",
  "FIGURES FINALIZED\n",
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
  "\nProportions:\n",
  file.path(
    composition_dir,
    "Yoshida_epithelial_cell_composition_publication.tiff"
  ),
  "\n"
)


cat(
  "\nCombined figure:\n",
  file.path(
    composition_dir,
    "Yoshida_UMAP_and_epithelial_composition_publication.tiff"
  ),
  "\n"
)
