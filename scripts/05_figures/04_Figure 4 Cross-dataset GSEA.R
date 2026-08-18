## =============================================================================
## FIGURE 4 - CROSS-COHORT GSEA
##
## Layout:
##
##                    COLUMN 1                       COLUMN 2
##
## ROW 1     Pediatric NES heatmap       Pediatric leading-edge
##
## ROW 2     Adult NES heatmap           Adult leading-edge
##
##
## Design:
##   - Pediatric first
##   - Adult second
##   - Full GEO accession names
##   - Wide rectangular heatmap cells
##   - Larger text throughout the entire figure
##   - Centered Pediatric / Adult titles
##   - No bold text
##   - No panel letters
##   - Shared NES color scale
##   - Shared leading-edge x scale
##
## =============================================================================


## =============================================================================
## 0. PACKAGES
## =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(grid)
library(stringr)


## =============================================================================
## 1. PATHS
## =============================================================================

TABLES_DIR <- here(
  "results",
  "gsea",
  "figures_final_v2",
  "tables"
)

LONG_FILE <- here(
  "results",
  "gsea",
  "GSEA_all_cohorts_long.csv"
)

GLOBAL_FILE <- file.path(
  TABLES_DIR,
  "Selected_global_concordant_pathways.csv"
)

LEADING_EDGE_FILE <- file.path(
  TABLES_DIR,
  "Leading_edge_selected_pathways.csv"
)

OUT_DIR <- here(
  "results",
  "gsea",
  "figures_final_v2"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


## =============================================================================
## 2. CHECK FILES
## =============================================================================

required_files <- c(
  LONG_FILE,
  GLOBAL_FILE,
  LEADING_EDGE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing required files:",
      paste(
        missing_files,
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}


## =============================================================================
## 3. COLORS
## =============================================================================

ADULT_COLOR <- "#3E76BC"

PEDIATRIC_COLOR <- "#FCCE24"

NEUTRAL_COLOR <- "#F7F7F7"


## =============================================================================
## 4. DATASET ORDER
## =============================================================================

dataset_order <- c(
  "GSE172274",
  "GSE179277",
  "GSE231409"
)


dataset_labels <- c(
  "GSE172274" = "GSE172274",
  "GSE179277" = "GSE179277",
  "GSE231409" = "GSE231409"
)


## Pediatric first

enrichment_order <- c(
  "Pediatric-enriched",
  "Adult-enriched"
)


## =============================================================================
## 5. LOAD DATA
## =============================================================================

gsea_long <- read.csv(
  LONG_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  
  mutate(
    
    Dataset = as.character(
      Dataset
    ),
    
    Collection = as.character(
      Collection
    ),
    
    pathway = as.character(
      pathway
    ),
    
    NES = as.numeric(
      NES
    ),
    
    padj = as.numeric(
      padj
    )
  )


selected_global <- read.csv(
  GLOBAL_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


leading_edge_selected <- read.csv(
  LEADING_EDGE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## =============================================================================
## 6. VALIDATE REQUIRED COLUMNS
## =============================================================================

required_gsea_columns <- c(
  "Dataset",
  "Collection",
  "pathway",
  "NES",
  "padj"
)

missing_columns <- setdiff(
  required_gsea_columns,
  colnames(gsea_long)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing columns in GSEA_all_cohorts_long.csv:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


required_global_columns <- c(
  "Collection",
  "pathway",
  "PathwayLabel",
  "EnrichmentGroup",
  "Mean_NES"
)

missing_columns <- setdiff(
  required_global_columns,
  colnames(selected_global)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing columns in Selected_global_concordant_pathways.csv:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


required_leading_columns <- c(
  "PathwayLabel",
  "EnrichmentGroup",
  "N_LE_shared_all3"
)

missing_columns <- setdiff(
  required_leading_columns,
  colnames(leading_edge_selected)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing columns in Leading_edge_selected_pathways.csv:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


## =============================================================================
## 7. PREPARE HEATMAP DATA
## =============================================================================

heatmap_data <- gsea_long %>%
  
  semi_join(
    
    selected_global %>%
      select(
        Collection,
        pathway
      ),
    
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  
  left_join(
    
    selected_global %>%
      select(
        Collection,
        pathway,
        PathwayLabel,
        EnrichmentGroup,
        Mean_NES
      ),
    
    by = c(
      "Collection",
      "pathway"
    )
  ) %>%
  
  mutate(
    
    Dataset = factor(
      Dataset,
      levels = dataset_order
    ),
    
    EnrichmentGroup = factor(
      EnrichmentGroup,
      levels = enrichment_order
    ),
    
    Significance = case_when(
      
      !is.na(padj) &
        padj < 0.01 ~ "**",
      
      !is.na(padj) &
        padj < 0.05 ~ "*",
      
      TRUE ~ ""
    ),
    
    ## Wrap long pathway names
    
    PathwayLabel_wrapped = str_wrap(
      PathwayLabel,
      width = 42
    )
  )


if (nrow(heatmap_data) == 0) {
  
  stop(
    "No pathways available for the heatmap."
  )
}


## =============================================================================
## 8. SPLIT HEATMAP DATA
## =============================================================================

heatmap_pediatric <- heatmap_data %>%
  
  filter(
    EnrichmentGroup == "Pediatric-enriched"
  )


heatmap_adult <- heatmap_data %>%
  
  filter(
    EnrichmentGroup == "Adult-enriched"
  )


## =============================================================================
## 9. ORDER PEDIATRIC PATHWAYS
##
## Highest positive Mean_NES toward the top
## =============================================================================

pediatric_order_original <- selected_global %>%
  
  filter(
    EnrichmentGroup == "Pediatric-enriched"
  ) %>%
  
  arrange(
    desc(Mean_NES)
  ) %>%
  
  pull(
    PathwayLabel
  ) %>%
  
  unique()


pediatric_order <- str_wrap(
  pediatric_order_original,
  width = 42
)


heatmap_pediatric$PathwayLabel_wrapped <- factor(
  
  heatmap_pediatric$PathwayLabel_wrapped,
  
  levels = rev(
    pediatric_order
  )
)


## =============================================================================
## 10. ORDER ADULT PATHWAYS
##
## Most negative Mean_NES toward the top
## =============================================================================

adult_order_original <- selected_global %>%
  
  filter(
    EnrichmentGroup == "Adult-enriched"
  ) %>%
  
  arrange(
    Mean_NES
  ) %>%
  
  pull(
    PathwayLabel
  ) %>%
  
  unique()


adult_order <- str_wrap(
  adult_order_original,
  width = 42
)


heatmap_adult$PathwayLabel_wrapped <- factor(
  
  heatmap_adult$PathwayLabel_wrapped,
  
  levels = rev(
    adult_order
  )
)


## =============================================================================
## 11. COMMON NES SCALE
## =============================================================================

valid_nes <- heatmap_data$NES[
  is.finite(
    heatmap_data$NES
  )
]


if (length(valid_nes) == 0) {
  
  stop(
    "No finite NES values available."
  )
}


nes_limit <- max(
  abs(valid_nes),
  na.rm = TRUE
)


if (
  !is.finite(nes_limit) ||
  nes_limit == 0
) {
  
  nes_limit <- 1
}


## Round upward to nearest 0.5

nes_limit_plot <- ceiling(
  nes_limit * 2
) / 2


nes_breaks <- pretty(
  c(
    -nes_limit_plot,
    nes_limit_plot
  ),
  n = 7
)


nes_breaks <- nes_breaks[
  nes_breaks >= -nes_limit_plot &
    nes_breaks <= nes_limit_plot
]


nes_scale <- scale_fill_gradient2(
  
  low = ADULT_COLOR,
  
  mid = NEUTRAL_COLOR,
  
  high = PEDIATRIC_COLOR,
  
  midpoint = 0,
  
  limits = c(
    -nes_limit_plot,
    nes_limit_plot
  ),
  
  breaks = nes_breaks,
  
  oob = scales::squish,
  
  name = "NES"
)


## =============================================================================
## 12. HEATMAP THEME
##
## LARGE TEXT THROUGHOUT
## CENTERED TITLES
## NO BOLD
## =============================================================================

heatmap_theme <- theme_minimal(
  base_size = 15
) +
  
  theme(
    
    panel.grid = element_blank(),
    
    axis.ticks = element_blank(),
    
    axis.line = element_blank(),
    
    axis.title = element_blank(),
    
    ## -------------------------------------------------------------------------
    ## GSE labels
    ## -------------------------------------------------------------------------
    
    axis.text.x = element_text(
      
      size = 16,
      
      colour = "#222222",
      
      angle = 0,
      
      hjust = 0.5,
      
      vjust = 0.5,
      
      margin = margin(
        b = 12
      )
    ),
    
    ## -------------------------------------------------------------------------
    ## Pathway labels
    ## -------------------------------------------------------------------------
    
    axis.text.y = element_text(
      
      size = 14,
      
      colour = "#222222",
      
      lineheight = 1.05,
      
      margin = margin(
        r = 12
      )
    ),
    
    ## -------------------------------------------------------------------------
    ## Pediatric / Adult titles
    ## -------------------------------------------------------------------------
    
    plot.title = element_text(
      
      size = 18,
      
      colour = "#222222",
      
      hjust = 0.5,
      
      margin = margin(
        b = 14
      )
    ),
    
    ## -------------------------------------------------------------------------
    ## NES legend
    ## -------------------------------------------------------------------------
    
    legend.title = element_text(
      size = 14
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.position = "bottom",
    
    plot.margin = margin(
      t = 10,
      r = 14,
      b = 10,
      l = 10
    )
  )


## =============================================================================
## 13. PEDIATRIC HEATMAP
##
## ratio = 0.28 -> wide rectangular cells
## =============================================================================

heatmap_plot_pediatric <- ggplot(
  
  heatmap_pediatric,
  
  aes(
    x = Dataset,
    y = PathwayLabel_wrapped,
    fill = NES
  )
) +
  
  geom_tile(
    
    width = 0.98,
    
    height = 0.76,
    
    colour = "white",
    
    linewidth = 0.5
  ) +
  
  ## Larger significance stars
  
  geom_text(
    
    aes(
      label = Significance
    ),
    
    size = 6,
    
    colour = "#222222"
  ) +
  
  nes_scale +
  
  scale_x_discrete(
    
    position = "top",
    
    labels = dataset_labels,
    
    expand = expansion(
      add = 0.06
    )
  ) +
  
  scale_y_discrete(
    
    expand = expansion(
      add = 0.06
    )
  ) +
  
  labs(
    
    title = "Pediatric-enriched pathways",
    
    x = NULL,
    
    y = NULL
  ) +
  
  coord_fixed(
    
    ratio = 0.28,
    
    clip = "off"
  ) +
  
  heatmap_theme


## =============================================================================
## 14. ADULT HEATMAP
## =============================================================================

heatmap_plot_adult <- ggplot(
  
  heatmap_adult,
  
  aes(
    x = Dataset,
    y = PathwayLabel_wrapped,
    fill = NES
  )
) +
  
  geom_tile(
    
    width = 0.98,
    
    height = 0.76,
    
    colour = "white",
    
    linewidth = 0.5
  ) +
  
  geom_text(
    
    aes(
      label = Significance
    ),
    
    size = 6,
    
    colour = "#222222"
  ) +
  
  nes_scale +
  
  scale_x_discrete(
    
    position = "top",
    
    labels = dataset_labels,
    
    expand = expansion(
      add = 0.06
    )
  ) +
  
  scale_y_discrete(
    
    expand = expansion(
      add = 0.06
    )
  ) +
  
  labs(
    
    title = "Adult-enriched pathways",
    
    x = NULL,
    
    y = NULL
  ) +
  
  coord_fixed(
    
    ratio = 0.28,
    
    clip = "off"
  ) +
  
  heatmap_theme


## =============================================================================
## 15. PREPARE LEADING-EDGE DATA
## =============================================================================

leading_edge_plot <- leading_edge_selected %>%
  
  mutate(
    
    N_LE_shared_all3 = as.numeric(
      N_LE_shared_all3
    ),
    
    EnrichmentGroup = factor(
      EnrichmentGroup,
      levels = enrichment_order
    ),
    
    PathwayLabel_wrapped = str_wrap(
      PathwayLabel,
      width = 42
    )
  ) %>%
  
  filter(
    
    !is.na(
      N_LE_shared_all3
    ),
    
    N_LE_shared_all3 > 0,
    
    !is.na(
      PathwayLabel
    ),
    
    PathwayLabel != ""
  )


if (nrow(leading_edge_plot) == 0) {
  
  stop(
    "No pathways with shared leading-edge genes."
  )
}


## =============================================================================
## 16. SPLIT LEADING-EDGE DATA
## =============================================================================

le_pediatric <- leading_edge_plot %>%
  
  filter(
    EnrichmentGroup == "Pediatric-enriched"
  ) %>%
  
  arrange(
    N_LE_shared_all3
  )


le_adult <- leading_edge_plot %>%
  
  filter(
    EnrichmentGroup == "Adult-enriched"
  ) %>%
  
  arrange(
    N_LE_shared_all3
  )


## =============================================================================
## 17. ORDER LEADING-EDGE PATHWAYS
## =============================================================================

le_pediatric$PathwayLabel_wrapped <- factor(
  
  le_pediatric$PathwayLabel_wrapped,
  
  levels = unique(
    le_pediatric$PathwayLabel_wrapped
  )
)


le_adult$PathwayLabel_wrapped <- factor(
  
  le_adult$PathwayLabel_wrapped,
  
  levels = unique(
    le_adult$PathwayLabel_wrapped
  )
)


## =============================================================================
## 18. COMMON LEADING-EDGE SCALE
## =============================================================================

max_le <- max(
  leading_edge_plot$N_LE_shared_all3,
  na.rm = TRUE
)


leading_x_max <- ceiling(
  max_le * 1.18
)


leading_breaks <- pretty(
  c(
    0,
    max_le
  ),
  n = 5
)


## =============================================================================
## 19. LEADING-EDGE THEME
##
## LARGE TEXT THROUGHOUT
## CENTERED TITLES
## =============================================================================

leading_theme <- theme_minimal(
  base_size = 15
) +
  
  theme(
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      
      colour = "#E5E5E5",
      
      linewidth = 0.45
    ),
    
    axis.title.y = element_blank(),
    
    ## -------------------------------------------------------------------------
    ## X-axis title
    ## -------------------------------------------------------------------------
    
    axis.title.x = element_text(
      
      size = 14,
      
      colour = "#222222",
      
      margin = margin(
        t = 12
      )
    ),
    
    ## -------------------------------------------------------------------------
    ## X-axis numbers
    ## -------------------------------------------------------------------------
    
    axis.text.x = element_text(
      
      size = 13,
      
      colour = "#444444"
    ),
    
    ## -------------------------------------------------------------------------
    ## Pathway names
    ## -------------------------------------------------------------------------
    
    axis.text.y = element_text(
      
      size = 14,
      
      colour = "#222222",
      
      lineheight = 1.05,
      
      margin = margin(
        r = 10
      )
    ),
    
    ## -------------------------------------------------------------------------
    ## Pediatric / Adult title
    ## -------------------------------------------------------------------------
    
    plot.title = element_text(
      
      size = 18,
      
      colour = "#222222",
      
      hjust = 0.5,
      
      margin = margin(
        b = 14
      )
    ),
    
    plot.margin = margin(
      t = 10,
      r = 30,
      b = 10,
      l = 16
    )
  )


## =============================================================================
## 20. PEDIATRIC LEADING EDGE
## =============================================================================

leading_pediatric <- ggplot(
  
  le_pediatric,
  
  aes(
    x = N_LE_shared_all3,
    y = PathwayLabel_wrapped
  )
) +
  
  geom_segment(
    
    aes(
      
      x = 0,
      
      xend = N_LE_shared_all3,
      
      yend = PathwayLabel_wrapped
    ),
    
    colour = PEDIATRIC_COLOR,
    
    linewidth = 1.4,
    
    lineend = "round"
  ) +
  
  geom_point(
    
    colour = PEDIATRIC_COLOR,
    
    size = 6
  ) +
  
  geom_text(
    
    aes(
      label = N_LE_shared_all3
    ),
    
    nudge_x = max_le * 0.03,
    
    size = 5,
    
    colour = "#444444"
  ) +
  
  scale_x_continuous(
    
    limits = c(
      0,
      leading_x_max
    ),
    
    breaks = leading_breaks,
    
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    
    title = "Pediatric-enriched pathways",
    
    x = "Leading-edge genes shared across all three cohorts",
    
    y = NULL
  ) +
  
  leading_theme +
  
  coord_cartesian(
    clip = "off"
  )


## =============================================================================
## 21. ADULT LEADING EDGE
## =============================================================================

leading_adult <- ggplot(
  
  le_adult,
  
  aes(
    x = N_LE_shared_all3,
    y = PathwayLabel_wrapped
  )
) +
  
  geom_segment(
    
    aes(
      
      x = 0,
      
      xend = N_LE_shared_all3,
      
      yend = PathwayLabel_wrapped
    ),
    
    colour = ADULT_COLOR,
    
    linewidth = 1.4,
    
    lineend = "round"
  ) +
  
  geom_point(
    
    colour = ADULT_COLOR,
    
    size = 6
  ) +
  
  geom_text(
    
    aes(
      label = N_LE_shared_all3
    ),
    
    nudge_x = max_le * 0.03,
    
    size = 5,
    
    colour = "#444444"
  ) +
  
  scale_x_continuous(
    
    limits = c(
      0,
      leading_x_max
    ),
    
    breaks = leading_breaks,
    
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    
    title = "Adult-enriched pathways",
    
    x = "Leading-edge genes shared across all three cohorts",
    
    y = NULL
  ) +
  
  leading_theme +
  
  coord_cartesian(
    clip = "off"
  )


## =============================================================================
## 22. NUMBER OF ROWS
## =============================================================================

n_heat_ped <- length(
  unique(
    heatmap_pediatric$PathwayLabel_wrapped
  )
)


n_heat_adult <- length(
  unique(
    heatmap_adult$PathwayLabel_wrapped
  )
)


n_le_ped <- nrow(
  le_pediatric
)


n_le_adult <- nrow(
  le_adult
)


## =============================================================================
## 23. ROW HEIGHTS
##
## Keeps Pediatric and Adult panels aligned between columns
## =============================================================================

row_height_pediatric <- max(
  n_heat_ped,
  n_le_ped,
  1
)


row_height_adult <- max(
  n_heat_adult,
  n_le_adult,
  1
)


## =============================================================================
## 24. TRUE 2 x 2 LAYOUT
##
##                COLUMN 1                    COLUMN 2
##
## ROW 1     Pediatric heatmap          Pediatric leading-edge
##
## ROW 2     Adult heatmap              Adult leading-edge
##
## =============================================================================

layout_design <- "
AB
CD
"


final_figure <- wrap_plots(
  
  A = heatmap_plot_pediatric,
  
  B = leading_pediatric,
  
  C = heatmap_plot_adult,
  
  D = leading_adult,
  
  design = layout_design,
  
  ## More room for the heatmaps
  widths = c(
    1.55,
    1
  ),
  
  heights = c(
    row_height_pediatric,
    row_height_adult
  ),
  
  guides = "collect"
) +
  
  plot_annotation(
    
    caption = paste0(
      
      "* FDR < 0.05; ** FDR < 0.01. ",
      
      "NES, normalized enrichment score."
    )
  ) &
  
  theme(
    
    legend.position = "bottom",
    
    ## Larger caption
    plot.caption = element_text(
      
      size = 12,
      
      colour = "#555555",
      
      hjust = 0,
      
      margin = margin(
        t = 12
      )
    )
  )


## =============================================================================
## 25. LARGER NES LEGEND
## =============================================================================

final_figure <- final_figure &
  
  guides(
    
    fill = guide_colorbar(
      
      title.position = "left",
      
      title.hjust = 0.5,
      
      barwidth = unit(
        10,
        "cm"
      ),
      
      barheight = unit(
        0.65,
        "cm"
      ),
      
      ticks = TRUE
    )
  )


## =============================================================================
## 26. DISPLAY
## =============================================================================

print(
  final_figure
)


## =============================================================================
## 27. EXPORT PLOTTING DATA
## =============================================================================

write.csv(
  
  heatmap_data,
  
  file.path(
    TABLES_DIR,
    "Figure4_heatmap_plotting_data_horizontal_v6.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  leading_edge_plot,
  
  file.path(
    TABLES_DIR,
    "Figure4_leading_edge_plotting_data_horizontal_v6.csv"
  ),
  
  row.names = FALSE
)


## =============================================================================
## 28. SAVE PNG
##
## Increased canvas size because text is larger
## =============================================================================

ggsave(
  
  filename = file.path(
    OUT_DIR,
    "Figure_4_GSEA_horizontal_v6.png"
  ),
  
  plot = final_figure,
  
  width = 20,
  
  height = 11,
  
  units = "in",
  
  dpi = 300,
  
  bg = "white",
  
  limitsize = FALSE
)


## =============================================================================
## 29. SAVE TIFF
## =============================================================================

ggsave(
  
  filename = file.path(
    OUT_DIR,
    "Figure_4_GSEA_horizontal_v6.tif"
  ),
  
  plot = final_figure,
  
  device = "tiff",
  
  width = 20,
  
  height = 11,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white",
  
  limitsize = FALSE
)


## =============================================================================
## 30. FINISHED
## =============================================================================

message(
  paste0(
    
    "\nFigure completed successfully.\n\n",
    
    "ROW 1:\n",
    "  Pediatric heatmap | Pediatric leading-edge\n\n",
    
    "ROW 2:\n",
    "  Adult heatmap | Adult leading-edge\n\n",
    
    "Typography:\n",
    "  Enrichment titles = 18\n",
    "  GSE labels = 16\n",
    "  Pathway labels = 14\n",
    "  Leading-edge X title = 14\n",
    "  Leading-edge X numbers = 13\n",
    "  Significance stars = 6\n",
    "  Leading-edge counts = 5\n",
    "  NES legend title = 14\n",
    "  NES legend values = 13\n",
    "  Caption = 12\n\n",
    
    "Output directory:\n",
    
    OUT_DIR,
    
    "\n"
  )
)