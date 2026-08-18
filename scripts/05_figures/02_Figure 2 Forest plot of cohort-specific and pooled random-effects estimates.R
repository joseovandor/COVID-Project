## ------------------------------------------------------------------------
## Figure 2 - Cohort-resolved forest plot
## XX-large text + separated statistics columns
##
## Features:
##   - Pediatric-enriched first
##   - Adult-enriched second
##   - Explicit numeric Y ordering
##   - Genes ordered by descending |RE_log2FC|
##   - Cohort-specific estimates + 95% CI
##   - Pooled random-effects estimate + 95% CI
##   - Clearly separated columns:
##       RE log2FC [95% CI]
##       FDR
##   - Extra-large text throughout
##   - Extra-large legend
##   - Legend on the right
##   - TIFF 600 dpi
##   - LZW compression
##   - Direct tiff() output
## ------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(grid)

## ============================================================
## PATHS
## ============================================================

TABLES_DIR <- here(
  "results",
  "meta_analysis",
  "tables"
)

OUT_DIR <- here(
  "results",
  "meta_analysis",
  "figures"
)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

## ============================================================
## SETTINGS
## ============================================================

COHORT_ORDER <- c(
  "GSE172274",
  "GSE179277",
  "GSE231409"
)

## Current palette
COHORT_COLOR <- c(
  GSE172274 = "#f0027f",
  GSE179277 = "#bf5b17",
  GSE231409 = "#7fc97f"
)

COHORT_SHAPE <- c(
  GSE172274 = 16,
  GSE179277 = 17,
  GSE231409 = 15
)

## Pooled random-effects colors
DIV_PED   <- "#FCCE24"
DIV_ADULT <- "#3E76BC"

## Vertical separation among cohort estimates
DODGE <- c(
  GSE172274 =  0.20,
  GSE179277 =  0.00,
  GSE231409 = -0.20
)

## ============================================================
## FUNCTIONS
## ============================================================

format_fdr <- function(x) {
  
  case_when(
    
    is.na(x) ~ "",
    
    x < 0.001 ~ formatC(
      x,
      format = "e",
      digits = 1
    ),
    
    TRUE ~ sprintf(
      "%.3f",
      x
    )
  )
}

## ============================================================
## LOAD DATA
## ============================================================

sig <- read.csv(
  file.path(
    TABLES_DIR,
    "Table_meta_significant_genes.csv"
  ),
  stringsAsFactors = FALSE
)

allg <- read.csv(
  file.path(
    TABLES_DIR,
    "Supplementary_meta_all_genes.csv"
  ),
  stringsAsFactors = FALSE
)

## ============================================================
## SELECT COHORT-SPECIFIC COLUMNS
## ============================================================

allg_sub <- allg %>%
  filter(
    SYMBOL %in% sig$SYMBOL
  ) %>%
  select(
    SYMBOL,
    matches(
      paste0(
        "^(",
        paste(
          COHORT_ORDER,
          collapse = "|"
        ),
        ")_(log2FC|SE)$"
      )
    )
  )

## ============================================================
## SELECT META-ANALYSIS COLUMNS
## ============================================================

sig_core <- sig %>%
  select(
    SYMBOL,
    k,
    RE_log2FC,
    RE_CI_lower,
    RE_CI_upper,
    RE_FDR,
    RobustnessClass
  )

## ============================================================
## MERGE
## ============================================================

df <- sig_core %>%
  left_join(
    allg_sub,
    by = "SYMBOL"
  ) %>%
  mutate(
    
    direction = case_when(
      RE_log2FC > 0 ~ "Pediatric-enriched",
      RE_log2FC < 0 ~ "Adult-enriched",
      TRUE          ~ NA_character_
    )
    
  ) %>%
  filter(
    !is.na(direction)
  )

## ============================================================
## PANEL ORDER
## ============================================================

df$direction <- factor(
  df$direction,
  levels = c(
    "Pediatric-enriched",
    "Adult-enriched"
  )
)

## ============================================================
## ORDER WITHIN EACH PANEL
##
## Largest |RE_log2FC| is displayed at the top.
## ============================================================

df <- df %>%
  group_by(
    direction
  ) %>%
  arrange(
    abs(RE_log2FC),
    RE_FDR,
    .by_group = TRUE
  ) %>%
  mutate(
    y = row_number()
  ) %>%
  ungroup()

## ============================================================
## CHECK ORDER
## ============================================================

order_check <- df %>%
  arrange(
    direction,
    desc(y)
  ) %>%
  select(
    direction,
    y,
    SYMBOL,
    RE_log2FC,
    RE_FDR
  )

print(
  order_check,
  n = Inf
)

## ============================================================
## TEXT LABELS
## ============================================================

df <- df %>%
  mutate(
    
    effect_ci = sprintf(
      "%.2f [%.2f, %.2f]",
      RE_log2FC,
      RE_CI_lower,
      RE_CI_upper
    ),
    
    fdr_text = format_fdr(
      RE_FDR
    )
  )

## ============================================================
## COHORT LONG FORMAT
## ============================================================

cohort_long <- df %>%
  select(
    SYMBOL,
    direction,
    y,
    all_of(
      paste0(
        COHORT_ORDER,
        "_log2FC"
      )
    ),
    all_of(
      paste0(
        COHORT_ORDER,
        "_SE"
      )
    )
  ) %>%
  
  pivot_longer(
    cols = -c(
      SYMBOL,
      direction,
      y
    ),
    
    names_to = c(
      "cohort",
      ".value"
    ),
    
    names_pattern =
      "(GSE\\d+)_(log2FC|SE)"
  ) %>%
  
  filter(
    !is.na(log2FC),
    !is.na(SE)
  ) %>%
  
  mutate(
    
    cohort = factor(
      cohort,
      levels = COHORT_ORDER
    ),
    
    lo = log2FC - 1.96 * SE,
    
    hi = log2FC + 1.96 * SE,
    
    y_dodge =
      y +
      DODGE[
        as.character(cohort)
      ]
  )

## ============================================================
## X RANGE
## ============================================================

x_data_abs <- max(
  abs(
    c(
      cohort_long$lo,
      cohort_long$hi,
      df$RE_CI_lower,
      df$RE_CI_upper
    )
  ),
  na.rm = TRUE
)

forest_abs <- ceiling(
  x_data_abs
)

## ============================================================
## COLUMN POSITIONS
##
## More separation between effect and FDR columns
## ============================================================

GENE_X <- -forest_abs - 0.85

STAT_EFFECT_X <- forest_abs + 0.95

STAT_FDR_X <- forest_abs + 5.90

TOTAL_X_MIN <- -forest_abs - 3.10

TOTAL_X_MAX <- forest_abs + 7.90

## ============================================================
## COLUMN HEADER DATA
## ============================================================

header_df <- df %>%
  group_by(
    direction
  ) %>%
  summarise(
    y_header = max(y) + 1.20,
    .groups = "drop"
  )

## ============================================================
## ALTERNATING ROWS
## ============================================================

row_background <- df %>%
  filter(
    y %% 2 == 0
  )

## ============================================================
## FACET LABELS
## ============================================================

facet_names <- c(
  
  "Pediatric-enriched" =
    "Pediatric-enriched genes",
  
  "Adult-enriched" =
    "Adult-enriched genes"
)

## ============================================================
## PLOT
## ============================================================

p <- ggplot() +
  
  ## ----------------------------------------------------------
## Alternating row backgrounds
## ----------------------------------------------------------

geom_rect(
  data = row_background,
  
  aes(
    xmin = -Inf,
    xmax = Inf,
    ymin = y - 0.47,
    ymax = y + 0.47
  ),
  
  inherit.aes = FALSE,
  
  fill = "#F7F7F7",
  color = NA
) +
  
  ## ----------------------------------------------------------
## Null line
## ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  color = "#8A8A8A",
  linewidth = 0.50,
  linetype = "22"
) +
  
  ## ----------------------------------------------------------
## Gene names
## ----------------------------------------------------------

geom_text(
  data = df,
  
  aes(
    x = GENE_X,
    y = y,
    label = SYMBOL
  ),
  
  inherit.aes = FALSE,
  
  hjust = 1,
  
  family = "sans",
  fontface = "italic",
  
  size = 6.00,
  
  color = "#111111"
) +
  
  ## ----------------------------------------------------------
## Cohort-specific confidence intervals
## ----------------------------------------------------------

geom_errorbarh(
  data = cohort_long,
  
  aes(
    xmin = lo,
    xmax = hi,
    y = y_dodge,
    color = cohort
  ),
  
  height = 0.07,
  
  linewidth = 0.65,
  
  alpha = 0.68
) +
  
  ## ----------------------------------------------------------
## Cohort-specific estimates
## ----------------------------------------------------------

geom_point(
  data = cohort_long,
  
  aes(
    x = log2FC,
    y = y_dodge,
    color = cohort,
    shape = cohort
  ),
  
  size = 2.95,
  
  alpha = 0.96
) +
  
  ## ----------------------------------------------------------
## Pooled random-effects CI
## ----------------------------------------------------------

geom_errorbarh(
  data = df,
  
  aes(
    xmin = RE_CI_lower,
    xmax = RE_CI_upper,
    y = y
  ),
  
  height = 0.16,
  
  linewidth = 1.05,
  
  color = "#151515"
) +
  
  ## ----------------------------------------------------------
## Pooled random-effects diamond
## ----------------------------------------------------------

geom_point(
  data = df,
  
  aes(
    x = RE_log2FC,
    y = y,
    fill = direction
  ),
  
  shape = 23,
  
  size = 5.40,
  
  stroke = 0.95,
  
  color = "#111111"
) +
  
  ## ----------------------------------------------------------
## Separator before statistics
## ----------------------------------------------------------

geom_vline(
  xintercept = forest_abs + 0.50,
  
  color = "#D0D0D0",
  
  linewidth = 0.48
) +
  
  ## ----------------------------------------------------------
## RE effect + CI values
## ----------------------------------------------------------

geom_text(
  data = df,
  
  aes(
    x = STAT_EFFECT_X,
    y = y,
    label = effect_ci
  ),
  
  inherit.aes = FALSE,
  
  hjust = 0,
  
  size = 5.55,
  
  color = "#222222"
) +
  
  ## ----------------------------------------------------------
## FDR values
## ----------------------------------------------------------

geom_text(
  data = df,
  
  aes(
    x = STAT_FDR_X,
    y = y,
    label = fdr_text
  ),
  
  inherit.aes = FALSE,
  
  hjust = 0,
  
  size = 5.55,
  
  color = "#222222"
) +
  
  ## ----------------------------------------------------------
## Header: RE effect
## ----------------------------------------------------------

geom_text(
  data = header_df,
  
  aes(
    x = STAT_EFFECT_X,
    y = y_header,
    label = "RE log2FC [95% CI]"
  ),
  
  inherit.aes = FALSE,
  
  hjust = 0,
  
  size = 5.95,
  
  fontface = "bold",
  
  color = "#111111"
) +
  
  ## ----------------------------------------------------------
## Header: FDR
## ----------------------------------------------------------

geom_text(
  data = header_df,
  
  aes(
    x = STAT_FDR_X,
    y = y_header,
    label = "FDR"
  ),
  
  inherit.aes = FALSE,
  
  hjust = 0,
  
  size = 5.95,
  
  fontface = "bold",
  
  color = "#111111"
) +
  
  ## ============================================================
## FACETS
## ============================================================

facet_wrap(
  ~ direction,
  
  ncol = 1,
  
  scales = "free_y",
  
  labeller = as_labeller(
    facet_names
  )
) +
  
  ## ============================================================
## COLORS / SHAPES
## ============================================================

scale_color_manual(
  name = "Cohort",
  values = COHORT_COLOR,
  breaks = COHORT_ORDER
) +
  
  scale_shape_manual(
    name = "Cohort",
    values = COHORT_SHAPE,
    breaks = COHORT_ORDER
  ) +
  
  scale_fill_manual(
    values = c(
      "Pediatric-enriched" = DIV_PED,
      "Adult-enriched" = DIV_ADULT
    ),
    
    guide = "none"
  ) +
  
  ## ============================================================
## Y AXIS
## ============================================================

scale_y_continuous(
  
  breaks = NULL,
  
  expand = expansion(
    add = c(
      0.65,
      1.75
    )
  )
) +
  
  ## ============================================================
## X AXIS
## ============================================================

scale_x_continuous(
  
  limits = c(
    TOTAL_X_MIN,
    TOTAL_X_MAX
  ),
  
  breaks = pretty(
    c(
      -forest_abs,
      forest_abs
    ),
    n = 7
  ),
  
  expand = expansion(
    mult = c(
      0,
      0
    )
  )
) +
  
  ## ============================================================
## LABELS
## ============================================================

labs(
  x = expression(
    log[2] *
      " fold change (Pediatric vs Adult)"
  ),
  
  y = NULL
) +
  
  ## ============================================================
## THEME
## ============================================================

theme_classic(
  base_size = 17,
  base_family = "sans"
) +
  
  theme(
    
    ## ----------------------------------------------------------
    ## Facet titles
    ## ----------------------------------------------------------
    
    strip.background = element_rect(
      fill = "#EFEFEF",
      color = NA
    ),
    
    strip.text = element_text(
      size = 19,
      face = "bold",
      color = "#222222",
      margin = margin(
        t = 10,
        b = 10
      )
    ),
    
    ## ----------------------------------------------------------
    ## Y axis hidden
    ## ----------------------------------------------------------
    
    axis.text.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.line.y = element_blank(),
    
    ## ----------------------------------------------------------
    ## X axis
    ## ----------------------------------------------------------
    
    axis.text.x = element_text(
      size = 15,
      color = "#333333"
    ),
    
    axis.title.x = element_text(
      size = 17,
      color = "#111111",
      margin = margin(
        t = 13
      )
    ),
    
    axis.ticks.x = element_line(
      linewidth = 0.52,
      color = "#333333"
    ),
    
    axis.line.x = element_line(
      linewidth = 0.62,
      color = "#333333"
    ),
    
    ## ----------------------------------------------------------
    ## LEGEND
    ## ----------------------------------------------------------
    
    legend.position = "right",
    
    legend.direction = "vertical",
    
    legend.title = element_text(
      face = "bold",
      size = 18
    ),
    
    legend.text = element_text(
      size = 16.5
    ),
    
    legend.key.height = unit(
      0.90,
      "cm"
    ),
    
    legend.key.width = unit(
      0.95,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.12,
      "cm"
    ),
    
    ## ----------------------------------------------------------
    ## Spacing
    ## ----------------------------------------------------------
    
    panel.spacing.y = unit(
      0.85,
      "cm"
    ),
    
    plot.margin = margin(
      t = 16,
      r = 20,
      b = 16,
      l = 16
    )
  ) +
  
  ## ============================================================
## LEGEND GUIDES
## ============================================================

guides(
  
  color = guide_legend(
    order = 1,
    
    override.aes = list(
      alpha = 1,
      size = 4.8
    )
  ),
  
  shape = guide_legend(
    order = 1,
    
    override.aes = list(
      alpha = 1,
      size = 4.8
    )
  )
)

## ============================================================
## DISPLAY
## ============================================================

p

## ============================================================
## SAVE TIFF DIRECTLY
## ============================================================

FIG_WIDTH_IN <- 16.5

FIG_HEIGHT_IN <- max(
  15.0,
  0.335 * nrow(df) + 4.0
)

tiff(
  filename = file.path(
    OUT_DIR,
    "Figure_2_Forest_Final_XXLargeText.tif"
  ),
  
  width = FIG_WIDTH_IN,
  
  height = FIG_HEIGHT_IN,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw",
  
  pointsize = 17
)

print(p)

dev.off()

## ============================================================
## FINAL MESSAGES
## ============================================================

message(
  "Figure completed with ",
  nrow(df),
  " genes."
)

message(
  "Pediatric-enriched displayed first."
)

message(
  "Adult-enriched displayed second."
)

message(
  "Genes ordered independently by descending |RE_log2FC|."
)

message(
  "Extra-large text and legend applied."
)

message(
  "Effect and FDR columns separated."
)

message(
  "Saved at: ",
  file.path(
    OUT_DIR,
    "Figure_2_Forest_Final_XXLargeText.tif"
  )
)