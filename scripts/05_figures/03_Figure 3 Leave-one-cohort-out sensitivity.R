## ------------------------------------------------------------------------
## Figure 3B - Leave-one-cohort-out cross-cohort reproducibility
##
## X = pooled log2FC from the remaining two cohorts
## Y = observed log2FC in the held-out cohort
##
## Statistics:
##   - Pearson r
##   - Lin's CCC
##   - Direction concordance
##
## Background:
##   - Yellow quadrant = Pediatric-enriched agreement
##   - Blue quadrant   = Adult-enriched agreement
##
## Only genes with k = 3 are included.
##
## Output:
##   Figure_3B_LOO_Reproducibility_Quadrants_Final.tif
##
## TIFF:
##   - 600 dpi
##   - LZW compression
##   - direct tiff(), no ggsave()
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

SENS_DIR <- here(
  "results",
  "meta_analysis",
  "sensitivity"
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

## Point colors
CONCORDANT_COLOR <- "#7fc97f"
DISCORDANT_COLOR <- "#D95F0E"

## Background quadrant colors
PED_BG   <- "#FCCE24"
ADULT_BG <- "#3E76BC"

## ============================================================
## LIN'S CONCORDANCE CORRELATION COEFFICIENT
## ============================================================

lin_ccc <- function(x, y) {
  
  ok <- complete.cases(x, y)
  
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) {
    return(NA_real_)
  }
  
  mx <- mean(x)
  my <- mean(y)
  
  vx <- var(x)
  vy <- var(y)
  
  cxy <- cov(x, y)
  
  ccc <- (2 * cxy) /
    (
      vx +
        vy +
        (mx - my)^2
    )
  
  return(ccc)
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

loo <- read.csv(
  file.path(
    SENS_DIR,
    "leave_one_out_results.csv"
  ),
  stringsAsFactors = FALSE
)

## ============================================================
## KEEP ONLY k = 3 GENES
## ============================================================

sig_loo <- sig %>%
  filter(
    k == 3
  ) %>%
  as_tibble()

cat(
  "\nLOO-eligible genes:",
  nrow(sig_loo),
  "/",
  nrow(sig),
  "\n"
)

## ============================================================
## HELD-OUT COHORT EFFECTS
## ============================================================

allg_idx <- allg %>%
  select(
    SYMBOL,
    all_of(
      paste0(
        COHORT_ORDER,
        "_log2FC"
      )
    )
  ) %>%
  distinct() %>%
  as_tibble()

## ============================================================
## BUILD LOO SCATTER DATA
## ============================================================

scatter_df <- lapply(
  COHORT_ORDER,
  function(omit) {
    
    held_col <- paste0(
      omit,
      "_log2FC"
    )
    
    loo %>%
      filter(
        Omitted_dataset == omit,
        SYMBOL %in% sig_loo$SYMBOL
      ) %>%
      
      left_join(
        allg_idx %>%
          select(
            SYMBOL,
            held_out = all_of(held_col)
          ),
        by = "SYMBOL"
      ) %>%
      
      filter(
        !is.na(Pooled_log2FC),
        !is.na(held_out)
      ) %>%
      
      transmute(
        
        SYMBOL,
        
        panel = omit,
        
        discovery = Pooled_log2FC,
        
        held_out = held_out,
        
        concordant =
          sign(Pooled_log2FC) ==
          sign(held_out)
      )
  }
) %>%
  bind_rows() %>%
  as_tibble()

scatter_df$panel <- factor(
  scatter_df$panel,
  levels = COHORT_ORDER
)

## ============================================================
## CHECK NUMBER OF GENES PER PANEL
## ============================================================

panel_n <- scatter_df %>%
  count(
    panel,
    name = "n_genes"
  )

print(
  panel_n
)

## ============================================================
## GLOBAL SYMMETRIC LIMIT
## ============================================================

global_lim <- max(
  abs(
    c(
      scatter_df$discovery,
      scatter_df$held_out
    )
  ),
  na.rm = TRUE
)

global_lim <- ceiling(
  global_lim * 1.10
)

## ============================================================
## PANEL STATISTICS
## ============================================================

panel_stats <- scatter_df %>%
  group_by(
    panel
  ) %>%
  summarise(
    
    pearson_r = cor(
      discovery,
      held_out,
      method = "pearson",
      use = "complete.obs"
    ),
    
    lin_ccc = lin_ccc(
      discovery,
      held_out
    ),
    
    direction_concordance =
      mean(
        concordant,
        na.rm = TRUE
      ) * 100,
    
    n_genes = n(),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    label = sprintf(
      "Pearson r = %.2f\nLin CCC = %.2f\nDirection concordance = %.0f%%",
      pearson_r,
      lin_ccc,
      direction_concordance
    ),
    
    stat_x = -global_lim * 0.92,
    
    stat_y = global_lim * 0.91
  ) %>%
  as_tibble()

print(
  panel_stats
)

## ============================================================
## DISCORDANT GENES
## ============================================================

discordant_df <- scatter_df %>%
  filter(
    !concordant
  ) %>%
  select(
    panel,
    SYMBOL,
    discovery,
    held_out
  ) %>%
  as_tibble()

cat(
  "\nDiscordant genes:\n"
)

print(
  discordant_df,
  n = Inf
)

## ============================================================
## FACET LABELS
## ============================================================

heldout_labels <- c(
  
  GSE172274 =
    "Held out: GSE172274",
  
  GSE179277 =
    "Held out: GSE179277",
  
  GSE231409 =
    "Held out: GSE231409"
)

## ============================================================
## FIGURE
## ============================================================

p_loo <- ggplot(
  scatter_df,
  aes(
    x = discovery,
    y = held_out
  )
) +
  
  ## ----------------------------------------------------------
## Pediatric-enriched concordant quadrant
## ----------------------------------------------------------

annotate(
  "rect",
  xmin = 0,
  xmax = global_lim,
  ymin = 0,
  ymax = global_lim,
  fill = PED_BG,
  alpha = 0.10
) +
  
  ## ----------------------------------------------------------
## Adult-enriched concordant quadrant
## ----------------------------------------------------------

annotate(
  "rect",
  xmin = -global_lim,
  xmax = 0,
  ymin = -global_lim,
  ymax = 0,
  fill = ADULT_BG,
  alpha = 0.10
) +
  
  ## ----------------------------------------------------------
## Pediatric-enriched quadrant label
## moved inward so text is fully visible
## ----------------------------------------------------------

annotate(
  "text",
  x = global_lim * 0.35,
  y = global_lim * 0.77,
  label = "Pediatric-enriched",
  size = 6.1,
  fontface = "bold",
  hjust = 0.5,
  color = "#8A6D00",
  alpha = 0.95
) +
  
  ## ----------------------------------------------------------
## Adult-enriched quadrant label
## ----------------------------------------------------------

annotate(
  "text",
  x = -global_lim * 0.50,
  y = -global_lim * 0.77,
  label = "Adult-enriched",
  size = 6.1,
  fontface = "bold",
  hjust = 0.5,
  color = "#244B7A",
  alpha = 0.95
) +
  
  ## ----------------------------------------------------------
## Horizontal zero line
## ----------------------------------------------------------

geom_hline(
  yintercept = 0,
  color = "#B8B8B8",
  linewidth = 0.75
) +
  
  ## ----------------------------------------------------------
## Vertical zero line
## ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  color = "#B8B8B8",
  linewidth = 0.75
) +
  
  ## ----------------------------------------------------------
## Identity line
## ----------------------------------------------------------

geom_abline(
  slope = 1,
  intercept = 0,
  color = "#686868",
  linewidth = 1.00,
  linetype = "22"
) +
  
  ## ----------------------------------------------------------
## Genes
## ----------------------------------------------------------

geom_point(
  aes(
    fill = concordant,
    shape = concordant
  ),
  size = 5.8,
  stroke = 0.75,
  color = "white",
  alpha = 0.95
) +
  
  ## ----------------------------------------------------------
## Label discordant genes
## ----------------------------------------------------------

geom_text(
  data = discordant_df,
  
  aes(
    x = discovery,
    y = held_out,
    label = SYMBOL
  ),
  
  inherit.aes = FALSE,
  
  nudge_x = 0.22,
  nudge_y = -0.18,
  
  hjust = 0,
  
  size = 5.2,
  
  fontface = "italic",
  
  color = "#8C3B00"
) +
  
  ## ----------------------------------------------------------
## Statistics
## ----------------------------------------------------------

geom_text(
  data = panel_stats,
  
  aes(
    x = stat_x,
    y = stat_y,
    label = label
  ),
  
  inherit.aes = FALSE,
  
  hjust = 0,
  vjust = 1,
  
  size = 5.8,
  
  lineheight = 1.12,
  
  color = "#222222"
) +
  
  ## ============================================================
## POINT COLOR
## ============================================================

scale_fill_manual(
  
  name = NULL,
  
  breaks = c(
    "TRUE",
    "FALSE"
  ),
  
  values = c(
    `TRUE` = CONCORDANT_COLOR,
    `FALSE` = DISCORDANT_COLOR
  ),
  
  labels = c(
    `TRUE` = "Concordant direction",
    `FALSE` = "Discordant direction"
  )
) +
  
  ## ============================================================
## POINT SHAPE
## ============================================================

scale_shape_manual(
  
  name = NULL,
  
  breaks = c(
    "TRUE",
    "FALSE"
  ),
  
  values = c(
    `TRUE` = 21,
    `FALSE` = 24
  ),
  
  labels = c(
    `TRUE` = "Concordant direction",
    `FALSE` = "Discordant direction"
  )
) +
  
  ## ============================================================
## FACETS
## ============================================================

facet_wrap(
  
  ~ panel,
  
  nrow = 1,
  
  scales = "fixed",
  
  labeller = as_labeller(
    heldout_labels
  )
) +
  
  ## ============================================================
## X SCALE
## ============================================================

scale_x_continuous(
  
  limits = c(
    -global_lim,
    global_lim
  ),
  
  breaks = pretty(
    c(
      -global_lim,
      global_lim
    ),
    n = 5
  ),
  
  expand = expansion(
    mult = c(
      0,
      0
    )
  )
) +
  
  ## ============================================================
## Y SCALE
## ============================================================

scale_y_continuous(
  
  limits = c(
    -global_lim,
    global_lim
  ),
  
  breaks = pretty(
    c(
      -global_lim,
      global_lim
    ),
    n = 5
  ),
  
  expand = expansion(
    mult = c(
      0,
      0
    )
  )
) +
  
  ## ============================================================
## AXIS LABELS
## ============================================================

labs(
  
  x = expression(
    "LOO pooled " *
      log[2] *
      "FC (remaining two cohorts)"
  ),
  
  y = expression(
    "Observed " *
      log[2] *
      "FC in held-out cohort"
  )
) +
  
  ## ============================================================
## THEME
## ============================================================

theme_classic(
  base_size = 18,
  base_family = "sans"
) +
  
  theme(
    
    ## Square panels
    aspect.ratio = 1,
    
    ## ----------------------------------------------------------
    ## Facet headers
    ## ----------------------------------------------------------
    
    strip.background = element_rect(
      fill = "#ECECEC",
      color = NA
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 20,
      color = "#222222",
      margin = margin(
        t = 10,
        b = 10
      )
    ),
    
    ## ----------------------------------------------------------
    ## Axis text
    ## ----------------------------------------------------------
    
    axis.text = element_text(
      size = 15.5,
      color = "#333333"
    ),
    
    ## ----------------------------------------------------------
    ## X title
    ## ----------------------------------------------------------
    
    axis.title.x = element_text(
      size = 18,
      color = "#111111",
      margin = margin(
        t = 12
      )
    ),
    
    ## ----------------------------------------------------------
    ## Y title
    ## ----------------------------------------------------------
    
    axis.title.y = element_text(
      size = 18,
      color = "#111111",
      margin = margin(
        r = 12
      )
    ),
    
    ## ----------------------------------------------------------
    ## Axes
    ## ----------------------------------------------------------
    
    axis.line = element_line(
      linewidth = 0.70,
      color = "#222222"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.60,
      color = "#222222"
    ),
    
    ## ----------------------------------------------------------
    ## LEGEND
    ## now below plot
    ## ----------------------------------------------------------
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.justification = "center",
    
    legend.box = "horizontal",
    
    legend.text = element_text(
      size = 16.5,
      color = "#222222"
    ),
    
    legend.key.height = unit(
      0.72,
      "cm"
    ),
    
    legend.key.width = unit(
      0.95,
      "cm"
    ),
    
    legend.spacing.x = unit(
      0.20,
      "cm"
    ),
    
    legend.margin = margin(
      t = 10,
      r = 0,
      b = 0,
      l = 0
    ),
    
    ## ----------------------------------------------------------
    ## Panel spacing
    ## ----------------------------------------------------------
    
    panel.spacing.x = unit(
      0.90,
      "cm"
    ),
    
    ## ----------------------------------------------------------
    ## Plot margins
    ## ----------------------------------------------------------
    
    plot.margin = margin(
      t = 14,
      r = 16,
      b = 14,
      l = 16
    )
  ) +
  
  ## ============================================================
## LEGEND GUIDES
##
## Concordant first, Discordant second
## ============================================================

guides(
  
  fill = guide_legend(
    order = 1,
    
    override.aes = list(
      size = 6.2,
      alpha = 1
    )
  ),
  
  shape = guide_legend(
    order = 1,
    
    override.aes = list(
      size = 6.2,
      alpha = 1
    )
  )
)

## ============================================================
## DISPLAY
## ============================================================

p_loo

## ============================================================
## SAVE TIFF DIRECTLY
## ============================================================

FIG_WIDTH_IN  <- 16.5
FIG_HEIGHT_IN <- 7.2

tiff(
  
  filename = file.path(
    OUT_DIR,
    "Figure_3B_LOO_Reproducibility_Quadrants_Final.tif"
  ),
  
  width = FIG_WIDTH_IN,
  
  height = FIG_HEIGHT_IN,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw",
  
  pointsize = 18
)

print(
  p_loo
)

dev.off()

## ============================================================
## FINAL MESSAGES
## ============================================================

message(
  "LOO reproducibility panel completed."
)

message(
  "Genes evaluated: ",
  nrow(sig_loo),
  " k=3 genes."
)

message(
  "Discordant observations: ",
  nrow(discordant_df)
)

message(
  "Yellow quadrant = Pediatric-enriched concordance."
)

message(
  "Blue quadrant = Adult-enriched concordance."
)

message(
  "Legend positioned below the panels."
)

message(
  "Saved at: ",
  file.path(
    OUT_DIR,
    "Figure_3B_LOO_Reproducibility_Quadrants_Final.tif"
  )
)