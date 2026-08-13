#############################################################
# Forest plot of significant genes from the REM meta-analysis
# - Gene symbols in italics
# - Colored by direction (same palette as Figure 1: Adult = #3E76BC,
#   Pediatric = #FCCE24)
# - Legend placed outside the plot, to the right
# - All plot text in English
#############################################################

if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
library(metafor)

# ---- 1. Load results ----
res <- read.csv("meta_analysis_REM_significant_log2FC_I2.csv", stringsAsFactors = FALSE)

# ---- 2. Sort by effect size ----
res <- res[order(res$estimate), ]

# ---- 3. Palette consistent with Figure 1 ----
COLOR_ADULT     <- "#3E76BC"
COLOR_PEDIATRIC <- "#FCCE24"

# positive = higher in Pediatric; negative = higher in Adult
res$point_col <- ifelse(res$estimate > 0, COLOR_PEDIATRIC, COLOR_ADULT)

# ---- 4. Italic gene symbols (plotmath expressions) ----
slab_expr <- do.call(expression, lapply(res$SYMBOL, function(g) bquote(italic(.(g)))))

# ---- 5. Forest plot ----
# [FIX] Passing col = res$point_col directly to forest() also colored the
# gene labels and the annotate text (log2FC [CI] values), not just the
# points. Instead: draw the whole forest plot in default black/white
# (labels, CI lines, annotate text all black), fix the row y-coordinates
# with "rows", then overlay colored solid points on top with points().
tiff("Forest_plot_REM_significant_genes_final.tiff", res = 600,
     width = 5500, height = 4400, compression = "lzw")

par(mar = c(5, 4, 2, 9))  # bottom, left, top, right (right enlarged for legend)

rows_use <- nrow(res):1  # top-to-bottom row order matching res's row order

forest(x = res$estimate,
       ci.lb = res$ci.lb,
       ci.ub = res$ci.ub,
       slab = slab_expr,
       rows = rows_use,
       psize = 1.3,
       col = "black",              # keep points/CI lines black for now...
       xlab = "Combined log2 Fold Change (REM)",
       refline = 0,
       header = c("Gene", "log2FC [95% CI]"),
       cex = 0.85)

# ...then overlay colored solid points on top -- labels/CI-lines/annotate
# text underneath stay black, only the point markers show color
points(res$estimate, rows_use, pch = 19, col = res$point_col, cex = 1.3)

# Legend outside the plot, to the right
par(xpd = TRUE)
usr <- par("usr")
legend(x = usr[2] + 0.02 * diff(usr[1:2]),
       y = mean(usr[3:4]),
       legend = c("Adult", "Pediatric"),
       pch = 19,
       col = c(COLOR_ADULT, COLOR_PEDIATRIC),
       bty = "n",
       cex = 0.85,
       xjust = 0, yjust = 0.5)
par(xpd = FALSE)

dev.off()

cat("Forest plot saved to Forest_plot_REM_significant_genes_final.tiff\n")

#############################################################
# NOTE: if the legend still overlaps the plot or gets clipped when you
# open the .tiff, increase the right margin (par(mar = c(5,4,2,X)),
# raise X) and/or increase the "width" of the tiff() call slightly.
#############################################################
