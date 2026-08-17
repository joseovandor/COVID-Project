# Author: Jose A Ovando-Ricardez
# Script for downloading, processing, and annotating GSE179277 dataset
# REVISED VERSION -- changes tied to Reviewer 1 comments are marked with [R#]

# Load necessary libraries
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Step 1: Set up directories for data storage
dir.create("R_GSE179277_Mick_Langelier/")
setwd("R_GSE179277_Mick_Langelier/")
dir.create("Samples/")

# Step 2: Download GEO supplemental files for GSE179277
getGEOSuppFiles("GSE179277", makeDirectory = TRUE, baseDir = "Samples/")
untar("Samples/GSE179277/GSE179277_RAW.tar", exdir = "Samples")

# Step 3: Load GEO dataset for metadata processing
gse <- getGEO("GSE179277", GSEMatrix = TRUE)
gse_data <- gse[[1]]
metadata <- pData(gse_data)

# Step 4: Filter and rename relevant metadata columns
# [R2] Sex is kept and will now be used as a covariate, not just descriptive metadata
filtered_metadata <- metadata[, c("geo_accession",
                                  "title",
                                  "infection:ch1",
                                  "age:ch1",
                                  "gender:ch1")]

colnames(filtered_metadata) <- c("GEO_ID",
                                 "Sample",
                                 "Status",
                                 "Age",
                                 "Sex")

# Step 5: Standardize the 'Status' column
filtered_metadata$Status <- ifelse(grepl("SC2", filtered_metadata$Status, ignore.case = TRUE),
                                   "SC2_Pos",
                                   ifelse(grepl("no_virus", filtered_metadata$Status, ignore.case = TRUE),
                                          "SC2_Neg", filtered_metadata$Status))

# Step 6: Process raw CEL files / raw counts
raw_data <- read.csv("Samples/GSE179277/GSE179277_all_adult_and_ped_samples_combined_counts_unfiltered.csv.gz")
rownames(raw_data) <- raw_data$X
raw_data$X <- NULL
raw_data$gene_name <- NULL

filtered_metadata <- filtered_metadata[filtered_metadata$Status %in% c("SC2_Pos", "SC2_Neg"), ]

filtered_metadata <- filtered_metadata %>%
  separate(Sample, into = c("Sample1", "Sample2", "Group"),
           sep = "_", extra = "merge",
           remove = FALSE)

filtered_metadata <- filtered_metadata %>%
  mutate(Sample = paste(Sample1, Sample2, sep = "_"))
filtered_metadata$Sample1 <- NULL
filtered_metadata$Sample2 <- NULL

filtered_metadata <- filtered_metadata %>%
  mutate(Status = paste(Status, Group, sep = "_"))

filtered_metadata <- filtered_metadata[filtered_metadata$Status %in% c("SC2_Pos_Peds", "SC2_Pos_Adult"), ]

raw_data <- raw_data[, colnames(raw_data) %in% filtered_metadata$Sample]

column_order <- filtered_metadata$Sample
raw_data <- raw_data[, column_order]

rownames(filtered_metadata) <- filtered_metadata$Sample

filtered_metadata$Status <- as.factor(filtered_metadata$Status)
filtered_metadata$Status <- relevel(filtered_metadata$Status, ref = "SC2_Pos_Adult")

# [R2] Standardize Sex as a factor for use as a covariate
filtered_metadata$Sex <- as.factor(filtered_metadata$Sex)

# ---------------------------------------------------------------------------
# [R4] Low-count filtering, explicitly reported.
# Keep genes with at least 10 counts in at least as many samples as the
# smallest group (standard DESeq2 vignette recommendation).
# ---------------------------------------------------------------------------
min_group_size <- min(table(filtered_metadata$Status))
keep <- rowSums(raw_data >= 10) >= min_group_size
raw_data_filtered <- raw_data[keep, ]
cat("Genes before filtering:", nrow(raw_data), "\n")
cat("Genes after low-count filtering:", nrow(raw_data_filtered), "\n")

# ---------------------------------------------------------------------------
# [R2] Design formula now adjusts for Sex as a covariate alongside Status.
# Age was not included as a covariate here because it is continuous and
# collinear with the Status grouping in this dataset (pediatric vs adult
# is itself an age-based split) -- report this reasoning in the methods text.
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = raw_data_filtered,
                               colData = filtered_metadata,
                               design = ~ Sex + Status)

# [R4] Normalization: DESeq2's median-of-ratios method is applied
# automatically inside DESeq() via estimateSizeFactors(). Reported explicitly
# here for the methods section instead of leaving it implicit.
dds <- DESeq(dds)

# [R4] Dispersion estimation diagnostic -- save for supplementary material
dir.create("Results/", showWarnings = FALSE)
tiff("Results/DispEsts_GSE179277.tiff", res = 600, width = 2500, height = 2000, compression = "lzw")
plotDispEsts(dds)
dev.off()

# ---------------------------------------------------------------------------
# PCA -- [R22] percentVar now reported on the plot axes, not just computed
# ---------------------------------------------------------------------------
rld <- vst(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = "Status", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

tiff("Results/PCA_R_GSE179277_Mick_Langelier_Eng.tiff", res = 600, width = 3700, height = 2000, compression = "lzw")
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(aes(colour = Status), size = 3) +
  ggtitle("GSE179277") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.4, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  ) +
  scale_color_manual(values = c("#3E76BC", "#FCCE24"), labels = c("SARS-CoV-2 + Adult", "SARS-CoV-2 + Pediatric")) +
  labs(color = "Status")
dev.off()

# ---------------------------------------------------------------------------
# [R3] Results with explicit FDR control (Benjamini-Hochberg, DESeq2 default)
# alpha = 0.05 sets the significance cutoff used for independent filtering,
# which optimizes power for the padj threshold you intend to report.
# [R4] lfcShrink applied (apeglm) for more reliable effect-size estimates,
# recommended over raw MLE log2FC for ranking/reporting.
# ---------------------------------------------------------------------------
res <- results(dds, contrast = c("Status", "SC2_Pos_Peds", "SC2_Pos_Adult"),
                alpha = 0.05)

# Coefficient name for apeglm shrinkage must match resultsNames(dds)
print(resultsNames(dds))
res_shrunk <- lfcShrink(dds, coef = "Status_SC2_Pos_Peds_vs_SC2_Pos_Adult",
                         type = "apeglm")

# ---------------------------------------------------------------------------
# [R4] Normalization -- report size factors used (median-of-ratios method)
# ---------------------------------------------------------------------------
size_factors <- sizeFactors(dds)
write.csv(data.frame(Sample = names(size_factors), SizeFactor = size_factors),
          "Results/sizeFactors_GSE179277.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# [R4] Dispersion estimation -- report which fit type was actually used
# (parametric is the DESeq2 default; it falls back automatically to "local"
# if the parametric fit fails, which changes how dispersions are shrunk)
# ---------------------------------------------------------------------------
cat("Dispersion fit type used:", dispersionFunction(dds)@fitType, "\n")

# ---------------------------------------------------------------------------
# [R4] Independent filtering -- report the FDR-optimized filtering threshold
# that results() applied (mean normalized count cutoff below which tests
# are excluded before BH correction, improving power for the chosen alpha)
# ---------------------------------------------------------------------------
filter_threshold <- metadata(res)$filterThreshold
cat("Independent filtering: mean count threshold =",
    round(as.numeric(filter_threshold), 3),
    "(", names(filter_threshold), "percentile )\n")

# ---------------------------------------------------------------------------
# [R4] Outlier handling -- DESeq2 flags outliers via Cook's distance and sets
# their padj to NA automatically (cooksCutoff, default TRUE). Report how many
# genes were affected so it can be stated explicitly in the methods text.
# ---------------------------------------------------------------------------
n_outlier_flagged <- sum(is.na(res$padj) & !is.na(res$pvalue))
cat("Genes flagged as outliers by Cook's distance (padj set to NA):",
    n_outlier_flagged, "\n")

table_CD <- data.frame(res_shrunk)
table_CD <- rownames_to_column(table_CD, "ENSEMBL")

# ---------------------------------------------------------------------------
# [R4] Gene annotation: AnnotationDbi::select can return multiple SYMBOL rows
# per ENSEMBL ID (one-to-many mapping), which silently duplicates DEGs
# downstream. Using mapIds() instead enforces a single best match per gene.
# ---------------------------------------------------------------------------
table_CD$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                          keys = table_CD$ENSEMBL,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")

table_CD <- drop_na(table_CD)

# [R3] DEGs should now be defined downstream using padj (BH-adjusted),
# e.g.: subset(table_CD, padj < 0.05 & abs(log2FoldChange) > 0.5)
# rather than the unadjusted pvalue used in the previous version.

write.csv(table_CD, "Results/DE_R_GSE179277_Mick_Langelier_revised.csv", row.names = FALSE)

# [R15] Save session info for reproducibility reporting
writeLines(capture.output(sessionInfo()), "Results/sessionInfo_GSE179277.txt")

# [R4] Save all workflow diagnostics together for direct citation in Methods
writeLines(c(
  paste("Genes before low-count filtering:", nrow(raw_data)),
  paste("Genes after low-count filtering:", nrow(raw_data_filtered)),
  paste("Design formula:", deparse(design(dds))),
  paste("Dispersion fit type:", dispersionFunction(dds)@fitType),
  paste("Independent filtering threshold (mean count):",
        round(as.numeric(filter_threshold), 3),
        "at", names(filter_threshold), "percentile"),
  paste("Genes flagged as outliers (Cook's distance, padj = NA):", n_outlier_flagged),
  paste("Log2FC shrinkage method: apeglm"),
  paste("FDR method: Benjamini-Hochberg (padj, DESeq2 default), alpha = 0.05")
), "Results/DESeq2_workflow_summary_GSE179277.txt")
