# Author: Jose A Ovando-Ricardez
# Script for downloading, processing, and annotating GSE231409 dataset
# REVISED VERSION -- changes tied to Reviewer 1 comments are marked with [R#]

# Load necessary libraries
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Step 1: Set up directories for data storage
dir.create("R_GSE231409_Hurst_Kelly/")
setwd("R_GSE231409_Hurst_Kelly/")
dir.create("Samples/")

# Step 2: Download GEO supplemental files for GSE231409
getGEOSuppFiles("GSE231409", makeDirectory = TRUE, baseDir = "Samples/")
untar("Samples/GSE231409/GSE231409_RAW.tar", exdir = "Samples")

# Step 3: Load GEO dataset for metadata processing
gse <- getGEO("GSE231409", GSEMatrix = TRUE)
gse_data <- gse[[1]]
metadata <- pData(gse_data)

# Step 4: Filter and rename relevant metadata columns
# [R2] Sex and Batch are kept and used as covariates, not just descriptive metadata
filtered_metadata <- metadata[, c("geo_accession",
                                  "title",
                                  "covid:ch1",
                                  "age:ch1",
                                  "Sex:ch1",
                                  "tissue:ch1",
                                  "batch:ch1")]

colnames(filtered_metadata) <- c("GEO_ID",
                                 "Sample",
                                 "Status",
                                 "Age",
                                 "Sex",
                                 "Tissue",
                                 "Batch")

# Step 5: Standardize the 'Status' column into SARS-CoV-2 positive/negative
filtered_metadata$Status <- ifelse(grepl("Positive", filtered_metadata$Status, ignore.case = TRUE),
                                   "SC2_Pos",
                                   ifelse(grepl("Negative", filtered_metadata$Status, ignore.case = TRUE),
                                          "SC2_Neg", filtered_metadata$Status))

# Step 6: Load raw gene-level count matrix (GEO-provided raw counts, not CEL files)
raw_data <- read.csv("Samples/GSE231409/GSE231409_BRAVE_RNASeq_counts.csv.gz")
rownames(raw_data) <- raw_data$target_id
raw_data$target_id <- NULL

filtered_metadata <- filtered_metadata[filtered_metadata$Status %in% c("SC2_Pos", "SC2_Neg"), ]
filtered_metadata <- filtered_metadata[filtered_metadata$Tissue %in% c("Nasal swab", "Nasopharyngeal swab"), ]

# [R19/R30] Group is labeled "Pediatric" (not "Peds") to match the
# terminology standardized throughout the manuscript text and figures.
filtered_metadata$Group <- ifelse(filtered_metadata$Age < 18, "Pediatric", "Adult")

filtered_metadata <- filtered_metadata %>%
  mutate(Status = paste(Status, Group, sep = "_"))

filtered_metadata <- filtered_metadata[filtered_metadata$Status %in% c("SC2_Pos_Pediatric", "SC2_Pos_Adult"), ]

raw_data <- raw_data[, colnames(raw_data) %in% filtered_metadata$Sample]

column_order <- filtered_metadata$Sample
raw_data <- raw_data[, column_order]

rownames(filtered_metadata) <- filtered_metadata$Sample

filtered_metadata$Status <- as.factor(filtered_metadata$Status)
filtered_metadata$Status <- relevel(filtered_metadata$Status, ref = "SC2_Pos_Adult")
filtered_metadata$Batch <- as.factor(filtered_metadata$Batch)

# [R2] Standardize Sex as a factor for use as a covariate
filtered_metadata$Sex <- as.factor(filtered_metadata$Sex)

# ---------------------------------------------------------------------------
# [R4] Batch handling -- REMOVED the previous ComBat() pre-correction step.
# ComBat() assumes continuous, roughly Gaussian data (its typical use case is
# microarray intensities); applying it to raw RNA-seq counts, rounding the
# output, and feeding that into DESeq2 breaks the negative-binomial
# assumptions DESeq2 relies on, and constitutes double-dipping (correcting
# and then testing on the same transformed data), which can inflate false
# positives. Batch is instead included directly as a covariate in the DESeq2
# design formula below -- the standard, defensible approach for count data.
# If batch-corrected counts are specifically needed elsewhere (e.g. for PCA
# visualization only, never for DE testing), use ComBat_seq() from the sva
# package, which is built for RNA-seq counts, not ComBat().
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# [R4] Low-count filtering, explicitly reported.
# ---------------------------------------------------------------------------
min_group_size <- min(table(filtered_metadata$Status))
keep <- rowSums(raw_data >= 10) >= min_group_size
raw_data_filtered <- raw_data[keep, ]
cat("Genes before filtering:", nrow(raw_data), "\n")
cat("Genes after low-count filtering:", nrow(raw_data_filtered), "\n")

# ---------------------------------------------------------------------------
# [R2/R4] Design formula adjusts for Batch and Sex as covariates alongside
# Status. Age is not included separately because Status/Group is itself an
# age-based split (collinear) -- report this reasoning in the methods text.
# Check resultsNames(dds) after fitting to confirm Batch has enough levels
# with sufficient samples per level; a batch level with very few samples can
# make the model rank-deficient.
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = raw_data_filtered,
                               colData = filtered_metadata,
                               design = ~ Batch + Sex + Status)

dds <- DESeq(dds)

dir.create("Results/", showWarnings = FALSE)
tiff("Results/DispEsts_GSE231409.tiff", res = 600, width = 2500, height = 2000, compression = "lzw")
plotDispEsts(dds)
dev.off()

# ---------------------------------------------------------------------------
# PCA -- [R22] percentVar reported on the plot axes.
# NOTE: this PCA is on the Batch+Sex+Status-fitted vst, so it still reflects
# raw batch structure visually; that is expected and fine for a diagnostic
# plot -- the correction happens in the DE test itself, not by pre-cleaning
# the counts used for visualization.
# ---------------------------------------------------------------------------
rld <- vst(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = c("Status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

tiff("Results/PCA_R_GSE231409_Hurst_Kelly_Eng.tiff", res = 600, width = 3700, height = 2000, compression = "lzw")
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(aes(colour = Status), size = 3) +
  ggtitle("GSE231409") +
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
# [R6] CONTRAST DIRECTION -- must match the convention used in GSE179277:
# positive log2FC = enriched in Pediatric. The previous version of this
# script had the contrast reversed (Adult vs Pediatric), which likely
# contributed to the sign/count mismatches flagged in reviewer comment 6.
# [R3] FDR control via alpha = 0.05 (BH-adjusted padj, DESeq2 default).
# [R4] lfcShrink applied (apeglm).
# ---------------------------------------------------------------------------
res <- results(dds, contrast = c("Status", "SC2_Pos_Pediatric", "SC2_Pos_Adult"),
                alpha = 0.05)

print(resultsNames(dds))
res_shrunk <- lfcShrink(dds, coef = "Status_SC2_Pos_Pediatric_vs_SC2_Pos_Adult",
                         type = "apeglm")

# ---------------------------------------------------------------------------
# [R4] Normalization -- report size factors used
# ---------------------------------------------------------------------------
size_factors <- sizeFactors(dds)
write.csv(data.frame(Sample = names(size_factors), SizeFactor = size_factors),
          "Results/sizeFactors_GSE231409.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# [R4] Dispersion estimation -- report which fit type was actually used
# ---------------------------------------------------------------------------
cat("Dispersion fit type used:", dispersionFunction(dds)@fitType, "\n")

# ---------------------------------------------------------------------------
# [R4] Independent filtering -- report the FDR-optimized filtering threshold
# ---------------------------------------------------------------------------
filter_threshold <- metadata(res)$filterThreshold
cat("Independent filtering: mean count threshold =",
    round(as.numeric(filter_threshold), 3),
    "(", names(filter_threshold), "percentile )\n")

# ---------------------------------------------------------------------------
# [R4] Outlier handling -- report genes flagged via Cook's distance
# ---------------------------------------------------------------------------
n_outlier_flagged <- sum(is.na(res$padj) & !is.na(res$pvalue))
cat("Genes flagged as outliers by Cook's distance (padj set to NA):",
    n_outlier_flagged, "\n")

table_CD <- data.frame(res_shrunk)
table_CD <- rownames_to_column(table_CD, "ENSEMBL")
table_CD$ENSEMBL <- sub("\\..*", "", table_CD$ENSEMBL)

# ---------------------------------------------------------------------------
# [R4] Gene annotation via mapIds() -- avoids silent duplication from
# one-to-many ENSEMBL-to-SYMBOL mappings that AnnotationDbi::select() can
# produce.
# ---------------------------------------------------------------------------
table_CD$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                          keys = table_CD$ENSEMBL,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")

table_CD <- drop_na(table_CD)

# [R3] DEGs should be defined downstream using padj (BH-adjusted):
# e.g. subset(table_CD, padj < 0.05 & abs(log2FoldChange) > 0.5)

write.csv(table_CD, "Results/DE_R_GSE231409_Hurst_Kelly_revised.csv", row.names = FALSE)

# [R15] Save session info for reproducibility reporting
writeLines(capture.output(sessionInfo()), "Results/sessionInfo_GSE231409.txt")

# [R4] Save all workflow diagnostics together for direct citation in Methods
writeLines(c(
  paste("Genes before low-count filtering:", nrow(raw_data)),
  paste("Genes after low-count filtering:", nrow(raw_data_filtered)),
  paste("Design formula:", deparse(design(dds))),
  paste("Batch handling: included as design covariate (not pre-corrected with ComBat)"),
  paste("Dispersion fit type:", dispersionFunction(dds)@fitType),
  paste("Independent filtering threshold (mean count):",
        round(as.numeric(filter_threshold), 3),
        "at", names(filter_threshold), "percentile"),
  paste("Genes flagged as outliers (Cook's distance, padj = NA):", n_outlier_flagged),
  paste("Log2FC shrinkage method: apeglm"),
  paste("FDR method: Benjamini-Hochberg (padj, DESeq2 default), alpha = 0.05"),
  paste("Contrast: Pediatric vs Adult (positive log2FC = Pediatric-enriched)")
), "Results/DESeq2_workflow_summary_GSE231409.txt")
