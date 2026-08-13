# Author: Jose A Ovando-Ricardez
# Script for downloading, processing, and annotating GSE172274 dataset
# REVISED VERSION -- changes tied to Reviewer 1 comments are marked with [R#]

# Load necessary libraries
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("R_GSE172274/")

# Step 3: Load GEO dataset for metadata processing
gse <- getGEO("GSE172274", GSEMatrix = TRUE)
gse_data <- gse[[1]]
metadata <- pData(gse_data)

# Step 4: Filter and rename relevant metadata columns
# [R2] Sex is kept and used as a covariate. NOTE: this dataset has no
# sequencing-batch variable in GEO metadata (unlike GSE231409), so no
# Batch covariate is included here.
filtered_metadata <- metadata[, c("geo_accession",
                                  "title",
                                  "diagnosis:ch1",
                                  "age:ch1",
                                  "Sex:ch1",
                                  "tissue:ch1")]

colnames(filtered_metadata) <- c("GEO_ID",
                                 "Sample",
                                 "Status",
                                 "Age",
                                 "Sex",
                                 "Tissue")

# Step 5: Standardize the 'Status' column into SARS-CoV-2 positive/negative
filtered_metadata$Status <- ifelse(grepl("Positive", filtered_metadata$Status, ignore.case = TRUE),
                                   "SC2_Pos",
                                   ifelse(grepl("Negative", filtered_metadata$Status, ignore.case = TRUE),
                                          "SC2_Neg", filtered_metadata$Status))

# Step 6: Load raw gene-level count matrix (Entrez-indexed raw counts)
raw_data <- read_tsv("GSE172274_raw_counts_GRCh38.p13_NCBI.tsv")
rownames(raw_data) <- raw_data$GeneID
raw_data$GeneID <- NULL

filtered_metadata$Age <- as.numeric(filtered_metadata$Age)

# [R19/R30] Group is labeled "Pediatric" (not "Peds") to match the
# terminology standardized throughout the manuscript text and figures.
# [R20] Age cutoff: <18 years = Pediatric, consistent with GSE231409.
# NOTE: some Adult rows are pooled samples (POOL 1-4) per the original
# study protocol -- see manuscript limitations text.
filtered_metadata$Group <- ifelse(filtered_metadata$Age < 18, "Pediatric", "Adult")

filtered_metadata <- filtered_metadata[filtered_metadata$GEO_ID %in% colnames(raw_data), ]

column_order <- filtered_metadata$GEO_ID
raw_data <- raw_data[, column_order]

rownames(filtered_metadata) <- filtered_metadata$GEO_ID

filtered_metadata$Status <- as.factor(filtered_metadata$Group)
filtered_metadata$Status <- relevel(filtered_metadata$Status, ref = "Adult")

# [R2] Standardize Sex as a factor for use as a covariate
filtered_metadata$Sex <- as.factor(filtered_metadata$Sex)

# ---------------------------------------------------------------------------
# [R4] Low-count filtering, explicitly reported.
# ---------------------------------------------------------------------------
min_group_size <- min(table(filtered_metadata$Status))
keep <- rowSums(raw_data >= 10) >= min_group_size
raw_data_filtered <- raw_data[keep, ]
cat("Genes before filtering:", nrow(raw_data), "\n")
cat("Genes after low-count filtering:", nrow(raw_data_filtered), "\n")

# ---------------------------------------------------------------------------
# [R2] Design formula adjusts for Sex as a covariate alongside Status.
# Age is not included separately because Status/Group is itself an
# age-based split (collinear) -- report this reasoning in the methods text.
# No Batch covariate: GSE172274 metadata does not include a batch variable.
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = raw_data_filtered,
                              colData = filtered_metadata,
                              design = ~ Sex + Status)

dds <- DESeq(dds)

dir.create("Results/", showWarnings = FALSE)
tiff("Results/DispEsts_GSE172274.tiff", res = 600, width = 2500, height = 2000, compression = "lzw")
plotDispEsts(dds)
dev.off()

# ---------------------------------------------------------------------------
# PCA -- [R22] percentVar reported on the plot axes.
# ---------------------------------------------------------------------------
rld <- vst(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = c("Status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

tiff("Results/PCA_R_GSE172274_Eng.tiff", res = 600, width = 3700, height = 2000, compression = "lzw")
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(aes(colour = Status), size = 3) +
  ggtitle("GSE172274") +
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
# [R6] CONTRAST DIRECTION -- consistent with GSE179277/GSE231409 convention:
# positive log2FC = enriched in Pediatric.
# [R3] FDR control via alpha = 0.05 (BH-adjusted padj, DESeq2 default).
# ---------------------------------------------------------------------------
res <- results(dds, contrast = c("Status", "Pediatric", "Adult"),
               alpha = 0.05)

# ---------------------------------------------------------------------------
# [Meta-analysis prep] Export UNSHRUNKEN log2FC + SE per gene, keyed by
# SYMBOL, for the cross-dataset REM meta-analysis (see REM_meta_analysis.R).
# Uses ENTREZID as keytype because GSE172274 counts are indexed by Entrez
# Gene ID, not ENSEMBL (unlike GSE179277 and GSE231409).
# ---------------------------------------------------------------------------
res_for_meta <- data.frame(res) %>%
  rownames_to_column("ENTREZID") %>%
  mutate(SYMBOL = AnnotationDbi::mapIds(org.Hs.eg.db,
                                        keys = ENTREZID,
                                        column = "SYMBOL",
                                        keytype = "ENTREZID",
                                        multiVals = "first")) %>%
  dplyr::select(SYMBOL, log2FoldChange, lfcSE) %>%
  drop_na(SYMBOL)

write.csv(res_for_meta, "Results/res_for_meta_GSE172274.csv", row.names = FALSE)

# [R4] lfcShrink applied (apeglm) -- for reporting/ranking within this
# dataset only; the meta-analysis above uses the unshrunken estimates.
print(resultsNames(dds))
res_shrunk <- lfcShrink(dds, coef = "Status_Pediatric_vs_Adult",
                        type = "apeglm")

# ---------------------------------------------------------------------------
# [R4] Normalization -- report size factors used
# ---------------------------------------------------------------------------
size_factors <- sizeFactors(dds)
write.csv(data.frame(Sample = names(size_factors), SizeFactor = size_factors),
          "Results/sizeFactors_GSE172274.csv", row.names = FALSE)

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
table_CD <- rownames_to_column(table_CD, "ENTREZID")

# ---------------------------------------------------------------------------
# [R4] Gene annotation via mapIds() -- avoids silent duplication from
# one-to-many ENTREZID-to-SYMBOL mappings that AnnotationDbi::select() can
# produce. NOTE: GSE172274 counts are indexed by Entrez Gene ID, not
# ENSEMBL (unlike GSE179277 and GSE231409) -- keytype changed accordingly.
# ---------------------------------------------------------------------------
table_CD$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                         keys = table_CD$ENTREZID,
                                         column = "SYMBOL",
                                         keytype = "ENTREZID",
                                         multiVals = "first")

table_CD <- drop_na(table_CD)

# [R3] DEGs should be defined downstream using padj (BH-adjusted):
# e.g. subset(table_CD, padj < 0.05 & abs(log2FoldChange) > 0.5)

write.csv(table_CD, "Results/DE_R_GSE172274_revised.csv", row.names = FALSE)

# [R15] Save session info for reproducibility reporting
writeLines(capture.output(sessionInfo()), "Results/sessionInfo_GSE172274.txt")

# [R4] Save all workflow diagnostics together for direct citation in Methods
writeLines(c(
  paste("Genes before low-count filtering:", nrow(raw_data)),
  paste("Genes after low-count filtering:", nrow(raw_data_filtered)),
  paste("Design formula:", deparse(design(dds))),
  paste("Batch handling: no batch variable available in GEO metadata for this dataset"),
  paste("Dispersion fit type:", dispersionFunction(dds)@fitType),
  paste("Independent filtering threshold (mean count):",
        round(as.numeric(filter_threshold), 3),
        "at", names(filter_threshold), "percentile"),
  paste("Genes flagged as outliers (Cook's distance, padj = NA):", n_outlier_flagged),
  paste("Log2FC shrinkage method: apeglm"),
  paste("FDR method: Benjamini-Hochberg (padj, DESeq2 default), alpha = 0.05"),
  paste("Contrast: Pediatric vs Adult (positive log2FC = Pediatric-enriched)")
), "Results/DESeq2_workflow_summary_GSE172274.txt")
