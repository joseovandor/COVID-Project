# Age-associated upper-airway transcriptomic signatures in SARS-CoV-2-positive individuals

Code for the analyses in *Cross-Dataset Identification of Age-Associated Upper-Airway
Transcriptomic Signatures in SARS-CoV-2-Positive Individuals*.

Three public bulk RNA-seq cohorts of SARS-CoV-2-positive upper-airway samples
(GSE172274, GSE179277 and GSE231409; 146 pediatric and 79 adult samples) are analysed
separately with DESeq2 and combined with a random-effects meta-analysis. The results are
then examined with sensitivity analyses, pre-ranked GSEA and a donor-level analysis of
the nasal single-cell dataset of Yoshida et al. (2022).

## Data

| Dataset | Type | Source |
|---|---|---|
| GSE172274 | bulk RNA-seq, nasopharyngeal swabs | GEO (Pierce et al., 2021) |
| GSE179277 | bulk RNA-seq, nasopharyngeal swabs | GEO (Mick et al., 2022) |
| GSE231409 | bulk RNA-seq, nasal/nasopharyngeal swabs | GEO (Hurst et al., 2025) |
| Yoshida et al. (2022) | scRNA-seq, nasal cavity | CELLxGENE (H5AD) |

The count matrices and the H5AD file are downloaded by the scripts into `data/` the
first time they are run (see `data/README.md`). Data and results are not tracked by git.

## Requirements

R 4.5.1 and Bioconductor 3.22. Main packages: DESeq2 1.50.2, metafor 5.2-1,
fgsea 1.36.2, msigdbr 26.1.1, edgeR 4.8.2, limma 3.66.0, SingleCellExperiment 1.32.0,
zellkonverter 1.20.1, rhdf5 2.54.1, scuttle 1.20.0, org.Hs.eg.db 3.22.0, apeglm 1.32.0,
GEOquery 2.78.0, tidyverse, here, patchwork and openxlsx. bayesmeta 3.5 is optional (used
only as a cross-check of the Bayesian model). The `sessionInfo()` of every step is saved
next to its results.

Open `sars-cov2-upper-airway-age-signature.Rproj` before running anything; all paths are
built with `here::here()` from the project root.

## Running the analysis

Run the scripts in this order. Each folder only needs the outputs of the previous ones.

| Step | Script | Main output |
|---|---|---|
| 1 | `scripts/01_differential_expression/01-03_deseq2_*.R` | `results/deseq2/<cohort>/` |
| 2 | `scripts/02_meta_analysis/01_harmonize_cohorts.R` | `results/meta_analysis/harmonization/` |
|   | `02_random_effects_model.R` | `results/meta_analysis/primary/` |
|   | `03_leave_one_out_robustness.R` | `results/meta_analysis/sensitivity/` |
|   | `04_meta_analysis_plots.R` | `results/meta_analysis/figures/` |
|   | `05_inference_methods.R` | `results/meta_analysis/inference_methods/` |
| 3 | `scripts/03_sensitivity_analyses/01_clinical_covariates.R` | `results/sensitivity/clinical_covariates/` |
|   | `02_continuous_age_sex_interaction.R` | `results/sensitivity/continuous_age_sex_interaction/` |
|   | `03_age_by_infection_GSE179277.R` | `results/sensitivity/age_by_infection/` |
| 4 | `scripts/04_gsea/01_gsea_by_cohort.R` | `results/gsea/by_cohort/` |
|   | `02_cross_cohort_concordance.R` | `results/gsea/cross_cohort/` |
|   | `03_gsea_plots.R` | `results/gsea/figures_final/` |
|   | `04_themes_partial_conjunction.R` | `results/gsea/themes/` |
| 5 | `scripts/05_single_cell/01_import_yoshida.R` | `results/scRNAseq_COVID_nasal_age/` |
|   | `02_umap_composition.R` | `.../cellular_composition/` |
|   | `03_pseudobulk_signature.R` | `.../pseudobulk_signature/` |
| 6 | `scripts/06_figures_tables/01_main_figures.R` | `results/manuscript/main_figures/` |
|   | `02_supplementary_tables.R` | `results/manuscript/supplementary_tables/` |
|   | `03_supplementary_figures.R` | `results/manuscript/supplementary_figures/` |

The three single-cell scripts share objects and must be run in the same R session
(`01` first). Step 5 reads the signature from step 2, and step 6 reads the outputs of
all previous steps.

The single-cell import is the slowest step and needs enough memory to read the H5AD
file.

## Main analysis choices

- Pediatric: age < 18 years; adult: age >= 18 years, assigned from the recorded age.
- DESeq2: `~ Sex + Group` (GSE172274, GSE179277) and `~ Batch + Sex + Group` (GSE231409);
  a gene is kept when its raw count is >= 10 in at least as many samples as the smaller
  age group.
- Meta-analysis: unshrunken log2 fold changes and standard errors, REML random-effects
  model with a z test, Benjamini-Hochberg FDR across all genes present in at least two
  cohorts, without any preselection.
- GSEA: genes ranked by the DESeq2 Wald statistic; MSigDB Hallmark, Reactome, KEGG and
  GO Biological Process.
- Single cell: pseudobulk of raw UMI counts per donor and cell type (edgeR
  quasi-likelihood) as the main analysis.

## Repository layout

```
config/       labels used for the GSEA themes in the figures
data/         downloaded input data (not tracked)
results/      all outputs (not tracked)
scripts/      analysis scripts, numbered in the order they are run
```

## Citation

If you use this code, please cite the article (reference to be added on publication).
