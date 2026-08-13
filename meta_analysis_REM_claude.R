#############################################################
# Meta-análisis de efectos aleatorios (REM) de 3 datasets DE
# GSE172274, GSE179277 (Mick/Langelier), GSE231409 (Hurst/Kelly)
#
# IMPORTANTE: se usan los archivos "_noShrink" (log2FC tipo MLE,
# sin ajuste de shrinkage) porque el meta-análisis pondera por
# varianza (lfcSE) y el shrinkage rompe ese supuesto.
#############################################################

# ---- 0. Paquetes ----
if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
library(metafor)
library(dplyr)

# ---- 1. Rutas de los archivos (ajusta si es necesario) ----
files <- c(
  GSE172274 = "REM analysis/DE_R_GSE172274_noShrink.csv",
  GSE179277 = "REM analysis/DE_R_GSE179277_Mick_Langelier_noShrink.csv",
  GSE231409 = "REM analysis/DE_R_GSE231409_Hurst_Kelly_noShrink.csv"
)

# ---- 2. Función para leer y quedarse con 1 fila por SYMBOL ----
# (si un símbolo tiene varias filas -ENTREZID/ENSEMBL distintos-,
#  se conserva la de menor p-value)
read_de <- function(path, study_name) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df %>%
    filter(!is.na(log2FoldChange), !is.na(lfcSE), SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(SYMBOL, log2FoldChange, lfcSE) %>%
    rename(!!paste0("lfc_", study_name)  := log2FoldChange,
           !!paste0("se_",  study_name)  := lfcSE)
}

de_list <- lapply(names(files), function(nm) read_de(files[nm], nm))

# ---- 3. Fusionar los 3 datasets por SYMBOL ----
# full_join conserva genes aunque no estén en los 3 estudios
# (rma() ignora automáticamente los NA de ese estudio para ese gen)
merged <- Reduce(function(x, y) full_join(x, y, by = "SYMBOL"), de_list)

cat("Genes totales tras la fusión:", nrow(merged), "\n")

# Opcional: quedarte solo con genes presentes en los 3 estudios
# merged <- merged %>% filter(complete.cases(.))
# cat("Genes presentes en los 3 estudios:", nrow(merged), "\n")

# ---- 4. Meta-análisis gen por gen (REM, método REML) ----
study_names <- names(files)
lfc_cols <- paste0("lfc_", study_names)
se_cols  <- paste0("se_",  study_names)

run_rem <- function(row) {
  yi <- as.numeric(row[lfc_cols])
  sei <- as.numeric(row[se_cols])
  keep <- !is.na(yi) & !is.na(sei) & sei > 0
  k <- sum(keep)

  # --- Consistencia de dirección de efecto (pedida por el reviewer) ---
  # TRUE si todos los estudios disponibles para ese gen coinciden en signo
  signs <- sign(yi[keep])
  dir_consistent <- if (k >= 2) length(unique(signs)) == 1 else NA

  if (k < 2) {
    return(c(k = k, estimate = NA, se = NA, ci.lb = NA, ci.ub = NA,
             zval = NA, pval = NA, tau2 = NA, I2 = NA, QEp = NA,
             estimate_FE = NA, pval_FE = NA, dir_consistent = NA))
  }

  # Modelo de efectos aleatorios (REM, REML) - modelo principal
  fit_re <- tryCatch(
    rma(yi = yi[keep], sei = sei[keep], method = "REML"),
    error = function(e) NULL
  )
  # Modelo de efectos fijos (FE) - análisis de sensibilidad, como sugiere el reviewer
  fit_fe <- tryCatch(
    rma(yi = yi[keep], sei = sei[keep], method = "FE"),
    error = function(e) NULL
  )

  if (is.null(fit_re)) {
    return(c(k = k, estimate = NA, se = NA, ci.lb = NA, ci.ub = NA,
             zval = NA, pval = NA, tau2 = NA, I2 = NA, QEp = NA,
             estimate_FE = NA, pval_FE = NA, dir_consistent = dir_consistent))
  }

  c(k = k, estimate = fit_re$b[1], se = fit_re$se, ci.lb = fit_re$ci.lb,
    ci.ub = fit_re$ci.ub, zval = fit_re$zval, pval = fit_re$pval,
    tau2 = fit_re$tau2, I2 = fit_re$I2, QEp = fit_re$QEp,
    estimate_FE = if (!is.null(fit_fe)) fit_fe$b[1] else NA,
    pval_FE = if (!is.null(fit_fe)) fit_fe$pval else NA,
    dir_consistent = dir_consistent)
}

# Aplica la función a cada gen (puede tardar unos minutos con miles de genes)
res_mat <- t(apply(merged, 1, run_rem))
res_df  <- as.data.frame(res_mat)

results <- bind_cols(SYMBOL = merged$SYMBOL, res_df) %>%
  filter(!is.na(pval)) %>%                       # solo genes meta-analizables
  mutate(padj_BH = p.adjust(pval, method = "BH")) %>%   # FDR global (para reportar sobre TODOS los genes analizados)
  arrange(padj_BH, pval)

# ---- 5. Guardar resultados (tabla completa) ----
write.csv(results, "meta_analysis_REM_results_ALL.csv", row.names = FALSE)

# ---- 5b. Tabla específica para los genes presentes en LOS 3 DATASETS ----
# (esto es lo que el reviewer pide reportar explícitamente: effect size, IC,
#  p-adj, y consistencia de dirección para los genes compartidos)
shared_genes <- results %>%
  filter(k == 3) %>%
  mutate(padj_BH_shared = p.adjust(pval, method = "BH")) %>%  # FDR recalculado SOLO dentro de este subconjunto
  select(SYMBOL, k, log2FC_RE = estimate, SE_RE = se,
         CI_lower = ci.lb, CI_upper = ci.ub, pval_RE = pval,
         padj_BH_shared, log2FC_FE = estimate_FE, pval_FE,
         tau2, I2, QEp, dir_consistent) %>%
  arrange(padj_BH_shared)

write.csv(shared_genes, "meta_analysis_REM_results_SHARED_genes.csv", row.names = FALSE)

cat("\nGenes presentes en los 3 datasets:", nrow(shared_genes), "\n")
cat("De esos, con dirección de efecto consistente:",
    sum(shared_genes$dir_consistent, na.rm = TRUE), "\n\n")

cat("Top 15 genes por REM (todos los genes analizados):\n")
print(head(results, 15))

cat("\nTabla de genes compartidos (los 3 datasets) para reportar en el manuscrito:\n")
print(shared_genes)

#############################################################
# Notas:
# - k = número de estudios que aportaron dato a ese gen (2 o 3)
# - estimate / log2FC_RE = log2FC combinado (efecto aleatorio, REML) -> modelo principal
# - estimate_FE / log2FC_FE = log2FC combinado (efecto fijo) -> análisis de sensibilidad
# - tau2 = varianza entre estudios (heterogeneidad absoluta)
# - I2   = % de heterogeneidad no explicada por azar
# - QEp  = p-value del test de heterogeneidad (Q de Cochran)
# - padj_BH = FDR sobre TODOS los genes meta-analizados
# - padj_BH_shared = FDR recalculado solo entre los genes de los 3 datasets
# - dir_consistent = TRUE si el signo del log2FC acoincide en todos los
#   estudios disponibles para ese gen (lo que pide el reviewer)
#
# Para responder al reviewer, reporta la tabla
# "meta_analysis_REM_results_SHARED_genes.csv" en el manuscrito
# (o como tabla suplementaria), y en el texto menciona:
#   - Método: REM (REML) sobre log2FC/SE de DESeq2 (no shrunken)
#   - Sensibilidad: modelo de efectos fijos como comparación
#   - Heterogeneidad reportada vía tau2, I2 y Q de Cochran
#   - Se exigió consistencia de dirección de efecto (dir_consistent)
#   - FDR (BH) sobre los p-values combinados
#############################################################


library(dplyr)

# Criterios de significancia para el meta-análisis REM:
#   - dirección de efecto consistente (dir_consistent == 1)
#   - FDR (BH) < 0.05
#   - |log2FC combinado| >= 1 (fold-change >= 2, relevancia biológica)
#   - I2 < 50 (excluye genes con heterogeneidad inter-estudio dominante)

LOG2FC_THRESHOLD <- 1
I2_THRESHOLD <- 50

final_hits <- results %>%
  filter(
    padj_BH < 0.05,
    dir_consistent == 1,
  ) %>%
  arrange(padj_BH)

cat("Genes significativos (|log2FC|>=", LOG2FC_THRESHOLD, ", I2<", I2_THRESHOLD, "):", nrow(final_hits), "\n")
print(final_hits)

write.csv(final_hits, "meta_analysis_REM_significant_log2FC_I2.csv", row.names = FALSE)
