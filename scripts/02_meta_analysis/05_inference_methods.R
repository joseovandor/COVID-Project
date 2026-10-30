# Sensitivity of the meta-analysis to the inference method
#
# With only two or three studies per gene, every eligible gene is re-estimated with:
#   REML_z         REML + Wald z (primary model, reproduced as a check)
#   REML_HKSJ      REML + Hartung-Knapp-Sidik-Jonkman (t, k - 1 df)
#   REML_HKSJ_mod  REML + HKSJ, ad hoc variant (metafor test = "adhoc")
#   DL_z / PM_z    DerSimonian-Laird / Paule-Mandel + z
#   PM_HKSJ        Paule-Mandel + HKSJ
#   FE_z           fixed effect
# BH-FDR is applied per method across all eligible genes.
#
# For the signature genes, a Bayesian normal-normal model (flat prior on mu;
# half-normal 0.5, half-normal 1.0 and half-Cauchy 0.5 priors on tau) is
# computed by numerical integration over a tau grid and, if installed, checked
# against bayesmeta. The Robust class is also recomputed over a grid of I^2 and
# leave-one-out thresholds, and each gene gets an evidence tier.
#
# Input : results/meta_analysis/harmonization/, primary/, sensitivity/
# Output: results/meta_analysis/inference_methods/

SCRIPT_BUILD <- "META_INFERENCE_METHODS_2026-10-26"
message("Running script build: ", SCRIPT_BUILD)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(metafor)
})

# Force dplyr verbs even if Bioconductor packages (AnnotationDbi, DESeq2,
# S4Vectors, ...) are attached in the session and mask them.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
count <- dplyr::count
first <- dplyr::first
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
distinct <- dplyr::distinct

ALPHA <- 0.05

out_dir <- here::here("results", "meta_analysis", "inference_methods")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# 01. Inputs ----

meta_long <- read.csv(
  here::here("results", "meta_analysis", "harmonization", "meta_input_combined_long.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE)
  )

primary <- read.csv(
  here::here("results", "meta_analysis", "primary", "meta_analysis_REM_all_genes.csv"),
  stringsAsFactors = FALSE
)

robust_tab <- read.csv(
  here::here("results", "meta_analysis", "sensitivity", "robustness_classification.csv"),
  stringsAsFactors = FALSE
)

loo <- read.csv(
  here::here("results", "meta_analysis", "sensitivity", "leave_one_out_results.csv"),
  stringsAsFactors = FALSE
)

signature <- primary %>%
  filter(RE_FDR < ALPHA) %>%
  arrange(RE_FDR)

cat("Primary universe:", nrow(primary), "genes | signature:", nrow(signature), "genes\n")

# Same universe as the primary analysis: genes with a valid primary fit.
meta_long <- meta_long %>%
  filter(SYMBOL %in% primary$SYMBOL)

gene_split <- split(
  meta_long[, c("Dataset", "log2FoldChange", "lfcSE")],
  meta_long$SYMBOL
)
gene_split <- gene_split[primary$SYMBOL]


# 02. Frequentist methods for every eligible gene ----

methods_spec <- list(
  REML_z        = list(method = "REML", test = "z"),
  REML_HKSJ     = list(method = "REML", test = "knha"),
  REML_HKSJ_mod = list(method = "REML", test = "adhoc"),
  DL_z          = list(method = "DL",   test = "z"),
  PM_z          = list(method = "PM",   test = "z"),
  PM_HKSJ       = list(method = "PM",   test = "knha"),
  FE_z          = list(method = "EE",   test = "z")
)

# metafor < 4.0 does not implement test = "adhoc"; fall back to HKSJ with the
# Knapp-Hartung truncation (q* = max(1, q)) computed by hand in that case.
adhoc_supported <- utils::packageVersion("metafor") >= "4.0.0"

fit_one <- function(yi, sei, method, test) {
  test_used <- if (test == "adhoc" && !adhoc_supported) "knha" else test
  # Same call as the primary analysis first; step-halving only if it fails.
  fit <- tryCatch(
    suppressWarnings(metafor::rma(yi = yi, sei = sei, method = method, test = test_used)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      suppressWarnings(metafor::rma(yi = yi, sei = sei, method = method, test = test_used,
                                    control = list(stepadj = 0.5, maxiter = 1000))),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    return(c(est = NA, se = NA, ci_lb = NA, ci_ub = NA, pval = NA, tau2 = NA))
  }
  if (test == "adhoc" && !adhoc_supported) {
    # Ad hoc correction: use the Wald SE whenever HKSJ would be smaller.
    fit_z <- metafor::rma(yi = yi, sei = sei, method = method, test = "z",
                          control = list(stepadj = 0.5, maxiter = 1000))
    se_adj <- max(fit$se, fit_z$se)
    df <- length(yi) - 1
    tval <- as.numeric(fit$b) / se_adj
    return(c(
      est = as.numeric(fit$b), se = se_adj,
      ci_lb = as.numeric(fit$b) - qt(0.975, df) * se_adj,
      ci_ub = as.numeric(fit$b) + qt(0.975, df) * se_adj,
      pval = 2 * pt(abs(tval), df, lower.tail = FALSE),
      tau2 = fit$tau2
    ))
  }
  c(
    est = as.numeric(fit$b), se = fit$se,
    ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
    pval = fit$pval, tau2 = fit$tau2
  )
}

message("Fitting ", length(methods_spec), " methods x ", length(gene_split), " genes...")
t0 <- Sys.time()

method_results <- list()

for (m in names(methods_spec)) {
  spec <- methods_spec[[m]]
  res_m <- t(vapply(
    gene_split,
    function(d) fit_one(d$log2FoldChange, d$lfcSE, spec$method, spec$test),
    numeric(6)
  ))
  res_m <- as.data.frame(res_m)
  res_m$SYMBOL <- names(gene_split)
  res_m$k <- vapply(gene_split, nrow, integer(1))
  res_m$method <- m
  res_m$FDR <- p.adjust(res_m$pval, method = "BH")
  method_results[[m]] <- res_m
  message("  ", m, " done (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
}

all_methods_long <- bind_rows(method_results)

# Sanity check: REML_z must reproduce the primary analysis.
check <- method_results$REML_z %>%
  inner_join(primary %>% select(SYMBOL, RE_pvalue, RE_log2FC), by = "SYMBOL")
max_p_diff <- max(abs(check$pval - check$RE_pvalue), na.rm = TRUE)
max_b_diff <- max(abs(check$est - check$RE_log2FC), na.rm = TRUE)
cat("Reproduction of primary REML_z: max |p diff| =", signif(max_p_diff, 3),
    "| max |log2FC diff| =", signif(max_b_diff, 3), "\n")
if (max_b_diff > 1e-3) warning("REML_z does not reproduce the primary estimates closely; inspect.")

write.csv(
  all_methods_long,
  file.path(out_dir, "meta_inference_all_genes_long.csv"),
  row.names = FALSE
)

# Summary: genome-wide number of significant genes per method.
genome_summary <- all_methods_long %>%
  group_by(method) %>%
  summarise(
    n_genes = n(),
    n_FDR05 = sum(FDR < ALPHA, na.rm = TRUE),
    n_FDR05_k3 = sum(FDR < ALPHA & k == 3, na.rm = TRUE),
    n_FDR05_k2 = sum(FDR < ALPHA & k == 2, na.rm = TRUE),
    n_signature_retained = sum(FDR < ALPHA & SYMBOL %in% signature$SYMBOL, na.rm = TRUE),
    min_p_genome = min(pval, na.rm = TRUE),
    # p value a single top gene would need to reach BH-FDR < 0.05
    BH_threshold_rank1 = ALPHA / sum(!is.na(pval)),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = names(methods_spec))) %>%
  arrange(method)

write.csv(genome_summary, file.path(out_dir, "meta_inference_genome_wide_summary.csv"), row.names = FALSE)
print(genome_summary)


# 03. Bayesian random-effects model for the signature genes
#
# Normal-normal hierarchical model:
#   y_i ~ N(theta_i, s_i^2), theta_i ~ N(mu, tau^2),
#   p(mu) flat, p(tau) = half-normal(scale) or half-Cauchy(scale).
# The marginal posterior of tau is evaluated on a fine grid; the posterior of
# mu is the corresponding mixture of normals (exact up to grid resolution).

bayes_grid <- function(y, s, prior = c("halfnormal", "halfcauchy"), scale = 0.5,
                       n_grid = 4000) {
  prior <- match.arg(prior)
  upper <- if (prior == "halfnormal") 6 * scale else 60 * scale
  tau <- seq(0, upper, length.out = n_grid)
  log_prior <- if (prior == "halfnormal") {
    dnorm(tau, 0, scale, log = TRUE)
  } else {
    dcauchy(tau, 0, scale, log = TRUE)
  }
  v <- outer(s^2, tau^2, "+")                       # k x grid
  w <- 1 / v
  V <- 1 / colSums(w)
  mu_hat <- colSums(w * y) * V
  resid2 <- colSums(w * (y - matrix(mu_hat, length(y), n_grid, byrow = TRUE))^2)
  log_lik <- 0.5 * log(V) - 0.5 * colSums(log(v)) - 0.5 * resid2
  log_post <- log_prior + log_lik
  wt <- exp(log_post - max(log_post))
  wt <- wt / sum(wt)
  cdf <- function(x) sum(wt * pnorm(x, mu_hat, sqrt(V)))
  lo <- min(mu_hat - 10 * sqrt(V)); hi <- max(mu_hat + 10 * sqrt(V))
  q <- function(p) uniroot(function(x) cdf(x) - p, c(lo, hi), tol = 1e-8)$root
  p_pos <- 1 - cdf(0)
  c(
    post_mean = sum(wt * mu_hat),
    post_median = q(0.5),
    cri_lb = q(0.025),
    cri_ub = q(0.975),
    p_mu_gt_0 = p_pos,
    p_direction = max(p_pos, 1 - p_pos),
    tau_post_median = tau[which(cumsum(wt) >= 0.5)[1]]
  )
}

bayes_priors <- list(
  HN0.5 = list(prior = "halfnormal", scale = 0.5),
  HN1.0 = list(prior = "halfnormal", scale = 1.0),
  HC0.5 = list(prior = "halfcauchy", scale = 0.5)
)

bayes_results <- list()

for (pr in names(bayes_priors)) {
  spec <- bayes_priors[[pr]]
  b <- t(vapply(
    signature$SYMBOL,
    function(g) {
      d <- gene_split[[g]]
      bayes_grid(d$log2FoldChange, d$lfcSE, spec$prior, spec$scale)
    },
    numeric(7)
  ))
  b <- as.data.frame(b)
  b$SYMBOL <- signature$SYMBOL
  b$prior <- pr
  b$cri_excludes_0 <- b$cri_lb > 0 | b$cri_ub < 0
  bayes_results[[pr]] <- b
}

bayes_long <- bind_rows(bayes_results)

# Optional cross-check with the bayesmeta package (if installed).
if (requireNamespace("bayesmeta", quietly = TRUE)) {
  bm_check <- bind_rows(lapply(signature$SYMBOL, function(g) {
    d <- gene_split[[g]]
    bm <- bayesmeta::bayesmeta(
      y = d$log2FoldChange, sigma = d$lfcSE,
      tau.prior = function(t) bayesmeta::dhalfnormal(t, scale = 0.5)
    )
    ci <- bm$post.interval(mu.level = 0.95, method = "central")
    data.frame(SYMBOL = g, bayesmeta_cri_lb = ci[1], bayesmeta_cri_ub = ci[2],
               bayesmeta_median = bm$summary["median", "mu"])
  }))
  bayes_long <- bayes_long %>%
    left_join(bm_check %>% mutate(prior = "HN0.5"), by = c("SYMBOL", "prior"))
  cat("bayesmeta cross-check (HN0.5): max |CrI difference| =",
      signif(max(abs(c(
        bm_check$bayesmeta_cri_lb - bayes_results$HN0.5$cri_lb,
        bm_check$bayesmeta_cri_ub - bayes_results$HN0.5$cri_ub
      ))), 3), "\n")
}

write.csv(bayes_long, file.path(out_dir, "bayesian_signature_genes.csv"), row.names = FALSE)


# 04. Signature-level table: every signature gene under every method ----

sig_methods_wide <- all_methods_long %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  select(SYMBOL, method, est, ci_lb, ci_ub, pval, FDR) %>%
  pivot_wider(
    names_from = method,
    values_from = c(est, ci_lb, ci_ub, pval, FDR),
    names_glue = "{method}_{.value}"
  )

bayes_wide <- bayes_long %>%
  select(SYMBOL, prior, post_median, cri_lb, cri_ub, p_direction, cri_excludes_0) %>%
  pivot_wider(
    names_from = prior,
    values_from = c(post_median, cri_lb, cri_ub, p_direction, cri_excludes_0),
    names_glue = "Bayes_{prior}_{.value}"
  )

signature_table <- signature %>%
  select(SYMBOL, k, RE_log2FC, RE_FDR, I2, tau2) %>%
  left_join(robust_tab %>% select(SYMBOL, RobustnessClass), by = "SYMBOL") %>%
  left_join(sig_methods_wide, by = "SYMBOL") %>%
  left_join(bayes_wide, by = "SYMBOL")

sig_flag_cols <- paste0(names(methods_spec), "_FDR")
signature_table$n_frequentist_methods_FDR05 <- rowSums(
  signature_table[, sig_flag_cols] < ALPHA, na.rm = TRUE
)
signature_table$retained_REML_HKSJ <- signature_table$REML_HKSJ_FDR < ALPHA
signature_table$retained_REML_HKSJ_mod <- signature_table$REML_HKSJ_mod_FDR < ALPHA
signature_table$retained_Bayes_HN0.5 <- signature_table$Bayes_HN0.5_cri_excludes_0

signature_table$REML_HKSJ_nominal_p05 <- signature_table$REML_HKSJ_pval < 0.05

# With k = 3 the HKSJ reference distribution is t(2): the smallest p value
# observed genome-wide is far above the BH threshold for ~18,000 tests, so
# genome-wide FDR control under HKSJ is unattainable by construction (see
# meta_inference_genome_wide_summary.csv). HKSJ is therefore used at the
# gene level (nominal p) for the signature genes.
# Proposed evidence tier (descriptive):
#   Tier 1: k = 3, HKSJ nominal p < 0.05 AND Bayesian HN(0.5) 95% CrI excludes 0
#   Tier 2: Bayesian HN(0.5) 95% CrI excludes 0 (not Tier 1)
#   Tier 3: supported only by Wald-type inference (exploratory)
signature_table <- signature_table %>%
  mutate(
    EvidenceTier = case_when(
      k == 3 & REML_HKSJ_nominal_p05 & retained_Bayes_HN0.5 ~ "Tier1_k3_HKSJ_and_Bayesian",
      retained_Bayes_HN0.5 ~ "Tier2_Bayesian",
      TRUE ~ "Tier3_exploratory"
    )
  ) %>%
  arrange(EvidenceTier, RE_FDR)

write.csv(signature_table, file.path(out_dir, "signature_genes_by_inference_method.csv"), row.names = FALSE)

nom_cols <- paste0(names(methods_spec), "_pval")

retention_summary <- data.frame(
  criterion = c(
    paste0(names(methods_spec), " FDR < 0.05"),
    paste0(names(methods_spec), " nominal p < 0.05"),
    paste0("Bayesian ", names(bayes_priors), " 95% CrI excludes 0")
  ),
  n_all = c(
    vapply(sig_flag_cols, function(cl) sum(signature_table[[cl]] < ALPHA, na.rm = TRUE), numeric(1)),
    vapply(nom_cols, function(cl) sum(signature_table[[cl]] < 0.05, na.rm = TRUE), numeric(1)),
    vapply(names(bayes_priors), function(pr) sum(signature_table[[paste0("Bayes_", pr, "_cri_excludes_0")]], na.rm = TRUE), numeric(1))
  ),
  n_k3 = c(
    vapply(sig_flag_cols, function(cl) sum(signature_table[[cl]] < ALPHA & signature_table$k == 3, na.rm = TRUE), numeric(1)),
    vapply(nom_cols, function(cl) sum(signature_table[[cl]] < 0.05 & signature_table$k == 3, na.rm = TRUE), numeric(1)),
    vapply(names(bayes_priors), function(pr) sum(signature_table[[paste0("Bayes_", pr, "_cri_excludes_0")]] & signature_table$k == 3, na.rm = TRUE), numeric(1))
  ),
  n_k2 = c(
    vapply(sig_flag_cols, function(cl) sum(signature_table[[cl]] < ALPHA & signature_table$k == 2, na.rm = TRUE), numeric(1)),
    vapply(nom_cols, function(cl) sum(signature_table[[cl]] < 0.05 & signature_table$k == 2, na.rm = TRUE), numeric(1)),
    vapply(names(bayes_priors), function(pr) sum(signature_table[[paste0("Bayes_", pr, "_cri_excludes_0")]] & signature_table$k == 2, na.rm = TRUE), numeric(1))
  ),
  total_signature = nrow(signature_table),
  total_k3 = sum(signature_table$k == 3),
  total_k2 = sum(signature_table$k == 2)
)

write.csv(retention_summary, file.path(out_dir, "signature_retention_by_method.csv"), row.names = FALSE)
print(retention_summary)
print(table(signature_table$EvidenceTier))


# 05. Alternative 'Robust' thresholds ----

loo_k3 <- loo %>%
  filter(SYMBOL %in% signature$SYMBOL[signature$k == 3]) %>%
  group_by(SYMBOL) %>%
  summarise(
    loo_all_same_direction = all(Same_direction_as_full_RE %in% TRUE),
    loo_max_p = max(pvalue, na.rm = TRUE),
    .groups = "drop"
  )

robust_base <- robust_tab %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  select(SYMBOL, k, I2, DirectionConsistent, RE_FE_same_direction, RobustnessClass) %>%
  left_join(loo_k3, by = "SYMBOL") %>%
  left_join(signature_table %>% select(SYMBOL, REML_HKSJ_FDR, REML_HKSJ_pval), by = "SYMBOL")

grid <- expand.grid(
  I2_max = c(75, 50, 25),
  loo_p_max = c(0.10, 0.05, 0.01),
  inference = c("REML_z", "REML_HKSJ_nominal"),
  stringsAsFactors = FALSE
)

robust_membership <- list()

grid$n_robust <- NA_integer_
grid$robust_genes <- NA_character_

for (i in seq_len(nrow(grid))) {
  is_robust <- with(
    robust_base,
    k == 3 &
      DirectionConsistent %in% TRUE &
      I2 < grid$I2_max[i] &
      RE_FE_same_direction %in% TRUE &
      loo_all_same_direction %in% TRUE &
      loo_max_p < grid$loo_p_max[i] &
      (grid$inference[i] == "REML_z" | REML_HKSJ_pval < 0.05)
  )
  is_robust[is.na(is_robust)] <- FALSE
  grid$n_robust[i] <- sum(is_robust)
  grid$robust_genes[i] <- paste(sort(robust_base$SYMBOL[is_robust]), collapse = ";")
  robust_membership[[i]] <- data.frame(
    SYMBOL = robust_base$SYMBOL,
    setting = paste0("I2<", grid$I2_max[i], " | LOO p<", grid$loo_p_max[i], " | ", grid$inference[i]),
    robust = is_robust
  )
}

primary_setting <- grid$I2_max == 75 & grid$loo_p_max == 0.05 & grid$inference == "REML_z"
cat("Primary Robust definition reproduced:", grid$n_robust[primary_setting],
    "(robustness_classification.csv:", sum(robust_base$RobustnessClass == "Robust"), ")\n")

primary_robust <- robust_base$SYMBOL[robust_base$RobustnessClass == "Robust"]
grid$n_primary_robust_retained <- vapply(
  strsplit(grid$robust_genes, ";"),
  function(x) sum(primary_robust %in% x),
  numeric(1)
)
grid$n_primary_robust <- length(primary_robust)

write.csv(grid, file.path(out_dir, "robust_threshold_sensitivity.csv"), row.names = FALSE)
write.csv(bind_rows(robust_membership), file.path(out_dir, "robust_threshold_membership_long.csv"), row.names = FALSE)
print(grid[, c("I2_max", "loo_p_max", "inference", "n_robust", "n_primary_robust_retained")])


# 06. Figures ----

method_labels <- c(
  REML_z = "REML + z\n(primary)", REML_HKSJ = "REML + HKSJ\n(nominal)", REML_HKSJ_mod = "REML + HKSJ\nad hoc (nominal)",
  DL_z = "DL + z", PM_z = "PM + z", PM_HKSJ = "PM + HKSJ\n(nominal)", FE_z = "Fixed effect",
  Bayes_HN0.5 = "Bayes\nHN(0.5)", Bayes_HN1.0 = "Bayes\nHN(1.0)", Bayes_HC0.5 = "Bayes\nHC(0.5)"
)

gene_order <- signature_table %>%
  arrange(k, RE_log2FC) %>%
  pull(SYMBOL)

tile_freq <- all_methods_long %>%
  filter(SYMBOL %in% signature$SYMBOL) %>%
  transmute(SYMBOL, method,
            supported = ifelse(grepl("HKSJ", method), pval < 0.05, FDR < ALPHA),
            value = -log10(FDR))

tile_bayes <- bayes_long %>%
  transmute(
    SYMBOL, method = paste0("Bayes_", prior), supported = cri_excludes_0,
    value = -log10(pmax(1 - p_direction, 1e-6))
  )

tile_df <- bind_rows(tile_freq, tile_bayes) %>%
  left_join(signature_table %>% select(SYMBOL, k), by = "SYMBOL") %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = gene_order),
    method = factor(method, levels = names(method_labels)),
    k_label = paste0("k = ", k),
    status = ifelse(supported, "Supported", "Not supported")
  )

p_tiles <- ggplot(tile_df, aes(x = method, y = SYMBOL, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(k_label ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = c("Supported" = "#2C6FB7", "Not supported" = "#E6E6E6")) +
  labs(x = NULL, y = NULL, fill = NULL,
       caption = "Supported: genome-wide FDR < 0.05 (z-based and fixed-effect methods); nominal p < 0.05 (HKSJ, t with k-1 df);\n95% credible interval excluding 0 (Bayesian). HN = half-normal, HC = half-Cauchy prior on tau.") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, size = 8),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

ggsave(file.path(fig_dir, "FigS_meta_inference_signature_tiles.pdf"), p_tiles, width = 10, height = 9.5)
ggsave(file.path(fig_dir, "FigS_meta_inference_signature_tiles.tiff"), p_tiles, width = 10, height = 9.5, dpi = 600, compression = "lzw")

ci_df <- bind_rows(
  all_methods_long %>%
    filter(SYMBOL %in% signature$SYMBOL, method %in% c("REML_z", "REML_HKSJ", "REML_HKSJ_mod")) %>%
    transmute(SYMBOL, method, est, lb = ci_lb, ub = ci_ub),
  bayes_long %>%
    filter(prior == "HN0.5") %>%
    transmute(SYMBOL, method = "Bayes_HN0.5", est = post_median, lb = cri_lb, ub = cri_ub)
) %>%
  left_join(signature_table %>% select(SYMBOL, k), by = "SYMBOL") %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = gene_order),
    method = factor(method, levels = c("REML_z", "REML_HKSJ", "REML_HKSJ_mod", "Bayes_HN0.5")),
    k_label = paste0("k = ", k),
    # truncate very wide HKSJ (t with 1 df) intervals for display
    lb_plot = pmax(lb, -6), ub_plot = pmin(ub, 6)
  )

p_ci <- ggplot(ci_df, aes(y = SYMBOL, x = est, xmin = lb_plot, xmax = ub_plot, colour = method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.15, linewidth = 0.4) +
  facet_grid(k_label ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(
    values = c(REML_z = "#222222", REML_HKSJ = "#D55E00", REML_HKSJ_mod = "#E69F00", Bayes_HN0.5 = "#0072B2"),
    labels = c("REML + z (primary)", "REML + HKSJ", "REML + HKSJ (ad hoc)", "Bayesian HN(0.5)")
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(x = "Pooled log2 fold change (Pediatric vs Adult), 95% interval\n(intervals truncated at +/-6 for display)", y = NULL, colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_meta_inference_intervals.pdf"), p_ci, width = 8, height = 10)
ggsave(file.path(fig_dir, "FigS_meta_inference_intervals.tiff"), p_ci, width = 8, height = 10, dpi = 600, compression = "lzw")

p_grid <- ggplot(grid, aes(x = factor(loo_p_max), y = factor(I2_max, levels = c(25, 50, 75)), fill = n_robust)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_robust), size = 4) +
  facet_wrap(~ inference, labeller = labeller(inference = c(REML_z = "Inference: REML + z (FDR < 0.05)", REML_HKSJ_nominal = "Additionally: REML + HKSJ p < 0.05"))) +
  scale_fill_gradient(low = "#F1F5FB", high = "#2C6FB7") +
  labs(x = "Leave-one-cohort-out: maximum p value", y = "Maximum I2 (%)", fill = "Robust genes") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "FigS_robust_threshold_grid.pdf"), p_grid, width = 8, height = 3.6)
ggsave(file.path(fig_dir, "FigS_robust_threshold_grid.tiff"), p_grid, width = 8, height = 3.6, dpi = 600, compression = "lzw")


# 07. Run summary and session information ----

writeLines(c(
  paste("Script build:", SCRIPT_BUILD),
  paste("Eligible genes:", nrow(primary), "| signature genes:", nrow(signature),
        "(k=3:", sum(signature$k == 3), "; k=2:", sum(signature$k == 2), ")"),
  paste("REML_z reproduction: max |log2FC diff| =", signif(max_b_diff, 3)),
  "",
  "Signature genes retained (FDR < 0.05 or 95% CrI excluding 0):",
  paste0("  ", retention_summary$criterion, ": ", retention_summary$n_all, "/", retention_summary$total_signature,
         " (k=3: ", retention_summary$n_k3, "/", retention_summary$total_k3,
         "; k=2: ", retention_summary$n_k2, "/", retention_summary$total_k2, ")"),
  "",
  "Evidence tiers:",
  paste0("  ", names(table(signature_table$EvidenceTier)), ": ", table(signature_table$EvidenceTier)),
  "",
  "Robust under alternative thresholds: see robust_threshold_sensitivity.csv",
  paste("metafor version:", utils::packageVersion("metafor"))
), file.path(out_dir, "run_summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done: ", out_dir)
