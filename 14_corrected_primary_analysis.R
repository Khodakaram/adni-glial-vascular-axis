# ==========================================================
# 14_corrected_primary_analysis.R
#
# Produces the corrected versions of Tables 1, 2 and 3 plus the
# new supplementary analyses requested by the reviewers.
#
# Two corrections are applied relative to the published models:
#   1. VEGF / sICAM-1 / sVCAM-1 are stored as log10(concentration)
#      in the ADNI QC multiplex release. The published models applied
#      a second natural log. Outcome is now ln(concentration).
#   2. `age` was built as baseline_age + years_from_baseline, making
#      it perfectly collinear with time within participants. The
#      longitudinal models now use baseline age, held fixed.
#
# Sections
#   A. Corrected outcome construction
#   B. Table 1  - descriptives in native units
#   C. Table 2  - baseline models, correct FDR family      (R1-4, R2-4)
#   D. Table 3  - longitudinal models, baseline age        (R2-4)
#   E. Random intercept vs random slope, LRT               (R2-10)
#   F. Residual diagnostics and collinearity               (R2-10)
#   G. Discrimination: AUC, Cohen's d, overlap             (R2-7, R2-9)
#   H. Correlations with bootstrap CIs                     (R2-8)
#   I. APOE e4 x diagnosis interaction                     (R2-11)
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
            "broom", "broom.mixed", "car", "pROC")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse); library(janitor); library(openxlsx)
library(lme4); library(lmerTest); library(broom); library(broom.mixed)
library(car); library(pROC)

# car masks dplyr::recode and pROC masks stats::var and stats::cov.
# Both are therefore called with explicit namespaces below. Do not
# drop the prefixes.

set.seed(20260829)

project_dir <- "R:/ADNI_Project"
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

out <- list()
parse_num <- function(x) readr::parse_number(as.character(x))

# ---------- Marker definitions ----------

marker_info <- tribble(
  ~marker,                ~label,           ~domain,    ~stored_scale, ~unit,
  "gfap_quanterix",       "Plasma GFAP",    "glial",    "raw",         "pg/mL",
  "strem2_msd_corrected", "Plasma sTREM2",  "glial",    "raw",         "pg/mL",
  "vegf_plasma_qc",       "Plasma VEGF",    "vascular", "log10",       "pg/mL",
  "sicam1_plasma_qc",     "Plasma sICAM-1", "vascular", "log10",       "ng/mL",
  "svcam1_plasma_qc",     "Plasma sVCAM-1", "vascular", "log10",       "ng/mL"
)
primary_markers <- marker_info$marker

# ---------- Load ----------

dat <- read_csv(file.path(clean_dir, "analysis_master_model_ready.csv"),
                show_col_types = FALSE, guess_max = 100000) %>%
  clean_names() %>%
  mutate(
    rid = as.character(rid), visit_key = as.character(visit_key),
    dx_label    = factor(dx_label,    levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender),
    age = parse_num(age), pteducat = parse_num(pteducat),
    apoe4 = parse_num(apoe4),
    years_from_baseline = parse_num(years_from_baseline)
  ) %>%
  mutate(baseline_age = age - years_from_baseline)   # time-invariant


# ==========================================================
# A. CORRECTED OUTCOME CONSTRUCTION
#
#   raw   markers: ln(value)
#   log10 markers: ln(concentration) = ln(10) * value
#
# Both give a coefficient interpretable as a proportional
# difference in concentration, so exp(beta) - 1 is a true
# percent change and GFAP is directly comparable to sVCAM-1.
# ==========================================================

for (i in seq_len(nrow(marker_info))) {
  m  <- marker_info$marker[i]
  sc <- marker_info$stored_scale[i]
  v  <- parse_num(dat[[m]])

  dat[[paste0("conc_", m)]] <- if (sc == "log10") 10^v else v          # native units
  dat[[paste0("ln_",   m)]] <- if (sc == "log10") log(10) * v else
                                 ifelse(!is.na(v) & v > 0, log(v), NA_real_)
}

baseline_dat <- dat %>% filter(visit_key == "bl", dx_label %in% c("CN", "MCI", "AD"))

covars_ok <- function(d) {
  d %>% filter(!is.na(age), !is.na(ptgender), !is.na(pteducat), !is.na(apoe4))
}

baseline_cc <- function(m) {
  baseline_dat %>%
    filter(!is.na(.data[[paste0("ln_", m)]]), !is.na(dx_label)) %>%
    covars_ok()
}

long_cc <- function(m) {
  dat %>%
    filter(!is.na(.data[[paste0("ln_", m)]]), !is.na(years_from_baseline),
           !is.na(baseline_dx), !is.na(baseline_age)) %>%
    covars_ok()
}

out$A_scale_key <- marker_info


# ==========================================================
# B. TABLE 1 - descriptives in native units
# ==========================================================

table1 <- map_dfr(primary_markers, function(m) {
  baseline_cc(m) %>%
    group_by(dx_label) %>%
    summarise(
      n      = n(),
      median = median(.data[[paste0("conc_", m)]]),
      q1     = quantile(.data[[paste0("conc_", m)]], 0.25),
      q3     = quantile(.data[[paste0("conc_", m)]], 0.75),
      .groups = "drop"
    ) %>%
    mutate(marker = m)
}) %>%
  left_join(marker_info, by = "marker") %>%
  mutate(median_iqr = sprintf("%.1f [%.1f-%.1f]", median, q1, q3)) %>%
  select(label, unit, dx_label, n, median_iqr, median, q1, q3)

out$B_table1 <- table1


# ==========================================================
# C. TABLE 2 - baseline adjusted models
#    FDR family = the 10 reported diagnosis contrasts.
# ==========================================================

baseline_fits <- map(set_names(primary_markers), function(m) {
  d <- baseline_cc(m)
  lm(as.formula(paste0("ln_", m, " ~ dx_label + age + ptgender + pteducat + apoe4")), data = d)
})

table2 <- imap_dfr(baseline_fits, function(fit, m) {
  tidy(fit) %>% mutate(marker = m, n = nobs(fit))
}) %>%
  filter(term %in% c("dx_labelMCI", "dx_labelAD")) %>%
  mutate(
    contrast = dplyr::recode(term,
                             dx_labelMCI = "MCI vs CN",
                             dx_labelAD  = "AD vs CN"),
    percent_change = (exp(estimate) - 1) * 100,
    pct_ci_low  = (exp(estimate - 1.96 * std.error) - 1) * 100,
    pct_ci_high = (exp(estimate + 1.96 * std.error) - 1) * 100,
    fdr_q = p.adjust(p.value, "BH")
  ) %>%
  left_join(marker_info, by = "marker") %>%
  select(label, contrast, estimate, std.error, percent_change,
         pct_ci_low, pct_ci_high, p.value, fdr_q, n)

out$C_table2 <- table2


# ==========================================================
# D. TABLE 3 - longitudinal models with baseline age
#    FDR family = the 25 reported longitudinal terms.
# ==========================================================

long_form <- function(m, re = "(1 | rid)") {
  as.formula(paste0("ln_", m,
    " ~ years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat + apoe4 + ", re))
}

long_fits_ri <- map(set_names(primary_markers), function(m) {
  lmer(long_form(m), data = long_cc(m), REML = FALSE)
})

reported_long <- c("years_from_baseline", "baseline_dxMCI", "baseline_dxAD",
                   "years_from_baseline:baseline_dxMCI", "years_from_baseline:baseline_dxAD")

table3 <- imap_dfr(long_fits_ri, function(fit, m) {
  tidy(fit, effects = "fixed") %>%
    mutate(marker = m, n_obs = nobs(fit), n_subj = ngrps(fit)[["rid"]])
}) %>%
  filter(term %in% reported_long) %>%
  mutate(
    effect_label = dplyr::recode(term,
      years_from_baseline                  = "Time, years (CN slope)",
      baseline_dxMCI                       = "MCI vs CN at baseline",
      baseline_dxAD                        = "AD vs CN at baseline",
      `years_from_baseline:baseline_dxMCI` = "Time x MCI",
      `years_from_baseline:baseline_dxAD`  = "Time x AD"),
    percent_change = (exp(estimate) - 1) * 100,
    pct_ci_low  = (exp(estimate - 1.96 * std.error) - 1) * 100,
    pct_ci_high = (exp(estimate + 1.96 * std.error) - 1) * 100,
    fdr_q = p.adjust(p.value, "BH")
  ) %>%
  left_join(marker_info, by = "marker") %>%
  select(label, effect_label, estimate, std.error, percent_change,
         pct_ci_low, pct_ci_high, p.value, fdr_q, n_obs, n_subj)

out$D_table3 <- table3

# Side-by-side with the published specification, for the response letter.
out$D_age_spec_comparison <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m)
  f_pub <- lmer(as.formula(paste0("ln_", m,
    " ~ years_from_baseline * baseline_dx + age + ptgender + pteducat + apoe4 + (1 | rid)")),
    data = d, REML = FALSE)
  bind_rows(
    tidy(f_pub, effects = "fixed") %>% mutate(spec = "published: time-varying age"),
    tidy(long_fits_ri[[m]], effects = "fixed") %>% mutate(spec = "corrected: baseline age")
  ) %>%
    filter(term %in% reported_long) %>%
    mutate(marker = m, percent_change = (exp(estimate) - 1) * 100)
})


# ==========================================================
# E. RANDOM INTERCEPT vs RANDOM SLOPE  (Reviewer 2, comment 10)
#
# Random slopes need within-person variation in time. The vascular
# panel has exactly two visits per participant, so a random slope
# is not identifiable there; the fit will be singular and that is
# the justification for keeping random intercepts.
# ==========================================================

re_comparison <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m)
  m_ri <- long_fits_ri[[m]]
  # Convergence warnings are expected here and are informative rather
  # than fatal, so only genuine errors abort the fit.
  m_rs <- tryCatch(
    suppressWarnings(
      lmer(long_form(m, "(1 + years_from_baseline | rid)"), data = d, REML = FALSE)),
    error = function(e) e
  )

  if (inherits(m_rs, "condition")) {
    return(tibble(marker = m, random_slope_status = paste("failed:", conditionMessage(m_rs)),
                  lrt_chisq = NA, lrt_df = NA, lrt_p = NA,
                  aic_ri = AIC(m_ri), aic_rs = NA, singular = NA,
                  preferred = "random intercept"))
  }

  lrt <- anova(m_ri, m_rs)
  sing <- isSingular(m_rs, tol = 1e-4)
  tibble(
    marker = m,
    random_slope_status = "converged",
    lrt_chisq = lrt$Chisq[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chisq)`[2],
    aic_ri = AIC(m_ri), aic_rs = AIC(m_rs),
    singular = sing,
    preferred = if (sing) "random intercept (slope model singular)"
                else if (!is.na(lrt$`Pr(>Chisq)`[2]) && lrt$`Pr(>Chisq)`[2] < 0.05)
                  "random slope" else "random intercept"
  )
})
out$E_random_effects <- re_comparison

# Where a random slope is preferred, report the fixed effects under it.
rs_markers <- re_comparison$marker[grepl("^random slope", re_comparison$preferred)]

out$E_random_slope_estimates <- if (length(rs_markers) == 0) {
  tibble(note = "No biomarker preferred a random slope over a random intercept.")
} else {
  map_dfr(rs_markers, function(m) {
    fit <- suppressWarnings(lmer(long_form(m, "(1 + years_from_baseline | rid)"),
                                 data = long_cc(m), REML = FALSE))
    tidy(fit, effects = "fixed") %>%
      filter(term %in% reported_long) %>%
      mutate(marker = m, percent_change = (exp(estimate) - 1) * 100)
  })
}


# ==========================================================
# F. RESIDUAL DIAGNOSTICS AND COLLINEARITY  (Reviewer 2, comment 10)
# ==========================================================

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4

resid_diag <- map_dfr(primary_markers, function(m) {
  rb <- residuals(baseline_fits[[m]])
  rl <- residuals(long_fits_ri[[m]])
  tibble(
    marker = m,
    baseline_resid_skew = skew(rb), baseline_resid_kurtosis = kurt(rb),
    baseline_ncv_p = tryCatch(ncvTest(baseline_fits[[m]])$p, error = function(e) NA_real_),
    longitudinal_resid_skew = skew(rl), longitudinal_resid_kurtosis = kurt(rl)
  )
})
out$F_residual_diagnostics <- resid_diag

vif_tbl <- map_dfr(primary_markers, function(m) {
  vb <- tryCatch(vif(baseline_fits[[m]]), error = function(e) NULL)
  vl <- tryCatch(vif(long_fits_ri[[m]]),  error = function(e) NULL)
  bind_rows(
    if (!is.null(vb)) tibble(marker = m, model = "baseline",
                             predictor = names(vb),
                             vif = if (is.matrix(vb)) vb[, 1] else as.numeric(vb)),
    if (!is.null(vl)) tibble(marker = m, model = "longitudinal",
                             predictor = names(vl),
                             vif = if (is.matrix(vl)) vl[, 1] else as.numeric(vl))
  )
})
out$F_vif <- vif_tbl

# QQ plots for the supplement
qq_df <- map_dfr(primary_markers, function(m) {
  r <- residuals(baseline_fits[[m]])
  tibble(marker = m, theoretical = qqnorm(r, plot.it = FALSE)$x,
         sample = qqnorm(r, plot.it = FALSE)$y)
}) %>% left_join(marker_info, by = "marker")

p_qq <- ggplot(qq_df, aes(theoretical, sample)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal QQ plots of baseline model residuals")

ggsave(file.path(results_dir, "FigureS3_residual_qq.tiff"), p_qq,
       width = 10, height = 6, dpi = 600, compression = "lzw")


# ==========================================================
# G. DISCRIMINATION  (Reviewer 2, comments 7 and 9)
#
# Percent change alone does not tell a reader whether a marker
# separates individuals. AUC, Cohen's d and the overlap
# coefficient do.
# ==========================================================

overlap_coef <- function(a, b) {
  rng <- range(c(a, b)); grid <- seq(rng[1], rng[2], length.out = 2048)
  da <- approx(density(a, from = rng[1], to = rng[2], n = 2048), xout = grid)$y
  db <- approx(density(b, from = rng[1], to = rng[2], n = 2048), xout = grid)$y
  sum(pmin(da, db)) * diff(grid)[1]
}

discrimination <- map_dfr(primary_markers, function(m) {
  d <- baseline_cc(m)
  y <- d[[paste0("ln_", m)]]

  # covariate-adjusted values, so AUC is comparable to the adjusted models
  adj <- residuals(lm(y ~ age + ptgender + pteducat + apoe4, data = d))

  map_dfr(list(c("AD", "CN"), c("MCI", "CN")), function(pair) {
    idx <- d$dx_label %in% pair
    if (sum(d$dx_label == pair[1]) < 10) return(NULL)
    lab <- factor(d$dx_label[idx], levels = rev(pair))
    r_raw <- roc(lab, y[idx],   quiet = TRUE, direction = "<")
    r_adj <- roc(lab, adj[idx], quiet = TRUE, direction = "<")
    ci_raw <- as.numeric(ci.auc(r_raw))
    g1 <- y[d$dx_label == pair[1]]; g2 <- y[d$dx_label == pair[2]]
    tibble(
      marker = m, contrast = paste(pair[1], "vs", pair[2]),
      n_group1 = length(g1), n_group2 = length(g2),
      auc_unadjusted = as.numeric(auc(r_raw)),
      auc_ci_low = ci_raw[1], auc_ci_high = ci_raw[3],
      auc_covariate_adjusted = as.numeric(auc(r_adj)),
      cohens_d = (mean(g1) - mean(g2)) /
        sqrt(((length(g1) - 1) * stats::var(g1) + (length(g2) - 1) * stats::var(g2)) /
               (length(g1) + length(g2) - 2)),
      overlap_coefficient = overlap_coef(g1, g2)
    )
  })
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

out$G_discrimination <- discrimination

# Observed power for the null baseline findings (Reviewer 2, comment 9):
# smallest AD-vs-CN difference detectable at 80% power in each panel.
power_tbl <- map_dfr(primary_markers, function(m) {
  d <- baseline_cc(m)
  n1 <- sum(d$dx_label == "AD"); n2 <- sum(d$dx_label == "CN")
  s  <- sd(d[[paste0("ln_", m)]])
  if (n1 < 5 || n2 < 5) return(NULL)
  d_detect <- (qnorm(0.975) + qnorm(0.80)) * sqrt(1 / n1 + 1 / n2)
  tibble(marker = m, n_AD = n1, n_CN = n2, residual_sd_log = s,
         min_detectable_log_diff = d_detect,
         min_detectable_percent = (exp(d_detect) - 1) * 100)
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$G_min_detectable_effect <- power_tbl


# ==========================================================
# H. CORRELATIONS WITH BOOTSTRAP CIs  (Reviewer 2, comment 8)
# ==========================================================

corr_vars <- paste0("ln_", primary_markers)
ci_input  <- baseline_dat %>% select(all_of(corr_vars))

boot_spearman_ci <- function(x, y, B = 2000) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(c(NA, NA, NA, n))
  rho <- cor(x, y, method = "spearman")
  bs <- replicate(B, { i <- sample.int(n, n, TRUE)
                       if (sd(x[i]) == 0 || sd(y[i]) == 0) NA_real_
                       else cor(x[i], y[i], method = "spearman") })
  c(rho, quantile(bs, 0.025, na.rm = TRUE), quantile(bs, 0.975, na.rm = TRUE), n)
}

pairs_grid <- t(combn(primary_markers, 2))
correlations <- map_dfr(seq_len(nrow(pairs_grid)), function(i) {
  a <- pairs_grid[i, 1]; b <- pairs_grid[i, 2]
  r <- boot_spearman_ci(ci_input[[paste0("ln_", a)]], ci_input[[paste0("ln_", b)]])
  ok <- complete.cases(ci_input[[paste0("ln_", a)]], ci_input[[paste0("ln_", b)]])
  tibble(marker_1 = a, marker_2 = b, spearman_rho = r[1],
         ci_low = r[2], ci_high = r[3], pairwise_n = r[4],
         p_value = suppressWarnings(cor.test(ci_input[[paste0("ln_", a)]][ok],
                                             ci_input[[paste0("ln_", b)]][ok],
                                             method = "spearman")$p.value))
}) %>%
  mutate(fdr_q = p.adjust(p_value, "BH"),
         ci_width = ci_high - ci_low,
         precision_flag = if_else(pairwise_n < 150,
                                  "underpowered, exploratory", "adequate")) %>%
  left_join(marker_info %>% select(marker, label_1 = label), by = c("marker_1" = "marker")) %>%
  left_join(marker_info %>% select(marker, label_2 = label), by = c("marker_2" = "marker")) %>%
  select(label_1, label_2, spearman_rho, ci_low, ci_high, pairwise_n,
         p_value, fdr_q, ci_width, precision_flag)

out$H_correlations <- correlations


# ==========================================================
# I. APOE e4 x DIAGNOSIS INTERACTION  (Reviewer 2, comment 11)
#    Sex x diagnosis added as well, given known GFAP sex effects.
# ==========================================================

interaction_tests <- map_dfr(primary_markers, function(m) {
  d <- baseline_cc(m)
  f0 <- lm(as.formula(paste0("ln_", m, " ~ dx_label + age + ptgender + pteducat + apoe4")),
           data = d)
  fa <- update(f0, . ~ . + dx_label:apoe4)
  fs <- update(f0, . ~ . + dx_label:ptgender)

  term_string <- function(fit, pattern) {
    tt <- tidy(fit) %>% filter(grepl(pattern, term))
    if (nrow(tt) == 0) return(NA_character_)
    paste(sprintf("%s: b=%.3f, p=%.3f", tt$term, tt$estimate, tt$p.value), collapse = "; ")
  }

  tibble(
    marker      = m,
    interaction = c("diagnosis x APOE e4", "diagnosis x sex"),
    lrt_p       = c(anova(f0, fa)$`Pr(>F)`[2], anova(f0, fs)$`Pr(>F)`[2]),
    terms       = c(term_string(fa, "dx_label.*apoe4"),
                    term_string(fs, "dx_label.*ptgender"))
  )
}) %>%
  mutate(fdr_q = p.adjust(lrt_p, "BH")) %>%
  left_join(marker_info, by = "marker") %>%
  select(label, interaction, lrt_p, fdr_q, terms)

out$I_interactions <- interaction_tests


# ==========================================================
# Save
# ==========================================================

# Any section that produced nothing is replaced with a note rather than
# an empty object, which openxlsx cannot write.
out <- imap(out, function(x, nm) {
  if (is.null(x) || (is.data.frame(x) && (nrow(x) == 0 || ncol(x) == 0)))
    tibble(note = paste("No rows produced for", nm)) else x
})

# Excel caps sheet names at 31 characters.
names(out) <- substr(make.unique(names(out)), 1, 31)

openxlsx::write.xlsx(out, file.path(results_dir, "corrected_primary_analysis.xlsx"),
                     overwrite = TRUE)

cat("\nDone.\n")
cat("Results: ", file.path(results_dir, "corrected_primary_analysis.xlsx"), "\n")
cat("Figure:  ", file.path(results_dir, "FigureS3_residual_qq.tiff"), "\n")
