# ==========================================================
# 17_fixes_and_risk_covariates.R
#
# Repairs three problems in the previous run and resolves the
# sICAM-1 discrepancy.
#
#   Part 1. Minimum detectable effect - formula corrected
#   Part 2. Random slopes - tested where they are identifiable
#   Part 3. Vascular assay QC flags, never used until now
#   Part 4. QC release vs ADMC release, visit by visit
#   Part 5. Vascular risk covariates from ADNIMERGE2 .rda tables
#
# Part 5 replaces script 16 entirely. Script 16 looked for CSVs;
# VITALS, MEDHIST and RECCMEDS live in the ADNIMERGE2 package as
# .rda files, so nothing loaded, and the creatinine fallback scan
# selected a 0/1 QC flag. The eGFR in vascular_risk_adjustment.xlsx
# is an artefact of age and sex only. Discard it.
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
            "broom", "broom.mixed")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse); library(janitor); library(openxlsx)
library(lme4); library(lmerTest); library(broom); library(broom.mixed)

project_dir <- "R:/ADNI_Project"
raw_dir     <- file.path(project_dir, "00_raw_data")
bio_dir     <- file.path(raw_dir, "biomarkers_excel")
rda_dir     <- file.path(raw_dir, "ADNIMERGE2", "data")
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

out <- list()
parse_num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))
first_nonmissing <- function(x) { x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else x[1] }

# ---------- Loader that handles both .rda and .csv ----------

load_adni <- function(stem) {
  f <- list.files(rda_dir, pattern = paste0("^", stem, "\\.rda$"),
                  full.names = TRUE, ignore.case = TRUE)
  if (length(f)) {
    env <- new.env(); nms <- load(f[1], envir = env)
    cat("Loaded (rda):", basename(f[1]), "\n")
    return(as_tibble(get(nms[1], envir = env)) %>% clean_names())
  }
  f <- list.files(raw_dir, pattern = paste0("^", stem, ".*\\.csv$"),
                  full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(f)) {
    cat("Loaded (csv):", basename(f[1]), "\n")
    return(read_csv(f[1], show_col_types = FALSE, guess_max = 100000) %>% clean_names())
  }
  cat("NOT FOUND:", stem, "\n"); NULL
}

# ---------- Analysis dataset with corrected outcomes ----------

marker_info <- tribble(
  ~marker,                ~label,           ~stored_scale,
  "gfap_quanterix",       "Plasma GFAP",    "raw",
  "strem2_msd_corrected", "Plasma sTREM2",  "raw",
  "vegf_plasma_qc",       "Plasma VEGF",    "log10",
  "sicam1_plasma_qc",     "Plasma sICAM-1", "log10",
  "svcam1_plasma_qc",     "Plasma sVCAM-1", "log10"
)
primary_markers  <- marker_info$marker
vascular_markers <- c("vegf_plasma_qc", "sicam1_plasma_qc", "svcam1_plasma_qc")

dat <- read_csv(file.path(clean_dir, "analysis_master_model_ready.csv"),
                show_col_types = FALSE, guess_max = 100000) %>%
  clean_names() %>%
  mutate(
    rid = as.character(rid), visit_key = as.character(visit_key),
    dx_label    = factor(dx_label,    levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender), age = parse_num(age),
    pteducat = parse_num(pteducat), apoe4 = parse_num(apoe4),
    years_from_baseline = parse_num(years_from_baseline)
  ) %>%
  mutate(baseline_age = age - years_from_baseline)

for (i in seq_len(nrow(marker_info))) {
  m <- marker_info$marker[i]; v <- parse_num(dat[[m]])
  dat[[paste0("ln_", m)]] <- if (marker_info$stored_scale[i] == "log10")
    log(10) * v else ifelse(!is.na(v) & v > 0, log(v), NA_real_)
}

baseline_cc <- function(m, d = dat) {
  d %>% filter(visit_key == "bl", !is.na(.data[[paste0("ln_", m)]]), !is.na(dx_label),
               !is.na(age), !is.na(ptgender), !is.na(pteducat), !is.na(apoe4))
}
long_cc <- function(m, d = dat) {
  d %>% filter(!is.na(.data[[paste0("ln_", m)]]), !is.na(years_from_baseline),
               !is.na(baseline_dx), !is.na(baseline_age), !is.na(ptgender),
               !is.na(pteducat), !is.na(apoe4))
}
reported_long <- c("years_from_baseline", "baseline_dxMCI", "baseline_dxAD",
                   "years_from_baseline:baseline_dxMCI", "years_from_baseline:baseline_dxAD")


# ==========================================================
# PART 1. MINIMUM DETECTABLE EFFECT, CORRECTED
#
# Previous run omitted the residual SD, so every panel returned
# the same number. Correct form:
#   delta_min = (z_0.975 + z_0.80) * s * sqrt(1/n1 + 1/n2)
# ==========================================================

cat("\n=== Part 1. Minimum detectable effect ===\n")

mde <- map_dfr(primary_markers, function(m) {
  d  <- baseline_cc(m)
  n1 <- sum(d$dx_label == "AD"); n2 <- sum(d$dx_label == "CN")
  if (n1 < 5 || n2 < 5) return(NULL)
  fit <- lm(as.formula(paste0("ln_", m, " ~ dx_label + age + ptgender + pteducat + apoe4")),
            data = d)
  s <- summary(fit)$sigma                      # residual SD, not raw SD
  delta <- (qnorm(0.975) + qnorm(0.80)) * s * sqrt(1 / n1 + 1 / n2)
  obs <- coef(fit)[["dx_labelAD"]]
  tibble(marker = m, n_AD = n1, n_CN = n2, residual_sd = s,
         min_detectable_log_diff = delta,
         min_detectable_percent  = (exp(delta) - 1) * 100,
         observed_percent        = (exp(obs) - 1) * 100,
         powered_for_observed    = abs(obs) >= delta)
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

print(as.data.frame(mde))
out$P1_min_detectable_effect <- mde


# ==========================================================
# PART 2. RANDOM SLOPES WHERE IDENTIFIABLE
#
# A random slope needs within-person variation in time. With 69%
# of GFAP participants contributing one observation, the full
# sample cannot support one. The test is therefore run on the
# subset with two or more measurements, and the inability to fit
# it on the full sample is reported as a finding.
# ==========================================================

cat("\n=== Part 2. Random slopes ===\n")

re_test <- map_dfr(primary_markers, function(m) {
  d_all <- long_cc(m)
  d2    <- d_all %>% group_by(rid) %>% filter(n() >= 2) %>% ungroup()
  n_rep <- n_distinct(d2$rid)

  # number of distinct time points per participant caps what is estimable
  max_times <- d2 %>% group_by(rid) %>% summarise(k = n_distinct(years_from_baseline),
                                                  .groups = "drop") %>% pull(k)
  n_3plus <- sum(max_times >= 3)

  base <- tibble(marker = m,
                 n_participants_full = n_distinct(d_all$rid),
                 n_participants_2plus = n_rep,
                 n_participants_3plus = n_3plus)

  if (n_rep < 30) return(base %>% mutate(
    result = "too few participants with repeat measures", lrt_p = NA_real_,
    conclusion = "random intercept retained"))

  f <- function(re) as.formula(paste0("ln_", m,
    " ~ years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat + apoe4 + ", re))

  m_ri <- lmer(f("(1 | rid)"), data = d2, REML = FALSE)
  m_rs <- tryCatch(suppressWarnings(
    lmer(f("(1 + years_from_baseline | rid)"), data = d2, REML = FALSE)),
    error = function(e) e)

  if (inherits(m_rs, "error")) return(base %>% mutate(
    result = paste("slope model not estimable:", conditionMessage(m_rs)),
    lrt_p = NA_real_, conclusion = "random intercept retained"))

  lrt  <- anova(m_ri, m_rs)
  sing <- isSingular(m_rs, tol = 1e-4)
  base %>% mutate(
    result = if (sing) "slope model singular" else "slope model converged",
    lrt_chisq = lrt$Chisq[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chisq)`[2],
    aic_intercept = AIC(m_ri), aic_slope = AIC(m_rs),
    conclusion = if (sing || is.na(lrt$`Pr(>Chisq)`[2]) || lrt$`Pr(>Chisq)`[2] >= 0.05)
      "random intercept retained" else "random slope preferred")
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

print(as.data.frame(re_test %>% select(label, n_participants_2plus,
                                       n_participants_3plus, result, lrt_p, conclusion)))
out$P2_random_slope_test <- re_test

# Fixed effects under a random slope, for any marker where it fits
out$P2_random_slope_fixed <- map_dfr(primary_markers, function(m) {
  d2 <- long_cc(m) %>% group_by(rid) %>% filter(n() >= 2) %>% ungroup()
  if (n_distinct(d2$rid) < 30) return(NULL)
  fit <- tryCatch(suppressWarnings(lmer(as.formula(paste0("ln_", m,
    " ~ years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat +",
    " apoe4 + (1 + years_from_baseline | rid)")), data = d2, REML = FALSE)),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tidy(fit, effects = "fixed") %>% filter(term %in% reported_long) %>%
    mutate(marker = m, percent_change = (exp(estimate) - 1) * 100,
           singular = isSingular(fit, tol = 1e-4))
})


# ==========================================================
# PART 3. VASCULAR ASSAY QC FLAGS
#
# The raw multiplex release carries below_ldd, read_low, read_hi
# and outler columns that the pipeline never used. If flag rates
# differ between the baseline and month-12 runs, that is a
# mechanism for the sICAM-1 shift.
# ==========================================================

cat("\n=== Part 3. Vascular QC flags ===\n")

raw_path <- file.path(bio_dir, "adni_plasma_raw_multiplex_11Nov2010.csv")

if (file.exists(raw_path)) {
  raw_mx <- read_csv(raw_path, show_col_types = FALSE, guess_max = 200000) %>% clean_names()

  analyte_map <- c(vegf_plasma_qc   = "vegf|vascular_endothelial_growth",
                   sicam1_plasma_qc = "icam|intercellular_adhesion",
                   svcam1_plasma_qc = "vcam|vascular_cell_adhesion")

  visit_col <- intersect(c("visit_code", "viscode2", "viscode"), names(raw_mx))[1]
  flag_cols <- intersect(c("below_ldd", "read_low", "read_hi", "outler", "outlier"),
                         names(raw_mx))
  cat("Flag columns present:", paste(flag_cols, collapse = ", "), "\n")

  raw_flags <- raw_mx %>%
    mutate(rid = as.character(rid),
           visit_key = as.character(.data[[visit_col]]),
           .a = tolower(as.character(analyte)),
           marker = case_when(
             str_detect(.a, analyte_map[["vegf_plasma_qc"]])   ~ "vegf_plasma_qc",
             str_detect(.a, analyte_map[["sicam1_plasma_qc"]]) ~ "sicam1_plasma_qc",
             str_detect(.a, analyte_map[["svcam1_plasma_qc"]]) ~ "svcam1_plasma_qc",
             TRUE ~ NA_character_)) %>%
    filter(!is.na(marker))

  # Flag rate by analyte and visit
  out$P3_flag_rates <- raw_flags %>%
    filter(visit_key %in% c("bl", "m12")) %>%
    group_by(marker, visit_key) %>%
    summarise(n = n(),
              across(all_of(flag_cols), ~ round(100 * mean(parse_num(.x) == 1, na.rm = TRUE), 1),
                     .names = "pct_{.col}"),
              .groups = "drop")
  print(as.data.frame(out$P3_flag_rates))

  # Any flag at all, per sample
  flag_any <- raw_flags %>%
    mutate(any_flag = as.integer(rowSums(
      across(all_of(flag_cols), ~ coalesce(parse_num(.x), 0) == 1)) > 0)) %>%
    group_by(rid, visit_key, marker) %>%
    summarise(any_flag = max(any_flag), .groups = "drop")

  out$P3_flag_by_group <- flag_any %>%
    left_join(dat %>% distinct(rid, baseline_dx), by = "rid") %>%
    filter(visit_key %in% c("bl", "m12"), !is.na(baseline_dx)) %>%
    count(marker, baseline_dx, visit_key, any_flag) %>%
    group_by(marker, baseline_dx, visit_key) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
    filter(any_flag == 1)

  # Refit the vascular models excluding flagged samples
  out$P3_models_excluding_flagged <- map_dfr(vascular_markers, function(m) {
    bad <- flag_any %>% filter(marker == m, any_flag == 1) %>%
      transmute(rid, visit_key, drop = TRUE)
    d <- long_cc(m) %>% left_join(bad, by = c("rid", "visit_key")) %>%
      filter(is.na(drop))
    if (n_distinct(d$rid) < 40) return(NULL)
    fit <- lmer(as.formula(paste0("ln_", m,
      " ~ years_from_baseline * baseline_dx + baseline_age + ptgender +",
      " pteducat + apoe4 + (1 | rid)")), data = d, REML = FALSE)
    tidy(fit, effects = "fixed") %>% filter(term %in% reported_long) %>%
      mutate(marker = m, analysis = "flagged samples excluded",
             percent_change = (exp(estimate) - 1) * 100,
             n_obs = nobs(fit), n_subj = ngrps(fit)[["rid"]])
  }) %>% left_join(marker_info, by = "marker") %>% relocate(label)

} else {
  cat("Raw multiplex file not found; skipping QC flag analysis.\n")
}


# ==========================================================
# PART 4. QC RELEASE vs ADMC RELEASE, VISIT BY VISIT
#
# The two releases disagreed about the 12-month change. This
# establishes whether they are the same measurement at all, and
# if so where they diverge.
# ==========================================================

cat("\n=== Part 4. QC vs ADMC ===\n")

admc_path <- file.path(bio_dir, "ADMC_CLINICALVARIABLES_16May2016.csv")
if (file.exists(admc_path)) {
  admc <- read_csv(admc_path, show_col_types = FALSE, guess_max = 100000) %>% clean_names()
  kv <- intersect(c("viscode2", "viscode"), names(admc))[1]

  admc_long <- admc %>%
    transmute(rid = as.character(rid), visit_key = as.character(.data[[kv]]),
              vegf_plasma_qc = parse_num(vegf), sicam1_plasma_qc = parse_num(icam),
              svcam1_plasma_qc = parse_num(vcam)) %>%
    group_by(rid, visit_key) %>%
    summarise(across(everything(), ~ first_nonmissing(.x)), .groups = "drop") %>%
    mutate(across(all_of(vascular_markers), as.numeric)) %>%
    pivot_longer(all_of(vascular_markers), names_to = "marker", values_to = "admc_value")

  qc_long <- dat %>%
    select(rid, visit_key, all_of(vascular_markers)) %>%
    pivot_longer(all_of(vascular_markers), names_to = "marker", values_to = "qc_log10")

  cmp <- inner_join(qc_long, admc_long, by = c("rid", "visit_key", "marker")) %>%
    filter(!is.na(qc_log10), !is.na(admc_value), admc_value > 0)

  # Is the ADMC release raw or already log-transformed?
  out$P4_admc_scale <- cmp %>% group_by(marker) %>%
    summarise(n = n(), median_admc = median(admc_value), median_qc_log10 = median(qc_log10),
              median_qc_backtransformed = 10^median(qc_log10),
              cor_with_qc      = cor(qc_log10, log10(admc_value)),
              cor_with_qc_raw  = cor(qc_log10, admc_value), .groups = "drop") %>%
    mutate(admc_appears_to_be = if_else(abs(median_admc - median_qc_backtransformed) /
                                          median_qc_backtransformed < 0.2,
                                        "raw concentration", "check manually"))
  print(as.data.frame(out$P4_admc_scale))

  # Agreement separately at baseline and month 12
  out$P4_agreement_by_visit <- cmp %>%
    filter(visit_key %in% c("bl", "m12")) %>%
    group_by(marker, visit_key) %>%
    summarise(n = n(),
              cor_log_scale = cor(qc_log10, log10(admc_value)),
              median_qc_backtransformed = 10^median(qc_log10),
              median_admc = median(admc_value),
              median_ratio_qc_over_admc = 10^median(qc_log10) / median(admc_value),
              .groups = "drop")
  print(as.data.frame(out$P4_agreement_by_visit))

  # Paired 12-month change computed identically in both releases
  out$P4_paired_change_both <- cmp %>%
    filter(visit_key %in% c("bl", "m12")) %>%
    left_join(dat %>% distinct(rid, baseline_dx), by = "rid") %>%
    mutate(ln_qc = log(10) * qc_log10, ln_admc = log(admc_value)) %>%
    select(rid, marker, baseline_dx, visit_key, ln_qc, ln_admc) %>%
    pivot_longer(c(ln_qc, ln_admc), names_to = "release", values_to = "ln_value") %>%
    pivot_wider(names_from = visit_key, values_from = ln_value) %>%
    filter(!is.na(bl), !is.na(m12)) %>%
    group_by(marker, release, baseline_dx) %>%
    summarise(n = n(), percent_change_12m = (exp(mean(m12 - bl)) - 1) * 100,
              t_p = t.test(m12 - bl)$p.value, .groups = "drop")
  print(as.data.frame(out$P4_paired_change_both))
}


# ==========================================================
# PART 5. VASCULAR RISK COVARIATES FROM ADNIMERGE2
# ==========================================================

cat("\n=== Part 5. Risk covariates ===\n")

# Inventory of everything available, so nothing is missed twice
rda_files <- list.files(rda_dir, pattern = "\\.rda$", ignore.case = TRUE)
out$P5_available_tables <- tibble(table = sort(sub("\\.rda$", "", rda_files)))
cat("ADNIMERGE2 tables available:", length(rda_files), "\n")

# --- VITALS: BMI and blood pressure ---
vitals <- load_adni("VITALS")
vitals_bl <- NULL
if (!is.null(vitals)) {
  kv <- intersect(c("viscode2", "viscode"), names(vitals))[1]
  g  <- function(n) if (n %in% names(vitals)) parse_num(vitals[[n]]) else NA_real_
  vitals_bl <- vitals %>%
    mutate(rid = as.character(rid), visit_key = as.character(.data[[kv]]),
           wt = g("vsweight"), wtu = g("vswtunit"),
           ht = g("vsheight"), htu = g("vshtunit"),
           sbp = g("vsbpsys"), dbp = g("vsbpdia")) %>%
    filter(visit_key %in% c("bl", "sc", "scmri")) %>%
    mutate(weight_kg = if_else(!is.na(wtu) & wtu == 1, wt * 0.45359237, wt),
           height_m  = if_else(!is.na(htu) & htu == 1, ht * 0.0254, ht / 100),
           bmi = weight_kg / height_m^2) %>%
    mutate(bmi = if_else(bmi > 12 & bmi < 70, bmi, NA_real_),
           sbp = if_else(sbp > 60 & sbp < 260, sbp, NA_real_),
           dbp = if_else(dbp > 30 & dbp < 160, dbp, NA_real_)) %>%
    group_by(rid) %>%
    summarise(bmi = first_nonmissing(bmi), sbp = first_nonmissing(sbp),
              dbp = first_nonmissing(dbp), .groups = "drop") %>%
    mutate(across(c(bmi, sbp, dbp), as.numeric))
  cat("  BMI:", sum(!is.na(vitals_bl$bmi)), " SBP:", sum(!is.na(vitals_bl$sbp)), "\n")
}

# --- MEDHIST: organ-system history ---
medhist <- load_adni("MEDHIST")
medhist_bl <- NULL
if (!is.null(medhist)) {
  g <- function(n) if (n %in% names(medhist)) parse_num(medhist[[n]]) else NA_real_
  medhist_bl <- medhist %>%
    mutate(rid = as.character(rid),
           mh_cardiovascular = g("mh4card"), mh_endocrine = g("mh9endo"),
           mh_renal = g("mh12rena"), mh_smoking = g("mh16smok")) %>%
    group_by(rid) %>%
    summarise(across(starts_with("mh_"), ~ first_nonmissing(.x)), .groups = "drop") %>%
    mutate(across(starts_with("mh_"), ~ as.integer(as.numeric(.x) == 1)))
  cat("  cardiovascular:", sum(medhist_bl$mh_cardiovascular, na.rm = TRUE),
      " endocrine:", sum(medhist_bl$mh_endocrine, na.rm = TRUE),
      " smoking:", sum(medhist_bl$mh_smoking, na.rm = TRUE), "\n")
}

# --- RECCMEDS: medication classes ---
meds <- load_adni("RECCMEDS")
meds_bl <- NULL
if (!is.null(meds)) {
  mc <- intersect(c("cmmed", "cmmeds", "medname", "cmtrt"), names(meds))[1]
  cat("  medication column:", mc, "\n")
  ahtn <- paste0("lisinopril|enalapril|ramipril|benazepril|captopril|quinapril|perindopril|",
    "losartan|valsartan|irbesartan|candesartan|olmesartan|telmisartan|amlodipine|nifedipine|",
    "diltiazem|verapamil|felodipine|metoprolol|atenolol|carvedilol|bisoprolol|propranolol|",
    "nebivolol|labetalol|hydrochlorothiazide|hctz|chlorthalidone|furosemide|indapamide|",
    "spironolactone|triamterene|clonidine|doxazosin|terazosin|hydralazine")
  adm <- paste0("metformin|glipizide|glyburide|glimepiride|glucophage|insulin|lantus|humalog|",
    "novolog|levemir|pioglitazone|rosiglitazone|actos|avandia|sitagliptin|januvia|saxagliptin|",
    "linagliptin|exenatide|liraglutide|byetta|victoza|empagliflozin|canagliflozin|",
    "dapagliflozin|acarbose|repaglinide|nateglinide")
  meds_bl <- meds %>%
    mutate(rid = as.character(rid), med = tolower(as.character(.data[[mc]]))) %>%
    group_by(rid) %>%
    summarise(med_antihypertensive = as.integer(any(str_detect(med, ahtn), na.rm = TRUE)),
              med_antidiabetic     = as.integer(any(str_detect(med, adm),  na.rm = TRUE)),
              .groups = "drop")
  cat("  antihypertensive:", sum(meds_bl$med_antihypertensive),
      " antidiabetic:", sum(meds_bl$med_antidiabetic), "\n")
}

# --- Creatinine: search every table, reject flag-like columns ---
cat("\nSearching all tables for creatinine...\n")
creat_bl <- NULL
creat_report <- list()

plausible_creatinine <- function(v) {
  v <- v[!is.na(v) & v > 0]
  length(v) > 100 && n_distinct(v) > 20 &&
    median(v) > 0.4 && median(v) < 2.0 && max(v) < 20
}

for (stem in sub("\\.rda$", "", rda_files)) {
  tb <- tryCatch(load_adni(stem), error = function(e) NULL)
  if (is.null(tb) || !"rid" %in% names(tb)) next

  # (a) a column named like creatinine
  hits <- grep("creat", names(tb), value = TRUE, ignore.case = TRUE)
  for (h in hits) {
    v <- parse_num(tb[[h]])
    creat_report[[length(creat_report) + 1]] <- tibble(
      table = stem, column = h, source = "name match", n = sum(!is.na(v) & v > 0),
      median = suppressWarnings(median(v[!is.na(v) & v > 0])),
      distinct = n_distinct(v), plausible = plausible_creatinine(v))
    if (is.null(creat_bl) && plausible_creatinine(v)) {
      kv <- intersect(c("viscode2", "viscode"), names(tb))[1]
      creat_bl <- tb %>% mutate(rid = as.character(rid), creatinine = parse_num(.data[[h]])) %>%
        { if (!is.na(kv)) filter(., as.character(.data[[kv]]) %in% c("bl", "sc", "scmri")) else . } %>%
        filter(!is.na(creatinine), creatinine > 0.1, creatinine < 20) %>%
        group_by(rid) %>% summarise(creatinine = first(creatinine), .groups = "drop")
      cat("  USING:", stem, "->", h, "(n =", nrow(creat_bl), ")\n")
    }
  }

  # (b) long-format labs with an analyte-name column
  nm_col <- intersect(c("testname", "test_name", "analyte", "biomarker", "loniuid"), names(tb))[1]
  vl_col <- intersect(c("value", "result", "testvalue", "avalue"), names(tb))[1]
  if (!is.na(nm_col) && !is.na(vl_col)) {
    sub <- tb %>% filter(str_detect(tolower(as.character(.data[[nm_col]])), "creatinine"))
    if (nrow(sub) > 50) {
      v <- parse_num(sub[[vl_col]])
      creat_report[[length(creat_report) + 1]] <- tibble(
        table = stem, column = paste0(nm_col, "/", vl_col), source = "long format",
        n = sum(!is.na(v) & v > 0), median = suppressWarnings(median(v[!is.na(v) & v > 0])),
        distinct = n_distinct(v), plausible = plausible_creatinine(v))
      if (is.null(creat_bl) && plausible_creatinine(v)) {
        creat_bl <- sub %>% mutate(rid = as.character(rid), creatinine = parse_num(.data[[vl_col]])) %>%
          filter(!is.na(creatinine), creatinine > 0.1, creatinine < 20) %>%
          group_by(rid) %>% summarise(creatinine = first(creatinine), .groups = "drop")
        cat("  USING:", stem, "long format (n =", nrow(creat_bl), ")\n")
      }
    }
  }
}

out$P5_creatinine_search <- if (length(creat_report))
  bind_rows(creat_report) else tibble(note = "no creatinine-like column found in any table")
print(as.data.frame(out$P5_creatinine_search))

# --- Assemble ---
risk <- dat %>% distinct(rid, baseline_age, ptgender)
for (tb in list(vitals_bl, medhist_bl, meds_bl, creat_bl))
  if (!is.null(tb)) risk <- left_join(risk, tb, by = "rid")

if ("creatinine" %in% names(risk)) {
  risk <- risk %>%
    mutate(fem = tolower(as.character(ptgender)) %in% c("female", "f", "2"),
           k = if_else(fem, 0.7, 0.9), a = if_else(fem, -0.241, -0.302),
           egfr = 142 * pmin(creatinine / k, 1)^a * pmax(creatinine / k, 1)^(-1.200) *
                  0.9938^baseline_age * if_else(fem, 1.012, 1)) %>%
    select(-fem, -k, -a)
}

# Composite flags, kept NA when no source contributed
comb <- function(...) {
  cols <- list(...)
  cols <- cols[!map_lgl(cols, ~ all(is.na(.x)))]
  if (!length(cols)) return(NA_integer_)
  reduce(cols, function(a, b) pmax(coalesce(a, 0L), coalesce(b, 0L)))
}
risk <- risk %>%
  mutate(
    hypertension = comb(
      if ("med_antihypertensive" %in% names(.)) med_antihypertensive else NA_integer_,
      if ("mh_cardiovascular" %in% names(.)) mh_cardiovascular else NA_integer_,
      if ("sbp" %in% names(.)) as.integer(sbp >= 140) else NA_integer_),
    diabetes = comb(
      if ("med_antidiabetic" %in% names(.)) med_antidiabetic else NA_integer_,
      if ("mh_endocrine" %in% names(.)) mh_endocrine else NA_integer_))

out$P5_covariate_availability <- risk %>%
  summarise(across(everything(), ~ sum(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_non_missing")
print(as.data.frame(out$P5_covariate_availability))
out$P5_risk_table <- risk

# --- Refit vascular models, unadjusted vs adjusted ---
dat_risk <- dat %>% left_join(risk %>% select(-baseline_age, -ptgender), by = "rid")

extra <- c("bmi", "hypertension", "diabetes", "mh_smoking", "egfr")
extra <- extra[extra %in% names(dat_risk)]
extra <- extra[map_lgl(extra, ~ sum(!is.na(dat_risk[[.x]])) > 200 &&
                             n_distinct(dat_risk[[.x]], na.rm = TRUE) > 1)]
cat("\nCovariates entering the models:", paste(extra, collapse = ", "), "\n")
out$P5_covariates_used <- tibble(covariate = if (length(extra)) extra else "none available")

if (length(extra)) {
  et <- paste("+", paste(extra, collapse = " + "))

  out$P5_baseline_adj <- map_dfr(vascular_markers, function(m) {
    d0 <- baseline_cc(m, dat_risk)
    d1 <- d0 %>% drop_na(all_of(extra))
    if (nrow(d1) < 50) return(NULL)
    f_un  <- lm(as.formula(paste0("ln_", m, " ~ dx_label + age + ptgender + pteducat + apoe4")),
                data = d1)
    f_adj <- lm(as.formula(paste0("ln_", m, " ~ dx_label + age + ptgender + pteducat + apoe4", et)),
                data = d1)
    bind_rows(tidy(f_un)  %>% mutate(model = "as published"),
              tidy(f_adj) %>% mutate(model = "plus vascular risk")) %>%
      filter(term %in% c("dx_labelMCI", "dx_labelAD")) %>%
      mutate(marker = m, percent_change = (exp(estimate) - 1) * 100,
             n = nrow(d1), n_lost_to_covariates = nrow(d0) - nrow(d1))
  }) %>% left_join(marker_info, by = "marker") %>% relocate(label)

  out$P5_longitudinal_adj <- map_dfr(vascular_markers, function(m) {
    d0 <- long_cc(m, dat_risk); d1 <- d0 %>% drop_na(all_of(extra))
    if (n_distinct(d1$rid) < 40) return(NULL)
    rhs <- "years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat + apoe4"
    f_un  <- lmer(as.formula(paste0("ln_", m, " ~ ", rhs, " + (1 | rid)")), data = d1, REML = FALSE)
    f_adj <- lmer(as.formula(paste0("ln_", m, " ~ ", rhs, et, " + (1 | rid)")), data = d1, REML = FALSE)
    bind_rows(tidy(f_un,  effects = "fixed") %>% mutate(model = "as published"),
              tidy(f_adj, effects = "fixed") %>% mutate(model = "plus vascular risk")) %>%
      filter(term %in% reported_long) %>%
      mutate(marker = m, percent_change = (exp(estimate) - 1) * 100,
             n_obs = nrow(d1), n_subj = n_distinct(d1$rid))
  }) %>% left_join(marker_info, by = "marker") %>% relocate(label)

  out$P5_attenuation <- out$P5_baseline_adj %>%
    select(label, term, model, estimate) %>%
    pivot_wider(names_from = model, values_from = estimate) %>%
    mutate(percent_attenuation = 100 * (1 - `plus vascular risk` / `as published`))

  out$P5_risk_by_diagnosis <- dat_risk %>%
    filter(visit_key == "bl", !is.na(dx_label)) %>% distinct(rid, .keep_all = TRUE) %>%
    group_by(dx_label) %>%
    summarise(n = n(), across(all_of(extra), ~ round(mean(.x, na.rm = TRUE), 2)), .groups = "drop")
  print(as.data.frame(out$P5_risk_by_diagnosis))
}


# ==========================================================
out <- imap(out, function(x, nm) {
  if (is.null(x) || (is.data.frame(x) && (nrow(x) == 0 || ncol(x) == 0)))
    tibble(note = paste("No rows produced for", nm)) else x
})
names(out) <- substr(make.unique(names(out)), 1, 31)

openxlsx::write.xlsx(out, file.path(results_dir, "fixes_and_risk_covariates.xlsx"),
                     overwrite = TRUE)
cat("\nDone:", file.path(results_dir, "fixes_and_risk_covariates.xlsx"), "\n")
