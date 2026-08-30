# ==========================================================
# 13_revision_diagnostics.R
#
# Revision step 1: diagnostics only. Nothing here changes the
# published models. It establishes the facts needed to decide
# what the corrected models should look like.
#
# Sections
#   A. Vascular biomarker scale verification  (R1-5, R2-7)
#   B. Cohort structure and follow-up          (R1-1)
#   C. Baseline vs longitudinal reconciliation (R1-2)
#   D. Biomarker-specific sample descriptions  (R1-3)
#   E. Included vs excluded comparison         (R2-5)
#   F. Age / time collinearity check           (new)
#   G. FDR family sensitivity                  (R2-4)
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
            "readr", "broom", "broom.mixed")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse)
library(janitor)
library(openxlsx)
library(lme4)
library(lmerTest)
library(readr)
library(broom)
library(broom.mixed)

project_dir <- "R:/ADNI_Project"

bio_dir     <- file.path(project_dir, "00_raw_data", "biomarkers_excel")
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

out <- list()   # everything gets written to one workbook at the end

# ---------- Helpers ----------

parse_num <- function(x) readr::parse_number(as.character(x))

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  x[1]
}

get_col <- function(df, col) {
  if (col %in% names(df)) df[[col]] else rep(NA, nrow(df))
}

coalesce_existing <- function(df, cols) {
  existing <- intersect(cols, names(df))
  if (length(existing) == 0) return(rep(NA_character_, nrow(df)))
  purrr::reduce(lapply(existing, function(x) as.character(df[[x]])), dplyr::coalesce)
}

standardize_keys <- function(df) {
  df <- clean_names(df)
  df %>% mutate(
    rid = as.character(get_col(., "rid")),
    visit_key = coalesce_existing(., c("viscode2", "viscode", "visit_code", "visit"))
  )
}

primary_markers <- c(
  "gfap_quanterix", "strem2_msd_corrected",
  "vegf_plasma_qc", "sicam1_plasma_qc", "svcam1_plasma_qc"
)
vascular_markers <- c("vegf_plasma_qc", "sicam1_plasma_qc", "svcam1_plasma_qc")

# ---------- Load the model-ready dataset ----------

dat <- read_csv(
  file.path(clean_dir, "analysis_master_model_ready.csv"),
  show_col_types = FALSE, guess_max = 100000
) %>%
  clean_names() %>%
  mutate(
    rid       = as.character(rid),
    visit_key = as.character(visit_key),
    dx_label    = factor(dx_label,    levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender  = factor(ptgender),
    age       = parse_num(age),
    pteducat  = parse_num(pteducat),
    apoe4     = parse_num(apoe4),
    years_from_baseline = parse_num(years_from_baseline)
  )

baseline_dat <- dat %>%
  filter(visit_key == "bl", dx_label %in% c("CN", "MCI", "AD"))

make_baseline_cc <- function(m) {
  baseline_dat %>%
    filter(!is.na(.data[[m]]), .data[[m]] > 0, !is.na(dx_label),
           !is.na(age), !is.na(ptgender), !is.na(pteducat), !is.na(apoe4))
}

make_long_cc <- function(m) {
  dat %>%
    filter(!is.na(.data[[m]]), .data[[m]] > 0, !is.na(years_from_baseline),
           !is.na(baseline_dx), !is.na(age), !is.na(ptgender),
           !is.na(pteducat), !is.na(apoe4))
}


# ==========================================================
# A. VASCULAR SCALE VERIFICATION
#
# Hypothesis: the QC multiplex file stores log10(concentration)
# even though the column names carry the original units. If so,
# safe_log() produced ln(log10(conc)) and every reported percent
# change for VEGF / sICAM-1 / sVCAM-1 is on the wrong scale.
# ==========================================================

cat("\n=== A. Vascular scale verification ===\n")

# A1. Observed range of the values actually used in the models.
#     Reference plasma ranges: VEGF ~50-800 pg/mL,
#     sICAM-1 ~100-400 ng/mL, sVCAM-1 ~400-1200 ng/mL.
#     If as_stored is ~2-3 and back_10 lands in those ranges,
#     the values are log10.

vascular_scale_check <- map_dfr(vascular_markers, function(m) {
  v <- dat[[m]]
  v <- v[!is.na(v)]
  tibble(
    marker          = m,
    n_values        = length(v),
    min_as_stored   = min(v),
    median_as_stored = median(v),
    max_as_stored   = max(v),
    median_back_10  = 10^median(v),   # if stored is log10
    median_back_e   = exp(median(v))  # if stored is natural log
  )
})
print(as.data.frame(vascular_scale_check))
out$A1_vascular_scale <- vascular_scale_check

# A2. Direct comparison against the raw (untransformed) multiplex file.
#     If the QC file is log10 of the raw file, then
#     cor(qc_value, log10(raw_value)) should be ~1 and the
#     regression slope of qc on log10(raw) should be ~1.

raw_path <- file.path(bio_dir, "adni_plasma_raw_multiplex_11Nov2010.csv")

if (file.exists(raw_path)) {

  raw_mx <- read_csv(raw_path, show_col_types = FALSE, guess_max = 100000) %>%
    standardize_keys()

  cat("\nRaw multiplex file columns:\n")
  print(names(raw_mx))

  # The raw release may be wide (one column per analyte) or long
  # (an analyte-name column plus a value column). Handle both.
  analyte_patterns <- c(
    vegf_plasma_qc   = "vegf|vascular_endothelial_growth",
    sicam1_plasma_qc = "icam|intercellular_adhesion",
    svcam1_plasma_qc = "vcam|vascular_cell_adhesion"
  )

  raw_long <- NULL

  wide_hits <- map(analyte_patterns, ~ grep(.x, names(raw_mx), value = TRUE, ignore.case = TRUE))

  if (all(lengths(wide_hits) > 0)) {
    raw_long <- imap_dfr(wide_hits, function(cols, marker) {
      tibble(
        rid = raw_mx$rid,
        visit_key = raw_mx$visit_key,
        marker = marker,
        raw_value = parse_num(raw_mx[[cols[1]]])
      )
    })
  } else {
    name_col  <- grep("analyte|test|protein|variable|item", names(raw_mx),
                      value = TRUE, ignore.case = TRUE)[1]
    value_col <- grep("value|result|conc", names(raw_mx),
                      value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(name_col) && !is.na(value_col)) {
      raw_long <- raw_mx %>%
        mutate(.nm = tolower(as.character(.data[[name_col]]))) %>%
        mutate(marker = case_when(
          grepl(analyte_patterns[["vegf_plasma_qc"]],   .nm) ~ "vegf_plasma_qc",
          grepl(analyte_patterns[["sicam1_plasma_qc"]], .nm) ~ "sicam1_plasma_qc",
          grepl(analyte_patterns[["svcam1_plasma_qc"]], .nm) ~ "svcam1_plasma_qc",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(marker)) %>%
        transmute(rid, visit_key, marker, raw_value = parse_num(.data[[value_col]]))
    } else {
      cat("\nCould not identify analyte/value columns in the raw file.\n",
          "Inspect names(raw_mx) above and set them manually.\n")
    }
  }

  if (!is.null(raw_long)) {

    qc_long <- dat %>%
      select(rid, visit_key, all_of(vascular_markers)) %>%
      pivot_longer(all_of(vascular_markers), names_to = "marker", values_to = "qc_value")

    scale_cmp <- raw_long %>%
      filter(!is.na(raw_value), raw_value > 0) %>%
      group_by(rid, visit_key, marker) %>%
      summarise(raw_value = first_nonmissing(raw_value), .groups = "drop") %>%
      mutate(raw_value = as.numeric(raw_value)) %>%
      inner_join(qc_long, by = c("rid", "visit_key", "marker")) %>%
      filter(!is.na(qc_value))

    scale_verdict <- scale_cmp %>%
      group_by(marker) %>%
      summarise(
        n_matched          = n(),
        median_raw         = median(raw_value),
        median_qc          = median(qc_value),
        cor_qc_vs_log10raw = cor(qc_value, log10(raw_value), use = "complete.obs"),
        cor_qc_vs_raw      = cor(qc_value, raw_value,        use = "complete.obs"),
        slope_qc_on_log10raw = coef(lm(qc_value ~ log10(raw_value)))[2],
        .groups = "drop"
      ) %>%
      mutate(
        verdict = if_else(
          cor_qc_vs_log10raw > 0.98 & abs(slope_qc_on_log10raw - 1) < 0.10,
          "QC file IS log10(raw) -> double log confirmed",
          "inconclusive, inspect manually"
        )
      )

    cat("\n")
    print(as.data.frame(scale_verdict))
    out$A2_scale_verdict <- scale_verdict
    out$A2_scale_pairs   <- head(scale_cmp, 2000)
  }

} else {
  cat("\nRaw multiplex file not found at:\n", raw_path,
      "\nDownload 'adni_plasma_raw_multiplex_11Nov2010.csv' from LONI,",
      "\nor rely on the A1 back-transformation check alone.\n")
}

# A3. What the corrected effect sizes would look like.
#     If the outcome should be ln(conc) = ln(10) * qc_value, then
#     refitting on ln(10)*qc_value gives beta on a true log-concentration
#     scale. This block previews the magnitude; section-by-section
#     refits come in the next script.

preview_corrected <- map_dfr(vascular_markers, function(m) {
  df <- make_baseline_cc(m)
  fit_old <- lm(log(df[[m]]) ~ dx_label + age + ptgender + pteducat + apoe4, data = df)
  fit_new <- lm(I(log(10) * df[[m]]) ~ dx_label + age + ptgender + pteducat + apoe4, data = df)
  bind_rows(
    broom::tidy(fit_old) %>% mutate(scale = "as published: ln(log10 conc)"),
    broom::tidy(fit_new) %>% mutate(scale = "corrected: ln(conc)")
  ) %>%
    filter(term %in% c("dx_labelMCI", "dx_labelAD")) %>%
    mutate(marker = m, percent_change = (exp(estimate) - 1) * 100)
})
cat("\nPreview of corrected vascular effect sizes:\n")
print(as.data.frame(preview_corrected %>%
  select(marker, scale, term, estimate, percent_change, p.value)))
out$A3_corrected_preview <- preview_corrected


# ==========================================================
# B. COHORT STRUCTURE AND FOLLOW-UP  (Reviewer 1, comment 1)
# ==========================================================

cat("\n=== B. Visits per participant and follow-up ===\n")

visits_per_participant <- map_dfr(primary_markers, function(m) {
  make_long_cc(m) %>%
    group_by(rid, baseline_dx) %>%
    summarise(
      n_visits = n(),
      followup_years = max(years_from_baseline) - min(years_from_baseline),
      .groups = "drop"
    ) %>%
    mutate(marker = m)
})

visit_count_distribution <- visits_per_participant %>%
  mutate(visit_group = case_when(
    n_visits == 1 ~ "1 measurement",
    n_visits == 2 ~ "2 measurements",
    TRUE          ~ "3 or more measurements"
  )) %>%
  count(marker, baseline_dx, visit_group) %>%
  group_by(marker, baseline_dx) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(marker, baseline_dx, visit_group)

followup_summary <- visits_per_participant %>%
  group_by(marker, baseline_dx) %>%
  summarise(
    n_participants   = n(),
    median_visits    = median(n_visits),
    mean_visits      = round(mean(n_visits), 2),
    median_followup  = round(median(followup_years), 2),
    q1_followup      = round(quantile(followup_years, 0.25), 2),
    q3_followup      = round(quantile(followup_years, 0.75), 2),
    max_followup     = round(max(followup_years), 2),
    pct_with_2plus   = round(100 * mean(n_visits >= 2), 1),
    .groups = "drop"
  )

print(as.data.frame(followup_summary))
out$B1_visit_distribution <- visit_count_distribution
out$B2_followup_summary   <- followup_summary
out$B3_visit_codes <- map_dfr(primary_markers, function(m) {
  make_long_cc(m) %>% count(marker = m, visit_key, name = "n_observations")
})

# Attrition by diagnosis and visit (Reviewer 2, comment 2)
attrition <- map_dfr(primary_markers, function(m) {
  make_long_cc(m) %>%
    count(marker = m, baseline_dx, visit_key, name = "n_observations")
}) %>%
  group_by(marker, baseline_dx) %>%
  mutate(pct_of_group_observations = round(100 * n_observations / sum(n_observations), 1)) %>%
  ungroup()
out$B4_attrition_by_visit <- attrition


# ==========================================================
# C. WHY LONGITUDINAL n > BASELINE n  (Reviewer 1, comment 2)
# ==========================================================

cat("\n=== C. Baseline vs longitudinal sample reconciliation ===\n")

sample_reconciliation <- map_dfr(primary_markers, function(m) {
  b <- make_baseline_cc(m)
  l <- make_long_cc(m)
  b_ids <- unique(b$rid)
  l_ids <- unique(l$rid)
  tibble(
    marker = m,
    baseline_participants          = length(b_ids),
    longitudinal_participants      = length(l_ids),
    longitudinal_observations      = nrow(l),
    in_long_with_baseline_value    = sum(l_ids %in% b_ids),
    in_long_without_baseline_value = sum(!l_ids %in% b_ids),
    pct_without_baseline_value     = round(100 * mean(!l_ids %in% b_ids), 1),
    earliest_visit_if_no_baseline  = {
      no_bl <- l %>% filter(!rid %in% b_ids)
      if (nrow(no_bl) == 0) NA_real_ else
        round(median(no_bl %>% group_by(rid) %>%
                       summarise(t = min(years_from_baseline), .groups = "drop") %>%
                       pull(t)), 2)
    }
  )
})

print(as.data.frame(sample_reconciliation))
out$C1_sample_reconciliation <- sample_reconciliation

# Diagnosis mix of participants entering the longitudinal model
# without a baseline biomarker value.
no_baseline_profile <- map_dfr(primary_markers, function(m) {
  b_ids <- unique(make_baseline_cc(m)$rid)
  make_long_cc(m) %>%
    distinct(rid, baseline_dx) %>%
    mutate(has_baseline_value = rid %in% b_ids, marker = m) %>%
    count(marker, baseline_dx, has_baseline_value)
})
out$C2_no_baseline_profile <- no_baseline_profile


# ==========================================================
# D. BIOMARKER-SPECIFIC SAMPLE DESCRIPTIONS  (Reviewer 1, comment 3)
#    Draft of a new supplementary table.
# ==========================================================

cat("\n=== D. Per-panel demographics ===\n")

describe_sample <- function(df, label) {
  df %>%
    group_by(dx_group = .data[[intersect(c("dx_label", "baseline_dx"), names(df))[1]]]) %>%
    summarise(
      n            = n_distinct(rid),
      age_mean     = round(mean(age, na.rm = TRUE), 1),
      age_sd       = round(sd(age, na.rm = TRUE), 1),
      female_pct   = round(100 * mean(tolower(as.character(ptgender)) %in% c("female", "f", "2")), 1),
      educ_mean    = round(mean(pteducat, na.rm = TRUE), 1),
      educ_sd      = round(sd(pteducat, na.rm = TRUE), 1),
      apoe4_pct    = round(100 * mean(apoe4 == 1, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(sample = label)
}

panel_demographics <- map_dfr(primary_markers, function(m) {
  bind_rows(
    describe_sample(make_baseline_cc(m) %>% distinct(rid, .keep_all = TRUE),
                    paste0(m, " (baseline model)")),
    describe_sample(make_long_cc(m) %>% distinct(rid, .keep_all = TRUE),
                    paste0(m, " (longitudinal model)"))
  )
}) %>%
  bind_rows(describe_sample(
    dat %>% filter(visit_key == "bl", !is.na(dx_label)) %>% distinct(rid, .keep_all = TRUE),
    "full covariate-complete cohort"
  )) %>%
  select(sample, dx_group, everything())

print(as.data.frame(panel_demographics))
out$D1_panel_demographics <- panel_demographics


# ==========================================================
# E. INCLUDED vs EXCLUDED  (Reviewer 2, comment 5)
#    Draft of a new supplementary table.
# ==========================================================

cat("\n=== E. Included vs excluded ===\n")

cohort_bl <- dat %>%
  filter(visit_key == "bl", dx_label %in% c("CN", "MCI", "AD")) %>%
  distinct(rid, .keep_all = TRUE)

included_vs_excluded <- map_dfr(primary_markers, function(m) {
  inc_ids <- unique(make_baseline_cc(m)$rid)
  cohort_bl %>%
    mutate(status = if_else(rid %in% inc_ids, "included", "excluded")) %>%
    group_by(status) %>%
    summarise(
      n          = n(),
      age_mean   = round(mean(age, na.rm = TRUE), 1),
      age_sd     = round(sd(age, na.rm = TRUE), 1),
      female_pct = round(100 * mean(tolower(as.character(ptgender)) %in% c("female", "f", "2")), 1),
      educ_mean  = round(mean(pteducat, na.rm = TRUE), 1),
      apoe4_pct  = round(100 * mean(apoe4 == 1, na.rm = TRUE), 1),
      mmse_mean  = if ("mmse" %in% names(cohort_bl)) round(mean(parse_num(mmse), na.rm = TRUE), 1) else NA_real_,
      adas13_mean = if ("adas13" %in% names(cohort_bl)) round(mean(parse_num(adas13), na.rm = TRUE), 1) else NA_real_,
      cn_pct     = round(100 * mean(dx_label == "CN"), 1),
      mci_pct    = round(100 * mean(dx_label == "MCI"), 1),
      ad_pct     = round(100 * mean(dx_label == "AD"), 1),
      .groups = "drop"
    ) %>%
    mutate(marker = m)
}) %>%
  select(marker, status, everything())

# Formal tests for the included/excluded contrast
inclusion_tests <- map_dfr(primary_markers, function(m) {
  inc_ids <- unique(make_baseline_cc(m)$rid)
  d <- cohort_bl %>% mutate(included = rid %in% inc_ids)
  tibble(
    marker      = m,
    p_age       = t.test(age ~ included, data = d)$p.value,
    p_education = t.test(pteducat ~ included, data = d)$p.value,
    p_sex       = chisq.test(table(d$ptgender, d$included))$p.value,
    p_apoe4     = chisq.test(table(d$apoe4, d$included))$p.value,
    p_diagnosis = chisq.test(table(d$dx_label, d$included))$p.value
  )
})

print(as.data.frame(included_vs_excluded))
out$E1_included_vs_excluded <- included_vs_excluded
out$E2_inclusion_tests      <- inclusion_tests


# ==========================================================
# F. AGE / TIME COLLINEARITY  (not raised by reviewers)
#
# age was built as baseline_age + years_from_baseline, so within
# any participant the two variables are identical. The reported
# "time" coefficient is therefore the rate of change net of the
# cross-sectional age effect, not the observed rate of change.
# ==========================================================

cat("\n=== F. Age / time collinearity ===\n")

collinearity_check <- map_dfr(primary_markers, function(m) {
  d <- make_long_cc(m)
  within_cor <- d %>%
    group_by(rid) %>%
    filter(n() >= 2) %>%
    mutate(age_c = age - mean(age), t_c = years_from_baseline - mean(years_from_baseline)) %>%
    ungroup() %>%
    summarise(r = cor(age_c, t_c, use = "complete.obs")) %>%
    pull(r)
  tibble(
    marker = m,
    overall_cor_age_time = cor(d$age, d$years_from_baseline, use = "complete.obs"),
    within_person_cor_age_time = within_cor,
    note = "within-person r near 1.0 means age adds no information beyond time"
  )
})

print(as.data.frame(collinearity_check))
out$F1_collinearity <- collinearity_check

# Side-by-side: time coefficient with time-varying age (as published)
# vs baseline age held fixed (the defensible specification).
age_spec_comparison <- map_dfr(primary_markers, function(m) {
  d <- make_long_cc(m) %>%
    group_by(rid) %>%
    mutate(baseline_age = age - years_from_baseline) %>%
    ungroup() %>%
    mutate(y = log(.data[[m]]))

  f_old <- lmer(y ~ years_from_baseline * baseline_dx + age + ptgender +
                  pteducat + apoe4 + (1 | rid), data = d, REML = FALSE)
  f_new <- lmer(y ~ years_from_baseline * baseline_dx + baseline_age + ptgender +
                  pteducat + apoe4 + (1 | rid), data = d, REML = FALSE)

  bind_rows(
    broom.mixed::tidy(f_old, effects = "fixed") %>% mutate(spec = "time-varying age (as published)"),
    broom.mixed::tidy(f_new, effects = "fixed") %>% mutate(spec = "baseline age (corrected)")
  ) %>%
    filter(grepl("years_from_baseline", term)) %>%
    mutate(marker = m, percent_change = (exp(estimate) - 1) * 100)
})
out$F2_age_specification <- age_spec_comparison


# ==========================================================
# G. FDR FAMILY SENSITIVITY  (Reviewer 2, comment 4)
#
# The published q-values came from p.adjust() applied to the full
# stacked coefficient table, so the family included five intercepts
# and all covariate terms: 35 tests at baseline, 50 longitudinally.
# This block reports q under four family definitions.
# ==========================================================

cat("\n=== G. FDR family sensitivity ===\n")

baseline_all_terms <- map_dfr(primary_markers, function(m) {
  d <- make_baseline_cc(m)
  fit <- lm(log(d[[m]]) ~ dx_label + age + ptgender + pteducat + apoe4, data = d)
  broom::tidy(fit) %>% mutate(marker = m, model = "baseline")
})

longitudinal_all_terms <- map_dfr(primary_markers, function(m) {
  d <- make_long_cc(m)
  fit <- lmer(log(d[[m]]) ~ years_from_baseline * baseline_dx + age + ptgender +
                pteducat + apoe4 + (1 | rid), data = d, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed") %>% mutate(marker = m, model = "longitudinal")
})

reported_baseline <- c("dx_labelMCI", "dx_labelAD")
reported_long <- c("years_from_baseline", "baseline_dxMCI", "baseline_dxAD",
                   "years_from_baseline:baseline_dxMCI", "years_from_baseline:baseline_dxAD")

all_terms <- bind_rows(baseline_all_terms, longitudinal_all_terms) %>%
  mutate(is_reported = (model == "baseline"     & term %in% reported_baseline) |
                       (model == "longitudinal" & term %in% reported_long))

# Family 1: as published (every coefficient, within model type)
fam_published <- all_terms %>%
  group_by(model) %>%
  mutate(q_as_published = p.adjust(p.value, "BH")) %>%
  ungroup()

# Family 2: reported terms only, within model type
fam_reported <- fam_published %>%
  filter(is_reported) %>%
  group_by(model) %>%
  mutate(q_reported_terms = p.adjust(p.value, "BH")) %>%
  ungroup()

# Family 3: per biomarker, reported terms only
fam_biomarker <- fam_reported %>%
  group_by(marker, model) %>%
  mutate(q_per_biomarker = p.adjust(p.value, "BH")) %>%
  ungroup()

# Family 4: global across all reported terms in the paper
fdr_sensitivity <- fam_biomarker %>%
  mutate(q_global = p.adjust(p.value, "BH")) %>%
  select(model, marker, term, estimate, std.error, p.value,
         q_as_published, q_reported_terms, q_per_biomarker, q_global) %>%
  mutate(across(starts_with("q_"), ~ round(.x, 6))) %>%
  arrange(model, marker, term)

# Does any conclusion move across families?
fdr_flips <- fdr_sensitivity %>%
  rowwise() %>%
  mutate(n_families_significant = sum(c(q_as_published, q_reported_terms,
                                        q_per_biomarker, q_global) < 0.05)) %>%
  ungroup() %>%
  filter(n_families_significant > 0 & n_families_significant < 4)

cat("\nTerms whose significance depends on the correction family:",
    nrow(fdr_flips), "\n")
if (nrow(fdr_flips) > 0) print(as.data.frame(fdr_flips))

out$G1_fdr_sensitivity <- fdr_sensitivity
out$G2_fdr_flips       <- fdr_flips
out$G3_family_sizes <- all_terms %>%
  count(model, name = "tests_in_published_family") %>%
  left_join(all_terms %>% filter(is_reported) %>%
              count(model, name = "tests_in_reported_family"), by = "model")


# ==========================================================
# Save
# ==========================================================

openxlsx::write.xlsx(out, file.path(results_dir, "revision_diagnostics.xlsx"),
                     overwrite = TRUE)

cat("\nDone. Diagnostics written to:\n",
    file.path(results_dir, "revision_diagnostics.xlsx"), "\n")
