# ==========================================================
# 15_sensitivity_analyses.R
#
# Sensitivity analyses requested by Reviewer 2 (comments 2, 3),
# plus a batch-effect investigation of the vascular longitudinal
# finding that neither reviewer raised but that determines how
# that finding should be framed.
#
# Sections
#   A. Vascular batch / assay-run investigation
#   B. Restriction to participants with 2+ measurements    (R2-2)
#   C. Restriction to 2+ years of follow-up                (R2-2)
#   D. MCI converters and time-varying diagnosis           (R2-3)
#   E. Informative dropout probe
#
# Run 14_corrected_primary_analysis.R first; this script repeats
# the corrected outcome construction so it can stand alone.
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
            "broom", "broom.mixed")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse); library(janitor); library(openxlsx)
library(lme4); library(lmerTest); library(broom); library(broom.mixed)

project_dir <- "R:/ADNI_Project"
bio_dir     <- file.path(project_dir, "00_raw_data", "biomarkers_excel")
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

out <- list()
parse_num <- function(x) readr::parse_number(as.character(x))
first_nonmissing <- function(x) { x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else x[1] }

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

long_cc <- function(m, d = dat) {
  d %>% filter(!is.na(.data[[paste0("ln_", m)]]), !is.na(years_from_baseline),
               !is.na(baseline_dx), !is.na(baseline_age), !is.na(ptgender),
               !is.na(pteducat), !is.na(apoe4))
}

reported_long <- c("years_from_baseline", "baseline_dxMCI", "baseline_dxAD",
                   "years_from_baseline:baseline_dxMCI", "years_from_baseline:baseline_dxAD")

fit_long <- function(m, d) {
  lmer(as.formula(paste0("ln_", m,
    " ~ years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat + apoe4 + (1 | rid)")),
    data = d, REML = FALSE)
}

extract_long <- function(fit, m, tag) {
  tidy(fit, effects = "fixed") %>%
    filter(term %in% reported_long) %>%
    mutate(marker = m, analysis = tag,
           percent_change = (exp(estimate) - 1) * 100,
           n_obs = nobs(fit), n_subj = ngrps(fit)[["rid"]])
}


# ==========================================================
# A. VASCULAR BATCH / ASSAY-RUN INVESTIGATION
#
# The vascular panel has exactly two visits, bl and m12, assayed
# by Myriad RBM as part of the Biomarkers Consortium project.
# The published sICAM-1 CN slope corresponds to roughly a 24%
# fall in concentration over 12 months, which is not a plausible
# biological rate. Three checks distinguish biology from a
# run-to-run shift.
# ==========================================================

cat("\n=== A. Vascular batch investigation ===\n")

# A1. Paired within-person change from bl to m12, by diagnosis.
#     A uniform shift across all three groups points to a batch
#     effect; a shift confined to one group points to biology.

paired_change <- map_dfr(vascular_markers, function(m) {
  ln <- paste0("ln_", m)
  dat %>%
    filter(visit_key %in% c("bl", "m12"), !is.na(.data[[ln]]), !is.na(baseline_dx)) %>%
    select(rid, baseline_dx, visit_key, value = all_of(ln)) %>%
    pivot_wider(names_from = visit_key, values_from = value) %>%
    filter(!is.na(bl), !is.na(m12)) %>%
    mutate(delta = m12 - bl) %>%
    group_by(baseline_dx) %>%
    summarise(
      n = n(),
      mean_delta_log = mean(delta),
      percent_change_12m = (exp(mean(delta)) - 1) * 100,
      sd_delta = sd(delta),
      t_p = t.test(delta)$p.value,
      pct_declining = round(100 * mean(delta < 0), 1),
      .groups = "drop"
    ) %>%
    mutate(marker = m)
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

print(as.data.frame(paired_change))
out$A1_paired_change <- paired_change

# A2. Is the bl-to-m12 shift the same in every group?
#     If the group differences vanish, the shift is common to all
#     samples and behaves like a run effect.
out$A2_shift_homogeneity <- map_dfr(vascular_markers, function(m) {
  ln <- paste0("ln_", m)
  d <- dat %>%
    filter(visit_key %in% c("bl", "m12"), !is.na(.data[[ln]]), !is.na(baseline_dx)) %>%
    select(rid, baseline_dx, visit_key, value = all_of(ln)) %>%
    pivot_wider(names_from = visit_key, values_from = value) %>%
    filter(!is.na(bl), !is.na(m12)) %>%
    mutate(delta = m12 - bl)
  a <- aov(delta ~ baseline_dx, data = d)
  tibble(marker = m,
         overall_mean_shift_pct = (exp(mean(d$delta)) - 1) * 100,
         anova_p_group_difference = summary(a)[[1]][["Pr(>F)"]][1],
         interpretation = if (summary(a)[[1]][["Pr(>F)"]][1] > 0.05)
           "shift is uniform across groups -> consistent with a run effect"
         else "shift differs by group -> not explained by a uniform run effect")
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

# A3. Independent source check. ADMC_CLINICALVARIABLES is a separate
#     processing of the same three analytes. If the 12-month drop
#     appears there too it is in the samples, not the QC pipeline.

admc_path <- file.path(bio_dir, "ADMC_CLINICALVARIABLES_16May2016.csv")
if (file.exists(admc_path)) {
  admc <- read_csv(admc_path, show_col_types = FALSE, guess_max = 100000) %>% clean_names()
  key_visit <- intersect(c("viscode2", "viscode"), names(admc))[1]

  admc_long <- admc %>%
    transmute(rid = as.character(rid), visit_key = as.character(.data[[key_visit]]),
              vegf_plasma_qc   = parse_num(vegf),
              sicam1_plasma_qc = parse_num(icam),
              svcam1_plasma_qc = parse_num(vcam)) %>%
    group_by(rid, visit_key) %>%
    summarise(across(everything(), first_nonmissing), .groups = "drop") %>%
    mutate(across(all_of(vascular_markers), as.numeric)) %>%
    pivot_longer(all_of(vascular_markers), names_to = "marker", values_to = "admc_value")

  out$A3_admc_check <- admc_long %>%
    filter(visit_key %in% c("bl", "m12"), !is.na(admc_value), admc_value > 0) %>%
    left_join(dat %>% distinct(rid, baseline_dx), by = "rid") %>%
    mutate(ln_admc = log(admc_value)) %>%
    select(rid, marker, visit_key, ln_admc, baseline_dx) %>%
    pivot_wider(names_from = visit_key, values_from = ln_admc) %>%
    filter(!is.na(bl), !is.na(m12)) %>%
    group_by(marker) %>%
    summarise(n = n(), mean_delta_log = mean(m12 - bl),
              percent_change_12m = (exp(mean(m12 - bl)) - 1) * 100,
              t_p = t.test(m12 - bl)$p.value, .groups = "drop") %>%
    mutate(note = "compare percent_change_12m with sheet A1; agreement means the shift is in the samples")
  print(as.data.frame(out$A3_admc_check))
} else {
  cat("ADMC file not found, skipping the independent source check.\n")
}

# A4. Does the raw multiplex release carry plate / run / batch fields?
raw_path <- file.path(bio_dir, "adni_plasma_raw_multiplex_11Nov2010.csv")
if (file.exists(raw_path)) {
  raw_mx <- read_csv(raw_path, show_col_types = FALSE, guess_max = 100000) %>% clean_names()
  batch_cols <- grep("plate|batch|run|assay|date|lot|kit", names(raw_mx),
                     value = TRUE, ignore.case = TRUE)
  out$A4_batch_fields <- tibble(
    candidate_batch_columns = if (length(batch_cols)) paste(batch_cols, collapse = ", ")
                              else "none found",
    all_columns = paste(names(raw_mx), collapse = " | ")
  )
  cat("Candidate batch fields in raw multiplex file:",
      if (length(batch_cols)) paste(batch_cols, collapse = ", ") else "none", "\n")
}


# ==========================================================
# B. RESTRICTION TO PARTICIPANTS WITH 2+ MEASUREMENTS
#
# Participants with a single observation contribute nothing to
# within-person slope estimation but do pull the pooled time
# coefficient toward a between-person contrast.
# ==========================================================

cat("\n=== B. Two-or-more measurements ===\n")

sens_2plus <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m) %>% group_by(rid) %>% filter(n() >= 2) %>% ungroup()
  if (n_distinct(d$rid) < 30 || n_distinct(d$baseline_dx) < 2) return(NULL)
  extract_long(fit_long(m, d), m, "restricted to 2+ measurements")
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$B_two_or_more_visits <- sens_2plus

# Within-between (Mundlak) decomposition, which uses the whole sample
# but separates within-person change from between-person differences.
mundlak <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m) %>%
    group_by(rid) %>%
    mutate(time_mean = mean(years_from_baseline),
           time_within = years_from_baseline - time_mean) %>%
    ungroup()
  if (sd(d$time_within, na.rm = TRUE) == 0) return(NULL)
  fit <- lmer(as.formula(paste0("ln_", m,
    " ~ time_within * baseline_dx + time_mean + baseline_age + ptgender +",
    " pteducat + apoe4 + (1 | rid)")), data = d, REML = FALSE)
  tidy(fit, effects = "fixed") %>%
    filter(grepl("time_within|time_mean", term)) %>%
    mutate(marker = m, percent_change = (exp(estimate) - 1) * 100,
           reading = if_else(grepl("time_mean", term),
                             "between-person component",
                             "within-person component"))
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$B_within_between <- mundlak


# ==========================================================
# C. RESTRICTION TO 2+ YEARS OF FOLLOW-UP
#    Not possible for the vascular panel, whose maximum follow-up
#    is about 1.6 years. That fact is itself part of the answer.
# ==========================================================

cat("\n=== C. Two-or-more years of follow-up ===\n")

followup_capacity <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m) %>% group_by(rid) %>%
    summarise(span = max(years_from_baseline) - min(years_from_baseline), .groups = "drop")
  tibble(marker = m, n_participants = nrow(d),
         n_with_2y_span = sum(d$span >= 2),
         pct_with_2y_span = round(100 * mean(d$span >= 2), 1),
         max_span = round(max(d$span), 2))
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$C_followup_capacity <- followup_capacity

sens_2y <- map_dfr(primary_markers, function(m) {
  ids <- long_cc(m) %>% group_by(rid) %>%
    summarise(span = max(years_from_baseline) - min(years_from_baseline), .groups = "drop") %>%
    filter(span >= 2) %>% pull(rid)
  d <- long_cc(m) %>% filter(rid %in% ids)
  if (n_distinct(d$rid) < 30 || n_distinct(d$baseline_dx) < 2) return(NULL)
  extract_long(fit_long(m, d), m, "restricted to 2+ years of follow-up")
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$C_two_year_followup <- sens_2y


# ==========================================================
# D. MCI CONVERTERS AND TIME-VARYING DIAGNOSIS  (Reviewer 2, comment 3)
#
# dx_label is already visit-level, so conversion status and a
# time-varying diagnosis model are both available without new data.
# ==========================================================

cat("\n=== D. Converters ===\n")

conversion <- dat %>%
  filter(!is.na(dx_label), !is.na(baseline_dx), !is.na(years_from_baseline)) %>%
  arrange(rid, years_from_baseline) %>%
  group_by(rid, baseline_dx) %>%
  summarise(
    n_visits_with_dx = n(),
    ever_ad  = any(dx_label == "AD"),
    ever_mci = any(dx_label == "MCI"),
    last_dx  = last(dx_label),
    time_to_first_ad = { i <- which(dx_label == "AD")
                         if (length(i)) years_from_baseline[i[1]] else NA_real_ },
    followup_span = max(years_from_baseline) - min(years_from_baseline),
    .groups = "drop"
  ) %>%
  mutate(converter = case_when(
    baseline_dx == "MCI" & ever_ad ~ "MCI to AD converter",
    baseline_dx == "MCI"           ~ "stable MCI",
    baseline_dx == "CN" & (ever_mci | ever_ad) ~ "CN progressor",
    baseline_dx == "CN"            ~ "stable CN",
    TRUE ~ "AD at baseline"
  ))

conversion_summary <- conversion %>%
  count(baseline_dx, converter) %>%
  group_by(baseline_dx) %>% mutate(percent = round(100 * n / sum(n), 1)) %>% ungroup()
print(as.data.frame(conversion_summary))
out$D1_conversion_summary <- conversion_summary

# Conversion rates within each biomarker analytic sample
out$D2_conversion_by_panel <- map_dfr(primary_markers, function(m) {
  ids <- unique(long_cc(m)$rid)
  conversion %>% filter(rid %in% ids) %>% count(marker = m, baseline_dx, converter) %>%
    group_by(marker, baseline_dx) %>% mutate(percent = round(100 * n / sum(n), 1)) %>% ungroup()
})

# D3. Time-varying diagnosis model
out$D3_time_varying_dx <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m) %>% filter(!is.na(dx_label))
  if (n_distinct(d$rid) < 30) return(NULL)
  fit <- lmer(as.formula(paste0("ln_", m,
    " ~ years_from_baseline * dx_label + baseline_age + ptgender + pteducat + apoe4 + (1 | rid)")),
    data = d, REML = FALSE)
  tidy(fit, effects = "fixed") %>%
    filter(grepl("years_from_baseline|dx_label", term)) %>%
    mutate(marker = m, analysis = "time-varying diagnosis",
           percent_change = (exp(estimate) - 1) * 100,
           n_obs = nobs(fit), n_subj = ngrps(fit)[["rid"]])
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)

# D4. Split the MCI arm into stable and converter
out$D4_mci_split <- map_dfr(primary_markers, function(m) {
  d <- long_cc(m) %>%
    left_join(conversion %>% select(rid, converter), by = "rid") %>%
    mutate(dx_group = factor(case_when(
      baseline_dx == "CN" ~ "CN",
      converter == "stable MCI" ~ "stable MCI",
      converter == "MCI to AD converter" ~ "MCI converter",
      baseline_dx == "AD" ~ "AD",
      TRUE ~ NA_character_),
      levels = c("CN", "stable MCI", "MCI converter", "AD"))) %>%
    filter(!is.na(dx_group))
  if (n_distinct(d$dx_group) < 3 || n_distinct(d$rid) < 40) return(NULL)
  fit <- lmer(as.formula(paste0("ln_", m,
    " ~ years_from_baseline * dx_group + baseline_age + ptgender + pteducat + apoe4 + (1 | rid)")),
    data = d, REML = FALSE)
  tidy(fit, effects = "fixed") %>%
    filter(grepl("years_from_baseline|dx_group", term)) %>%
    mutate(marker = m, analysis = "MCI split by conversion",
           percent_change = (exp(estimate) - 1) * 100)
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)


# ==========================================================
# E. INFORMATIVE DROPOUT PROBE
#
# Does the baseline biomarker value predict whether a participant
# returns for a follow-up measurement? If it does, dropout is
# informative and the slopes need that caveat.
# ==========================================================

cat("\n=== E. Informative dropout ===\n")

dropout <- map_dfr(primary_markers, function(m) {
  ln <- paste0("ln_", m)
  d <- long_cc(m) %>% group_by(rid) %>%
    mutate(returned = as.integer(n() >= 2)) %>%
    filter(years_from_baseline == min(years_from_baseline)) %>%
    slice(1) %>% ungroup()
  if (n_distinct(d$returned) < 2) return(NULL)
  fit <- glm(returned ~ get(ln) + baseline_dx + baseline_age + ptgender + pteducat + apoe4,
             data = d, family = binomial)
  tidy(fit) %>% filter(grepl("get\\(ln\\)|baseline_dx", term)) %>%
    mutate(marker = m, odds_ratio = exp(estimate),
           term = dplyr::recode(term, `get(ln)` = "baseline biomarker level (log)"))
}) %>% left_join(marker_info, by = "marker") %>% relocate(label)
out$E_informative_dropout <- dropout


out <- imap(out, function(x, nm) {
  if (is.null(x) || (is.data.frame(x) && (nrow(x) == 0 || ncol(x) == 0)))
    tibble(note = paste("No rows produced for", nm)) else x
})

# Excel caps sheet names at 31 characters.
names(out) <- substr(make.unique(names(out)), 1, 31)

openxlsx::write.xlsx(out, file.path(results_dir, "sensitivity_analyses.xlsx"),
                     overwrite = TRUE)
cat("\nDone. Results:", file.path(results_dir, "sensitivity_analyses.xlsx"), "\n")
