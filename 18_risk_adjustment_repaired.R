# ==========================================================
# 18_risk_adjustment_repaired.R
#
# Replaces Part 5 of script 17, which failed three ways:
#   - BMI unit conversion produced impossible values (mean 37,
#     75th percentile 52) because the ADNIMERGE2 unit fields are
#     not the numeric 1/2 codes the CSV release uses.
#   - MEDHIST organ-system fields came back empty because
#     ADNIMERGE2 names them differently.
#   - Requiring BMI-complete cases collapsed the vascular
#     longitudinal model from 566 participants to 90.
#
# Fixes here:
#   Part A. Print what VITALS and MEDHIST actually contain
#   Part B. Derive BMI from the values themselves, not unit codes
#   Part C. Locate creatinine, including a LONI CSV if present
#   Part D. Tiered adjustment, so sample size is never silently lost
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
            "broom", "broom.mixed")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse); library(janitor); library(openxlsx)
library(lme4); library(lmerTest); library(broom); library(broom.mixed)

project_dir <- "R:/ADNI_Project"
raw_dir     <- file.path(project_dir, "00_raw_data")
rda_dir     <- file.path(raw_dir, "ADNIMERGE2", "data")
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

out <- list()
parse_num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))
first_nonmissing <- function(x) { x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else x[1] }

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

# ---------- Analysis dataset ----------

marker_info <- tribble(
  ~marker,            ~label,           ~stored_scale,
  "vegf_plasma_qc",   "Plasma VEGF",    "log10",
  "sicam1_plasma_qc", "Plasma sICAM-1", "log10",
  "svcam1_plasma_qc", "Plasma sVCAM-1", "log10"
)
vascular_markers <- marker_info$marker

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

for (m in vascular_markers) dat[[paste0("ln_", m)]] <- log(10) * parse_num(dat[[m]])

reported_long <- c("years_from_baseline", "baseline_dxMCI", "baseline_dxAD",
                   "years_from_baseline:baseline_dxMCI", "years_from_baseline:baseline_dxAD")


# ==========================================================
# PART A. WHAT DO VITALS AND MEDHIST ACTUALLY CONTAIN?
# ==========================================================

cat("\n=== Part A. Source inspection ===\n")

vitals  <- load_adni("VITALS")
medhist <- load_adni("MEDHIST")

describe_cols <- function(tb, tag) {
  if (is.null(tb)) return(tibble(source = tag, column = NA, note = "table not found"))
  map_dfr(names(tb), function(cc) {
    v <- tb[[cc]]
    nv <- parse_num(v)
    tibble(source = tag, column = cc,
           class = class(v)[1],
           n_non_missing = sum(!is.na(v)),
           n_distinct = dplyr::n_distinct(v),
           example_values = paste(head(unique(as.character(v[!is.na(v)])), 5), collapse = " | "),
           numeric_median = if (sum(!is.na(nv)) > 0) median(nv, na.rm = TRUE) else NA_real_)
  })
}

out$A_vitals_columns  <- describe_cols(vitals,  "VITALS")
out$A_medhist_columns <- describe_cols(medhist, "MEDHIST")

# The unit fields specifically: what values do they hold?
if (!is.null(vitals)) {
  unit_cols <- grep("unit|unt", names(vitals), value = TRUE, ignore.case = TRUE)
  cat("VITALS unit-like columns:", paste(unit_cols, collapse = ", "), "\n")
  out$A_vitals_unit_values <- if (length(unit_cols)) {
    map_dfr(unit_cols, function(cc)
      tibble(column = cc, value = as.character(vitals[[cc]])) %>%
        count(column, value, name = "n"))
  } else tibble(note = "no unit columns found")
  print(as.data.frame(out$A_vitals_unit_values))
}

# MEDHIST binary history fields, whatever they are called
if (!is.null(medhist)) {
  mh_bin <- names(medhist)[map_lgl(medhist, function(v) {
    u <- unique(na.omit(as.character(v)))
    length(u) > 1 && length(u) <= 4 && all(u %in% c("0", "1", "2", "-4", "Yes", "No"))
  })]
  cat("MEDHIST binary-looking fields:", paste(head(mh_bin, 40), collapse = ", "), "\n")
  out$A_medhist_binary_fields <- tibble(field = mh_bin)
}


# ==========================================================
# PART B. BMI FROM THE VALUES, NOT THE UNIT CODES
#
# Unit codes are unreliable here, but the magnitudes are not.
# Adult weight above 130 is pounds; height above 120 is
# centimetres and height between 45 and 90 is inches.
# ==========================================================

cat("\n=== Part B. BMI ===\n")

vitals_bl <- NULL
if (!is.null(vitals)) {
  kv <- intersect(c("viscode2", "viscode"), names(vitals))[1]
  g  <- function(n) if (n %in% names(vitals)) parse_num(vitals[[n]]) else NA_real_

  v0 <- vitals %>%
    mutate(rid = as.character(rid),
           visit_key = if (!is.na(kv)) as.character(.data[[kv]]) else NA_character_,
           wt_raw = g("vsweight"), ht_raw = g("vsheight"),
           sbp = g("vsbpsys"), dbp = g("vsbpdia"))

  cat("Raw weight  median:", median(v0$wt_raw, na.rm = TRUE),
      " range:", paste(round(range(v0$wt_raw, na.rm = TRUE), 1), collapse = "-"), "\n")
  cat("Raw height  median:", median(v0$ht_raw, na.rm = TRUE),
      " range:", paste(round(range(v0$ht_raw, na.rm = TRUE), 1), collapse = "-"), "\n")

  v0 <- v0 %>%
    mutate(
      weight_kg = case_when(
        is.na(wt_raw)               ~ NA_real_,
        wt_raw >= 130 & wt_raw < 500 ~ wt_raw * 0.45359237,   # pounds
        wt_raw >= 30  & wt_raw < 130 ~ wt_raw,                # kilograms
        TRUE                        ~ NA_real_),
      height_m = case_when(
        is.na(ht_raw)              ~ NA_real_,
        ht_raw >= 120 & ht_raw < 220 ~ ht_raw / 100,          # centimetres
        ht_raw >= 45  & ht_raw < 90  ~ ht_raw * 0.0254,       # inches
        ht_raw >= 1.2 & ht_raw < 2.2 ~ ht_raw,                # already metres
        TRUE                       ~ NA_real_),
      bmi = weight_kg / height_m^2,
      bmi = if_else(bmi >= 14 & bmi <= 60, bmi, NA_real_),
      sbp = if_else(sbp > 60 & sbp < 260, sbp, NA_real_),
      dbp = if_else(dbp > 30 & dbp < 160, dbp, NA_real_))

  cat("Derived BMI  n:", sum(!is.na(v0$bmi)),
      " median:", round(median(v0$bmi, na.rm = TRUE), 1),
      " IQR:", paste(round(quantile(v0$bmi, c(.25, .75), na.rm = TRUE), 1), collapse = "-"), "\n")
  cat("  (expect a median near 26-27; if not, inspect A_vitals_columns)\n")

  out$B_bmi_check <- tibble(
    n_weight = sum(!is.na(v0$wt_raw)), n_height = sum(!is.na(v0$ht_raw)),
    n_bmi = sum(!is.na(v0$bmi)),
    bmi_median = median(v0$bmi, na.rm = TRUE),
    bmi_q1 = quantile(v0$bmi, .25, na.rm = TRUE),
    bmi_q3 = quantile(v0$bmi, .75, na.rm = TRUE),
    pct_weight_pounds = round(100 * mean(v0$wt_raw >= 130, na.rm = TRUE), 1),
    pct_height_cm     = round(100 * mean(v0$ht_raw >= 120, na.rm = TRUE), 1))

  # Baseline-ish value per participant, falling back to the earliest available
  vitals_bl <- v0 %>%
    mutate(is_bl = visit_key %in% c("bl", "sc", "scmri")) %>%
    arrange(rid, desc(is_bl)) %>%
    group_by(rid) %>%
    summarise(bmi = first_nonmissing(bmi), sbp = first_nonmissing(sbp),
              dbp = first_nonmissing(dbp), .groups = "drop") %>%
    mutate(across(c(bmi, sbp, dbp), as.numeric))
  cat("Participants with BMI:", sum(!is.na(vitals_bl$bmi)),
      " SBP:", sum(!is.na(vitals_bl$sbp)), "\n")
}


# ==========================================================
# PART C. MEDICATIONS, MEDICAL HISTORY, CREATININE
# ==========================================================

cat("\n=== Part C. Medications, history, creatinine ===\n")

meds <- load_adni("RECCMEDS")
meds_bl <- NULL
if (!is.null(meds)) {
  mc <- intersect(c("cmmed", "cmmeds", "medname", "cmtrt"), names(meds))[1]
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
  cat("Antihypertensive:", sum(meds_bl$med_antihypertensive),
      " antidiabetic:", sum(meds_bl$med_antidiabetic),
      " (ever prescribed, across all visits)\n")
}

# Medical history: match by pattern rather than a fixed field name
medhist_bl <- NULL
if (!is.null(medhist)) {
  find_mh <- function(pat) {
    cc <- grep(pat, names(medhist), value = TRUE, ignore.case = TRUE)
    if (!length(cc)) return(NA_character_)
    cc[1]
  }
  fields <- c(mh_cardiovascular = find_mh("card"),
              mh_endocrine      = find_mh("endo|metab"),
              mh_renal          = find_mh("rena|renl|genito"),
              mh_smoking        = find_mh("smok"))
  cat("MEDHIST fields matched:\n"); print(fields)
  out$C_medhist_fields_used <- tibble(concept = names(fields), column = unname(fields))

  keep <- fields[!is.na(fields)]
  if (length(keep)) {
    medhist_bl <- medhist %>%
      mutate(rid = as.character(rid)) %>%
      select(rid, all_of(unname(keep))) %>%
      rename(!!!setNames(unname(keep), names(keep))) %>%
      group_by(rid) %>%
      summarise(across(everything(), ~ first_nonmissing(.x)), .groups = "drop") %>%
      mutate(across(-rid, ~ as.integer(as.character(.x) %in% c("1", "Yes", "yes"))))
    cat("Positive history counts:\n")
    print(colSums(medhist_bl[-1], na.rm = TRUE))
  }
}

# Creatinine: prefer the file exported by 19_extract_creatinine_from_adnimerge.R,
# which resolves ADNI's coded lab column names through the data dictionary.
cat("\nLooking for creatinine...\n")
creat_bl <- NULL

creat_file <- file.path(clean_dir, "creatinine_baseline.csv")
if (file.exists(creat_file)) {
  creat_bl <- read_csv(creat_file, show_col_types = FALSE) %>%
    mutate(rid = as.character(rid), creatinine = as.numeric(creatinine)) %>%
    select(rid, creatinine)
  cat("Using creatinine_baseline.csv from script 19 - n =", nrow(creat_bl),
      "| median", round(median(creat_bl$creatinine), 2), "mg/dL\n")
}

lab_candidates <- c("LABDATA", "labdata", "ADNI_LAB", "LABORATORY", "BLOODCHEM", "CHEM")
if (is.null(creat_bl)) for (stem in lab_candidates) {
  labs <- load_adni(stem)
  if (is.null(labs)) next

  cc <- grep("creat", names(labs), value = TRUE, ignore.case = TRUE)
  if (length(cc)) {
    v <- parse_num(labs[[cc[1]]])
    if (sum(!is.na(v) & v > 0) > 100 && dplyr::n_distinct(v) > 20) {
      kv <- intersect(c("viscode2", "viscode"), names(labs))[1]
      creat_bl <- labs %>%
        mutate(rid = as.character(rid), creatinine = parse_num(.data[[cc[1]]])) %>%
        { if (!is.na(kv)) filter(., as.character(.data[[kv]]) %in% c("bl", "sc", "scmri")) else . } %>%
        filter(!is.na(creatinine), creatinine > 0.2, creatinine < 15) %>%
        group_by(rid) %>% summarise(creatinine = first(creatinine), .groups = "drop")
      cat("Creatinine found in", stem, "column", cc[1], "- n =", nrow(creat_bl), "\n")
      break
    }
  }
}
if (is.null(creat_bl)) {
  cat("\nNO CREATININE AVAILABLE.\n",
      "Run 19_extract_creatinine_from_adnimerge.R first. ADNI's lab table uses\n",
      "coded column names, so a plain name search will not find creatinine.\n",
      "Tier 4 will be reported as not fitted.\n")
  out$C_creatinine_status <- tibble(
    status = "not available", action = "run script 19 first")
} else {
  out$C_creatinine_status <- tibble(status = "available", n = nrow(creat_bl))
}


# ==========================================================
# PART D. TIERED ADJUSTMENT
#
# Adjusting on a covariate that is missing for most participants
# silently destroys the sample. Each tier is therefore fitted on
# its own complete-case set, and every tier reports its n, so the
# trade-off between adjustment and power is visible.
# ==========================================================

cat("\n=== Part D. Tiered adjustment ===\n")

# One row per participant. baseline_age is computed per row as
# age - years_from_baseline, and because age is integer-valued it
# varies slightly between visits. distinct() on it would return
# several rows per participant and fan out the join below.
risk <- dat %>%
  arrange(rid, years_from_baseline) %>%
  group_by(rid) %>%
  summarise(baseline_age = first(baseline_age[!is.na(baseline_age)]),
            ptgender     = first(as.character(ptgender)[!is.na(ptgender)]),
            .groups = "drop")

stopifnot(!any(duplicated(risk$rid)))

for (tb in list(vitals_bl, meds_bl, medhist_bl, creat_bl)) {
  if (is.null(tb)) next
  if (any(duplicated(tb$rid))) tb <- tb %>% distinct(rid, .keep_all = TRUE)
  risk <- left_join(risk, tb, by = "rid")
}
stopifnot(!any(duplicated(risk$rid)))

if ("creatinine" %in% names(risk)) {
  risk <- risk %>%
    mutate(fem = tolower(as.character(ptgender)) %in% c("female", "f", "2"),
           k = if_else(fem, 0.7, 0.9), a = if_else(fem, -0.241, -0.302),
           egfr = 142 * pmin(creatinine / k, 1)^a * pmax(creatinine / k, 1)^(-1.200) *
                  0.9938^baseline_age * if_else(fem, 1.012, 1)) %>%
    select(-fem, -k, -a)
}

has <- function(v) v %in% names(risk) && sum(!is.na(risk[[v]])) > 200 &&
  dplyr::n_distinct(risk[[v]], na.rm = TRUE) > 1

tiers <- list(
  "1: as published"                = character(0),
  "2: + blood pressure and meds"   = c("sbp", "med_antihypertensive", "med_antidiabetic"),
  "3: + BMI and smoking"           = c("sbp", "med_antihypertensive", "med_antidiabetic",
                                       "bmi", "mh_smoking"),
  "4: + renal function"            = c("sbp", "med_antihypertensive", "med_antidiabetic",
                                       "bmi", "mh_smoking", "egfr")
)
tiers <- map(tiers, ~ .x[map_lgl(.x, has)])
tiers <- tiers[!duplicated(map_chr(tiers, ~ paste(sort(.x), collapse = "+")))]
cat("Tiers to be fitted:\n"); print(tiers)
out$D_tiers <- tibble(tier = names(tiers),
                      covariates = map_chr(tiers, ~ if (length(.x)) paste(.x, collapse = ", ") else "none"))

n_before <- nrow(dat)
dat_risk <- dat %>% left_join(risk %>% select(-baseline_age, -ptgender), by = "rid")
if (nrow(dat_risk) != n_before)
  stop("Covariate join duplicated rows: ", n_before, " -> ", nrow(dat_risk),
       ". A covariate table has more than one row per participant.")
cat("Covariate join preserved row count:", nrow(dat_risk), "rows\n")

out$D_baseline_tiers <- map_dfr(names(tiers), function(tn) {
  ex <- tiers[[tn]]
  map_dfr(vascular_markers, function(m) {
    d <- dat_risk %>%
      filter(visit_key == "bl", !is.na(.data[[paste0("ln_", m)]]), !is.na(dx_label),
             !is.na(age), !is.na(ptgender), !is.na(pteducat), !is.na(apoe4))
    if (length(ex)) d <- d %>% drop_na(all_of(ex))
    if (nrow(d) < 60) return(tibble(label = marker_info$label[marker_info$marker == m],
                                    tier = tn, term = "NOT FITTED",
                                    note = paste("only", nrow(d), "participants")))
    rhs <- paste("dx_label + age + ptgender + pteducat + apoe4",
                 if (length(ex)) paste("+", paste(ex, collapse = " + ")) else "")
    fit <- lm(as.formula(paste0("ln_", m, " ~ ", rhs)), data = d)
    tidy(fit) %>% filter(term %in% c("dx_labelMCI", "dx_labelAD")) %>%
      mutate(label = marker_info$label[marker_info$marker == m], tier = tn,
             percent_change = (exp(estimate) - 1) * 100,
             pct_ci_low  = (exp(estimate - 1.96 * std.error) - 1) * 100,
             pct_ci_high = (exp(estimate + 1.96 * std.error) - 1) * 100,
             n = nrow(d))
  })
}) %>% relocate(label, tier, term)
print(as.data.frame(out$D_baseline_tiers %>% select(any_of(
  c("label", "tier", "term", "percent_change", "p.value", "n", "note")))))

out$D_longitudinal_tiers <- map_dfr(names(tiers), function(tn) {
  ex <- tiers[[tn]]
  map_dfr(vascular_markers, function(m) {
    d <- dat_risk %>%
      filter(!is.na(.data[[paste0("ln_", m)]]), !is.na(years_from_baseline),
             !is.na(baseline_dx), !is.na(baseline_age), !is.na(ptgender),
             !is.na(pteducat), !is.na(apoe4))
    if (length(ex)) d <- d %>% drop_na(all_of(ex))
    if (dplyr::n_distinct(d$rid) < 150)
      return(tibble(label = marker_info$label[marker_info$marker == m], tier = tn,
                    term = "NOT FITTED",
                    note = paste("only", dplyr::n_distinct(d$rid), "participants;",
                                 "too few to interpret")))
    rhs <- paste("years_from_baseline * baseline_dx + baseline_age + ptgender + pteducat + apoe4",
                 if (length(ex)) paste("+", paste(ex, collapse = " + ")) else "")
    fit <- lmer(as.formula(paste0("ln_", m, " ~ ", rhs, " + (1 | rid)")),
                data = d, REML = FALSE)
    tidy(fit, effects = "fixed") %>% filter(term %in% reported_long) %>%
      mutate(label = marker_info$label[marker_info$marker == m], tier = tn,
             percent_change = (exp(estimate) - 1) * 100,
             n_obs = nobs(fit), n_subj = ngrps(fit)[["rid"]])
  })
}) %>% relocate(label, tier, term)
print(as.data.frame(out$D_longitudinal_tiers %>% select(any_of(
  c("label", "tier", "term", "percent_change", "p.value", "n_subj", "note")))))

# Attenuation of the AD contrast across tiers
out$D_attenuation <- out$D_baseline_tiers %>%
  filter(term == "dx_labelAD") %>%
  select(label, tier, estimate, n) %>%
  group_by(label) %>%
  mutate(percent_attenuation_vs_tier1 =
           100 * (1 - estimate / estimate[tier == "1: as published"])) %>%
  ungroup()

# Are the risk factors distributed differently by diagnosis?
risk_vars <- unique(unlist(tiers))
if (length(risk_vars)) {
  out$D_risk_by_diagnosis <- dat_risk %>%
    filter(visit_key == "bl", !is.na(dx_label)) %>% distinct(rid, .keep_all = TRUE) %>%
    group_by(dx_label) %>%
    summarise(n = n(), across(all_of(risk_vars),
                              list(mean = ~ round(mean(.x, na.rm = TRUE), 2),
                                   n    = ~ sum(!is.na(.x))),
                              .names = "{.col}_{.fn}"), .groups = "drop")
  print(as.data.frame(out$D_risk_by_diagnosis))
}

out$D_risk_table <- risk


# ==========================================================
out <- imap(out, function(x, nm) {
  if (is.null(x) || (is.data.frame(x) && (nrow(x) == 0 || ncol(x) == 0)))
    tibble(note = paste("No rows produced for", nm)) else x
})
names(out) <- substr(make.unique(names(out)), 1, 31)

openxlsx::write.xlsx(out, file.path(results_dir, "risk_adjustment_repaired.xlsx"),
                     overwrite = TRUE)
cat("\nDone:", file.path(results_dir, "risk_adjustment_repaired.xlsx"), "\n")
