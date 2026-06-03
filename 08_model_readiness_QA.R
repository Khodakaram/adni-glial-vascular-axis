# ==========================================================
# 08_model_readiness_QA.R
# Purpose: Diagnose model readiness and fix baseline_dx / age problems
# ==========================================================

library(tidyverse)
library(janitor)
library(openxlsx)
library(lubridate)
library(readr)

project_dir <- "R:/ADNI_Project"

clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------

safe_log <- function(x) {
  x <- parse_number(as.character(x))
  ifelse(!is.na(x) & x > 0, log(x), NA_real_)
}

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  
  x_chr <- as.character(x)
  x_chr[x_chr == "" | x_chr == "NA"] <- NA_character_
  
  out <- suppressWarnings(ymd(x_chr))
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(mdy(x_chr[idx]))
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(dmy(x_chr[idx]))
  
  return(out)
}

# ---------- 1. Load latest analysis file ----------

analysis_file_1 <- file.path(clean_dir, "analysis_master_with_log_variables.csv")
analysis_file_2 <- file.path(clean_dir, "analysis_master_biomarkers_longitudinal.csv")

if (file.exists(analysis_file_1)) {
  dat <- read_csv(analysis_file_1, show_col_types = FALSE) %>% clean_names()
} else {
  dat <- read_csv(analysis_file_2, show_col_types = FALSE) %>% clean_names()
}

# ---------- 2. Standardize important variables ----------

dat <- dat %>%
  select(-any_of(c("baseline_dx", "baseline_dx.x", "baseline_dx.y"))) %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key),
    dx_label = as.character(dx_label),
    examdate = parse_date_safe(examdate),
    ptdobyy = parse_number(as.character(ptdobyy)),
    pteducat = parse_number(as.character(pteducat)),
    apoe4 = parse_number(as.character(apoe4)),
    years_from_baseline = parse_number(as.character(years_from_baseline))
  )

# ---------- 3. Recalculate age ----------

dat <- dat %>%
  mutate(
    age_recalculated = if_else(
      !is.na(examdate) & !is.na(ptdobyy),
      as.numeric(year(examdate) - ptdobyy),
      NA_real_
    )
  )

# If age exists but is empty, replace it
if ("age" %in% names(dat)) {
  dat <- dat %>%
    mutate(
      age = parse_number(as.character(age)),
      age = if_else(is.na(age), age_recalculated, age)
    )
} else {
  dat <- dat %>%
    mutate(age = age_recalculated)
}

# ---------- 4. Create clean baseline diagnosis table ----------

baseline_dx_tbl <- dat %>%
  filter(dx_label %in% c("CN", "MCI", "AD")) %>%
  arrange(
    rid,
    case_when(
      visit_key == "bl" ~ 1,
      visit_key == "sc" ~ 2,
      TRUE ~ 3
    ),
    years_from_baseline
  ) %>%
  group_by(rid) %>%
  summarise(
    baseline_dx = first(dx_label),
    .groups = "drop"
  )

dat <- dat %>%
  left_join(baseline_dx_tbl, by = "rid") %>%
  mutate(
    dx_label = factor(dx_label, levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender)
  )

# ---------- 5. Define primary markers ----------

primary_markers <- c(
  "gfap_quanterix",
  "strem2_msd_corrected",
  "vegf_plasma_qc",
  "sicam1_plasma_qc",
  "svcam1_plasma_qc"
)

primary_markers <- primary_markers[primary_markers %in% names(dat)]

for (m in primary_markers) {
  log_m <- paste0("log_", m)
  if (!log_m %in% names(dat)) {
    dat[[log_m]] <- safe_log(dat[[m]])
  }
}

# ---------- 6. Overall covariate missingness ----------

covariate_missingness <- dat %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    non_missing_dx_label = sum(!is.na(dx_label) & dx_label %in% c("CN", "MCI", "AD")),
    non_missing_baseline_dx = sum(!is.na(baseline_dx) & baseline_dx %in% c("CN", "MCI", "AD")),
    non_missing_age = sum(!is.na(age)),
    non_missing_examdate = sum(!is.na(examdate)),
    non_missing_ptdobyy = sum(!is.na(ptdobyy)),
    non_missing_ptgender = sum(!is.na(ptgender)),
    non_missing_pteducat = sum(!is.na(pteducat)),
    non_missing_apoe4 = sum(!is.na(apoe4)),
    non_missing_years_from_baseline = sum(!is.na(years_from_baseline))
  )

# ---------- 7. Age summary ----------

age_values <- dat$age[!is.na(dat$age)]

if (length(age_values) > 0) {
  age_summary <- tibble(
    n_age = length(age_values),
    min_age = min(age_values),
    q1_age = as.numeric(quantile(age_values, 0.25)),
    median_age = median(age_values),
    q3_age = as.numeric(quantile(age_values, 0.75)),
    max_age = max(age_values)
  )
} else {
  age_summary <- tibble(
    n_age = 0,
    min_age = NA_real_,
    q1_age = NA_real_,
    median_age = NA_real_,
    q3_age = NA_real_,
    max_age = NA_real_
  )
}library(tidyverse)
library(janitor)
library(lubridate)
library(readr)

dat_check <- read_csv(
  "R:/ADNI_Project/02_clean_data/analysis_master_biomarkers_longitudinal.csv",
  show_col_types = FALSE
) %>%
  clean_names()

names(dat_check)

dat_check %>%
  summarise(
    n_rows = n(),
    non_missing_examdate = sum(!is.na(examdate)),
    non_missing_ptdobyy = sum(!is.na(ptdobyy)),
    sample_examdate = paste(head(unique(examdate), 10), collapse = " | "),
    sample_ptdobyy = paste(head(unique(ptdobyy), 10), collapse = " | ")
  )
# ---------- 8. Baseline model readiness ----------

baseline_dat <- dat %>%
  filter(visit_key == "bl") %>%
  filter(dx_label %in% c("CN", "MCI", "AD"))

baseline_model_readiness <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  step0 <- baseline_dat
  step1 <- step0 %>% filter(!is.na(.data[[log_m]]))
  step2 <- step1 %>% filter(!is.na(.data$dx_label))
  step3 <- step2 %>% filter(!is.na(.data$age))
  step4 <- step3 %>% filter(!is.na(.data$ptgender))
  step5 <- step4 %>% filter(!is.na(.data$pteducat))
  step6 <- step5 %>% filter(!is.na(.data$apoe4))
  
  tibble(
    marker = m,
    baseline_all_rows = nrow(step0),
    marker_nonmissing_rows = nrow(step1),
    after_dx_rows = nrow(step2),
    after_age_rows = nrow(step3),
    after_ptgender_rows = nrow(step4),
    after_pteducat_rows = nrow(step5),
    after_apoe4_rows = nrow(step6),
    final_subjects = n_distinct(step6$rid),
    final_dx_groups = n_distinct(step6$dx_label)
  )
})

# ---------- 9. Longitudinal model readiness ----------

longitudinal_model_readiness <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  step0 <- dat
  step1 <- step0 %>% filter(!is.na(.data[[log_m]]))
  step2 <- step1 %>% filter(!is.na(.data$years_from_baseline))
  step3 <- step2 %>% filter(!is.na(.data$baseline_dx))
  step4 <- step3 %>% filter(!is.na(.data$age))
  step5 <- step4 %>% filter(!is.na(.data$ptgender))
  step6 <- step5 %>% filter(!is.na(.data$pteducat))
  step7 <- step6 %>% filter(!is.na(.data$apoe4))
  
  tibble(
    marker = m,
    longitudinal_all_rows = nrow(step0),
    marker_nonmissing_rows = nrow(step1),
    after_years_rows = nrow(step2),
    after_baseline_dx_rows = nrow(step3),
    after_age_rows = nrow(step4),
    after_ptgender_rows = nrow(step5),
    after_pteducat_rows = nrow(step6),
    after_apoe4_rows = nrow(step7),
    final_subjects = n_distinct(step7$rid),
    final_dx_groups = n_distinct(step7$baseline_dx)
  )
})

# ---------- 10. Diagnosis counts after full filters ----------

baseline_final_dx_counts <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  baseline_dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(.data$dx_label),
      !is.na(.data$age),
      !is.na(.data$ptgender),
      !is.na(.data$pteducat),
      !is.na(.data$apoe4)
    ) %>%
    count(marker = m, dx_label, name = "n")
})

longitudinal_final_dx_counts <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(.data$years_from_baseline),
      !is.na(.data$baseline_dx),
      !is.na(.data$age),
      !is.na(.data$ptgender),
      !is.na(.data$pteducat),
      !is.na(.data$apoe4)
    ) %>%
    count(marker = m, baseline_dx, name = "n")
})

# ---------- 11. Save QA ----------

openxlsx::write.xlsx(
  list(
    covariate_missingness = covariate_missingness,
    age_summary = age_summary,
    baseline_model_readiness = baseline_model_readiness,
    longitudinal_model_readiness = longitudinal_model_readiness,
    baseline_final_dx_counts = baseline_final_dx_counts,
    longitudinal_final_dx_counts = longitudinal_final_dx_counts
  ),
  file.path(results_dir, "model_readiness_QA.xlsx"),
  overwrite = TRUE
)

write_csv(
  dat,
  file.path(clean_dir, "analysis_master_model_ready_check.csv")
)

cat("Done.\n")
cat("Check this file:\n")
cat(file.path(results_dir, "model_readiness_QA.xlsx"), "\n")