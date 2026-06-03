# ==========================================================
# 11_repair_age_from_ptdemog.R
# Purpose: Repair age using PTDEMOG.rda VISDATE + PTDOBYY
# ==========================================================

library(tidyverse)
library(janitor)
library(lubridate)
library(readr)
library(openxlsx)

project_dir <- "R:/ADNI_Project"

data_dir <- file.path(project_dir, "00_raw_data", "ADNIMERGE2", "data")
clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------

load_rda_df <- function(file_stem) {
  
  file <- list.files(
    data_dir,
    pattern = paste0("^", file_stem, "\\.rda$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(file) == 0) {
    stop(paste("File not found:", file_stem))
  }
  
  env <- new.env()
  object_names <- load(file[1], envir = env)
  
  df <- get(object_names[1], envir = env) %>%
    as_tibble() %>%
    clean_names()
  
  return(df)
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  return(x[1])
}

parse_date_robust <- function(x) {
  
  if (inherits(x, "Date")) return(x)
  
  x_chr <- as.character(x)
  x_chr[x_chr == "" | x_chr == "NA" | x_chr == "NaN"] <- NA_character_
  
  out <- suppressWarnings(ymd(x_chr))
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(mdy(x_chr[idx]))
  
  idx <- is.na(out)
  out[idx] <- suppressWarnings(dmy(x_chr[idx]))
  
  numeric_x <- suppressWarnings(as.numeric(x_chr))
  idx <- is.na(out) & !is.na(numeric_x)
  
  if (any(idx)) {
    
    candidate_r <- as.Date(numeric_x[idx], origin = "1970-01-01")
    candidate_excel <- as.Date(numeric_x[idx], origin = "1899-12-30")
    
    plausible_r <- candidate_r >= as.Date("2000-01-01") &
      candidate_r <= as.Date("2035-12-31")
    
    candidate <- candidate_r
    candidate[!plausible_r] <- candidate_excel[!plausible_r]
    
    out[idx] <- candidate
  }
  
  return(out)
}

coalesce_date_existing <- function(df, cols) {
  
  existing <- intersect(cols, names(df))
  
  if (length(existing) == 0) {
    return(as.Date(rep(NA_character_, nrow(df))))
  }
  
  date_list <- lapply(existing, function(col) {
    parse_date_robust(df[[col]])
  })
  
  purrr::reduce(date_list, dplyr::coalesce)
}

safe_log <- function(x) {
  x <- parse_number(as.character(x))
  ifelse(!is.na(x) & x > 0, log(x), NA_real_)
}

# ---------- 1. Load current merged biomarker dataset ----------

analysis_file <- file.path(clean_dir, "analysis_master_biomarkers_longitudinal.csv")

dat <- read_csv(analysis_file, show_col_types = FALSE, guess_max = 100000) %>%
  clean_names() %>%
  select(-any_of(c("age", "baseline_dx", "baseline_dx_x", "baseline_dx_y"))) %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key),
    dx_label = as.character(dx_label),
    ptdobyy = parse_number(as.character(ptdobyy)),
    pteducat = parse_number(as.character(pteducat)),
    apoe4 = parse_number(as.character(apoe4)),
    years_from_baseline = parse_number(as.character(years_from_baseline))
  )

# ---------- 2. Load PTDEMOG ----------

ptdemog <- load_rda_df("PTDEMOG")

# ---------- 3. Create age source from PTDEMOG ----------

ptdemog_clean <- ptdemog %>%
  mutate(
    rid = as.character(rid),
    ptdobyy_demo = parse_number(as.character(ptdobyy))
  )

ptdemog_clean$ptdemog_date <- coalesce_date_existing(
  ptdemog_clean,
  c("visdate", "examdate", "userdate", "userdate2", "update_stamp")
)

ptdemog_age <- ptdemog_clean %>%
  filter(!is.na(rid)) %>%
  group_by(rid) %>%
  summarise(
    ptdobyy_demo = first_nonmissing(ptdobyy_demo),
    first_ptdemog_date = suppressWarnings(min(ptdemog_date, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_ptdemog_date = if_else(
      is.infinite(as.numeric(first_ptdemog_date)),
      as.Date(NA),
      first_ptdemog_date
    ),
    baseline_age_ptdemog = if_else(
      !is.na(first_ptdemog_date) & !is.na(ptdobyy_demo),
      as.numeric(year(first_ptdemog_date) - ptdobyy_demo),
      NA_real_
    )
  )

# ---------- 4. Join age into analysis dataset ----------

dat_ready <- dat %>%
  left_join(ptdemog_age, by = "rid") %>%
  mutate(
    ptdobyy_final = coalesce(ptdobyy, ptdobyy_demo),
    
    age = case_when(
      !is.na(baseline_age_ptdemog) & !is.na(years_from_baseline) ~
        baseline_age_ptdemog + years_from_baseline,
      
      !is.na(baseline_age_ptdemog) ~
        baseline_age_ptdemog,
      
      TRUE ~ NA_real_
    )
  )

# ---------- 5. Recreate baseline diagnosis ----------

baseline_dx_tbl <- dat_ready %>%
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

dat_ready <- dat_ready %>%
  left_join(baseline_dx_tbl, by = "rid") %>%
  mutate(
    dx_label = factor(dx_label, levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender)
  )

# ---------- 6. Recreate log biomarkers ----------

primary_markers <- c(
  "gfap_quanterix",
  "strem2_msd_corrected",
  "vegf_plasma_qc",
  "sicam1_plasma_qc",
  "svcam1_plasma_qc"
)

primary_markers <- primary_markers[primary_markers %in% names(dat_ready)]

for (m in primary_markers) {
  dat_ready[[paste0("log_", m)]] <- safe_log(dat_ready[[m]])
}

# ---------- 7. QA ----------

covariate_check <- dat_ready %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    non_missing_ptdemog_date = sum(!is.na(first_ptdemog_date)),
    non_missing_ptdobyy_final = sum(!is.na(ptdobyy_final)),
    non_missing_age = sum(!is.na(age)),
    non_missing_ptgender = sum(!is.na(ptgender)),
    non_missing_pteducat = sum(!is.na(pteducat)),
    non_missing_apoe4 = sum(!is.na(apoe4)),
    non_missing_years_from_baseline = sum(!is.na(years_from_baseline))
  )

age_check <- dat_ready %>%
  filter(!is.na(age)) %>%
  summarise(
    n_age = n(),
    min_age = min(age),
    q1_age = quantile(age, 0.25),
    median_age = median(age),
    q3_age = quantile(age, 0.75),
    max_age = max(age)
  )

baseline_dat <- dat_ready %>%
  filter(visit_key == "bl") %>%
  filter(dx_label %in% c("CN", "MCI", "AD"))

baseline_model_readiness <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  df <- baseline_dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(dx_label),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
  
  tibble(
    marker = m,
    final_rows = nrow(df),
    final_subjects = n_distinct(df$rid),
    final_dx_groups = n_distinct(df$dx_label)
  )
})

longitudinal_model_readiness <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  df <- dat_ready %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
  
  tibble(
    marker = m,
    final_rows = nrow(df),
    final_subjects = n_distinct(df$rid),
    final_dx_groups = n_distinct(df$baseline_dx)
  )
})

ptdemog_date_diagnostics <- tibble(
  ptdemog_columns = paste(names(ptdemog_clean), collapse = " | "),
  n_ptdemog_rows = nrow(ptdemog_clean),
  n_ptdemog_subjects = n_distinct(ptdemog_clean$rid),
  non_missing_ptdemog_date = sum(!is.na(ptdemog_clean$ptdemog_date)),
  non_missing_ptdobyy_demo = sum(!is.na(ptdemog_clean$ptdobyy_demo))
)

# ---------- 8. Save ----------

write_csv(
  dat_ready,
  file.path(clean_dir, "analysis_master_model_ready.csv")
)

openxlsx::write.xlsx(
  list(
    covariate_check = covariate_check,
    age_check = age_check,
    baseline_model_readiness = baseline_model_readiness,
    longitudinal_model_readiness = longitudinal_model_readiness,
    ptdemog_date_diagnostics = ptdemog_date_diagnostics
  ),
  file.path(results_dir, "age_repair_from_ptdemog_QA.xlsx"),
  overwrite = TRUE
)

cat("Done.\n")
cat("Model-ready dataset:\n")
cat(file.path(clean_dir, "analysis_master_model_ready.csv"), "\n\n")
cat("QA file:\n")
cat(file.path(results_dir, "age_repair_from_ptdemog_QA.xlsx"), "\n")