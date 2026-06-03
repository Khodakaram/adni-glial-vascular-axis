# ==========================================================
# 04_build_clinical_core_dataset.R
# Purpose: Build clean longitudinal clinical core from ADNIMERGE2 R data files
# ==========================================================

library(tidyverse)
library(janitor)
library(lubridate)
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

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  
  x <- as.character(x)
  
  out <- suppressWarnings(ymd(x))
  missing_idx <- is.na(out)
  
  out[missing_idx] <- suppressWarnings(mdy(x[missing_idx]))
  
  return(out)
}

make_visit_key <- function(df) {
  df %>%
    mutate(
      viscode = if ("viscode" %in% names(.)) as.character(viscode) else NA_character_,
      viscode2 = if ("viscode2" %in% names(.)) as.character(viscode2) else NA_character_,
      visit_key = coalesce(
        na_if(viscode2, ""),
        na_if(viscode, "")
      )
    )
}

# ---------- Load selected ADNIMERGE2 tables ----------

dxsum <- load_rda_df("DXSUM")
ptdemog <- load_rda_df("PTDEMOG")
adas <- load_rda_df("ADAS")
mmse <- load_rda_df("MMSE")
apoe <- load_rda_df("APOERES")
wmh <- load_rda_df("UCD_WMH")

# ---------- Clean diagnosis table ----------

dx_clean <- dxsum %>%
  make_visit_key() %>%
  transmute(
    rid,
    ptid,
    viscode,
    viscode2,
    visit_key,
    examdate = parse_date_safe(examdate),
    diagnosis
  ) %>%
  distinct()

# ---------- Clean demographics ----------

demo_clean <- ptdemog %>%
  make_visit_key() %>%
  mutate(visdate = parse_date_safe(visdate)) %>%
  arrange(rid, visdate) %>%
  group_by(rid) %>%
  summarise(
    ptid_demo = first_nonmissing(ptid),
    ptgender = first_nonmissing(ptgender),
    pteducat = first_nonmissing(pteducat),
    ptdobyy = first_nonmissing(ptdobyy),
    .groups = "drop"
  )

# ---------- Clean ADAS ----------

adas_clean <- adas %>%
  make_visit_key() %>%
  transmute(
    rid,
    visit_key,
    adas_visdate = parse_date_safe(visdate),
    adas11 = totscore,
    adas13 = total13
  ) %>%
  group_by(rid, visit_key) %>%
  summarise(
    adas_visdate = first_nonmissing(adas_visdate),
    adas11 = first_nonmissing(adas11),
    adas13 = first_nonmissing(adas13),
    .groups = "drop"
  )

# ---------- Clean MMSE ----------

mmse_clean <- mmse %>%
  make_visit_key() %>%
  transmute(
    rid,
    visit_key,
    mmse_visdate = parse_date_safe(visdate),
    mmse = mmscore
  ) %>%
  group_by(rid, visit_key) %>%
  summarise(
    mmse_visdate = first_nonmissing(mmse_visdate),
    mmse = first_nonmissing(mmse),
    .groups = "drop"
  )

# ---------- Clean APOE ----------

apoe_clean <- apoe %>%
  group_by(rid) %>%
  summarise(
    apoe_genotype = first_nonmissing(genotype),
    .groups = "drop"
  ) %>%
  mutate(
    apoe4 = if_else(str_detect(as.character(apoe_genotype), "4"), 1, 0)
  )

# ---------- Clean WMH ----------

wmh_clean <- wmh %>%
  make_visit_key() %>%
  transmute(
    rid,
    visit_key,
    wmh_examdate = parse_date_safe(examdate),
    total_wmh
  ) %>%
  group_by(rid, visit_key) %>%
  summarise(
    wmh_examdate = first_nonmissing(wmh_examdate),
    total_wmh = first_nonmissing(total_wmh),
    .groups = "drop"
  )

# ---------- Merge clinical core ----------

clinical_core <- dx_clean %>%
  left_join(demo_clean, by = "rid") %>%
  left_join(apoe_clean, by = "rid") %>%
  left_join(adas_clean, by = c("rid", "visit_key")) %>%
  left_join(mmse_clean, by = c("rid", "visit_key")) %>%
  left_join(wmh_clean, by = c("rid", "visit_key")) %>%
  mutate(
    diagnosis_chr = as.character(diagnosis),
    dx_label = case_when(
      diagnosis_chr %in% c("1", "NL", "CN", "Normal") ~ "CN",
      diagnosis_chr %in% c("2", "MCI") ~ "MCI",
      diagnosis_chr %in% c("3", "AD", "Dementia") ~ "AD",
      TRUE ~ diagnosis_chr
    )
  ) %>%
  group_by(rid) %>%
  mutate(
    baseline_date = if_else(
      all(is.na(examdate)),
      as.Date(NA),
      min(examdate, na.rm = TRUE)
    ),
    years_from_baseline = as.numeric(examdate - baseline_date) / 365.25
  ) %>%
  ungroup()

# ---------- Quality checks ----------

table_summary <- tibble(
  table = c("DXSUM", "PTDEMOG", "ADAS", "MMSE", "APOERES", "UCD_WMH", "clinical_core"),
  n_rows = c(
    nrow(dxsum),
    nrow(ptdemog),
    nrow(adas),
    nrow(mmse),
    nrow(apoe),
    nrow(wmh),
    nrow(clinical_core)
  ),
  n_subjects = c(
    n_distinct(dxsum$rid),
    n_distinct(ptdemog$rid),
    n_distinct(adas$rid),
    n_distinct(mmse$rid),
    n_distinct(apoe$rid),
    n_distinct(wmh$rid),
    n_distinct(clinical_core$rid)
  )
)

diagnosis_counts <- clinical_core %>%
  count(dx_label, sort = TRUE)

visit_counts <- clinical_core %>%
  count(visit_key, sort = TRUE)

missingness <- clinical_core %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    non_missing_dx = sum(!is.na(dx_label)),
    non_missing_ptgender = sum(!is.na(ptgender)),
    non_missing_pteducat = sum(!is.na(pteducat)),
    non_missing_apoe4 = sum(!is.na(apoe4)),
    non_missing_adas13 = sum(!is.na(adas13)),
    non_missing_mmse = sum(!is.na(mmse)),
    non_missing_total_wmh = sum(!is.na(total_wmh))
  )

duplicate_keys <- clinical_core %>%
  count(rid, visit_key) %>%
  filter(n > 1)

# ---------- Save outputs ----------

write_csv(
  clinical_core,
  file.path(clean_dir, "clinical_core_longitudinal.csv")
)

openxlsx::write.xlsx(
  list(
    table_summary = table_summary,
    diagnosis_counts = diagnosis_counts,
    visit_counts = visit_counts,
    missingness = missingness,
    duplicate_keys = duplicate_keys
  ),
  file.path(results_dir, "clinical_core_QA.xlsx"),
  overwrite = TRUE
)

cat("Done.\n")
cat("Clinical core file:\n")
cat(file.path(clean_dir, "clinical_core_longitudinal.csv"), "\n\n")
cat("QA file:\n")
cat(file.path(results_dir, "clinical_core_QA.xlsx"), "\n")