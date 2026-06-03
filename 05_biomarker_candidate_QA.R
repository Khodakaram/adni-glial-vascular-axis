# ==========================================================
# 05_biomarker_candidate_QA.R
# Purpose: Inspect candidate biomarker files before merging
# ==========================================================

library(tidyverse)
library(readxl)
library(janitor)
library(openxlsx)
library(stringr)

project_dir <- "R:/ADNI_Project"

bio_dir <- file.path(project_dir, "00_raw_data", "biomarkers_excel")
clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. Load clinical core ----------

clinical_unique_path <- file.path(clean_dir, "clinical_core_longitudinal_unique.csv")
clinical_path <- file.path(clean_dir, "clinical_core_longitudinal.csv")

if (file.exists(clinical_unique_path)) {
  clinical_core <- read_csv(clinical_unique_path, show_col_types = FALSE)
} else {
  clinical_core <- read_csv(clinical_path, show_col_types = FALSE) %>%
    group_by(rid, visit_key) %>%
    summarise(across(everything(), ~ .x[which(!is.na(.x))[1]]), .groups = "drop")
}

clinical_keys <- clinical_core %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key)
  ) %>%
  distinct(rid, visit_key)

# ---------- 2. Candidate files from inventory ----------

candidate_files <- tibble::tribble(
  ~marker_group, ~priority, ~file_name, ~notes,
  
  "sTREM2", "Primary",
  "ADNI_HAASS_WASHU_LAB_01Jun2026.csv",
  "Best candidate for soluble TREM2; use corrected variables after QA.",
  
  "sTREM2", "Duplicate older/date-near",
  "ADNI_HAASS_WASHU_LAB_31May2026.csv",
  "Do not use together with 01Jun2026 unless needed for comparison.",
  
  "GFAP plasma", "Primary",
  "UPENN_PLASMA_FUJIREBIO_QUANTERIX_01Jun2026.csv",
  "Contains GFAP_Q and GFAP_F.",
  
  "GFAP plasma", "Duplicate older/date-near",
  "UPENN_PLASMA_FUJIREBIO_QUANTERIX_31May2026.csv",
  "Do not use together with 01Jun2026 unless needed for comparison.",
  
  "GFAP plasma", "Secondary",
  "PLASMA_ABETA_PROJECT_ADX_VUMC_31May2026.csv",
  "Contains GFAP; useful as secondary/sensitivity dataset.",
  
  "GFAP CSF / YKL40", "Primary CSF proteomics",
  "CSFMRM_31May2026.csv",
  "Contains GFAP_ALAAELNQLR and CH3L1 peptide variables.",
  
  "GFAP CSF / YKL40", "Older duplicate/proteomics",
  "Biomarkers Consortium CSF Proteomics MRM Data.csv",
  "Similar CSF MRM data; use only if needed.",
  
  "YKL40 / TREM2 proteomics", "Secondary",
  "adni_csfproteomics.csv",
  "Contains CH3L1 and TREM2 peptide variables with HUMAN-style names.",
  
  "VEGF / ICAM / VCAM", "Primary wide candidate",
  "ADMC_CLINICALVARIABLES_16May2016.csv",
  "Contains VEGF, ICAM, VCAM in wide format.",
  
  "VEGF / ICAM / VCAM plasma", "Secondary wide QC",
  "adni_plasma_qc_multiplex_11Nov2010.csv",
  "Contains VEGF, ICAM, VCAM with long column names.",
  
  "VEGF / ICAM / VCAM plasma", "Raw long format",
  "adni_plasma_raw_multiplex_11Nov2010.csv",
  "May contain biomarkers in analyte rows, not in column names.",
  
  "VEGF / ICAM / VCAM CSF", "Secondary CSF QC",
  "Biomarkers Consortium ADNI CSF QC Multiplex data.csv",
  "Contains VEGF, ICAM, VCAM but merge key must be checked.",
  
  "VEGF / ICAM / VCAM CSF", "Raw long format",
  "Biomarkers Consortium ADNI CSF Multiplex Raw Data.csv",
  "May contain biomarkers in analyte rows, not in column names.",
  
  "WMH", "Primary",
  "UCD_WMH_01Jun2026.csv",
  "Contains TOTAL_WMH.",
  
  "WMH", "Duplicate older/date-near",
  "UCD_WMH_31May2026.csv",
  "Do not use together with 01Jun2026 unless needed for comparison."
) %>%
  mutate(file_path = file.path(bio_dir, file_name))

# ---------- 3. Helper functions ----------

read_any_file <- function(path) {
  ext <- tools::file_ext(path) %>% tolower()
  
  if (ext == "csv") {
    read_csv(path, show_col_types = FALSE, guess_max = 100000)
  } else if (ext %in% c("xlsx", "xls")) {
    read_excel(path, sheet = 1)
  } else {
    stop("Unsupported file type: ", path)
  }
}

coalesce_existing <- function(df, cols) {
  existing <- intersect(cols, names(df))
  
  if (length(existing) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  
  values <- lapply(existing, function(x) as.character(df[[x]]))
  purrr::reduce(values, dplyr::coalesce)
}

make_keys <- function(df) {
  df <- df %>% clean_names()
  
  df$rid <- if ("rid" %in% names(df)) as.character(df$rid) else NA_character_
  
  df$visit_key <- coalesce_existing(
    df,
    c("viscode2", "viscode", "visit_code", "visit")
  )
  
  df
}

marker_patterns <- list(
  gfap = "GFAP",
  strem2 = "STREM2|TREM2",
  ykl40_ch3l1 = "YKL|CHI3L1|CH3L1",
  vegf = "VEGF|VASCULAR ENDOTHELIAL",
  icam = "ICAM|INTERCELLULAR ADHESION",
  vcam = "VCAM|VASCULAR CELL ADHESION",
  wmh = "TOTAL_WMH|WMH"
)

scan_one_file <- function(file_path, file_name, marker_group, priority, notes) {
  
  if (!file.exists(file_path)) {
    return(tibble(
      file_name = file_name,
      marker_group = marker_group,
      priority = priority,
      file_exists = FALSE,
      n_rows = NA_integer_,
      n_subjects = NA_integer_,
      n_visit_keys = NA_integer_,
      clinical_key_matches = NA_integer_,
      marker = NA_character_,
      matched_columns = NA_character_,
      non_missing_in_matched_columns = NA_integer_,
      analyte_row_hits = NA_integer_,
      analyte_examples = NA_character_,
      notes = notes
    ))
  }
  
  raw <- read_any_file(file_path)
  df <- make_keys(raw)
  
  n_rows <- nrow(df)
  n_subjects <- if ("rid" %in% names(df)) n_distinct(df$rid, na.rm = TRUE) else NA_integer_
  n_visit_keys <- if ("visit_key" %in% names(df)) n_distinct(df$visit_key, na.rm = TRUE) else NA_integer_
  
  file_keys <- df %>%
    filter(!is.na(rid), !is.na(visit_key)) %>%
    distinct(rid, visit_key)
  
  clinical_key_matches <- inner_join(
    file_keys,
    clinical_keys,
    by = c("rid", "visit_key")
  ) %>%
    nrow()
  
  analyte_cols <- names(df)[str_detect(names(df), "analyte|ana_unit|biomarker|test")]
  
  map_dfr(names(marker_patterns), function(marker_name) {
    
    pattern <- marker_patterns[[marker_name]]
    
    matched_cols <- names(df)[str_detect(toupper(names(df)), pattern)]
    
    non_missing_col_values <- if (length(matched_cols) > 0) {
      sum(map_int(matched_cols, function(col) {
        sum(!is.na(df[[col]]) & as.character(df[[col]]) != "")
      }))
    } else {
      0
    }
    
    analyte_hit_rows <- 0
    analyte_examples <- NA_character_
    
    if (length(analyte_cols) > 0) {
      analyte_text <- df %>%
        select(all_of(analyte_cols)) %>%
        mutate(across(everything(), as.character))
      
      row_hit_logical <- apply(analyte_text, 1, function(x) {
        any(str_detect(toupper(x), pattern), na.rm = TRUE)
      })
      
      analyte_hit_rows <- sum(row_hit_logical, na.rm = TRUE)
      
      if (analyte_hit_rows > 0) {
        examples <- analyte_text[row_hit_logical, , drop = FALSE] %>%
          unlist(use.names = FALSE) %>%
          unique()
        
        analyte_examples <- paste(head(examples, 5), collapse = " | ")
      }
    }
    
    tibble(
      file_name = file_name,
      marker_group = marker_group,
      priority = priority,
      file_exists = TRUE,
      n_rows = n_rows,
      n_subjects = n_subjects,
      n_visit_keys = n_visit_keys,
      clinical_key_matches = clinical_key_matches,
      marker = marker_name,
      matched_columns = paste(matched_cols, collapse = " | "),
      non_missing_in_matched_columns = non_missing_col_values,
      analyte_row_hits = analyte_hit_rows,
      analyte_examples = analyte_examples,
      notes = notes
    )
  })
}

# ---------- 4. Run scan ----------

biomarker_QA <- pmap_dfr(
  candidate_files,
  function(marker_group, priority, file_name, notes, file_path) {
    scan_one_file(file_path, file_name, marker_group, priority, notes)
  }
)

file_level_summary <- biomarker_QA %>%
  group_by(file_name, marker_group, priority, file_exists, n_rows, n_subjects, n_visit_keys, clinical_key_matches, notes) %>%
  summarise(
    detected_markers = paste(marker[
      matched_columns != "" | analyte_row_hits > 0
    ], collapse = ", "),
    .groups = "drop"
  )

useful_marker_hits <- biomarker_QA %>%
  filter(
    matched_columns != "" |
      analyte_row_hits > 0
  ) %>%
  arrange(marker, desc(clinical_key_matches), file_name)

# ---------- 5. Save output ----------

openxlsx::write.xlsx(
  list(
    file_level_summary = file_level_summary,
    useful_marker_hits = useful_marker_hits,
    full_biomarker_QA = biomarker_QA,
    candidate_files = candidate_files
  ),
  file.path(results_dir, "biomarker_candidate_QA.xlsx"),
  overwrite = TRUE
)

cat("Done. Check this file:\n")
cat(file.path(results_dir, "biomarker_candidate_QA.xlsx"), "\n")