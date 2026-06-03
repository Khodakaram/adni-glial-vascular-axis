# ==========================================================
# 06_merge_selected_biomarkers.R
# Purpose: Merge selected ADNI biomarkers with clinical core
# ==========================================================

library(tidyverse)
library(janitor)
library(openxlsx)
library(readr)

project_dir <- "R:/ADNI_Project"

bio_dir <- file.path(project_dir, "00_raw_data", "biomarkers_excel")
clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  return(x[1])
}

parse_num <- function(x) {
  readr::parse_number(as.character(x))
}

get_col <- function(df, col) {
  if (col %in% names(df)) {
    return(df[[col]])
  } else {
    return(rep(NA, nrow(df)))
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

standardize_keys <- function(df) {
  df <- df %>% clean_names()
  
  df %>%
    mutate(
      rid = as.character(get_col(., "rid")),
      visit_key = coalesce_existing(
        .,
        c("viscode2", "viscode", "visit_code", "visit")
      )
    )
}

collapse_by_visit <- function(df) {
  df %>%
    filter(!is.na(rid), !is.na(visit_key)) %>%
    group_by(rid, visit_key) %>%
    summarise(
      across(everything(), first_nonmissing),
      .groups = "drop"
    )
}

row_mean_safe <- function(df, cols) {
  existing <- intersect(cols, names(df))
  
  if (length(existing) == 0) {
    return(rep(NA_real_, nrow(df)))
  }
  
  mat <- df %>%
    select(all_of(existing)) %>%
    mutate(across(everything(), parse_num)) %>%
    as.matrix()
  
  out <- rowMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  return(out)
}

read_bio_csv <- function(file_name) {
  path <- file.path(bio_dir, file_name)
  
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  
  read_csv(path, show_col_types = FALSE, guess_max = 100000) %>%
    standardize_keys()
}

# ---------- 1. Load clinical core ----------

clinical_unique_path <- file.path(clean_dir, "clinical_core_longitudinal_unique.csv")
clinical_path <- file.path(clean_dir, "clinical_core_longitudinal.csv")

if (file.exists(clinical_unique_path)) {
  clinical_core <- read_csv(clinical_unique_path, show_col_types = FALSE)
} else {
  clinical_core <- read_csv(clinical_path, show_col_types = FALSE) %>%
    mutate(
      rid = as.character(rid),
      visit_key = as.character(visit_key)
    ) %>%
    group_by(rid, visit_key) %>%
    summarise(
      across(everything(), first_nonmissing),
      .groups = "drop"
    )
}

clinical_core <- clinical_core %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key)
  )

# ---------- 2. sTREM2: ADNI_HAASS_WASHU_LAB_01Jun2026 ----------

strem2_raw <- read_bio_csv("ADNI_HAASS_WASHU_LAB_01Jun2026.csv")

strem2_clean <- strem2_raw %>%
  transmute(
    rid,
    visit_key,
    strem2_wu = parse_num(get_col(., "wu_strem2")),
    strem2_wu_cv = parse_num(get_col(., "wu_strem2_cv")),
    strem2_wu_corrected = parse_num(get_col(., "wu_strem2corrected")),
    strem2_msd = parse_num(get_col(., "msd_strem2")),
    strem2_msd_cv = parse_num(get_col(., "msd_strem2_cv")),
    strem2_msd_corrected = parse_num(get_col(., "msd_strem2corrected")),
    strem2_outlier = as.character(get_col(., "trem2outlier"))
  ) %>%
  collapse_by_visit()

# ---------- 3. Plasma GFAP: UPENN plasma Fujirebio / Quanterix ----------

gfap_plasma_raw <- read_bio_csv("UPENN_PLASMA_FUJIREBIO_QUANTERIX_01Jun2026.csv")

gfap_plasma_clean <- gfap_plasma_raw %>%
  transmute(
    rid,
    visit_key,
    gfap_quanterix = parse_num(get_col(., "gfap_q")),
    gfap_fujirebio = parse_num(get_col(., "gfap_f"))
  ) %>%
  collapse_by_visit()

# ---------- 4. Vascular plasma markers: QC multiplex ----------

vascular_qc_raw <- read_bio_csv("adni_plasma_qc_multiplex_11Nov2010.csv")

vascular_qc_clean <- vascular_qc_raw %>%
  transmute(
    rid,
    visit_key,
    vegf_plasma_qc = parse_num(get_col(., "vascular_endothelial_growth_factor_vegf_pg_m_l")),
    sicam1_plasma_qc = parse_num(get_col(., "intercellular_adhesion_molecule_1_icam_ng_m_l")),
    svcam1_plasma_qc = parse_num(get_col(., "vascular_cell_adhesion_molecule_1_vcam_ng_m_l"))
  ) %>%
  collapse_by_visit()

# ---------- 5. Vascular ADMC source, kept for checking/sensitivity ----------

admc_raw <- read_bio_csv("ADMC_CLINICALVARIABLES_16May2016.csv")

admc_clean <- admc_raw %>%
  transmute(
    rid,
    visit_key,
    vegf_admc = parse_num(get_col(., "vegf")),
    sicam1_admc = parse_num(get_col(., "icam")),
    svcam1_admc = parse_num(get_col(., "vcam"))
  ) %>%
  collapse_by_visit()

# ---------- 6. CSF MRM: GFAP and CH3L1 / YKL-40 ----------

csfmrm_raw <- read_bio_csv("CSFMRM_31May2026.csv")

csfmrm_clean <- csfmrm_raw %>%
  mutate(
    gfap_csf_mrm = parse_num(get_col(., "gfap_alaaelnqlr")),
    ch3l1_ilgqqvpyatk = parse_num(get_col(., "ch3l1_ilgqqvpyatk")),
    ch3l1_sftlassetgvgapisgpgipgr = parse_num(get_col(., "ch3l1_sftlassetgvgapisgpgipgr")),
    ch3l1_vtidssydiak = parse_num(get_col(., "ch3l1_vtidssydiak"))
  ) %>%
  mutate(
    ykl40_csf_mrm_mean = row_mean_safe(
      .,
      c(
        "ch3l1_ilgqqvpyatk",
        "ch3l1_sftlassetgvgapisgpgipgr",
        "ch3l1_vtidssydiak"
      )
    )
  ) %>%
  select(
    rid,
    visit_key,
    gfap_csf_mrm,
    ch3l1_ilgqqvpyatk,
    ch3l1_sftlassetgvgapisgpgipgr,
    ch3l1_vtidssydiak,
    ykl40_csf_mrm_mean
  ) %>%
  collapse_by_visit()

# ---------- 7. CSF proteomics longitudinal: CH3L1/YKL-40 and TREM2 peptides ----------

csfprot_raw <- read_bio_csv("adni_csfproteomics.csv")

csfprot_clean <- csfprot_raw %>%
  mutate(
    trem2_prot_vlvevladpldhr_y2 = parse_num(get_col(., "trem2_human_vlvevladpldhr_y2")),
    trem2_prot_vlvevladpldhr_y5 = parse_num(get_col(., "trem2_human_vlvevladpldhr_y5")),
    trem2_prot_vvsthnlwllsflr_y4 = parse_num(get_col(., "trem2_human_vvsthnlwllsflr_y4")),
    trem2_prot_vvsthnlwllsflr_y5 = parse_num(get_col(., "trem2_human_vvsthnlwllsflr_y5")),
    
    ch3l1_prot_ilgqqvpyatk_y5 = parse_num(get_col(., "ch3l1_human_ilgqqvpyatk_y5")),
    ch3l1_prot_ilgqqvpyatk_y9 = parse_num(get_col(., "ch3l1_human_ilgqqvpyatk_y9")),
    ch3l1_prot_vtidssydiak_y7 = parse_num(get_col(., "ch3l1_human_vtidssydiak_y7")),
    ch3l1_prot_vtidssydiak_y8 = parse_num(get_col(., "ch3l1_human_vtidssydiak_y8"))
  ) %>%
  mutate(
    trem2_csfprot_mean = row_mean_safe(
      .,
      c(
        "trem2_prot_vlvevladpldhr_y2",
        "trem2_prot_vlvevladpldhr_y5",
        "trem2_prot_vvsthnlwllsflr_y4",
        "trem2_prot_vvsthnlwllsflr_y5"
      )
    ),
    ykl40_csfprot_mean = row_mean_safe(
      .,
      c(
        "ch3l1_prot_ilgqqvpyatk_y5",
        "ch3l1_prot_ilgqqvpyatk_y9",
        "ch3l1_prot_vtidssydiak_y7",
        "ch3l1_prot_vtidssydiak_y8"
      )
    )
  ) %>%
  select(
    rid,
    visit_key,
    trem2_prot_vlvevladpldhr_y2,
    trem2_prot_vlvevladpldhr_y5,
    trem2_prot_vvsthnlwllsflr_y4,
    trem2_prot_vvsthnlwllsflr_y5,
    trem2_csfprot_mean,
    ch3l1_prot_ilgqqvpyatk_y5,
    ch3l1_prot_ilgqqvpyatk_y9,
    ch3l1_prot_vtidssydiak_y7,
    ch3l1_prot_vtidssydiak_y8,
    ykl40_csfprot_mean
  ) %>%
  collapse_by_visit()

# ---------- 8. WMH CSV source ----------

wmh_raw <- read_bio_csv("UCD_WMH_01Jun2026.csv")

wmh_clean <- wmh_raw %>%
  transmute(
    rid,
    visit_key,
    total_wmh_csv = parse_num(get_col(., "total_wmh"))
  ) %>%
  collapse_by_visit()

# ---------- 9. Merge all selected biomarkers with clinical core ----------

analysis_master <- clinical_core %>%
  left_join(strem2_clean, by = c("rid", "visit_key")) %>%
  left_join(gfap_plasma_clean, by = c("rid", "visit_key")) %>%
  left_join(vascular_qc_clean, by = c("rid", "visit_key")) %>%
  left_join(admc_clean, by = c("rid", "visit_key")) %>%
  left_join(csfmrm_clean, by = c("rid", "visit_key")) %>%
  left_join(csfprot_clean, by = c("rid", "visit_key")) %>%
  left_join(wmh_clean, by = c("rid", "visit_key"))

# ---------- 10. Create QA outputs ----------

biomarker_columns <- c(
  "strem2_wu_corrected",
  "strem2_msd_corrected",
  "gfap_quanterix",
  "gfap_fujirebio",
  "vegf_plasma_qc",
  "sicam1_plasma_qc",
  "svcam1_plasma_qc",
  "vegf_admc",
  "sicam1_admc",
  "svcam1_admc",
  "gfap_csf_mrm",
  "ykl40_csf_mrm_mean",
  "trem2_csfprot_mean",
  "ykl40_csfprot_mean",
  "total_wmh_csv"
)

biomarker_missingness <- analysis_master %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    across(
      any_of(biomarker_columns),
      ~ sum(!is.na(.)),
      .names = "non_missing_{.col}"
    )
  )

biomarker_subject_counts <- analysis_master %>%
  summarise(
    across(
      any_of(biomarker_columns),
      ~ n_distinct(rid[!is.na(.)]),
      .names = "subjects_{.col}"
    )
  )

diagnosis_by_marker <- map_dfr(
  intersect(biomarker_columns, names(analysis_master)),
  function(marker) {
    analysis_master %>%
      filter(!is.na(.data[[marker]])) %>%
      count(marker = marker, dx_label, name = "n_rows")
  }
)

duplicate_keys_after_merge <- analysis_master %>%
  count(rid, visit_key) %>%
  filter(n > 1)

table_summary <- tibble(
  table = c(
    "clinical_core",
    "strem2_clean",
    "gfap_plasma_clean",
    "vascular_qc_clean",
    "admc_clean",
    "csfmrm_clean",
    "csfprot_clean",
    "wmh_clean",
    "analysis_master"
  ),
  n_rows = c(
    nrow(clinical_core),
    nrow(strem2_clean),
    nrow(gfap_plasma_clean),
    nrow(vascular_qc_clean),
    nrow(admc_clean),
    nrow(csfmrm_clean),
    nrow(csfprot_clean),
    nrow(wmh_clean),
    nrow(analysis_master)
  ),
  n_subjects = c(
    n_distinct(clinical_core$rid),
    n_distinct(strem2_clean$rid),
    n_distinct(gfap_plasma_clean$rid),
    n_distinct(vascular_qc_clean$rid),
    n_distinct(admc_clean$rid),
    n_distinct(csfmrm_clean$rid),
    n_distinct(csfprot_clean$rid),
    n_distinct(wmh_clean$rid),
    n_distinct(analysis_master$rid)
  )
)

# ---------- 11. Save outputs ----------

write_csv(
  analysis_master,
  file.path(clean_dir, "analysis_master_biomarkers_longitudinal.csv")
)

openxlsx::write.xlsx(
  list(
    table_summary = table_summary,
    biomarker_missingness = biomarker_missingness,
    biomarker_subject_counts = biomarker_subject_counts,
    diagnosis_by_marker = diagnosis_by_marker,
    duplicate_keys_after_merge = duplicate_keys_after_merge
  ),
  file.path(results_dir, "analysis_master_biomarkers_QA.xlsx"),
  overwrite = TRUE
)

cat("Done.\n")
cat("Merged dataset:\n")
cat(file.path(clean_dir, "analysis_master_biomarkers_longitudinal.csv"), "\n\n")
cat("QA file:\n")
cat(file.path(results_dir, "analysis_master_biomarkers_QA.xlsx"), "\n")