# ==========================================================
# 03_select_adnimerge2_core_candidates.R
# Purpose: Find ADNIMERGE2 tables useful for the clinical core dataset
# ==========================================================

library(tidyverse)
library(readxl)
library(openxlsx)
library(janitor)
library(stringr)

project_dir <- "R:/ADNI_Project"
inventory_file <- file.path(project_dir, "01_inventory", "ADNIMERGE2_table_inventory.xlsx")
out_dir <- file.path(project_dir, "01_inventory")

# ---- 1. Read inventory ----
adni_inv <- read_excel(inventory_file, sheet = "important_tables") %>%
  clean_names()

names(adni_inv)

# ---- 2. Define variable groups needed for your project ----
id_terms <- "RID|PTID|VISCODE|VISCODE2|EXAMDATE|VISIT"
diagnosis_terms <- "DX|DIAGNOSIS|DXCHANGE|PHC|DEMODX"
demo_terms <- "AGE|SEX|PTGENDER|GENDER|PTEDUCAT|EDUCATION"
apoe_terms <- "APOE|APOE4|GENOTYPE"
adas_terms <- "ADAS|ADAS11|ADAS13|TOTAL13"
mmse_terms <- "MMSE|MMSCORE"
wmh_terms <- "WMH|WHITE|HYPERINTENSITY"

# ---- 3. Score each table based on useful variables ----
candidate_tables <- adni_inv %>%
  mutate(
    has_id = str_detect(toupper(columns), id_terms),
    has_diagnosis = str_detect(toupper(columns), diagnosis_terms),
    has_demo = str_detect(toupper(columns), demo_terms),
    has_apoe = str_detect(toupper(columns), apoe_terms),
    has_adas = str_detect(toupper(columns), adas_terms),
    has_mmse = str_detect(toupper(columns), mmse_terms),
    has_wmh = str_detect(toupper(columns), wmh_terms),
    
    clinical_core_score =
      as.integer(has_id) +
      as.integer(has_diagnosis) +
      as.integer(has_demo) +
      as.integer(has_apoe) +
      as.integer(has_adas) +
      as.integer(has_mmse) +
      as.integer(has_wmh)
  ) %>%
  arrange(desc(clinical_core_score), file_name)

# ---- 4. Create separate candidate lists ----
diagnosis_demo_candidates <- candidate_tables %>%
  filter(has_id & (has_diagnosis | has_demo)) %>%
  select(file_name, object_name, n_rows, n_cols, clinical_core_score, columns)

apoe_candidates <- candidate_tables %>%
  filter(has_id & has_apoe) %>%
  select(file_name, object_name, n_rows, n_cols, clinical_core_score, columns)

adas_candidates <- candidate_tables %>%
  filter(has_id & has_adas) %>%
  select(file_name, object_name, n_rows, n_cols, clinical_core_score, columns)

mmse_candidates <- candidate_tables %>%
  filter(has_id & has_mmse) %>%
  select(file_name, object_name, n_rows, n_cols, clinical_core_score, columns)

wmh_candidates <- candidate_tables %>%
  filter(has_id & has_wmh) %>%
  select(file_name, object_name, n_rows, n_cols, clinical_core_score, columns)

# ---- 5. Save output ----
openxlsx::write.xlsx(
  list(
    all_scored_tables = candidate_tables,
    diagnosis_demo_candidates = diagnosis_demo_candidates,
    apoe_candidates = apoe_candidates,
    adas_candidates = adas_candidates,
    mmse_candidates = mmse_candidates,
    wmh_candidates = wmh_candidates
  ),
  file.path(out_dir, "ADNIMERGE2_core_candidate_tables.xlsx"),
  overwrite = TRUE
)

cat("Done. Check this file:\n")
cat(file.path(out_dir, "ADNIMERGE2_core_candidate_tables.xlsx"))