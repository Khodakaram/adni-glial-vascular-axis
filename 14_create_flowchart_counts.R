# ==========================================================
# 14_create_flowchart_counts.R
# Purpose: Generate participant-selection counts for manuscript flowchart
# ==========================================================

library(tidyverse)
library(janitor)
library(openxlsx)
library(readr)

project_dir <- "R:/ADNI_Project"

clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. Load model-ready dataset ----------

dat <- read_csv(
  file.path(clean_dir, "analysis_master_model_ready.csv"),
  show_col_types = FALSE,
  guess_max = 100000
) %>%
  clean_names() %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key),
    dx_label = as.character(dx_label),
    baseline_dx = as.character(baseline_dx),
    age = readr::parse_number(as.character(age)),
    pteducat = readr::parse_number(as.character(pteducat)),
    apoe4 = readr::parse_number(as.character(apoe4)),
    years_from_baseline = readr::parse_number(as.character(years_from_baseline))
  )

# ---------- 2. Define baseline dataset ----------

baseline_dat <- dat %>%
  filter(visit_key == "bl")

# ---------- 3. General cohort flow ----------

n_all <- n_distinct(dat$rid)

n_baseline <- baseline_dat %>%
  summarise(n = n_distinct(rid)) %>%
  pull(n)

baseline_dx_complete <- baseline_dat %>%
  filter(dx_label %in% c("CN", "MCI", "AD"))

n_baseline_dx <- n_distinct(baseline_dx_complete$rid)

baseline_covariate_complete <- baseline_dx_complete %>%
  filter(
    !is.na(age),
    !is.na(ptgender),
    !is.na(pteducat),
    !is.na(apoe4)
  )

n_baseline_covariate <- n_distinct(baseline_covariate_complete$rid)

baseline_cognition_complete <- baseline_covariate_complete %>%
  filter(
    !is.na(adas13),
    !is.na(mmse)
  )

n_baseline_cognition <- n_distinct(baseline_cognition_complete$rid)

general_flow <- tibble(
  step = c(
    "ADNI participants in clinical core",
    "Participants with baseline visit",
    "Participants with baseline CN/MCI/AD diagnosis",
    "Participants with complete baseline covariates",
    "Participants with complete baseline cognitive data"
  ),
  n_participants = c(
    n_all,
    n_baseline,
    n_baseline_dx,
    n_baseline_covariate,
    n_baseline_cognition
  ),
  excluded_from_previous_step = c(
    NA,
    n_all - n_baseline,
    n_baseline - n_baseline_dx,
    n_baseline_dx - n_baseline_covariate,
    n_baseline_covariate - n_baseline_cognition
  )
)

# ---------- 4. Missingness reasons at baseline ----------

baseline_missing_reasons <- tibble(
  variable = c(
    "Diagnosis",
    "Age",
    "Sex",
    "Education",
    "APOE4",
    "ADAS13",
    "MMSE",
    "GFAP",
    "sTREM2",
    "VEGF",
    "sICAM1",
    "sVCAM1"
  ),
  missing_n_participants = c(
    n_distinct(baseline_dat$rid[!baseline_dat$dx_label %in% c("CN", "MCI", "AD") | is.na(baseline_dat$dx_label)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$age)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$ptgender)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$pteducat)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$apoe4)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$adas13)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$mmse)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$gfap_quanterix)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$strem2_msd_corrected)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$vegf_plasma_qc)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$sicam1_plasma_qc)]),
    n_distinct(baseline_dat$rid[is.na(baseline_dat$svcam1_plasma_qc)])
  )
)

# ---------- 5. Baseline biomarker-specific analysis cohorts ----------

baseline_biomarker_cohorts <- tibble(
  analysis_group = c(
    "GFAP baseline model",
    "sTREM2 baseline model",
    "VEGF baseline model",
    "sICAM1 baseline model",
    "sVCAM1 baseline model",
    "Complete vascular panel baseline",
    "Complete all five primary biomarkers baseline"
  ),
  n_participants = c(
    baseline_covariate_complete %>%
      filter(!is.na(gfap_quanterix)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(!is.na(strem2_msd_corrected)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(!is.na(vegf_plasma_qc)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(!is.na(sicam1_plasma_qc)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(!is.na(svcam1_plasma_qc)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(
        !is.na(vegf_plasma_qc),
        !is.na(sicam1_plasma_qc),
        !is.na(svcam1_plasma_qc)
      ) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    
    baseline_covariate_complete %>%
      filter(
        !is.na(gfap_quanterix),
        !is.na(strem2_msd_corrected),
        !is.na(vegf_plasma_qc),
        !is.na(sicam1_plasma_qc),
        !is.na(svcam1_plasma_qc)
      ) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n)
  )
)

# ---------- 6. Longitudinal biomarker-specific analysis cohorts ----------

longitudinal_covariate_complete <- dat %>%
  filter(
    baseline_dx %in% c("CN", "MCI", "AD"),
    !is.na(years_from_baseline),
    !is.na(age),
    !is.na(ptgender),
    !is.na(pteducat),
    !is.na(apoe4)
  )

longitudinal_biomarker_cohorts <- tibble(
  analysis_group = c(
    "GFAP longitudinal model",
    "sTREM2 longitudinal model",
    "VEGF longitudinal model",
    "sICAM1 longitudinal model",
    "sVCAM1 longitudinal model",
    "Complete vascular panel longitudinal",
    "Complete all five primary biomarkers longitudinal"
  ),
  n_rows = c(
    longitudinal_covariate_complete %>% filter(!is.na(gfap_quanterix)) %>% nrow(),
    longitudinal_covariate_complete %>% filter(!is.na(strem2_msd_corrected)) %>% nrow(),
    longitudinal_covariate_complete %>% filter(!is.na(vegf_plasma_qc)) %>% nrow(),
    longitudinal_covariate_complete %>% filter(!is.na(sicam1_plasma_qc)) %>% nrow(),
    longitudinal_covariate_complete %>% filter(!is.na(svcam1_plasma_qc)) %>% nrow(),
    longitudinal_covariate_complete %>%
      filter(!is.na(vegf_plasma_qc), !is.na(sicam1_plasma_qc), !is.na(svcam1_plasma_qc)) %>%
      nrow(),
    longitudinal_covariate_complete %>%
      filter(
        !is.na(gfap_quanterix),
        !is.na(strem2_msd_corrected),
        !is.na(vegf_plasma_qc),
        !is.na(sicam1_plasma_qc),
        !is.na(svcam1_plasma_qc)
      ) %>%
      nrow()
  ),
  n_participants = c(
    longitudinal_covariate_complete %>% filter(!is.na(gfap_quanterix)) %>% summarise(n = n_distinct(rid)) %>% pull(n),
    longitudinal_covariate_complete %>% filter(!is.na(strem2_msd_corrected)) %>% summarise(n = n_distinct(rid)) %>% pull(n),
    longitudinal_covariate_complete %>% filter(!is.na(vegf_plasma_qc)) %>% summarise(n = n_distinct(rid)) %>% pull(n),
    longitudinal_covariate_complete %>% filter(!is.na(sicam1_plasma_qc)) %>% summarise(n = n_distinct(rid)) %>% pull(n),
    longitudinal_covariate_complete %>% filter(!is.na(svcam1_plasma_qc)) %>% summarise(n = n_distinct(rid)) %>% pull(n),
    longitudinal_covariate_complete %>%
      filter(!is.na(vegf_plasma_qc), !is.na(sicam1_plasma_qc), !is.na(svcam1_plasma_qc)) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n),
    longitudinal_covariate_complete %>%
      filter(
        !is.na(gfap_quanterix),
        !is.na(strem2_msd_corrected),
        !is.na(vegf_plasma_qc),
        !is.na(sicam1_plasma_qc),
        !is.na(svcam1_plasma_qc)
      ) %>%
      summarise(n = n_distinct(rid)) %>%
      pull(n)
  )
)

# ---------- 7. Diagnosis counts for baseline flow ----------

baseline_dx_counts <- baseline_covariate_complete %>%
  count(dx_label, name = "n_participants")

# ---------- 8. Save outputs ----------

openxlsx::write.xlsx(
  list(
    general_flow = general_flow,
    baseline_missing_reasons = baseline_missing_reasons,
    baseline_biomarker_cohorts = baseline_biomarker_cohorts,
    longitudinal_biomarker_cohorts = longitudinal_biomarker_cohorts,
    baseline_dx_counts = baseline_dx_counts
  ),
  file.path(results_dir, "flowchart_counts.xlsx"),
  overwrite = TRUE
)

cat("Done. Check:\n")
cat(file.path(results_dir, "flowchart_counts.xlsx"), "\n")