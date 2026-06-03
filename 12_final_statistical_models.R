# ==========================================================
# 12_final_statistical_models.R
# Purpose: Run final baseline and longitudinal biomarker models
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest", "broom")

missing <- needed[!needed %in% rownames(installed.packages())]

if (length(missing) > 0) {
  install.packages(missing)
}

library(tidyverse)
library(janitor)
library(openxlsx)
library(lme4)
library(lmerTest)
library(broom)

project_dir <- "R:/ADNI_Project"

clean_dir <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")
figures_dir <- file.path(project_dir, "05_figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------

safe_log <- function(x) {
  x <- readr::parse_number(as.character(x))
  ifelse(!is.na(x) & x > 0, log(x), NA_real_)
}

coef_table <- function(model, model_name) {
  as.data.frame(coef(summary(model))) %>%
    rownames_to_column("term") %>%
    mutate(model = model_name) %>%
    relocate(model, term)
}

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
    dx_label = factor(dx_label, levels = c("CN", "MCI", "AD")),
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender),
    age = readr::parse_number(as.character(age)),
    pteducat = readr::parse_number(as.character(pteducat)),
    apoe4 = readr::parse_number(as.character(apoe4)),
    years_from_baseline = readr::parse_number(as.character(years_from_baseline))
  )

# ---------- 2. Define main biomarkers ----------

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

# ---------- 3. Baseline dataset ----------

baseline_dat <- dat %>%
  filter(visit_key == "bl") %>%
  filter(dx_label %in% c("CN", "MCI", "AD"))

# ---------- 4. Overall and marker availability ----------

overall_summary <- dat %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    baseline_subjects = n_distinct(rid[visit_key == "bl"]),
    non_missing_age = sum(!is.na(age)),
    non_missing_apoe4 = sum(!is.na(apoe4)),
    non_missing_pteducat = sum(!is.na(pteducat))
  )

marker_availability <- map_dfr(primary_markers, function(m) {
  tibble(
    marker = m,
    total_non_missing_rows = sum(!is.na(dat[[m]])),
    total_non_missing_subjects = n_distinct(dat$rid[!is.na(dat[[m]])]),
    baseline_non_missing_rows = sum(!is.na(baseline_dat[[m]])),
    baseline_non_missing_subjects = n_distinct(baseline_dat$rid[!is.na(baseline_dat[[m]])])
  )
})

# ---------- 5. Baseline descriptive statistics ----------

baseline_descriptives <- map_dfr(primary_markers, function(m) {
  baseline_dat %>%
    group_by(dx_label) %>%
    summarise(
      marker = m,
      n = sum(!is.na(.data[[m]])),
      mean = mean(.data[[m]], na.rm = TRUE),
      sd = sd(.data[[m]], na.rm = TRUE),
      median = median(.data[[m]], na.rm = TRUE),
      q1 = as.numeric(quantile(.data[[m]], 0.25, na.rm = TRUE)),
      q3 = as.numeric(quantile(.data[[m]], 0.75, na.rm = TRUE)),
      .groups = "drop"
    )
})


# ---------- 6. Final baseline adjusted linear models ----------

baseline_lm_results <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  model_df <- baseline_dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(dx_label),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
  
  n_rows_model <- nrow(model_df)
  n_subjects_model <- n_distinct(model_df$rid)
  
  formula <- as.formula(
    paste0(log_m, " ~ dx_label + age + ptgender + pteducat + apoe4")
  )
  
  model <- lm(formula, data = model_df)
  
  coef_table(model, paste0("baseline_lm_", m)) %>%
    mutate(
      marker = m,
      n_rows = n_rows_model,
      n_subjects = n_subjects_model
    )
})

# ---------- 7. Final longitudinal mixed-effects models ----------

longitudinal_lmer_results <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  model_df <- dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
  
  n_rows_model <- nrow(model_df)
  n_subjects_model <- n_distinct(model_df$rid)
  
  formula <- as.formula(
    paste0(
      log_m,
      " ~ years_from_baseline * baseline_dx + age + ptgender + pteducat + apoe4 + (1 | rid)"
    )
  )
  
  model <- lmer(formula, data = model_df, REML = FALSE)
  
  coef_table(model, paste0("lmer_", m)) %>%
    mutate(
      marker = m,
      n_rows = n_rows_model,
      n_subjects = n_subjects_model
    )
})

# ---------- 8. Baseline glial-vascular correlations ----------

corr_vars <- paste0("log_", primary_markers)
corr_vars <- corr_vars[corr_vars %in% names(dat)]

corr_input <- baseline_dat %>%
  select(all_of(corr_vars))

cor_matrix <- cor(
  corr_input,
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor_matrix_df <- as.data.frame(cor_matrix) %>%
  rownames_to_column("marker")

cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(marker_1 = Var1, marker_2 = Var2, spearman_r = Freq)

# ---------- 9. Save baseline plots ----------

for (m in primary_markers) {
  
  log_m <- paste0("log_", m)
  
  plot_df <- baseline_dat %>%
    filter(!is.na(.data[[log_m]]), dx_label %in% c("CN", "MCI", "AD"))
  
  p <- ggplot(plot_df, aes(x = dx_label, y = .data[[log_m]])) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Baseline", m, "by diagnosis"),
      x = "Diagnosis",
      y = paste("log", m)
    )
  
  ggsave(
    filename = file.path(figures_dir, paste0("final_baseline_", m, "_by_diagnosis.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---------- 10. Save longitudinal plots ----------

for (m in primary_markers) {
  
  log_m <- paste0("log_", m)
  
  plot_df <- dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx)
    )
  
  p <- ggplot(plot_df, aes(x = years_from_baseline, y = .data[[log_m]])) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(method = "loess", se = TRUE) +
    facet_wrap(~ baseline_dx) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Longitudinal trajectory of", m),
      x = "Years from baseline",
      y = paste("log", m)
    )
  
  ggsave(
    filename = file.path(figures_dir, paste0("final_longitudinal_", m, "_trajectory.png")),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ---------- 11. Save correlation heatmap ----------

p_corr <- ggplot(cor_long, aes(x = marker_1, y = marker_2, fill = spearman_r)) +
  geom_tile() +
  geom_text(aes(label = round(spearman_r, 2)), size = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Baseline glial-vascular biomarker Spearman correlation matrix",
    x = "",
    y = "",
    fill = "Spearman r"
  )

ggsave(
  filename = file.path(figures_dir, "final_baseline_glial_vascular_correlation_heatmap.png"),
  plot = p_corr,
  width = 8,
  height = 6,
  dpi = 300
)

# ---------- 12. Save all outputs ----------

openxlsx::write.xlsx(
  list(
    overall_summary = overall_summary,
    marker_availability = marker_availability,
    baseline_descriptives = baseline_descriptives,
    baseline_lm_results = baseline_lm_results,
    longitudinal_lmer_results = longitudinal_lmer_results,
    correlation_matrix = cor_matrix_df,
    correlation_long = cor_long
  ),
  file.path(results_dir, "final_statistical_models_results.xlsx"),
  overwrite = TRUE
)

cat("Done.\n")
cat("Final results file:\n")
cat(file.path(results_dir, "final_statistical_models_results.xlsx"), "\n\n")
cat("Figures saved in:\n")
cat(figures_dir, "\n")