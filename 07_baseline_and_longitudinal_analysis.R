# ==========================================================
# 07_baseline_and_longitudinal_analysis.R
# Purpose: First real ADNI analysis after biomarker merge
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx", "lme4", "lmerTest", "lubridate")

missing <- needed[!needed %in% rownames(installed.packages())]

if (length(missing) > 0) {
  install.packages(missing)
}

library(tidyverse)
library(janitor)
library(openxlsx)
library(lme4)
library(lmerTest)
library(lubridate)

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

safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_q1 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, 0.25))
}

safe_q3 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, 0.75))
}

coef_table <- function(model, model_name) {
  as.data.frame(coef(summary(model))) %>%
    rownames_to_column("term") %>%
    mutate(model = model_name) %>%
    relocate(model, term)
}

# ---------- 1. Load merged biomarker dataset ----------

analysis_file <- file.path(clean_dir, "analysis_master_biomarkers_longitudinal.csv")

dat <- read_csv(analysis_file, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    rid = as.character(rid),
    visit_key = as.character(visit_key),
    dx_label = as.character(dx_label),
    examdate = suppressWarnings(ymd(examdate)),
    ptdobyy = readr::parse_number(as.character(ptdobyy)),
    pteducat = readr::parse_number(as.character(pteducat)),
    apoe4 = readr::parse_number(as.character(apoe4)),
    years_from_baseline = readr::parse_number(as.character(years_from_baseline))
  )

# ---------- 2. Create age if it does not already exist ----------

if (!"age" %in% names(dat)) {
  dat <- dat %>%
    mutate(
      age = if_else(
        !is.na(examdate) & !is.na(ptdobyy),
        as.numeric(year(examdate) - ptdobyy),
        NA_real_
      )
    )
}

# ---------- 3. Define marker variables ----------

primary_markers <- c(
  "gfap_quanterix",
  "strem2_msd_corrected",
  "vegf_plasma_qc",
  "sicam1_plasma_qc",
  "svcam1_plasma_qc"
)

secondary_markers <- c(
  "gfap_fujirebio",
  "strem2_wu_corrected",
  "gfap_csf_mrm",
  "ykl40_csf_mrm_mean",
  "trem2_csfprot_mean",
  "ykl40_csfprot_mean",
  "total_wmh_csv"
)

all_markers <- c(primary_markers, secondary_markers)
all_markers <- all_markers[all_markers %in% names(dat)]

# ---------- 4. Log-transform biomarkers ----------

for (m in all_markers) {
  dat[[paste0("log_", m)]] <- safe_log(dat[[m]])
}

# ---------- 5. Create baseline diagnosis per subject ----------

baseline_dx <- dat %>%
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
  left_join(baseline_dx, by = "rid") %>%
  mutate(
    baseline_dx = factor(baseline_dx, levels = c("CN", "MCI", "AD")),
    dx_label = factor(dx_label, levels = c("CN", "MCI", "AD")),
    ptgender = factor(ptgender)
  )

# ---------- 6. Create baseline dataset ----------

baseline_dat <- dat %>%
  filter(visit_key == "bl") %>%
  filter(dx_label %in% c("CN", "MCI", "AD"))

# ---------- 7. Basic sample summaries ----------

overall_summary <- dat %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(rid),
    n_baseline_subjects = n_distinct(rid[visit_key == "bl"]),
    n_cn = n_distinct(rid[baseline_dx == "CN"]),
    n_mci = n_distinct(rid[baseline_dx == "MCI"]),
    n_ad = n_distinct(rid[baseline_dx == "AD"])
  )

baseline_diagnosis_counts <- baseline_dat %>%
  count(dx_label, name = "n_rows")

marker_availability <- map_dfr(all_markers, function(m) {
  dat %>%
    summarise(
      marker = m,
      non_missing_rows = sum(!is.na(.data[[m]])),
      non_missing_subjects = n_distinct(rid[!is.na(.data[[m]])]),
      baseline_non_missing_rows = sum(visit_key == "bl" & !is.na(.data[[m]])),
      baseline_non_missing_subjects = n_distinct(rid[visit_key == "bl" & !is.na(.data[[m]])])
    )
})

# ---------- 8. Baseline descriptive statistics by diagnosis ----------

baseline_marker_summary <- map_dfr(all_markers, function(m) {
  baseline_dat %>%
    group_by(dx_label) %>%
    summarise(
      marker = m,
      n = sum(!is.na(.data[[m]])),
      mean = mean(.data[[m]], na.rm = TRUE),
      sd = sd(.data[[m]], na.rm = TRUE),
      median = safe_median(.data[[m]]),
      q1 = safe_q1(.data[[m]]),
      q3 = safe_q3(.data[[m]]),
      .groups = "drop"
    )
})

# ---------- 9. Baseline adjusted linear models ----------

baseline_lm_results <- map_dfr(primary_markers, function(m) {
  
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
  
  if (nrow(df) < 50 || n_distinct(df$dx_label) < 2) {
    return(tibble(
      model = paste0("baseline_lm_", m),
      term = "MODEL_NOT_RUN",
      note = "Too few rows or diagnosis groups"
    ))
  }
  
  formula <- as.formula(
    paste0(log_m, " ~ dx_label + age + ptgender + pteducat + apoe4")
  )
  
  model <- lm(formula, data = df)
  
  coef_table(model, paste0("baseline_lm_", m)) %>%
    mutate(
      n_rows = nrow(df),
      n_subjects = n_distinct(df$rid)
    )
})

# ---------- 10. Longitudinal mixed-effects models ----------

longitudinal_lmer_results <- map_dfr(primary_markers, function(m) {
  
  log_m <- paste0("log_", m)
  
  df <- dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
  
  if (nrow(df) < 80 || n_distinct(df$rid) < 30 || n_distinct(df$baseline_dx) < 2) {
    return(tibble(
      model = paste0("lmer_", m),
      term = "MODEL_NOT_RUN",
      note = "Too few rows, subjects, or diagnosis groups"
    ))
  }
  
  formula <- as.formula(
    paste0(log_m, " ~ years_from_baseline * baseline_dx + age + ptgender + pteducat + apoe4 + (1 | rid)")
  )
  
  model <- tryCatch(
    lmer(formula, data = df, REML = FALSE),
    error = function(e) e
  )
  
  if (inherits(model, "error")) {
    return(tibble(
      model = paste0("lmer_", m),
      term = "MODEL_ERROR",
      note = model$message
    ))
  }
  
  coef_table(model, paste0("lmer_", m)) %>%
    mutate(
      n_rows = nrow(df),
      n_subjects = n_distinct(df$rid)
    )
})

# ---------- 11. Correlation matrix for glial-vascular markers ----------

corr_vars <- paste0("log_", primary_markers)
corr_vars <- corr_vars[corr_vars %in% names(dat)]

corr_input <- dat %>%
  filter(visit_key == "bl") %>%
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

# ---------- 12. Save baseline plots ----------

for (m in primary_markers) {
  
  log_m <- paste0("log_", m)
  
  if (!log_m %in% names(baseline_dat)) next
  
  plot_df <- baseline_dat %>%
    filter(!is.na(.data[[log_m]]), dx_label %in% c("CN", "MCI", "AD"))
  
  if (nrow(plot_df) < 20) next
  
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
    filename = file.path(figures_dir, paste0("baseline_", m, "_by_diagnosis.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---------- 13. Save longitudinal plots ----------

for (m in primary_markers) {
  
  log_m <- paste0("log_", m)
  
  if (!log_m %in% names(dat)) next
  
  plot_df <- dat %>%
    filter(
      !is.na(.data[[log_m]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx)
    )
  
  if (nrow(plot_df) < 50) next
  
  p <- ggplot(plot_df, aes(x = years_from_baseline, y = .data[[log_m]])) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(aes(group = baseline_dx), method = "loess", se = TRUE) +
    facet_wrap(~ baseline_dx) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Longitudinal trajectory of", m),
      x = "Years from baseline",
      y = paste("log", m)
    )
  
  ggsave(
    filename = file.path(figures_dir, paste0("longitudinal_", m, "_trajectory.png")),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ---------- 14. Save correlation heatmap ----------

p_corr <- ggplot(cor_long, aes(x = marker_1, y = marker_2, fill = spearman_r)) +
  geom_tile() +
  geom_text(aes(label = round(spearman_r, 2)), size = 3) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Baseline glial-vascular biomarker Spearman correlation matrix",
    x = "",
    y = "",
    fill = "Spearman r"
  )

ggsave(
  filename = file.path(figures_dir, "baseline_glial_vascular_correlation_heatmap.png"),
  plot = p_corr,
  width = 8,
  height = 6,
  dpi = 300
)

# ---------- 15. Save all results ----------

openxlsx::write.xlsx(
  list(
    overall_summary = overall_summary,
    baseline_diagnosis_counts = baseline_diagnosis_counts,
    marker_availability = marker_availability,
    baseline_marker_summary = baseline_marker_summary,
    baseline_lm_results = baseline_lm_results,
    longitudinal_lmer_results = longitudinal_lmer_results,
    correlation_matrix = cor_matrix_df,
    correlation_long = cor_long
  ),
  file.path(results_dir, "first_statistical_analysis_results.xlsx"),
  overwrite = TRUE
)

write_csv(
  dat,
  file.path(clean_dir, "analysis_master_with_log_variables.csv")
)

cat("Done.\n")
cat("Results file:\n")
cat(file.path(results_dir, "first_statistical_analysis_results.xlsx"), "\n\n")
cat("Updated dataset:\n")
cat(file.path(clean_dir, "analysis_master_with_log_variables.csv"), "\n\n")
cat("Figures saved in:\n")
cat(figures_dir, "\n")