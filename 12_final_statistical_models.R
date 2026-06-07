# ==========================================================
# 12_final_statistical_models.R
# Purpose:
# Run final baseline and longitudinal biomarker models.
#
# IMPORTANT REVISION:
# Baseline descriptives and baseline Figure 2 now use the SAME
# biomarker-specific complete-case analytic samples as the adjusted
# baseline regression models. This prevents mismatches such as:
#   descriptive GFAP n != adjusted GFAP model n.
# ==========================================================

needed <- c(
  "tidyverse", "janitor", "openxlsx", "lme4", "lmerTest",
  "broom", "readr"
)

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
library(readr)

project_dir <- "R:/ADNI_Project"

clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")
figures_dir <- file.path(project_dir, "05_figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------

safe_log <- function(x) {
  x_num <- readr::parse_number(as.character(x))
  ifelse(!is.na(x_num) & x_num > 0, log(x_num), NA_real_)
}

coef_table <- function(model, model_name) {
  as.data.frame(coef(summary(model))) %>%
    rownames_to_column("term") %>%
    mutate(model = model_name) %>%
    relocate(model, term)
}

format_q <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = 2),
      formatC(signif(x, 3), format = "fg", digits = 3)
    )
  )
}

clean_model_table <- function(df) {

  if (!"Estimate" %in% names(df)) {
    return(df)
  }

  p_col <- intersect(c("Pr(>|t|)", "p.value", "p_value"), names(df))[1]

  if (is.na(p_col)) {
    stop("No p-value column found in model table.")
  }

  df %>%
    mutate(
      percent_change = (exp(Estimate) - 1) * 100,
      p_value = .data[[p_col]],
      fdr_p = p.adjust(p_value, method = "BH"),
      significance = case_when(
        fdr_p < 0.001 ~ "***",
        fdr_p < 0.01 ~ "**",
        fdr_p < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
}

summarise_filter_step <- function(df, marker_name, step_name) {
  df %>%
    summarise(
      marker = marker_name,
      step = step_name,
      n_rows = n(),
      n_subjects = n_distinct(rid),
      CN = sum(dx_label == "CN", na.rm = TRUE),
      MCI = sum(dx_label == "MCI", na.rm = TRUE),
      AD = sum(dx_label == "AD", na.rm = TRUE),
      .groups = "drop"
    )
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

marker_info <- tibble::tribble(
  ~marker,                 ~marker_label,         ~domain,
  "gfap_quanterix",         "Plasma GFAP",         "glial",
  "strem2_msd_corrected",   "Plasma sTREM2",       "glial",
  "vegf_plasma_qc",         "Plasma VEGF",         "vascular",
  "sicam1_plasma_qc",       "Plasma sICAM-1",      "vascular",
  "svcam1_plasma_qc",       "Plasma sVCAM-1",      "vascular"
) %>%
  filter(marker %in% names(dat))

primary_markers <- marker_info$marker

# Create log variables if missing.
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

# ---------- 4. Biomarker-specific complete-case function ----------

make_baseline_complete_case <- function(marker_name) {

  log_marker <- paste0("log_", marker_name)

  baseline_dat %>%
    filter(
      !is.na(.data[[marker_name]]),
      !is.na(.data[[log_marker]]),
      !is.na(dx_label),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
}

make_longitudinal_complete_case <- function(marker_name) {

  log_marker <- paste0("log_", marker_name)

  dat %>%
    filter(
      !is.na(.data[[marker_name]]),
      !is.na(.data[[log_marker]]),
      !is.na(years_from_baseline),
      !is.na(baseline_dx),
      !is.na(age),
      !is.na(ptgender),
      !is.na(pteducat),
      !is.na(apoe4)
    )
}

# ---------- 5. Overall and marker availability ----------

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

  baseline_cc <- make_baseline_complete_case(m)
  long_cc <- make_longitudinal_complete_case(m)

  tibble(
    marker = m,
    total_non_missing_rows = sum(!is.na(dat[[m]])),
    total_non_missing_subjects = n_distinct(dat$rid[!is.na(dat[[m]])]),
    baseline_non_missing_rows = sum(!is.na(baseline_dat[[m]])),
    baseline_non_missing_subjects = n_distinct(baseline_dat$rid[!is.na(baseline_dat[[m]])]),
    baseline_complete_case_rows = nrow(baseline_cc),
    baseline_complete_case_subjects = n_distinct(baseline_cc$rid),
    longitudinal_complete_case_rows = nrow(long_cc),
    longitudinal_complete_case_subjects = n_distinct(long_cc$rid)
  )
}) %>%
  left_join(marker_info, by = "marker")

# ---------- 6. Baseline complete-case filter audit ----------

baseline_filter_audit <- map_dfr(primary_markers, function(m) {

  log_m <- paste0("log_", m)

  step0 <- baseline_dat
  step1 <- step0 %>% filter(!is.na(.data[[m]]))
  step2 <- step1 %>% filter(!is.na(.data[[log_m]]))
  step3 <- step2 %>% filter(!is.na(dx_label))
  step4 <- step3 %>% filter(!is.na(age))
  step5 <- step4 %>% filter(!is.na(ptgender))
  step6 <- step5 %>% filter(!is.na(pteducat))
  step7 <- step6 %>% filter(!is.na(apoe4))

  bind_rows(
    summarise_filter_step(step0, m, "0_baseline_CN_MCI_AD"),
    summarise_filter_step(step1, m, "1_nonmissing_raw_marker"),
    summarise_filter_step(step2, m, "2_positive_log_marker"),
    summarise_filter_step(step3, m, "3_nonmissing_diagnosis"),
    summarise_filter_step(step4, m, "4_nonmissing_age"),
    summarise_filter_step(step5, m, "5_nonmissing_sex"),
    summarise_filter_step(step6, m, "6_nonmissing_education"),
    summarise_filter_step(step7, m, "7_final_complete_case")
  )
}) %>%
  left_join(marker_info, by = "marker")

gfap_filter_audit <- baseline_filter_audit %>%
  filter(marker == "gfap_quanterix")

# ---------- 7. Corrected baseline descriptives ----------
# These are complete-case descriptives. They should be used for Table 1.

baseline_descriptives <- map_dfr(primary_markers, function(m) {

  df <- make_baseline_complete_case(m)

  df %>%
    group_by(dx_label) %>%
    summarise(
      marker = m,
      n = n(),
      mean = mean(.data[[m]], na.rm = TRUE),
      sd = sd(.data[[m]], na.rm = TRUE),
      median = median(.data[[m]], na.rm = TRUE),
      q1 = as.numeric(quantile(.data[[m]], 0.25, na.rm = TRUE)),
      q3 = as.numeric(quantile(.data[[m]], 0.75, na.rm = TRUE)),
      median_iqr = paste0(
        round(median, 2),
        " [",
        round(q1, 2),
        "-",
        round(q3, 2),
        "]"
      ),
      sample_definition = "baseline biomarker-specific complete case",
      .groups = "drop"
    )
}) %>%
  left_join(marker_info, by = "marker") %>%
  select(
    marker, marker_label, domain, dx_label, n,
    mean, sd, median, q1, q3, median_iqr, sample_definition
  )

baseline_log_descriptives <- map_dfr(primary_markers, function(m) {

  log_m <- paste0("log_", m)
  df <- make_baseline_complete_case(m)

  df %>%
    group_by(dx_label) %>%
    summarise(
      marker = m,
      n = n(),
      mean_log = mean(.data[[log_m]], na.rm = TRUE),
      sd_log = sd(.data[[log_m]], na.rm = TRUE),
      median_log = median(.data[[log_m]], na.rm = TRUE),
      q1_log = as.numeric(quantile(.data[[log_m]], 0.25, na.rm = TRUE)),
      q3_log = as.numeric(quantile(.data[[log_m]], 0.75, na.rm = TRUE)),
      .groups = "drop"
    )
}) %>%
  left_join(marker_info, by = "marker")

# ---------- 8. Final baseline adjusted linear models ----------

baseline_lm_results <- map_dfr(primary_markers, function(m) {

  log_m <- paste0("log_", m)
  model_df <- make_baseline_complete_case(m)

  if (nrow(model_df) < 50 || n_distinct(model_df$dx_label) < 2) {
    return(tibble(
      marker = m,
      model = paste0("baseline_lm_", m),
      term = "MODEL_NOT_RUN",
      note = "Too few rows or diagnosis groups",
      n_rows = nrow(model_df),
      n_subjects = n_distinct(model_df$rid)
    ))
  }

  formula <- as.formula(
    paste0(log_m, " ~ dx_label + age + ptgender + pteducat + apoe4")
  )

  model <- lm(formula, data = model_df)

  coef_table(model, paste0("baseline_lm_", m)) %>%
    mutate(
      marker = m,
      n_rows = nrow(model_df),
      n_subjects = n_distinct(model_df$rid)
    )
})

baseline_lm_clean <- clean_model_table(baseline_lm_results) %>%
  left_join(marker_info, by = "marker")

baseline_dx_effects <- baseline_lm_clean %>%
  filter(term %in% c("dx_labelMCI", "dx_labelAD")) %>%
  mutate(
    contrast = case_when(
      term == "dx_labelMCI" ~ "MCI vs CN",
      term == "dx_labelAD" ~ "AD vs CN",
      TRUE ~ term
    )
  ) %>%
  select(
    marker, marker_label, contrast, Estimate, `Std. Error`,
    percent_change, p_value, fdr_p, significance, n_rows, n_subjects
  )

# ---------- 9. Final longitudinal mixed-effects models ----------

longitudinal_lmer_results <- map_dfr(primary_markers, function(m) {

  log_m <- paste0("log_", m)
  model_df <- make_longitudinal_complete_case(m)

  if (nrow(model_df) < 80 ||
      n_distinct(model_df$rid) < 30 ||
      n_distinct(model_df$baseline_dx) < 2) {
    return(tibble(
      marker = m,
      model = paste0("lmer_", m),
      term = "MODEL_NOT_RUN",
      note = "Too few rows, subjects, or diagnosis groups",
      n_rows = nrow(model_df),
      n_subjects = n_distinct(model_df$rid)
    ))
  }

  formula <- as.formula(
    paste0(
      log_m,
      " ~ years_from_baseline * baseline_dx + age + ptgender + pteducat + apoe4 + (1 | rid)"
    )
  )

  model <- tryCatch(
    lmer(formula, data = model_df, REML = FALSE),
    error = function(e) e
  )

  if (inherits(model, "error")) {
    return(tibble(
      marker = m,
      model = paste0("lmer_", m),
      term = "MODEL_ERROR",
      note = model$message,
      n_rows = nrow(model_df),
      n_subjects = n_distinct(model_df$rid)
    ))
  }

  coef_table(model, paste0("lmer_", m)) %>%
    mutate(
      marker = m,
      n_rows = nrow(model_df),
      n_subjects = n_distinct(model_df$rid)
    )
})

longitudinal_lmer_clean <- clean_model_table(longitudinal_lmer_results) %>%
  left_join(marker_info, by = "marker")

longitudinal_main_effects <- longitudinal_lmer_clean %>%
  filter(term %in% c(
    "years_from_baseline",
    "baseline_dxMCI",
    "baseline_dxAD",
    "years_from_baseline:baseline_dxMCI",
    "years_from_baseline:baseline_dxAD"
  ))

# ---------- 10. Baseline glial-vascular correlations ----------
# Pairwise complete observations are used for correlations, as described
# in the manuscript.

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

# Pairwise n for each correlation cell.
cor_pairwise_n <- expand_grid(
  marker_1 = corr_vars,
  marker_2 = corr_vars
) %>%
  mutate(
    pairwise_n = map2_int(marker_1, marker_2, function(a, b) {
      sum(!is.na(corr_input[[a]]) & !is.na(corr_input[[b]]))
    })
  )

cor_long <- cor_long %>%
  left_join(cor_pairwise_n, by = c("marker_1", "marker_2"))

# ---------- 11. Corrected Figure 2: complete-case baseline distributions ----------

figure2_plot_df <- map_dfr(primary_markers, function(m) {

  log_m <- paste0("log_", m)

  make_baseline_complete_case(m) %>%
    transmute(
      rid,
      dx_label,
      marker = m,
      log_value = .data[[log_m]]
    )
}) %>%
  left_join(marker_info, by = "marker") %>%
  mutate(
    marker_label = factor(
      marker_label,
      levels = c(
        "Plasma GFAP",
        "Plasma sICAM-1",
        "Plasma sTREM2",
        "Plasma sVCAM-1",
        "Plasma VEGF"
      )
    ),
    dx_label = factor(dx_label, levels = c("CN", "MCI", "AD"))
  )

panel_ranges <- figure2_plot_df %>%
  group_by(marker_label) %>%
  summarise(
    y_min = min(log_value, na.rm = TRUE),
    y_max = max(log_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_n = y_min + 0.10 * (y_max - y_min),
    y_q = y_max - 0.08 * (y_max - y_min)
  )

n_labels <- figure2_plot_df %>%
  count(marker_label, dx_label, name = "n") %>%
  left_join(panel_ranges, by = "marker_label") %>%
  mutate(label = paste0("n=", n))

q_labels <- baseline_dx_effects %>%
  mutate(
    marker_label = factor(
      marker_label,
      levels = c(
        "Plasma GFAP",
        "Plasma sICAM-1",
        "Plasma sTREM2",
        "Plasma sVCAM-1",
        "Plasma VEGF"
      )
    ),
    contrast_short = case_when(
      contrast == "MCI vs CN" ~ "mci_q",
      contrast == "AD vs CN" ~ "ad_q",
      TRUE ~ NA_character_
    )
  ) %>%
  select(marker_label, contrast_short, fdr_p) %>%
  filter(!is.na(contrast_short)) %>%
  pivot_wider(names_from = contrast_short, values_from = fdr_p) %>%
  mutate(
    dx_label = factor("MCI", levels = c("CN", "MCI", "AD")),
    label = paste0(
      "MCI vs CN: q = ", format_q(mci_q),
      "\nAD vs CN: q = ", format_q(ad_q)
    )
  ) %>%
  left_join(panel_ranges, by = "marker_label")

p_figure2 <- ggplot(figure2_plot_df, aes(x = dx_label, y = log_value)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 0.7) +
  geom_text(
    data = n_labels,
    aes(x = dx_label, y = y_n, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_text(
    data = q_labels,
    aes(x = dx_label, y = y_q, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_wrap(~ marker_label, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Baseline biomarker distributions across diagnostic groups",
    x = "Diagnosis",
    y = "log-transformed biomarker level"
  )

ggsave(
  filename = file.path(figures_dir, "Figure2_baseline_biomarker_distributions_complete_case.png"),
  plot = p_figure2,
  width = 14,
  height = 8,
  dpi = 300
)

# ---------- 12. Optional individual baseline plots ----------

for (m in primary_markers) {

  log_m <- paste0("log_", m)

  plot_df <- make_baseline_complete_case(m)

  if (nrow(plot_df) < 20) next

  p <- ggplot(plot_df, aes(x = dx_label, y = .data[[log_m]])) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Baseline", m, "by diagnosis, complete-case sample"),
      x = "Diagnosis",
      y = paste("log", m)
    )

  ggsave(
    filename = file.path(figures_dir, paste0("final_baseline_", m, "_by_diagnosis_complete_case.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---------- 13. Optional longitudinal exploratory plots ----------

for (m in primary_markers) {

  log_m <- paste0("log_", m)

  plot_df <- make_longitudinal_complete_case(m)

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
    filename = file.path(figures_dir, paste0("final_longitudinal_", m, "_trajectory.png")),
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

# ---------- 15. Save all outputs ----------
# Sheet names are kept <=31 characters for Excel compatibility.

openxlsx::write.xlsx(
  list(
    overall_summary = overall_summary,
    marker_availability = marker_availability,
    baseline_filter_audit = baseline_filter_audit,
    gfap_filter_audit = gfap_filter_audit,
    baseline_descriptives = baseline_descriptives,
    baseline_log_desc = baseline_log_descriptives,
    baseline_lm_results = baseline_lm_results,
    baseline_lm_clean = baseline_lm_clean,
    baseline_dx_effects = baseline_dx_effects,
    longitudinal_lmer_results = longitudinal_lmer_results,
    longitudinal_lmer_clean = longitudinal_lmer_clean,
    longitudinal_main_fx = longitudinal_main_effects,
    correlation_matrix = cor_matrix_df,
    correlation_long = cor_long
  ),
  file.path(results_dir, "final_statistical_models_results.xlsx"),
  overwrite = TRUE
)

cat("Done.\n")
cat("Final results file:\n")
cat(file.path(results_dir, "final_statistical_models_results.xlsx"), "\n\n")
cat("Corrected Figure 2:\n")
cat(file.path(figures_dir, "Figure2_baseline_biomarker_distributions_complete_case.png"), "\n\n")
cat("Figures saved in:\n")
cat(figures_dir, "\n")
