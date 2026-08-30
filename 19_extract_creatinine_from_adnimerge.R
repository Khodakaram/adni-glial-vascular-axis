# ==========================================================
# 19_extract_creatinine_from_adnimerge.R
#
# ADNI's laboratory table is wide with coded column names such as
# RCT392, and the readable descriptions live in the data
# dictionary rather than in the column names. Script 17 searched
# column names, so it found nothing. This script resolves the
# codes through four routes and exports baseline creatinine.
#
#   Route 1. datadic, the ADNI data dictionary table
#   Route 2. Hmisc-style label attributes on the columns
#   Route 3. the package .Rd documentation
#   Route 4. a distribution scan, as a last resort
#
# Output: 02_clean_data/creatinine_baseline.csv
#         04_results/creatinine_extraction.xlsx  (audit trail)
#
# Run this before 18_risk_adjustment_repaired.R.
# ==========================================================

needed <- c("tidyverse", "janitor", "openxlsx")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(tidyverse); library(janitor); library(openxlsx)

project_dir <- "R:/ADNI_Project"
raw_dir     <- file.path(project_dir, "00_raw_data")
clean_dir   <- file.path(project_dir, "02_clean_data")
results_dir <- file.path(project_dir, "04_results")

out <- list()
parse_num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

# ==========================================================
# 1. Find the ADNIMERGE data, wherever it lives
# ==========================================================

cat("=== Locating ADNIMERGE data ===\n")

data_dirs <- character(0)

for (pkg in c("ADNIMERGE", "ADNIMERGE2")) {
  p <- tryCatch(find.package(pkg), error = function(e) NULL)
  if (!is.null(p)) {
    cat("Installed package found:", pkg, "at", p, "\n")
    data_dirs <- c(data_dirs, file.path(p, "data"))
  }
}

# Also any unpacked copies sitting in the project tree
loose <- list.dirs(raw_dir, recursive = TRUE, full.names = TRUE)
loose <- loose[grepl("ADNIMERGE.*[/\\\\]data$", loose, ignore.case = TRUE)]
data_dirs <- unique(c(data_dirs, loose))

if (!length(data_dirs)) stop("No ADNIMERGE data directory found.")
cat("Directories to search:\n"); print(data_dirs)

# Load an object by name from any of those directories, or from
# the installed package's lazy-load database.
load_object <- function(stem) {
  for (d in data_dirs) {
    f <- list.files(d, pattern = paste0("^", stem, "\\.rda$"),
                    full.names = TRUE, ignore.case = TRUE)
    if (length(f)) {
      env <- new.env(); nms <- load(f[1], envir = env)
      cat("  loaded", stem, "from", basename(dirname(dirname(f[1]))), "\n")
      return(get(nms[1], envir = env))
    }
  }
  for (pkg in c("ADNIMERGE", "ADNIMERGE2")) {
    obj <- tryCatch({
      e <- new.env(); utils::data(list = stem, package = pkg, envir = e); get(stem, envir = e)
    }, error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(obj)) { cat("  loaded", stem, "from package", pkg, "\n"); return(obj) }
  }
  NULL
}

# Inventory, so we can see what is actually there
all_rda <- unique(basename(unlist(lapply(data_dirs, list.files, pattern = "\\.rda$"))))
out$available_tables <- tibble(table = sort(sub("\\.rda$", "", all_rda)))
cat("Tables available:", nrow(out$available_tables), "\n")

labs <- NULL
for (stem in c("labdata", "LABDATA", "Labdata")) {
  labs <- load_object(stem)
  if (!is.null(labs)) break
}
if (is.null(labs)) stop("Could not load labdata. Check the available_tables sheet.")

labs <- as_tibble(labs)
cat("labdata:", nrow(labs), "rows,", ncol(labs), "columns\n")


# ==========================================================
# 2. Build a code-to-description map from every available route
# ==========================================================

cat("\n=== Resolving column codes to descriptions ===\n")

label_map <- tibble(column = character(), description = character(), route = character())

# Route 1: the data dictionary
datadic <- load_object("datadic")
if (!is.null(datadic)) {
  dd <- as_tibble(datadic) %>% clean_names()
  fld <- intersect(c("fldname", "field_name", "variable"), names(dd))[1]
  txt <- intersect(c("text", "description", "label"), names(dd))[1]
  tbl <- intersect(c("tblname", "crfname", "table_name", "phase"), names(dd))[1]
  if (!is.na(fld) && !is.na(txt)) {
    d1 <- dd %>%
      { if (!is.na(tbl)) filter(., grepl("lab", as.character(.data[[tbl]]), ignore.case = TRUE)) |>
          bind_rows(dd) else . } %>%
      transmute(column = as.character(.data[[fld]]),
                description = as.character(.data[[txt]]), route = "datadic") %>%
      distinct(column, .keep_all = TRUE)
    label_map <- bind_rows(label_map, d1)
    cat("  datadic contributed", nrow(d1), "descriptions\n")
  }
}

# Route 2: label attributes carried on the columns themselves
attr_lab <- map_dfr(names(labs), function(cc) {
  l <- attr(labs[[cc]], "label")
  if (is.null(l) || !nzchar(as.character(l)[1])) return(NULL)
  tibble(column = cc, description = as.character(l)[1], route = "column label")
})
label_map <- bind_rows(label_map, attr_lab)
cat("  column attributes contributed", nrow(attr_lab), "descriptions\n")

# Route 3: the package documentation
rd_hits <- tryCatch({
  db <- tools::Rd_db("ADNIMERGE")
  key <- grep("labdata", names(db), ignore.case = TRUE, value = TRUE)
  if (!length(key)) NULL else {
    txt <- paste(capture.output(db[[key[1]]]), collapse = "\n")
    lines <- unlist(strsplit(txt, "\n"))
    hits <- grep("creatinine", lines, ignore.case = TRUE, value = TRUE)
    if (length(hits)) tibble(doc_line = hits) else NULL
  }
}, error = function(e) NULL)
out$rd_creatinine_lines <- if (is.null(rd_hits))
  tibble(note = "no creatinine lines found in package documentation") else rd_hits
if (!is.null(rd_hits)) { cat("  documentation mentions creatinine:\n"); print(as.data.frame(rd_hits)) }

label_map <- label_map %>%
  filter(column %in% names(labs)) %>%
  distinct(column, route, .keep_all = TRUE)
out$label_map <- label_map


# ==========================================================
# 3. Identify creatinine
# ==========================================================

cat("\n=== Creatinine candidates ===\n")

# mg/dL sits near 0.9; umol/L sits near 80
classify_units <- function(v) {
  v <- v[!is.na(v) & v > 0]
  if (!length(v)) return(NA_character_)
  m <- median(v)
  if (m > 0.3 && m < 3)   return("mg/dL")
  if (m > 30 && m < 300)  return("umol/L")
  NA_character_
}

described_hits <- label_map %>%
  filter(grepl("creatinin", description, ignore.case = TRUE))

candidates <- bind_rows(
  described_hits %>% mutate(source = "description match"),
  tibble(column = grep("creat|rct392", names(labs), value = TRUE, ignore.case = TRUE),
         description = NA_character_, route = "name match", source = "name match")
) %>%
  distinct(column, .keep_all = TRUE) %>%
  rowwise() %>%
  mutate(
    n_values    = sum(!is.na(parse_num(labs[[column]])) & parse_num(labs[[column]]) > 0),
    n_distinct  = dplyr::n_distinct(parse_num(labs[[column]])),
    median_value = suppressWarnings(median(parse_num(labs[[column]])[
      !is.na(parse_num(labs[[column]])) & parse_num(labs[[column]]) > 0])),
    inferred_units = classify_units(parse_num(labs[[column]]))
  ) %>%
  ungroup() %>%
  mutate(usable = n_values > 100 & n_distinct > 20 & !is.na(inferred_units)) %>%
  arrange(desc(usable), desc(n_values))

# Route 4: distribution scan, only if nothing above worked
if (!any(candidates$usable)) {
  cat("  No labelled match; falling back to a distribution scan.\n")
  scan <- map_dfr(names(labs), function(cc) {
    v <- parse_num(labs[[cc]]); v <- v[!is.na(v) & v > 0]
    if (length(v) < 200 || dplyr::n_distinct(v) < 25) return(NULL)
    u <- classify_units(v)
    if (is.na(u)) return(NULL)
    tibble(column = cc, description = NA_character_, route = "distribution scan",
           source = "distribution scan", n_values = length(v),
           n_distinct = dplyr::n_distinct(v), median_value = median(v),
           inferred_units = u, usable = TRUE)
  })
  candidates <- bind_rows(candidates, scan) %>% arrange(desc(usable), desc(n_values))
  cat("  NOTE: a scan match is a guess. Confirm the chosen column before publishing.\n")
}

print(as.data.frame(candidates %>% select(column, description, route, n_values,
                                          median_value, inferred_units, usable)))
out$creatinine_candidates <- candidates

chosen <- candidates %>% filter(usable) %>% slice(1)
if (nrow(chosen) == 0) stop("No usable creatinine column. Inspect creatinine_candidates.")

cat("\nCHOSEN:", chosen$column, "|", chosen$description %||% "(no description)",
    "| units:", chosen$inferred_units, "\n")


# ==========================================================
# 4. Extract baseline creatinine
# ==========================================================

labs2 <- labs %>% clean_names()
col   <- make_clean_names(chosen$column)
kv    <- intersect(c("viscode2", "viscode"), names(labs2))[1]

creatinine_baseline <- labs2 %>%
  mutate(rid = as.character(rid),
         creatinine_raw = parse_num(.data[[col]]),
         creatinine = if (chosen$inferred_units == "umol/L")
           creatinine_raw / 88.4 else creatinine_raw,   # to mg/dL
         visit_key = if (!is.na(kv)) as.character(.data[[kv]]) else NA_character_) %>%
  filter(!is.na(creatinine), creatinine > 0.2, creatinine < 15) %>%
  mutate(is_bl = is.na(visit_key) | visit_key %in% c("bl", "sc", "scmri", "screen")) %>%
  arrange(rid, desc(is_bl)) %>%
  group_by(rid) %>%
  summarise(creatinine = first(creatinine),
            visit_used = first(visit_key), .groups = "drop")

cat("\nBaseline creatinine for", nrow(creatinine_baseline), "participants\n")
cat("  median", round(median(creatinine_baseline$creatinine), 2), "mg/dL",
    " IQR", paste(round(quantile(creatinine_baseline$creatinine, c(.25, .75)), 2),
                  collapse = "-"), "\n")
cat("  (expect a median near 0.8-1.0 mg/dL in this age group)\n")

out$creatinine_summary <- tibble(
  chosen_column = chosen$column, description = chosen$description,
  route = chosen$route, inferred_units = chosen$inferred_units,
  n_participants = nrow(creatinine_baseline),
  median_mgdl = median(creatinine_baseline$creatinine),
  q1 = quantile(creatinine_baseline$creatinine, .25),
  q3 = quantile(creatinine_baseline$creatinine, .75))

out$creatinine_baseline <- creatinine_baseline

write_csv(creatinine_baseline, file.path(clean_dir, "creatinine_baseline.csv"))

# How much of each biomarker panel does this actually cover?
model_ready <- tryCatch(
  read_csv(file.path(clean_dir, "analysis_master_model_ready.csv"),
           show_col_types = FALSE, guess_max = 100000) %>% clean_names(),
  error = function(e) NULL)

if (!is.null(model_ready)) {
  out$coverage <- map_dfr(c("gfap_quanterix", "strem2_msd_corrected",
                            "vegf_plasma_qc", "sicam1_plasma_qc", "svcam1_plasma_qc"),
    function(m) {
      ids <- model_ready %>% filter(!is.na(.data[[m]])) %>%
        pull(rid) %>% as.character() %>% unique()
      tibble(marker = m, panel_n = length(ids),
             with_creatinine = sum(ids %in% creatinine_baseline$rid),
             pct_covered = round(100 * mean(ids %in% creatinine_baseline$rid), 1))
    })
  cat("\nCreatinine coverage by panel:\n"); print(as.data.frame(out$coverage))
  cat("\nIf vascular coverage is below about 60%, adjusting for eGFR will cost\n",
      "more power than it buys. The tiered models in script 18 will show this.\n")
}

names(out) <- substr(make.unique(names(out)), 1, 31)
openxlsx::write.xlsx(out, file.path(results_dir, "creatinine_extraction.xlsx"),
                     overwrite = TRUE)

cat("\nDone.\n")
cat("Creatinine file:", file.path(clean_dir, "creatinine_baseline.csv"), "\n")
cat("Audit trail:    ", file.path(results_dir, "creatinine_extraction.xlsx"), "\n")
