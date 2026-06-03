# ==========================================================
# 01_scan_adni_files.R
# Purpose: Scan all ADNI CSV/XLS/XLSX files and find relevant variables
# ==========================================================

library(tidyverse)
library(readxl)
library(openxlsx)
library(janitor)
library(stringr)

# ---- 1. Set your project path ----
project_dir <- "R:/ADNI_Project"   # CHANGE THIS to your own folder

raw_dir <- file.path(project_dir, "00_raw_data")
bio_dir <- file.path(raw_dir, "biomarkers_excel")
out_dir <- file.path(project_dir, "01_inventory")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. List all files ----
all_files <- list.files(
  bio_dir,
  pattern = "\\.(csv|xlsx|xls)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

length(all_files)
print(all_files)

# ---- 3. Keywords based on your proposal ----
keywords <- c(
  "GFAP", "TREM2", "sTREM2", "YKL", "YKL40", "CHI3L1",
  "VEGF", "ICAM", "sICAM", "VCAM", "sVCAM",
  "ADAS", "ADAS13", "MMSE",
  "APOE", "APOE4",
  "WMH", "WHITE", "HYPERINTENSITY"
)

# ---- 4. Function to safely read file columns ----
read_file_preview <- function(file) {
  
  ext <- tools::file_ext(file) |> tolower()
  
  result <- tryCatch({
    
    if (ext == "csv") {
      df <- readr::read_csv(file, n_max = 20, show_col_types = FALSE)
      sheet_name <- NA_character_
    } else {
      sheets <- readxl::excel_sheets(file)
      df <- readxl::read_excel(file, sheet = sheets[1], n_max = 20)
      sheet_name <- sheets[1]
    }
    
    tibble(
      file_path = file,
      file_name = basename(file),
      sheet = sheet_name,
      columns = paste(names(df), collapse = " | ")
    )
    
  }, error = function(e) {
    tibble(
      file_path = file,
      file_name = basename(file),
      sheet = NA_character_,
      columns = paste("ERROR:", e$message)
    )
  })
  
  return(result)
}

# ---- 5. Scan all files ----
file_inventory <- map_dfr(all_files, read_file_preview)

# ---- 6. Search columns for proposal keywords ----
inventory_hits <- file_inventory %>%
  mutate(
    matched_keywords = map_chr(columns, function(x) {
      hits <- keywords[str_detect(toupper(x), toupper(keywords))]
      paste(unique(hits), collapse = ", ")
    }),
    has_relevant_keyword = matched_keywords != ""
  )

# ---- 7. Save results ----
write_csv(file_inventory, file.path(out_dir, "all_adni_file_inventory.csv"))
write_csv(inventory_hits, file.path(out_dir, "relevant_keyword_hits.csv"))

openxlsx::write.xlsx(
  list(
    all_files = file_inventory,
    relevant_hits = inventory_hits
  ),
  file.path(out_dir, "ADNI_file_inventory.xlsx"),
  overwrite = TRUE
)

cat("Done. Check this file:\n")
cat(file.path(out_dir, "ADNI_file_inventory.xlsx"))