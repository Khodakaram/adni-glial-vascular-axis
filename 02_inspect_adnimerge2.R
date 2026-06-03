# ==========================================================
# 02_inspect_adnimerge2.R
# Purpose: Inspect ADNIMERGE2 R Data files and identify useful tables
# ==========================================================

library(tidyverse)
library(janitor)
library(openxlsx)
library(stringr)

project_dir <- "R:/ADNI_Project"

adnimerge2_dir <- file.path(project_dir, "00_raw_data", "ADNIMERGE2")
data_dir <- file.path(adnimerge2_dir, "data")
out_dir <- file.path(project_dir, "01_inventory")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Check folders ----
cat("ADNIMERGE2 folder exists:", dir.exists(adnimerge2_dir), "\n")
cat("ADNIMERGE2/data folder exists:", dir.exists(data_dir), "\n")

# ---- 2. List all R data files ----
rdata_files <- list.files(
  data_dir,
  pattern = "\\.(rda|RData)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("Number of R Data files found:", length(rdata_files), "\n")
print(rdata_files)

if (length(rdata_files) == 0) {
  stop("No .rda or .RData files found in ADNIMERGE2/data.")
}

# ---- 3. Keywords needed for your project ----
keywords <- c(
  "RID", "PTID", "VISCODE", "VISCODE2",
  "DX", "DIAGNOSIS", "DXCHANGE",
  "AGE", "PTGENDER", "SEX", "PTEDUCAT",
  "APOE", "APOE4",
  "ADAS", "ADAS13",
  "MMSE",
  "EXAMDATE", "USERDATE",
  "WMH", "ICV"
)

keyword_regex <- paste(keywords, collapse = "|")

# ---- 4. Function to inspect one RData file ----
inspect_rdata_file <- function(file) {
  
  env <- new.env()
  
  result <- tryCatch({
    
    object_names <- load(file, envir = env)
    
    map_dfr(object_names, function(obj_name) {
      
      obj <- get(obj_name, envir = env)
      
      if (inherits(obj, "data.frame")) {
        
        cols <- names(obj)
        matched_cols <- cols[str_detect(toupper(cols), keyword_regex)]
        
        tibble(
          file_name = basename(file),
          file_path = file,
          object_name = obj_name,
          object_class = paste(class(obj), collapse = " | "),
          n_rows = nrow(obj),
          n_cols = ncol(obj),
          matched_columns = paste(matched_cols, collapse = " | "),
          all_columns = paste(cols, collapse = " | ")
        )
        
      } else {
        
        tibble(
          file_name = basename(file),
          file_path = file,
          object_name = obj_name,
          object_class = paste(class(obj), collapse = " | "),
          n_rows = NA_integer_,
          n_cols = NA_integer_,
          matched_columns = NA_character_,
          all_columns = NA_character_
        )
      }
    })
    
  }, error = function(e) {
    
    tibble(
      file_name = basename(file),
      file_path = file,
      object_name = NA_character_,
      object_class = "ERROR",
      n_rows = NA_integer_,
      n_cols = NA_integer_,
      matched_columns = e$message,
      all_columns = NA_character_
    )
  })
  
  return(result)
}

# ---- 5. Inspect all R data files ----
adnimerge2_inventory <- map_dfr(rdata_files, inspect_rdata_file)

# ---- 6. Find candidate useful tables ----
candidate_tables <- adnimerge2_inventory %>%
  filter(!is.na(all_columns)) %>%
  filter(
    str_detect(toupper(all_columns), "RID|PTID") |
      str_detect(toupper(all_columns), "VISCODE|VISIT")
  ) %>%
  arrange(desc(n_rows))

# ---- 7. Find files likely relevant to your proposal ----
proposal_related_tables <- adnimerge2_inventory %>%
  filter(
    str_detect(toupper(file_name), "ADAS|MMSE|APOE|DX|DEMO|PTDEMOG|VISIT|WMH|MERGE") |
      str_detect(toupper(matched_columns), "ADAS|MMSE|APOE|DX|AGE|PTGENDER|VISCODE|WMH")
  ) %>%
  arrange(file_name)

# ---- 8. Save output ----
openxlsx::write.xlsx(
  list(
    all_rdata_files = adnimerge2_inventory,
    candidate_tables = candidate_tables,
    proposal_related_tables = proposal_related_tables
  ),
  file.path(out_dir, "ADNIMERGE2_rdata_inventory.xlsx"),
  overwrite = TRUE
)

write_csv(
  adnimerge2_inventory,
  file.path(out_dir, "ADNIMERGE2_rdata_inventory.csv")
)

write_csv(
  proposal_related_tables,
  file.path(out_dir, "ADNIMERGE2_proposal_related_tables.csv")
)

cat("Done. Check:\n")
cat(file.path(out_dir, "ADNIMERGE2_rdata_inventory.xlsx"), "\n")