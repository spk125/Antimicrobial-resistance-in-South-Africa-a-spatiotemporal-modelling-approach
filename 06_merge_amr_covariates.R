# Step 5: Read imputed socioeconomic data

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

step2_path <- file.path("data", "step2_imputation_gam_output.xlsx")

if (!file.exists(step2_path)) {
  stop("Expected file not found: ", step2_path, ". Run step2_imputation_gam first.")
}

imputed_data <- read_xlsx(step2_path)
message("Loaded imputed data: ", nrow(imputed_data), " rows, ", ncol(imputed_data), " columns")
print(head(imputed_data, 5))
message("Column names:")
print(colnames(imputed_data))

bug_drug_dir <- file.path("data", "amr_data_bug_drug")
crosswalk_dir <- file.path(bug_drug_dir, "for_crosswalk")
# Use existing output directory casing to avoid creating duplicates
output_dir <- file.path("data", "amr_with_socio_economic")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

normalize_imputed_keys <- function(df) {
  df %>%
    mutate(
      district_code = trimws(as.character(district_code)),
      year = as.integer(year)
    )
}

normalize_bug_keys <- function(df) {
  df %>%
    mutate(
      district_code = trimws(as.character(district_code)),
      year_taken = as.integer(year_taken),
      year = year_taken
    )
}

imputed_for_join <- normalize_imputed_keys(imputed_data)

merge_bug_drug_dir <- function(dir_path, out_dir, name_suffix = NULL) {
  bug_drug_files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(bug_drug_files) == 0) {
    message("No bug–drug CSV files found in ", dir_path)
    return(invisible(NULL))
  }

  for (csv_path in bug_drug_files) {
    bd_df <- read.csv(csv_path, stringsAsFactors = FALSE)
    bd_df <- normalize_bug_keys(bd_df)

    missing_cols <- setdiff(c("district_code", "year"), colnames(bd_df))
    if (length(missing_cols) > 0) {
      warning("Skipping ", csv_path, " (missing keys: ", paste(missing_cols, collapse = ", "), ")")
      next
    }

    merged <- bd_df %>%
      full_join(imputed_for_join, by = c("district_code", "year"))

    merged <- merged %>%
      filter(
        !is.na(district_code),
        district_code != "",
        district_code != "NaN"
      ) %>%
      filter(year >= 2013, year <= 2023)

    base_name <- tools::file_path_sans_ext(basename(csv_path))
    if (!is.null(name_suffix) && nzchar(name_suffix)) {
      base_name <- paste0(base_name, "_", name_suffix)
    }
    out_name <- paste0(base_name, "_merged.csv")
    out_path <- file.path(out_dir, out_name)
    write.csv(merged, out_path, row.names = FALSE)
    message("Merged and wrote: ", out_path, " (", nrow(merged), " rows)")
  }
}

merge_bug_drug_dir(bug_drug_dir, output_dir)
merge_bug_drug_dir(crosswalk_dir, output_dir, "withcrosswalk")

# ------------------------------------------------------------------
# Summarize abaumanii_carbapenem_level1_merged counts by district_code
# ------------------------------------------------------------------
target_file <- file.path(output_dir, "abaumanii_carbapenem_level1_merged.csv")

if (file.exists(target_file)) {
  counts_under_11 <- read.csv(target_file, stringsAsFactors = FALSE) %>%
    group_by(district_code) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n < 11) %>%
    arrange(district_code)

  message("District codes with < 11 records in ", basename(target_file), ":")
  print(counts_under_11)
} else {
  message("Target merged file not found: ", target_file)
}
