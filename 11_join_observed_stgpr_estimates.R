#!/usr/bin/env Rscript
# For scenario_one *_prediction.csv files:
# 1) Keep rows where above.30 == "no"
# 2) Left join STGPR predictions (stgpr_pred) by district_code + year
# 3) Write one output CSV per bug-drug

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)

scenario_one_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("data", "data_for_modelling", "scenario_one")
}

crosswalk_model_dir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path("data", "models", "scenario_crosswalkeddata")
}

out_dir <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  file.path("data", "data_for_modelling", "scenario_one_above30_no_with_stgpr")
}

fs::dir_create(out_dir)

prediction_files <- list.files(
  scenario_one_dir,
  pattern = "_alllevels_merged_prediction\\.csv$",
  full.names = TRUE
)

if (length(prediction_files) == 0) {
  message("No *_alllevels_merged_prediction.csv files found in: ", scenario_one_dir)
  quit(status = 0)
}

for (pf in prediction_files) {
  bug_core <- sub("_alllevels_merged_prediction\\.csv$", "", basename(pf))

  scen_df <- tryCatch(readr::read_csv(pf, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(scen_df)) {
    message("Could not read: ", pf, " -- skipping")
    next
  }

  required_scen <- c("district_code", "year", "above.30")
  if (!all(required_scen %in% names(scen_df))) {
    message("Missing one of district_code/year/above.30 in ", basename(pf), " -- skipping")
    next
  }

  scen_no <- scen_df %>%
    mutate(
      district_code = as.character(district_code),
      year = as.integer(year),
      above.30 = tolower(trimws(as.character(.data[["above.30"]])))
    ) %>%
    filter(above.30 == "no")

  stgpr_fp <- file.path(
    crosswalk_model_dir,
    paste0(bug_core, "_final_withcrosswalk_merged"),
    "stgpr",
    paste0(bug_core, "_final_withcrosswalk_merged_prediction_with_stgpr.csv")
  )

  if (!file.exists(stgpr_fp)) {
    message("STGPR file not found for ", bug_core, ": ", stgpr_fp, " -- writing filtered rows without stgpr_pred")
    out_fp <- file.path(out_dir, paste0(bug_core, "_alllevels_prediction_above30no_with_stgpr.csv"))
    readr::write_csv(scen_no, out_fp)
    next
  }

  stgpr_df <- tryCatch(readr::read_csv(stgpr_fp, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(stgpr_df)) {
    message("Could not read STGPR file: ", stgpr_fp, " -- writing filtered rows without stgpr_pred")
    out_fp <- file.path(out_dir, paste0(bug_core, "_alllevels_prediction_above30no_with_stgpr.csv"))
    readr::write_csv(scen_no, out_fp)
    next
  }

  stgpr_year_col <- if ("year" %in% names(stgpr_df)) {
    "year"
  } else if ("year_taken" %in% names(stgpr_df)) {
    "year_taken"
  } else {
    NA_character_
  }

  if (!("district_code" %in% names(stgpr_df)) || is.na(stgpr_year_col) || !("stgpr_pred" %in% names(stgpr_df))) {
    message("STGPR file missing district/year/stgpr_pred for ", bug_core, " -- writing filtered rows without join")
    out_fp <- file.path(out_dir, paste0(bug_core, "_alllevels_prediction_above30no_with_stgpr.csv"))
    readr::write_csv(scen_no, out_fp)
    next
  }

  stgpr_small <- stgpr_df %>%
    transmute(
      district_code = as.character(district_code),
      year = as.integer(.data[[stgpr_year_col]]),
      stgpr_pred = as.numeric(stgpr_pred)
    ) %>%
    distinct(district_code, year, .keep_all = TRUE)

  joined <- scen_no %>%
    left_join(stgpr_small, by = c("district_code", "year"))

  core_cols <- c(
    "district_code",
    "year_taken",
    paste0(bug_core, "_resistant"),
    paste0(bug_core, "_sensitive"),
    paste0(bug_core, "_tested"),
    "above.30",
    paste0(bug_core, "_resistance_per"),
    "stgpr_pred"
  )

  joined <- joined %>%
    select(any_of(core_cols))

  out_fp <- file.path(out_dir, paste0(bug_core, "_alllevels_prediction_above30no_with_stgpr.csv"))
  readr::write_csv(joined, out_fp)
  message("Wrote: ", out_fp, " (rows=", nrow(joined), ")")
}

message("Done. Outputs in: ", out_dir)
