#!/usr/bin/env Rscript
# Step 7: Create different modelling scenarios
# Scenario one: keep only modelling rows with above.30 == "yes" and move above.30 == "no" into prediction.
# Scenario two: collect all merged crosswalk files into a single folder.

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
})

src_dir <- file.path("data", "amr_with_socio_economic_scaled_modelling_prediction")
scenario_one_dir <- file.path("data", "data_for_modelling", "scenario_one")
scenario_two_dir <- file.path("data", "data_for_modelling", "scenario_two")

if (!dir.exists(src_dir)) {
  stop("Source directory does not exist: ", src_dir)
}
if (!dir.exists(scenario_one_dir)) {
  dir.create(scenario_one_dir, recursive = TRUE)
}
if (!dir.exists(scenario_two_dir)) {
  dir.create(scenario_two_dir, recursive = TRUE)
}
# Scenario one: filter modelling on above.30 == "yes"; prediction = original prediction + modelling rows where above.30 == "no"
model_files <- list.files(src_dir, pattern = "_alllevels_merged_modelling\\.csv$", full.names = TRUE)

if (length(model_files) == 0) {
  message("No modelling files found for scenario one in ", src_dir)
} else {
  for (mod_path in model_files) {
    pred_path <- sub("_modelling\\.csv$", "_prediction.csv", mod_path)
    if (!file.exists(pred_path)) {
      warning("Prediction file missing for ", basename(mod_path), "; skipping")
      next
    }

    mod_df <- read.csv(mod_path, stringsAsFactors = FALSE, check.names = FALSE)
    pred_df <- read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE)

    if (!"above.30" %in% colnames(mod_df)) {
      warning("Column 'above.30' not found in ", basename(mod_path), "; skipping")
      next
    }

    mod_yes <- mod_df %>% filter(above.30 == "yes")
    mod_no <- mod_df %>% filter(above.30 == "no")
    pred_out <- bind_rows(pred_df, mod_no)

    mod_out_path <- file.path(scenario_one_dir, basename(mod_path))
    pred_out_path <- file.path(scenario_one_dir, basename(pred_path))

    write.csv(mod_yes, mod_out_path, row.names = FALSE)
    write.csv(pred_out, pred_out_path, row.names = FALSE)

    message("Scenario one written: ", basename(mod_out_path), " (", nrow(mod_yes), " rows), ",
            basename(pred_out_path), " (", nrow(pred_out), " rows)")
  }
}

# Scenario two: collect all crosswalk merged CSVs into a single folder.
crosswalk_files <- list.files(src_dir, pattern = "_withcrosswalk_merged_.*\\.csv$", full.names = TRUE)

if (length(crosswalk_files) == 0) {
  message("No crosswalk files found for scenario two in ", src_dir)
} else {
  copied2 <- file.copy(crosswalk_files, scenario_two_dir, overwrite = TRUE)
  message("Scenario two: copied ", sum(copied2), " files to ", scenario_two_dir)
  if (any(!copied2)) {
    warning("Scenario two: some files were not copied: ",
            paste(basename(crosswalk_files[!copied2]), collapse = ", "))
  }
}
