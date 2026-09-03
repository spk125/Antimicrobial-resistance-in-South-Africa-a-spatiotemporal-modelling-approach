# Step 6: Prepare data for modelling

suppressPackageStartupMessages({
  library(dplyr)
})

input_dir <- file.path("data", "amr_with_socio_economic")
output_dir <- file.path("data", "amr_with_socio_econmic_scaled")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
mod_pred_dir <- file.path("data", "amr_with_socio_economic_scaled_modelling_prediction")
if (!dir.exists(mod_pred_dir)) dir.create(mod_pred_dir, recursive = TRUE)

socio_cols <- c(
  "Hospital beds per 10 000 target population (rescaled)",
  "Antenatal client initiated on ART rate",
  "Diarrhoea case fatality under 5 years rate",
  "Hospital beds per 10 000 target population",
  "Pneumonia case fatality under 5 years rate",
  "Severe acute malnutrition case fatality under 5 years rate",
  "Adult ART Total",
  "Medical practitioners per 100 000 population",
  "No of households",
  "Percentage of households with access to improved sanitation (index)",
  "Pharmacists per 100 000 population",
  "Population.y",
  "Clients remaining on ART rate",
  "Delivery by Caesarean section rate (district hospitals)",
  "Proportion of health facilities with essential medicines",
  "Tracer items stock-out rate (fixed clinic/CHC/CDC)",
  "Antenatal client HIV 1st test positive",
  "HIV test positive 15 years and older (excl ANC)",
  "Death under 5 years against live birth rate",
  "Access to services (Sewerage and sanitation)",
  "Access to services (Water)",
  "Acute malnutrition under five",
  "Adult literacy",
  "GDPR per capita",
  "Gini coefficient",
  "Human Development Index",
  "Immunisation rate",
  "Infant mortality",
  "Low Birth Weight",
  "Mean years of schooling",
  "Patients remaining on ART",
  "Patients starting ART treatment",
  "Under five mortality",
  "Expected years of schooling",
  "HIV test positive client 15 years and older rate (incl ANC)",
  "Maternal mortality"
)

scale_numeric <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) return(x_num)
  as.numeric(scale(x_num))
}

files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  message("No CSV files found in ", input_dir)
} else {
  for (csv_path in files) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

    cols_present <- intersect(socio_cols, colnames(df))
    cols_missing <- setdiff(socio_cols, colnames(df))
    if (length(cols_missing) > 0) {
      message("Warning: skipping missing columns in ", basename(csv_path), ": ", paste(cols_missing, collapse = ", "))
    }

    df_scaled <- df %>%
      mutate(across(all_of(cols_present), scale_numeric))

    out_path <- file.path(output_dir, basename(csv_path))
    write.csv(df_scaled, out_path, row.names = FALSE)
    message("Scaled and wrote: ", out_path)

    # Split into modelling (year_taken present) and prediction (year_taken missing)
    if (!"year_taken" %in% colnames(df_scaled)) {
      message("Skipping modelling/prediction split for ", basename(csv_path), " (no year_taken column)")
    } else {
      base_name <- tools::file_path_sans_ext(basename(csv_path))
      modelling <- df_scaled %>% filter(!is.na(year_taken))
      prediction <- df_scaled %>% filter(is.na(year_taken))

      out_mod <- file.path(mod_pred_dir, paste0(base_name, "_modelling.csv"))
      out_pred <- file.path(mod_pred_dir, paste0(base_name, "_prediction.csv"))

      write.csv(modelling, out_mod, row.names = FALSE)
      write.csv(prediction, out_pred, row.names = FALSE)

      message("  Wrote modelling: ", out_mod, " (", nrow(modelling), " rows)")
      message("  Wrote prediction: ", out_pred, " (", nrow(prediction), " rows)")
    }
  }
}
