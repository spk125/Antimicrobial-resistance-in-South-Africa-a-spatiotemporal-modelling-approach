#!/usr/bin/env Rscript
# Step 8: Create training and testing datasets
# Splits each *_modelling.csv using a district holdout based on an ordered list,
# targeting ~20% of districts for testing.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

set.seed(123)  # for reproducible splits; adjust/remove if fully random each run is preferred

scenario_dirs <- file.path(
  "data",
  "data_for_modelling",
  c("scenario_one", "scenario_two", "scenario_three")
)

for (scenario_dir in scenario_dirs) {
  if (!dir.exists(scenario_dir)) {
    dir.create(scenario_dir, recursive = TRUE)
  }
  if (!dir.exists(scenario_dir)) {
    message("Scenario directory does not exist, skipping: ", scenario_dir)
    next
  }

  mod_files <- list.files(scenario_dir, pattern = "_modelling\\.csv$", full.names = TRUE)

  if (length(mod_files) == 0) {
    message("No modelling files found in ", scenario_dir)
    next
  }

  if (basename(scenario_dir) == "scenario_three") {
    message("Scenario three: leaving modelling files unchanged in ", scenario_dir)
    next
  }

  for (f in mod_files) {
    df <- read_csv(f, show_col_types = FALSE)
    n <- nrow(df)

    if (n == 0) {
      message("Empty file, skipping split: ", basename(f))
      next
    }

    # District-based split: hold out whole districts to avoid overlap
    district_col <- if ("district_code" %in% colnames(df)) "district_code" else if ("district_id" %in% colnames(df)) "district_id" else NA
    if (is.na(district_col)) {
      warning("No district_code or district_id column found in ", basename(f), "; skipping")
      next
    }

    custom_order <- c(
      "DC31", "DC18", "EKU", "DC40", "TSH", "JHB", "DC28", "CPT",
      "MAN", "DC35", "DC9", "DC21", "DC42", "ETH", "DC22", "NMA",
      "DC15", "DC29", "BUF", "DC37", "DC32", "DC48", "DC13", "DC30"
    )

    base <- tools::file_path_sans_ext(basename(f))
    order_metrics_path <- file.path(
      "data", "models", "order",
      paste0(base, "_training_districtsplit_stgpr_loo_district_metrics.csv")
    )

    district_levels <- unique(df[[district_col]])

    if (file.exists(order_metrics_path)) {
      order_df <- read_csv(order_metrics_path, show_col_types = FALSE)
      if (all(c("district", "rsq") %in% names(order_df))) {
        ordered_districts <- order_df %>%
          arrange(rsq) %>%
          pull(district)
      } else {
        warning("Order file missing district/rsq columns: ", order_metrics_path, "; using custom order")
        remaining_levels <- setdiff(district_levels, custom_order)
        ordered_districts <- c(custom_order, remaining_levels)
      }
    } else {
      remaining_levels <- setdiff(district_levels, custom_order)
      ordered_districts <- c(custom_order, remaining_levels)
    }

    ordered_districts <- ordered_districts[ordered_districts %in% district_levels]

    n_test_districts <- max(1, ceiling(length(ordered_districts) * 0.20))
    test_districts <- tail(ordered_districts, n_test_districts)

    if (length(test_districts) == 0 || length(test_districts) == length(ordered_districts)) {
      warning("Unable to create district-based split for ", basename(f), "; skipping")
      next
    }

    train_df_d <- df[!(df[[district_col]] %in% test_districts), , drop = FALSE]
    test_df_d  <- df[df[[district_col]] %in% test_districts, , drop = FALSE]

    train_path_d <- file.path(scenario_dir, paste0(base, "_training.csv"))
    test_path_d  <- file.path(scenario_dir, paste0(base, "_testing.csv"))

    write_csv(train_df_d, train_path_d)
    write_csv(test_df_d, test_path_d)

    message("District split ", basename(f), " -> training: ", nrow(train_df_d), " rows; testing: ", nrow(test_df_d), " rows")
  }
}
