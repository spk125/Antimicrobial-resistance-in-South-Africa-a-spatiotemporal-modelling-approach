#!/usr/bin/env Rscript
# calc_annualized_rate_by_bugdrug_stgpr2023_only.R
# Re-run the annualized-rate calculations and plots but restrict to districts that
# have STGPR model predictions for the END_YEAR (default 2023). Outputs are
# written to outputs/annualized_rates/stgpr2023_only/ with filenames suffixed accordingly.

required_pkgs <- c("dplyr", "ggplot2", "readr", "fs", "broom")
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

library(dplyr)
library(ggplot2)
library(readr)
library(fs)
library(broom)

input_dir <- "data/observed_stgpr_joined"
output_dir_base <- "outputs/annualized_rates/stgpr2023_only"
fs::dir_create(output_dir_base)

# Years to compare
START_YEAR <- 2014
END_YEAR <- 2023
TOLERANCE_YEARS <- 1

csv_files <- dir(input_dir, pattern = "\\.csv$", full.names = TRUE)

all_summaries <- list()

process_file <- function(file_path) {
  message("Processing (STGPR2023-only): ", file_path)
  df <- readr::read_csv(file_path, show_col_types = FALSE)

  # ensure year numeric
  df <- df %>% mutate(year = as.integer(year))

  # find districts with a non-missing stgpr_pred at END_YEAR
  districts_with_stgpr_end <- df %>% filter(year == END_YEAR & !is.na(stgpr_pred)) %>% pull(district_code) %>% unique()
  if (length(districts_with_stgpr_end) == 0) {
    message("  No STGPR predictions at ", END_YEAR, " in this file — skipping.")
    return(NULL)
  }

  # restrict dataset to rows for those districts and within the year window (+/- tolerance)
  df_sub <- df %>%
    filter(district_code %in% districts_with_stgpr_end) %>%
    filter(!is.na(year) & year >= (START_YEAR - TOLERANCE_YEARS) & year <= (END_YEAR + TOLERANCE_YEARS)) %>%
    mutate(value = dplyr::coalesce(resistance_per, stgpr_pred))

  # compute per-district summary comparing START_YEAR and END_YEAR (same logic as original script)
  summary <- df_sub %>%
    group_by(district_code) %>%
    group_modify(~{
      d <- .x %>% arrange(year)
      n_records <- sum(!is.na(d$value))

      # start candidate
      start_cands <- which(!is.na(d$value) & abs(d$year - START_YEAR) <= TOLERANCE_YEARS)
      if (length(start_cands) > 0) {
        sel <- start_cands[which.min(abs(d$year[start_cands] - START_YEAR))]
        start_year_actual <- d$year[sel]
        start_value <- d$value[sel]
      } else {
        start_year_actual <- NA_integer_
        start_value <- NA_real_
      }

      # end candidate
      end_cands <- which(!is.na(d$value) & abs(d$year - END_YEAR) <= TOLERANCE_YEARS)
      if (length(end_cands) > 0) {
        sel2 <- end_cands[which.min(abs(d$year[end_cands] - END_YEAR))]
        end_year_actual <- d$year[sel2]
        end_value <- d$value[sel2]
      } else {
        end_year_actual <- NA_integer_
        end_value <- NA_real_
      }

      tibble(
        intended_start_year = START_YEAR,
        intended_end_year = END_YEAR,
        start_year_actual = start_year_actual,
        start_value = start_value,
        end_year_actual = end_year_actual,
        end_value = end_value,
        n_records = n_records
      )
    }) %>% ungroup()

  # compute annualized changes
  summary <- summary %>% mutate(
    n_years = ifelse(!is.na(start_year_actual) & !is.na(end_year_actual), end_year_actual - start_year_actual, NA_real_),
    annual_absolute_change = ifelse(!is.na(start_value) & !is.na(end_value) & n_years > 0, (end_value - start_value) / n_years, NA_real_),
    annual_relative_change = ifelse(!is.na(start_value) & !is.na(end_value) & n_years > 0 & start_value > 0, ((end_value / start_value)^(1 / n_years) - 1), NA_real_)
  ) %>% mutate(
    annual_absolute_change_pct = ifelse(!is.na(annual_absolute_change), annual_absolute_change * 100, NA_real_),
    annual_relative_change_pct = ifelse(!is.na(annual_relative_change), annual_relative_change * 100, NA_real_)
  )

  # compute slope_per_year where possible (using df_sub)
  slopes <- df_sub %>% filter(!is.na(value)) %>% group_by(district_code) %>% filter(n() >= 2) %>%
    do(broom::tidy(lm(value ~ year, data = .)) %>% filter(term == "year") %>% select(estimate)) %>% rename(slope_per_year = estimate)

  summary <- summary %>% left_join(slopes, by = "district_code")

  # file-safe base and write
  base <- tools::file_path_sans_ext(basename(file_path))
  out_base <- paste0(base, "_stgpr2023_only")
  readr::write_csv(summary, file.path(output_dir_base, paste0(out_base, "_annualized_rates.csv")))

  # timeseries plot for the selected districts only
  p_ts <- df_sub %>% filter(!is.na(value)) %>% ggplot(aes(x = year, y = value, group = district_code)) +
    geom_line(alpha = 0.6, linewidth = 0.7, colour = "grey40") +
    geom_point(alpha = 0.8, size = 1, colour = "grey20") +
    labs(title = paste(base, "— value over time (districts with STGPR at", END_YEAR, ")"), y = "Value", x = "Year") +
    theme_minimal()
  ggsave(filename = file.path(output_dir_base, paste0(out_base, "_timeseries.png")), p_ts, width = 10, height = 6, dpi = 150)

  # bar plot of annual absolute change (pct)
  p_bar <- summary %>% arrange(annual_absolute_change_pct) %>% mutate(district_code = factor(district_code, levels = unique(district_code))) %>%
    ggplot(aes(x = district_code, y = annual_absolute_change_pct)) + geom_col(fill = "steelblue") + coord_flip() +
    labs(title = paste(base, "— annual absolute change (STGPR2023 districts)"), x = "District", y = "Annual absolute change (%)") + theme_minimal() +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1, suffix = "%"))
  ggsave(filename = file.path(output_dir_base, paste0(out_base, "_annual_absolute_change_bar.png")), p_bar, width = 10, height = 8, dpi = 150)

  summary <- summary %>% mutate(bugdrug = base)
  return(summary)
}

if (length(csv_files) == 0) stop("No CSV files found in ", input_dir)

for (f in csv_files) {
  out <- tryCatch(process_file(f), error = function(e) { message("Error processing ", f, ": ", e$message); NULL })
  if (!is.null(out)) all_summaries[[length(all_summaries) + 1]] <- out
}

if (length(all_summaries) > 0) {
  master <- bind_rows(all_summaries)
  readr::write_csv(master, file.path(output_dir_base, "all_bugdrug_annualized_rates_stgpr2023_only.csv"))
  message("Wrote combined summary to ", file.path(output_dir_base, "all_bugdrug_annualized_rates_stgpr2023_only.csv"))
} else {
  message("No summaries were produced for STGPR2023-only processing.")
}

message("Done. Per-file summaries and plots are in: ", output_dir_base)
