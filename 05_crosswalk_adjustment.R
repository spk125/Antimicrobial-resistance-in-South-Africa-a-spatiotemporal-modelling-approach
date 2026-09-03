#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
})

output_dir <- file.path("data", "amr_data_bug_drug")
if (!dir.exists(output_dir)) {
  stop("Expected directory not found: ", output_dir)
}

amr_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(amr_files) == 0) {
  stop("No files found in: ", output_dir)
}

base_names <- basename(amr_files)
level1_files <- amr_files[grepl("_level1\\.csv$", base_names)]
without_files <- amr_files[grepl("_without_level1\\.csv$", base_names)]
base_level1 <- sub("_level1\\.csv$", "", basename(level1_files))
base_without <- sub("_without_level1\\.csv$", "", basename(without_files))
bases <- intersect(base_level1, base_without)
if (length(bases) == 0) {
  stop("No matching level1/without_level1 file pairs found in: ", output_dir)
}

district_province <- readr::read_csv(
  "District_Code,District_Name,Province\nBUF,Buffalo City Metropolitan Municipality,Eastern Cape\nNMA,Nelson Mandela Bay Metropolitan Municipality,Eastern Cape\nDC10,Sarah Baartman District Municipality,Eastern Cape\nDC12,Amathole District Municipality,Eastern Cape\nDC13,Chris Hani District Municipality,Eastern Cape\nDC14,Joe Gqabi District Municipality,Eastern Cape\nDC15,OR Tambo District Municipality,Eastern Cape\nDC44,Alfred Nzo District Municipality,Eastern Cape\nMAN,Mangaung Metropolitan Municipality,Free State\nDC16,Xhariep District Municipality,Free State\nDC18,Lejweleputswa District Municipality,Free State\nDC19,Thabo Mofutsanyana District Municipality,Free State\nDC20,Fezile Dabi District Municipality,Free State\nEKU,City of Ekurhuleni Metropolitan Municipality,Gauteng\nJHB,City of Johannesburg Metropolitan Municipality,Gauteng\nTSH,City of Tshwane Metropolitan Municipality,Gauteng\nDC42,Sedibeng District Municipality,Gauteng\nDC48,West Rand District Municipality,Gauteng\nETH,eThekwini Metropolitan Municipality,KwaZulu-Natal\nDC21,Ugu District Municipality,KwaZulu-Natal\nDC22,uMgungundlovu District Municipality,KwaZulu-Natal\nDC23,uThukela District Municipality,KwaZulu-Natal\nDC24,uMzinyathi District Municipality,KwaZulu-Natal\nDC25,Amajuba District Municipality,KwaZulu-Natal\nDC26,Zululand District Municipality,KwaZulu-Natal\nDC27,uMkhanyakude District Municipality,KwaZulu-Natal\nDC28,King Cetshwayo District Municipality,KwaZulu-Natal\nDC29,iLembe District Municipality,KwaZulu-Natal\nDC43,Harry Gwala District Municipality,KwaZulu-Natal\nDC33,Mopani District Municipality,Limpopo\nDC34,Vhembe District Municipality,Limpopo\nDC35,Capricorn District Municipality,Limpopo\nDC36,Waterberg District Municipality,Limpopo\nDC47,Sekhukhune District Municipality,Limpopo\nDC30,Gert Sibande District Municipality,Mpumalanga\nDC31,Nkangala District Municipality,Mpumalanga\nDC32,Ehlanzeni District Municipality,Mpumalanga\nDC37,Bojanala Platinum District Municipality,North West\nDC38,Ngaka Modiri Molema District Municipality,North West\nDC39,Dr Ruth Segomotsi Mompati District Municipality,North West\nDC40,Dr Kenneth Kaunda District Municipality,North West\nDC6,Namakwa District Municipality,Northern Cape\nDC7,Pixley ka Seme District Municipality,Northern Cape\nDC8,ZF Mgcawu District Municipality,Northern Cape\nDC9,Frances Baard District Municipality,Northern Cape\nDC45,John Taolo Gaetsewe District Municipality,Northern Cape\nCPT,City of Cape Town Metropolitan Municipality,Western Cape\nDC1,West Coast District Municipality,Western Cape\nDC2,Cape Winelands District Municipality,Western Cape\nDC3,Overberg District Municipality,Western Cape\nDC4,Garden Route District Municipality,Western Cape\nDC5,Central Karoo District Municipality,Western Cape\n",
  show_col_types = FALSE
)

offset <- 0.001
logit <- function(p) {
  log(p / (1 - p))
}
adjust_prob <- function(p, offset) {
  pmin(pmax(p, offset), 1 - offset)
}

output_crosswalk_dir <- file.path("data", "amr_data_bug_drug", "for_crosswalk")
if (!dir.exists(output_crosswalk_dir)) {
  dir.create(output_crosswalk_dir, recursive = TRUE)
}

process_bug_drug <- function(base_name) {
  level1_path <- file.path(output_dir, paste0(base_name, "_level1.csv"))
  without_path <- file.path(output_dir, paste0(base_name, "_without_level1.csv"))

  if (!file.exists(level1_path) || !file.exists(without_path)) {
    warning("Missing input files for: ", base_name)
    return(invisible(NULL))
  }

  df_level1 <- read_csv(level1_path, show_col_types = FALSE)
  df_without <- read_csv(without_path, show_col_types = FALSE)

  merged_df <- merge(
    df_level1,
    df_without,
    by = c("district_code", "year_taken"),
    all = TRUE,
    suffixes = c("_level1", "_without_level1")
  )

  base_cols <- c(
    "district_code",
    "year_taken",
    paste0(base_name, "_tested_level1"),
    paste0(base_name, "_resistant_level1"),
    paste0(base_name, "_resistance_per_level1"),
    paste0(base_name, "_tested_without_level1"),
    paste0(base_name, "_resistant_without_level1"),
    paste0(base_name, "_resistance_per_without_level1")
  )
  missing_cols <- setdiff(base_cols, names(merged_df))
  if (length(missing_cols) > 0) {
    warning("Missing expected columns for: ", base_name, " (", paste(missing_cols, collapse = ", "), ")")
    return(invisible(NULL))
  }

  out_df <- merged_df[, base_cols]
  names(out_df) <- c(
    "district_code",
    "year_taken",
    "level1_tested",
    "level1_resistant",
    "level1_resistance_per",
    "without_level1_tested",
    "without_level1_resistant",
    "without_level1_resistance_per"
  )

  out_df <- merge(
    out_df,
    district_province,
    by.x = "district_code",
    by.y = "District_Code",
    all.x = TRUE
  )

  out_df[is.na(out_df)] <- 0

  province_year_summary <- aggregate(
    cbind(
      level1_tested,
      level1_resistant,
      without_level1_tested,
      without_level1_resistant
    ) ~ Province + year_taken,
    data = out_df,
    FUN = sum
  )
  province_year_summary$level1_resistance_per <- with(
    province_year_summary,
    ifelse(level1_tested > 0, level1_resistant / level1_tested, 0)
  )

  province_keep <- province_year_summary[, c(
    "Province",
    "year_taken",
    "level1_tested",
    "level1_resistance_per"
  )]
  names(province_keep) <- c(
    "Province",
    "year_taken",
    "province_level1_tested",
    "province_level1_resistance_per"
  )
  out_df <- merge(
    out_df,
    province_keep,
    by = c("Province", "year_taken"),
    all.x = TRUE
  )

  national_year_summary <- aggregate(
    cbind(level1_tested, level1_resistant) ~ year_taken,
    data = out_df,
    FUN = sum
  )
  names(national_year_summary) <- c(
    "year_taken",
    "national_level1_tested",
    "national_level1_resistant"
  )
  out_df <- merge(
    out_df,
    national_year_summary,
    by = "year_taken",
    all.x = TRUE
  )
  out_df$national_level1_resistance_per <- with(
    out_df,
    ifelse(national_level1_tested > 0, national_level1_resistant / national_level1_tested, 0)
  )

  p_l1_adj <- adjust_prob(out_df$level1_resistance_per, offset)
  p_non_l1_adj <- adjust_prob(out_df$without_level1_resistance_per, offset)
  p_prov_adj <- adjust_prob(out_df$province_level1_resistance_per, offset)
  p_nat_adj <- adjust_prob(out_df$national_level1_resistance_per, offset)

  logit_l1 <- logit(p_l1_adj)
  logit_non_l1 <- logit(p_non_l1_adj)
  logit_prov <- logit(p_prov_adj)
  logit_nat <- logit(p_nat_adj)

  use_l1 <- out_df$level1_tested >= 30
  need_fallback <- !use_l1

  use_prov <- need_fallback & out_df$province_level1_tested >= 30
  use_nat <- need_fallback & !use_prov & out_df$national_level1_tested >= 30

  logit_corrected <- rep(NA_real_, nrow(out_df))
  source_corrected <- rep(NA_character_, nrow(out_df))
  r_corrected <- rep(NA_real_, nrow(out_df))
  n_corrected <- rep(NA_real_, nrow(out_df))

  delta_l1 <- logit_non_l1 - logit_l1
  delta_prov <- logit_non_l1 - logit_prov
  delta_nat <- logit_non_l1 - logit_nat

  logit_corrected[use_l1] <- logit_non_l1[use_l1] - delta_l1[use_l1]
  source_corrected[use_l1] <- "level1_delta"

  logit_corrected[use_prov] <- logit_non_l1[use_prov] - delta_prov[use_prov]
  source_corrected[use_prov] <- "province_level1_delta"

  logit_corrected[use_nat] <- logit_non_l1[use_nat] - delta_nat[use_nat]
  source_corrected[use_nat] <- "national_level1_delta"

  p_corrected <- ifelse(is.na(logit_corrected), NA_real_, 1 / (1 + exp(-logit_corrected)))
  use_corrected <- use_l1 | use_prov | use_nat
  r_corrected[use_corrected] <- p_corrected[use_corrected] * out_df$without_level1_tested[use_corrected]
  n_corrected[use_corrected] <- out_df$without_level1_tested[use_corrected]

  out_df$corrected_resistance_per <- p_corrected
  out_df$corrected_resistant <- r_corrected
  out_df$corrected_tested <- n_corrected
  out_df$source_corrected <- source_corrected

  total_tested <- rowSums(
    cbind(out_df$level1_tested, out_df$without_level1_tested),
    na.rm = TRUE
  )
  out_df <- out_df[total_tested >= 30, , drop = FALSE]

  final_cols <- c(
    "year_taken",
    "district_code",
    "corrected_resistance_per",
    "corrected_resistant",
    "corrected_tested"
  )
  out_df <- out_df[, final_cols, drop = FALSE]
  names(out_df) <- c(
    "year_taken",
    "district_code",
    paste0(base_name, "_resistance_per"),
    paste0(base_name, "_resistant"),
    paste0(base_name, "_tested")
  )

  output_path <- file.path(output_crosswalk_dir, paste0(base_name, "_final.csv"))
  write_csv(out_df, output_path)
  message("Wrote: ", output_path)
  print(head(out_df, 5))
}

invisible(lapply(bases, process_bug_drug))
