# Step 10: Evaluate Mahalanobis distance for applicability

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)
})

config <- list(
  prune_low_variance = TRUE,
  variance_thresh = 1e-8,
  prune_high_corr = TRUE,
  corr_thresh = 0.98,
  use_pca = TRUE,
  pca_var = 0.95,
  use_robust_cov = TRUE,
  stratify_by_year = TRUE,
  min_n_per_stratum = 10,
  plot_xmax = suppressWarnings(as.numeric(Sys.getenv("MAHALANOBIS_PLOT_XMAX", unset = NA))),
  plot_bin_width = suppressWarnings(as.numeric(Sys.getenv("MAHALANOBIS_PLOT_BIN_WIDTH", unset = 1)))
)

scenario_dir <- file.path("data", "data_for_modelling", "scenario_two")
eval_dir <- file.path("data", "evals")
if (!dir.exists(eval_dir)) dir.create(eval_dir, recursive = TRUE)

modelling_files <- list.files(
  scenario_dir,
  pattern = "_modelling\\.csv$",
  full.names = TRUE
)

if (length(modelling_files) == 0) {
  stop("No modelling files found in ", scenario_dir)
}

run_one <- function(modelling_path) {
  base_name <- sub("_modelling\\.csv$", "", basename(modelling_path))
  prediction_path <- file.path(scenario_dir, paste0(base_name, "_prediction.csv"))

  if (!file.exists(prediction_path)) {
    warning("Missing prediction file for ", base_name, "; skipping")
    return(invisible(NULL))
  }

  df_train <- read.csv(modelling_path, stringsAsFactors = FALSE, check.names = FALSE)
  df_pred <- read.csv(prediction_path, stringsAsFactors = FALSE, check.names = FALSE)

  message("Loaded ", base_name, ": train=", nrow(df_train), " pred=", nrow(df_pred))

  # -----------------------------
  # 1) Choose covariates (X only)
  # -----------------------------
  drop_cols <- c(
    "district_code",
    "year_taken",
    "year",
    "lon",
    "lat",
    "abaumanii_carbapenem_resistance_per",
    "abaumanii_carbapenem_resistant",
    "abaumanii_carbapenem_tested",
    # Removed due to low overlap between train and prediction
    "Medical practitioners per 100 000 population",
    "Human Development Index",
    "Delivery by Caesarean section rate (district hospitals)",
    "Gini coefficient",
    "No of households",
    "Mean years of schooling",
    "Adult ART Total",
    "Death under 5 years against live birth rate",
    "Hospital beds per 10 000 target population",
    "Population.y",
    "Antenatal client HIV 1st test positive",
    "Hospital beds per 10 000 target population (rescaled)",
    # Additional low-overlap features (top 15 after removals)
    "Acute malnutrition under five",
    "Maternal mortality",
    "Antenatal client initiated on ART rate",
    "HIV test positive 15 years and older (excl ANC)",
    "Pharmacists per 100 000 population",
    "Adult literacy",
    "Percentage of households with access to improved sanitation (index)",
    "ART Adult client viral load done (VLD)"
  )

  x_cols <- setdiff(names(df_train), drop_cols)

  # Sanity: ensure same covariates exist in both
  stopifnot(all(x_cols %in% names(df_pred)))

  # -----------------------------
  # 2) Build matrices
  # -----------------------------
  X_train <- as.matrix(df_train[, x_cols])
  X_pred <- as.matrix(df_pred[, x_cols])

  # -----------------------------
  # 2a) Prune low-variance + high-correlation features (train-based)
  # -----------------------------
  if (config$prune_low_variance) {
    v <- apply(X_train, 2, function(z) var(z, na.rm = TRUE))
    keep <- !(is.na(v) | v <= config$variance_thresh)
    X_train <- X_train[, keep, drop = FALSE]
    X_pred <- X_pred[, keep, drop = FALSE]
    x_cols <- x_cols[keep]
  }

  if (config$prune_high_corr && ncol(X_train) > 1) {
    cmat <- suppressWarnings(cor(X_train, use = "pairwise.complete.obs"))
    if (any(is.na(cmat))) cmat[is.na(cmat)] <- 0
    drop <- rep(FALSE, ncol(cmat))
    for (i in seq_len(ncol(cmat) - 1)) {
      if (drop[i]) next
      for (j in (i + 1):ncol(cmat)) {
        if (!drop[j] && abs(cmat[i, j]) >= config$corr_thresh) {
          drop[j] <- TRUE
        }
      }
    }
    keep <- !drop
    X_train <- X_train[, keep, drop = FALSE]
    X_pred <- X_pred[, keep, drop = FALSE]
    x_cols <- x_cols[keep]
  }

  # -----------------------------
  # 3) Impute missing values (TRAIN medians)
  # -----------------------------
  train_medians <- apply(X_train, 2, function(z) median(z, na.rm = TRUE))

  impute_with <- function(X, med) {
    for (j in seq_len(ncol(X))) {
      nas <- is.na(X[, j])
      if (any(nas)) X[nas, j] <- med[j]
    }
    X
  }

  X_train <- impute_with(X_train, train_medians)
  X_pred <- impute_with(X_pred, train_medians)

  # -----------------------------
  # 4) Standardize using TRAIN mean/sd
  # -----------------------------
  train_mean <- colMeans(X_train)
  train_sd <- apply(X_train, 2, sd)
  train_sd[train_sd == 0] <- 1

  scale_with <- function(X, mu, s) {
    sweep(sweep(X, 2, mu, "-"), 2, s, "/")
  }

  X_train_s <- scale_with(X_train, train_mean, train_sd)
  X_pred_s <- scale_with(X_pred, train_mean, train_sd)

  # -----------------------------
  # 4a) Optional PCA (train-based)
  # -----------------------------
  if (config$use_pca && ncol(X_train_s) > 1) {
    pca <- prcomp(X_train_s, center = FALSE, scale. = FALSE)
    var_prop <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
    k <- which(var_prop >= config$pca_var)[1]
    X_train_s <- pca$x[, seq_len(k), drop = FALSE]
    X_pred_s <- predict(pca, newdata = X_pred_s)[, seq_len(k), drop = FALSE]
  }

  # -----------------------------
  # 5) Covariance inverse (regularized) + Mahalanobis distances
  # -----------------------------
  robust_cov_safe <- function(X) {
    n_min <- 2 * ncol(X) + 1
    if (nrow(X) < n_min) {
      return(cov(X))
    }
    out <- tryCatch(cov.rob(X)$cov, error = function(e) cov(X))
    out
  }

  S <- if (config$use_robust_cov) robust_cov_safe(X_train_s) else cov(X_train_s)

  lambda <- 1e-6
  S_reg <- S + diag(lambda, ncol(S))

  S_inv <- ginv(S_reg)
  mu <- colMeans(X_train_s)

  mahal_dist <- function(X, mu, S_inv) {
    apply(X, 1, function(x) {
      sqrt(t(x - mu) %*% S_inv %*% (x - mu))
    })
  }

  df_pred$mahalanobis_dist <- as.numeric(mahal_dist(X_pred_s, mu, S_inv))
  df_train$mahalanobis_dist <- as.numeric(mahal_dist(X_train_s, mu, S_inv))

  # -----------------------------
  # 5a) Optional stratified distances by year
  # -----------------------------
  year_var <- if ("year_taken" %in% names(df_train)) "year_taken" else "year"
  if (config$stratify_by_year && year_var %in% names(df_train) && year_var %in% names(df_pred)) {
    years <- sort(unique(df_train[[year_var]]))
    df_train$mahalanobis_dist_strat <- NA_real_
    df_pred$mahalanobis_dist_strat <- NA_real_

    for (yr in years) {
      idx_train <- which(df_train[[year_var]] == yr)
      if (length(idx_train) < config$min_n_per_stratum) next

      Xtr <- X_train_s[idx_train, , drop = FALSE]
      S_yr <- if (config$use_robust_cov) robust_cov_safe(Xtr) else cov(Xtr)
      S_yr <- S_yr + diag(lambda, ncol(S_yr))
      S_yr_inv <- ginv(S_yr)
      mu_yr <- colMeans(Xtr)

      df_train$mahalanobis_dist_strat[idx_train] <- as.numeric(mahal_dist(Xtr, mu_yr, S_yr_inv))

      idx_pred <- which(df_pred[[year_var]] == yr)
      if (length(idx_pred) > 0) {
        Xpr <- X_pred_s[idx_pred, , drop = FALSE]
        df_pred$mahalanobis_dist_strat[idx_pred] <- as.numeric(mahal_dist(Xpr, mu_yr, S_yr_inv))
      }
    }

    df_train$mahalanobis_dist_strat[is.na(df_train$mahalanobis_dist_strat)] <- df_train$mahalanobis_dist[is.na(df_train$mahalanobis_dist_strat)]
    df_pred$mahalanobis_dist_strat[is.na(df_pred$mahalanobis_dist_strat)] <- df_pred$mahalanobis_dist[is.na(df_pred$mahalanobis_dist_strat)]
  }

  # -----------------------------
  # 6) Thresholds based on training distribution (optional)
  # -----------------------------
  cut90 <- as.numeric(quantile(df_train$mahalanobis_dist, 0.90, na.rm = TRUE))
  cut95 <- as.numeric(quantile(df_train$mahalanobis_dist, 0.95, na.rm = TRUE))

  df_pred$applicability_flag <- dplyr::case_when(
    df_pred$mahalanobis_dist <= cut90 ~ "Within training support",
    df_pred$mahalanobis_dist <= cut95 ~ "Near edge",
    TRUE ~ "Out-of-support (extrapolation)"
  )

  if ("mahalanobis_dist_strat" %in% names(df_pred)) {
    cut90_s <- as.numeric(quantile(df_train$mahalanobis_dist_strat, 0.90, na.rm = TRUE))
    cut95_s <- as.numeric(quantile(df_train$mahalanobis_dist_strat, 0.95, na.rm = TRUE))
    df_pred$applicability_flag_strat <- dplyr::case_when(
      df_pred$mahalanobis_dist_strat <= cut90_s ~ "Within training support",
      df_pred$mahalanobis_dist_strat <= cut95_s ~ "Near edge",
      TRUE ~ "Out-of-support (extrapolation)"
    )
  }

  # -----------------------------
  # 7) District-level summaries
  # -----------------------------
  district_applicability <- df_pred %>%
    group_by(district_code) %>%
    summarise(
      md_median = median(mahalanobis_dist, na.rm = TRUE),
      md_p90 = quantile(mahalanobis_dist, 0.90, na.rm = TRUE),
      md_max = max(mahalanobis_dist, na.rm = TRUE),
      worst_flag = applicability_flag[which.max(mahalanobis_dist)],
      .groups = "drop"
    ) %>%
    arrange(desc(md_median))

  bug_eval_dir <- file.path(eval_dir, base_name)
  if (!dir.exists(bug_eval_dir)) dir.create(bug_eval_dir, recursive = TRUE)

  write.csv(
    df_pred,
    file.path(bug_eval_dir, paste0(base_name, "_prediction_with_mahalanobis.csv")),
    row.names = FALSE
  )
  write.csv(
    district_applicability,
    file.path(bug_eval_dir, paste0(base_name, "_district_applicability.csv")),
    row.names = FALSE
  )

  # -----------------------------
  # 8) Mahalanobis distance chart
  # -----------------------------
  chart_path <- file.path(bug_eval_dir, paste0(base_name, "_mahalanobis_distance.png"))
  png(chart_path, width = 1200, height = 800, res = 150)
  par(mar = c(5, 5, 4, 2) + 0.1)

  bug_drug_label <- base_name
  bug_drug_label <- gsub("_final_withcrosswalk_merged$", "", bug_drug_label)
  bug_drug_label <- gsub("_final_withcrosswalk$", "", bug_drug_label)
  bug_drug_label <- gsub("_merged$", "", bug_drug_label)
  bug_drug_label <- gsub("_", " ", bug_drug_label)
  bug_drug_label <- paste0(toupper(substr(bug_drug_label, 1, 1)), substr(bug_drug_label, 2, nchar(bug_drug_label)))

  if (!is.na(config$plot_xmax) && config$plot_xmax > 0) {
    break_max <- max(config$plot_xmax, ceiling(max(df_pred$mahalanobis_dist, na.rm = TRUE)))
    hist_breaks <- seq(0, break_max + config$plot_bin_width, by = config$plot_bin_width)
    hist_xlim <- c(0, config$plot_xmax)
  } else {
    hist_breaks <- 30
    hist_xlim <- NULL
  }

  hist_args <- list(
    x = df_pred$mahalanobis_dist,
    breaks = hist_breaks,
    col = "#4C78A8",
    border = "white",
    main = "",
    xlab = "Mahalanobis distance",
    ylab = "Number of modelled district-years"
  )
  if (!is.null(hist_xlim)) hist_args$xlim <- hist_xlim
  do.call(hist, hist_args)
  axis(1)
  axis(2)
  box()
  abline(v = cut95, col = "#E45756", lwd = 2, lty = 2)
  legend(
    "topright",
    legend = "95th percentile (training data)",
    col = "#E45756",
    lty = 2,
    lwd = 2,
    bty = "n"
  )
  dev.off()
  message("Saved chart to ", chart_path)

  message("Saved outputs for ", base_name, " to ", bug_eval_dir)
}

for (modelling_path in modelling_files) {
  run_one(modelling_path)
}
