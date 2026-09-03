#!/usr/bin/env Rscript
# Step 9.1: Train models for scenario_one (all levels above threshold)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(caret)
  library(randomForest)
  library(xgboost)
  library(gbm)
  library(Cubist)
  library(nnet)
  library(kernlab)
  library(mgcv)
  library(ggplot2)
  library(stringr)
  library(sf)
  library(cowplot)
  library(geodata)
})

set.seed(2025)

args <- commandArgs(trailingOnly = TRUE)
default_scenario_dir <- file.path("data", "data_for_modelling", "scenario_two")
default_model_root <- file.path("data", "models", "scenario_crosswalkeddata")

scenario_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else if (nzchar(Sys.getenv("SCENARIO_DIR"))) {
  Sys.getenv("SCENARIO_DIR")
} else {
  default_scenario_dir
}

model_root <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else if (nzchar(Sys.getenv("MODEL_ROOT"))) {
  Sys.getenv("MODEL_ROOT")
} else {
  default_model_root
}

if (!dir.exists(scenario_dir)) {
  stop("Scenario directory does not exist: ", scenario_dir)
}
if (!dir.exists(model_root)) {
  dir.create(model_root, recursive = TRUE)
}

features <- c(
  "year", "lon", "lat",
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

train_pattern <- "_modelling_training\\.csv$"

train_files <- list.files(
  scenario_dir,
  pattern = train_pattern,
  full.names = TRUE
)

if (length(train_files) == 0) {
  message("No training files found in ", scenario_dir)
  quit(status = 0)
}

ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)

for (train_path in train_files) {
  test_path <- sub("_training\\.csv$", "_testing.csv", train_path)
  pred_path <- sub("_modelling_training(_districtsplit)?\\.csv$", "_prediction.csv", train_path)
  if (!file.exists(test_path)) {
    warning("Testing file missing for ", basename(train_path), "; skipping")
    next
  }

  train_df <- read.csv(train_path, stringsAsFactors = FALSE, check.names = FALSE)
  test_df  <- read.csv(test_path, stringsAsFactors = FALSE, check.names = FALSE)
  pred_df  <- if (file.exists(pred_path)) read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE) else NULL

  target_candidates <- grep("_resistance_per$", colnames(train_df), value = TRUE)
  if (length(target_candidates) == 0) {
    warning("No *_resistance_per column in ", basename(train_path), "; skipping")
    next
  }
  target <- target_candidates[1]

  missing_cols <- setdiff(c(target, features), colnames(train_df))
  if (length(missing_cols) > 0) {
    warning("Missing columns in training data for ", basename(train_path), ": ", paste(missing_cols, collapse = ", "))
    next
  }

  train_df <- train_df[, c(target, features), drop = FALSE] %>%
    filter(!is.na(.data[[target]]))

  test_df <- test_df[, c(target, features), drop = FALSE]
  if (!is.null(pred_df)) {
    pred_df <- pred_df[, c("district_code", target, features), drop = FALSE]
  }

  metrics_rows <- list()

  fm <- as.formula(
    paste0("`", target, "` ~ ", paste(paste0("`", features, "`"), collapse = " + "))
  )

  # Random Forest
  rf_mod <- train(
    fm,
    data       = train_df,
    method     = "rf",
    trControl  = ctrl,
    tuneLength = 5
  )

  train_pred  <- predict(rf_mod, train_df)
  in_metrics  <- postResample(pred = train_pred, obs = train_df[[target]])

  test_pred   <- predict(rf_mod, test_df)
  out_metrics <- postResample(pred = test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> RF: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    in_metrics["Rsquared"], in_metrics["RMSE"],
    out_metrics["Rsquared"], out_metrics["RMSE"]
  ))

  bug_name <- sub("_modelling_training$", "", tools::file_path_sans_ext(basename(train_path)))
  bug_dir <- file.path(model_root, bug_name, "randomforest")
  if (!dir.exists(bug_dir)) dir.create(bug_dir, recursive = TRUE)

  model_path <- file.path(bug_dir, paste0(bug_name, "_rf_model.rds"))
  metrics_path <- file.path(bug_dir, paste0(bug_name, "_rf_metrics.csv"))
  varimp_path <- file.path(bug_dir, paste0(bug_name, "_rf_var_importance.csv"))

  saveRDS(rf_mod, model_path)

  metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(in_metrics["Rsquared"], out_metrics["Rsquared"]),
    rmse = c(in_metrics["RMSE"], out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(metrics_df, metrics_path, row.names = FALSE)
  rf_imp <- tryCatch(caret::varImp(rf_mod), error = function(e) NULL)
  if (!is.null(rf_imp)) {
    rf_imp_df <- as.data.frame(rf_imp$importance)
    rf_imp_df <- data.frame(variable = rownames(rf_imp_df), rf_imp_df, row.names = NULL)
    write.csv(rf_imp_df, varimp_path, row.names = FALSE)
  }
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "rf",
    dataset = c("train", "test"),
    rsq = c(in_metrics["Rsquared"], out_metrics["Rsquared"]),
    rmse = c(in_metrics["RMSE"], out_metrics["RMSE"])
  )

  # XGBoost
  xgb_mod <- suppressWarnings(
    train(
      fm,
      data       = train_df,
      method     = "xgbTree",
      trControl  = ctrl,
      tuneLength = 5,
      verbose    = 0
    )
  )

  xgb_train_pred  <- suppressWarnings(predict(xgb_mod, train_df))
  xgb_in_metrics  <- postResample(pred = xgb_train_pred, obs = train_df[[target]])

  xgb_test_pred   <- suppressWarnings(predict(xgb_mod, test_df))
  xgb_out_metrics <- postResample(pred = xgb_test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> XGB: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    xgb_in_metrics["Rsquared"], xgb_in_metrics["RMSE"],
    xgb_out_metrics["Rsquared"], xgb_out_metrics["RMSE"]
  ))

  xgb_dir <- file.path(model_root, bug_name, "xgboost")
  if (!dir.exists(xgb_dir)) dir.create(xgb_dir, recursive = TRUE)

  xgb_model_path <- file.path(xgb_dir, paste0(bug_name, "_xgb_model.rds"))
  xgb_metrics_path <- file.path(xgb_dir, paste0(bug_name, "_xgb_metrics.csv"))

  saveRDS(xgb_mod, xgb_model_path)

  xgb_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(xgb_in_metrics["Rsquared"], xgb_out_metrics["Rsquared"]),
    rmse = c(xgb_in_metrics["RMSE"], xgb_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(xgb_metrics_df, xgb_metrics_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "xgb",
    dataset = c("train", "test"),
    rsq = c(xgb_in_metrics["Rsquared"], xgb_out_metrics["Rsquared"]),
    rmse = c(xgb_in_metrics["RMSE"], xgb_out_metrics["RMSE"])
  )

  # GBM
  gbm_mod <- train(
    fm,
    data       = train_df,
    method     = "gbm",
    trControl  = ctrl,
    tuneLength = 5,
    verbose    = FALSE
  )

  gbm_train_pred  <- predict(gbm_mod, train_df)
  gbm_in_metrics  <- postResample(pred = gbm_train_pred, obs = train_df[[target]])

  gbm_test_pred   <- predict(gbm_mod, test_df)
  gbm_out_metrics <- postResample(pred = gbm_test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> GBM: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    gbm_in_metrics["Rsquared"], gbm_in_metrics["RMSE"],
    gbm_out_metrics["Rsquared"], gbm_out_metrics["RMSE"]
  ))

  gbm_dir <- file.path(model_root, bug_name, "gbm")
  if (!dir.exists(gbm_dir)) dir.create(gbm_dir, recursive = TRUE)

  gbm_model_path <- file.path(gbm_dir, paste0(bug_name, "_gbm_model.rds"))
  gbm_metrics_path <- file.path(gbm_dir, paste0(bug_name, "_gbm_metrics.csv"))

  saveRDS(gbm_mod, gbm_model_path)

  gbm_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(gbm_in_metrics["Rsquared"], gbm_out_metrics["Rsquared"]),
    rmse = c(gbm_in_metrics["RMSE"], gbm_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(gbm_metrics_df, gbm_metrics_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "gbm",
    dataset = c("train", "test"),
    rsq = c(gbm_in_metrics["Rsquared"], gbm_out_metrics["Rsquared"]),
    rmse = c(gbm_in_metrics["RMSE"], gbm_out_metrics["RMSE"])
  )

  # Cubist
  cubist_mod <- train(
    fm,
    data       = train_df,
    method     = "cubist",
    trControl  = ctrl,
    tuneLength = 5
  )

  cubist_train_pred  <- predict(cubist_mod, train_df)
  cubist_in_metrics  <- postResample(pred = cubist_train_pred, obs = train_df[[target]])

  cubist_test_pred   <- predict(cubist_mod, test_df)
  cubist_out_metrics <- postResample(pred = cubist_test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> Cubist: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    cubist_in_metrics["Rsquared"], cubist_in_metrics["RMSE"],
    cubist_out_metrics["Rsquared"], cubist_out_metrics["RMSE"]
  ))

  cubist_dir <- file.path(model_root, bug_name, "cubist")
  if (!dir.exists(cubist_dir)) dir.create(cubist_dir, recursive = TRUE)

  cubist_model_path <- file.path(cubist_dir, paste0(bug_name, "_cubist_model.rds"))
  cubist_metrics_path <- file.path(cubist_dir, paste0(bug_name, "_cubist_metrics.csv"))

  saveRDS(cubist_mod, cubist_model_path)

  cubist_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(cubist_in_metrics["Rsquared"], cubist_out_metrics["Rsquared"]),
    rmse = c(cubist_in_metrics["RMSE"], cubist_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(cubist_metrics_df, cubist_metrics_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "cubist",
    dataset = c("train", "test"),
    rsq = c(cubist_in_metrics["Rsquared"], cubist_out_metrics["Rsquared"]),
    rmse = c(cubist_in_metrics["RMSE"], cubist_out_metrics["RMSE"])
  )

  # Neural network (single hidden layer via nnet)
  nnet_mod <- train(
    fm,
    data       = train_df,
    method     = "nnet",
    trControl  = ctrl,
    tuneLength = 5,
    trace      = FALSE
  )

  nnet_train_pred  <- predict(nnet_mod, train_df)
  nnet_in_metrics  <- postResample(pred = nnet_train_pred, obs = train_df[[target]])

  nnet_test_pred   <- predict(nnet_mod, test_df)
  nnet_out_metrics <- postResample(pred = nnet_test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> NNET: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    nnet_in_metrics["Rsquared"], nnet_in_metrics["RMSE"],
    nnet_out_metrics["Rsquared"], nnet_out_metrics["RMSE"]
  ))

  nnet_dir <- file.path(model_root, bug_name, "nnet")
  if (!dir.exists(nnet_dir)) dir.create(nnet_dir, recursive = TRUE)

  nnet_model_path <- file.path(nnet_dir, paste0(bug_name, "_nnet_model.rds"))
  nnet_metrics_path <- file.path(nnet_dir, paste0(bug_name, "_nnet_metrics.csv"))

  saveRDS(nnet_mod, nnet_model_path)

  nnet_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(nnet_in_metrics["Rsquared"], nnet_out_metrics["Rsquared"]),
    rmse = c(nnet_in_metrics["RMSE"], nnet_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(nnet_metrics_df, nnet_metrics_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "nnet",
    dataset = c("train", "test"),
    rsq = c(nnet_in_metrics["Rsquared"], nnet_out_metrics["Rsquared"]),
    rmse = c(nnet_in_metrics["RMSE"], nnet_out_metrics["RMSE"])
  )

  # SVM (radial kernel)
  svm_mod <- train(
    fm,
    data       = train_df,
    method     = "svmRadial",
    trControl  = ctrl,
    tuneLength = 5
  )

  svm_train_pred  <- predict(svm_mod, train_df)
  svm_in_metrics  <- postResample(pred = svm_train_pred, obs = train_df[[target]])

  svm_test_pred   <- predict(svm_mod, test_df)
  svm_out_metrics <- postResample(pred = svm_test_pred, obs = test_df[[target]])

  cat(sprintf(
    "%s -> SVM: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    svm_in_metrics["Rsquared"], svm_in_metrics["RMSE"],
    svm_out_metrics["Rsquared"], svm_out_metrics["RMSE"]
  ))

  svm_dir <- file.path(model_root, bug_name, "svm")
  if (!dir.exists(svm_dir)) dir.create(svm_dir, recursive = TRUE)

  svm_model_path <- file.path(svm_dir, paste0(bug_name, "_svm_model.rds"))
  svm_metrics_path <- file.path(svm_dir, paste0(bug_name, "_svm_metrics.csv"))

  saveRDS(svm_mod, svm_model_path)

  svm_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(svm_in_metrics["Rsquared"], svm_out_metrics["Rsquared"]),
    rmse = c(svm_in_metrics["RMSE"], svm_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(svm_metrics_df, svm_metrics_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "svm",
    dataset = c("train", "test"),
    rsq = c(svm_in_metrics["Rsquared"], svm_out_metrics["Rsquared"]),
    rmse = c(svm_in_metrics["RMSE"], svm_out_metrics["RMSE"])
  )

  # GAM (spatial + temporal smooths)
  # Clean names for mgcv: make syntactic, then replace non-alnum with underscores
  clean_names <- function(nms) {
    nms <- make.names(nms, unique = TRUE)
    gsub("[^A-Za-z0-9_]", "_", nms)
  }

  train_gam <- train_df
  test_gam  <- test_df
  names(train_gam) <- clean_names(names(train_gam))
  names(test_gam)  <- clean_names(names(test_gam))

  # Align target/feature names in cleaned data
  target_gam   <- clean_names(target)
  feature_gam  <- clean_names(features)

  preds <- setdiff(names(train_gam), target_gam)
  linear_preds <- setdiff(preds, c("lon", "lat", "year"))

  n_spatial <- nrow(distinct(train_gam, lon, lat))
  n_time    <- length(unique(train_gam$year))
  k_space   <- max(3, min(10, floor(n_spatial / 2)))
  k_time    <- max(3, min(6,  floor(n_time    / 2)))

  gam_formula <- as.formula(
    paste0(
      target_gam, " ~ ",
      sprintf("te(lon, lat, bs = c('tp','tp'), k = c(%d,%d)) + ", k_space, k_space),
      sprintf("s(year, bs = 'cr', k = %d)", k_time),
      if (length(linear_preds) > 0) paste0(" + ", paste(linear_preds, collapse = " + ")) else ""
    )
  )

  set.seed(2025)
  gam_mod <- bam(
    gam_formula,
    data     = train_gam,
    method   = "fREML",
    discrete = TRUE
  )

  gam_train_pred  <- predict(gam_mod, newdata = train_gam)
  gam_in_metrics  <- postResample(pred = gam_train_pred, obs = train_gam[[target_gam]])

  gam_test_pred   <- predict(gam_mod, newdata = test_gam)
  gam_out_metrics <- postResample(pred = gam_test_pred, obs = test_gam[[target_gam]])

  cat(sprintf(
    "%s -> GAM: In-sample R²=%.3f RMSE=%.3f | Out-of-sample R²=%.3f RMSE=%.3f\n",
    basename(train_path),
    gam_in_metrics["Rsquared"], gam_in_metrics["RMSE"],
    gam_out_metrics["Rsquared"], gam_out_metrics["RMSE"]
  ))

  gam_dir <- file.path(model_root, bug_name, "gam")
  if (!dir.exists(gam_dir)) dir.create(gam_dir, recursive = TRUE)

  gam_model_path <- file.path(gam_dir, paste0(bug_name, "_gam_model.rds"))
  gam_metrics_path <- file.path(gam_dir, paste0(bug_name, "_gam_metrics.csv"))
  gam_varimp_path <- file.path(gam_dir, paste0(bug_name, "_gam_var_importance.csv"))

  saveRDS(gam_mod, gam_model_path)

  gam_metrics_df <- data.frame(
    dataset = c("train", "test"),
    rsq = c(gam_in_metrics["Rsquared"], gam_out_metrics["Rsquared"]),
    rmse = c(gam_in_metrics["RMSE"], gam_out_metrics["RMSE"]),
    row.names = NULL
  )
  write.csv(gam_metrics_df, gam_metrics_path, row.names = FALSE)
  gam_smry <- summary(gam_mod)
  gam_smooth <- if (!is.null(gam_smry$s.table)) {
    data.frame(
      term = rownames(gam_smry$s.table),
      edf = gam_smry$s.table[, "edf"],
      ref_df = gam_smry$s.table[, "Ref.df"],
      F = gam_smry$s.table[, "F"],
      p_value = gam_smry$s.table[, "p-value"],
      row.names = NULL
    )
  } else {
    data.frame()
  }
  gam_param <- if (!is.null(gam_smry$p.table)) {
    data.frame(
      term = rownames(gam_smry$p.table),
      estimate = gam_smry$p.table[, "Estimate"],
      std_error = gam_smry$p.table[, "Std. Error"],
      t_value = gam_smry$p.table[, "t value"],
      p_value = gam_smry$p.table[, "Pr(>|t|)"],
      row.names = NULL
    )
  } else {
    data.frame()
  }
gam_imp <- dplyr::bind_rows(
  dplyr::mutate(gam_smooth, type = "smooth"),
  dplyr::mutate(gam_param,  type = "param")
)
  write.csv(gam_imp, gam_varimp_path, row.names = FALSE)
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    model = "gam",
    dataset = c("train", "test"),
    rsq = c(gam_in_metrics["Rsquared"], gam_out_metrics["Rsquared"]),
    rmse = c(gam_in_metrics["RMSE"], gam_out_metrics["RMSE"])
  )

  # ------------------------------------------------------------------
  # Stacked meta-data (model predictions + coordinates/time + target)
  # ------------------------------------------------------------------
  year_var <- if ("year_taken" %in% names(train_df)) "year_taken" else "year"

  meta_train <- data.frame(
    rf     = predict(rf_mod,    train_df),
    xgb    = predict(xgb_mod,   train_df),
    gbm    = predict(gbm_mod,   train_df),
    cubist = predict(cubist_mod, train_df),
    svm    = predict(svm_mod,   train_df),
    gam    = predict(gam_mod,   newdata = train_gam),
    lon    = train_df$lon,
    lat    = train_df$lat,
    year   = train_df[[year_var]],
    y      = train_df[[target]]
  )

  meta_test <- data.frame(
    rf     = predict(rf_mod,    test_df),
    xgb    = predict(xgb_mod,   test_df),
    gbm    = predict(gbm_mod,   test_df),
    cubist = predict(cubist_mod, test_df),
    svm    = predict(svm_mod,   test_df),
    gam    = predict(gam_mod,   newdata = test_gam),
    lon    = test_df$lon,
    lat    = test_df$lat,
    year   = test_df[[year_var]],
    y      = test_df[[target]]
  )

  stack_dir <- file.path(model_root, bug_name, "stack")
  if (!dir.exists(stack_dir)) dir.create(stack_dir, recursive = TRUE)
  write.csv(meta_train, file.path(stack_dir, paste0(bug_name, "_meta_train.csv")), row.names = FALSE)
  write.csv(meta_test,  file.path(stack_dir, paste0(bug_name, "_meta_test.csv")),  row.names = FALSE)
  meta_pred <- NULL
  if (!is.null(pred_df)) {
    meta_pred <- data.frame(
      rf     = predict(rf_mod,    pred_df),
      xgb    = predict(xgb_mod,   pred_df),
      gbm    = predict(gbm_mod,   pred_df),
      cubist = predict(cubist_mod, pred_df),
      svm    = predict(svm_mod,   pred_df),
      gam    = predict(gam_mod,   newdata = {
        tmp <- pred_df
        names(tmp) <- clean_names(names(tmp))
        tmp
      }),
      lon    = pred_df$lon,
      lat    = pred_df$lat,
      year   = pred_df[[year_var]],
      y      = pred_df[[target]]
    )
    write.csv(meta_pred, file.path(stack_dir, paste0(bug_name, "_meta_pred.csv")), row.names = FALSE)
  }

  combined_metrics <- bind_rows(metrics_rows)
  combined_metrics$model <- factor(combined_metrics$model, levels = unique(combined_metrics$model))
  write.csv(combined_metrics, file.path(model_root, bug_name, paste0(bug_name, "_all_models_metrics.csv")), row.names = FALSE)

  test_metrics <- combined_metrics %>%
    filter(dataset == "test")

  pareto_front <- function(df) {
    keep <- rep(TRUE, nrow(df))
    for (i in seq_len(nrow(df))) {
      for (j in seq_len(nrow(df))) {
        if (i == j) next
        if (df$rsq[j] >= df$rsq[i] &&
            df$rmse[j] <= df$rmse[i] &&
            (df$rsq[j] > df$rsq[i] || df$rmse[j] < df$rmse[i])) {
          keep[i] <- FALSE
          break
        }
      }
    }
    df[keep, , drop = FALSE]
  }

  pareto_models <- pareto_front(test_metrics)

  top3 <- pareto_models %>%
    arrange(desc(rsq), rmse) %>%
    slice_head(n = 3)

  write.csv(top3, file.path(model_root, bug_name, paste0(bug_name, "_top3_models_by_test_rsq.csv")), row.names = FALSE)

  # ------------------------------------------------------------------
  # ST-GPR using top 3 models from stacked data (if DiceKriging available)
  # ------------------------------------------------------------------
  if (requireNamespace("DiceKriging", quietly = TRUE)) {
    top_models <- as.character(top3$model)
    # Ensure the columns exist in the meta data
    top_models <- top_models[top_models %in% colnames(meta_train)]

    if (length(top_models) > 0) {
      inputs <- c(as.character(top_models), "lon", "lat", if ("year_taken" %in% names(train_df)) "year_taken" else "year")

      Xtr <- as.matrix(meta_train[, inputs, drop = FALSE])
      ytr <- meta_train$y
      Xte <- as.matrix(meta_test[, inputs, drop = FALSE])
      yte <- meta_test$y

      set.seed(2025)
      stgpr_mod <- DiceKriging::km(
        formula   = ~ .,
        design    = Xtr,
        response  = ytr,
        covtype   = "matern5_2",
        control   = list(parinit = "ARD")
      )

      pred <- predict(stgpr_mod, newdata = Xte, type = "SK")
      yhat <- pred$mean

      stgpr_metrics <- postResample(pred = yhat, obs = yte)
      cat(sprintf(
        "%s -> ST-GPR on (%s): R²=%.3f RMSE=%.3f\n",
        basename(train_path),
        paste(top_models, collapse = ", "),
        stgpr_metrics["Rsquared"], stgpr_metrics["RMSE"]
      ))

      stgpr_dir <- file.path(model_root, bug_name, "stgpr")
      if (!dir.exists(stgpr_dir)) dir.create(stgpr_dir, recursive = TRUE)

      saveRDS(stgpr_mod, file.path(stgpr_dir, paste0(bug_name, "_stgpr_model.rds")))

      stgpr_metrics_df <- data.frame(
        model = paste(top_models, collapse = ", "),
        rsq = stgpr_metrics["Rsquared"],
        rmse = stgpr_metrics["RMSE"],
        row.names = NULL
      )
      write.csv(stgpr_metrics_df, file.path(stgpr_dir, paste0(bug_name, "_stgpr_metrics.csv")), row.names = FALSE)

      # Plot predicted vs actual on test set
      plot_df <- data.frame(actual = yte, predicted = yhat)
      parts <- strsplit(bug_name, "_")[[1]]
      core_name <- paste(parts[1:min(2, length(parts))], collapse = " ")
      bug_title <- str_to_title(core_name)
      p <- ggplot(plot_df, aes(x = actual, y = predicted)) +
        geom_point(alpha = 0.6, color = "blue") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        labs(
          title = bug_title,
          x = "Observed prevalence",
          y = "Modelled prevalence"
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        theme_bw() +
        theme(panel.grid = element_blank())

      ggsave(
        filename = file.path(stgpr_dir, paste0(bug_name, "_stgpr_pred_vs_actual.png")),
        plot = p, width = 6, height = 5, dpi = 300, device = "png"
      )
      message("Saved ST-GPR predicted vs actual plot to ", file.path(stgpr_dir, paste0(bug_name, "_stgpr_pred_vs_actual.png")))
    } else {
      warning("Top model columns not found in meta data for ", basename(train_path), "; skipping ST-GPR")
    }
  } else {
    warning("DiceKriging not installed; skipping ST-GPR for ", basename(train_path))
  }

  # ST-GPR prediction on full prediction set if available and model was fitted
  if (exists("stgpr_mod") && !is.null(pred_df) && !is.null(meta_pred)) {
    if (exists("top_models") && length(top_models) > 0) {
      pred_inputs <- c(as.character(top_models), "lon", "lat", year_var)
      Xpred <- as.matrix(meta_pred[, pred_inputs, drop = FALSE])
      stgpr_pred <- predict(stgpr_mod, newdata = Xpred, type = "SK")$mean
      pred_out <- pred_df
      pred_out$stgpr_pred <- stgpr_pred
      if ("district_code" %in% names(pred_df)) {
        pred_out$district_code <- pred_df$district_code
      }

      pred_save_dir <- file.path(model_root, bug_name, "stgpr")
      if (!dir.exists(pred_save_dir)) dir.create(pred_save_dir, recursive = TRUE)
      pred_save_path <- file.path(pred_save_dir, paste0(bug_name, "_prediction_with_stgpr.csv"))
      write.csv(pred_out, pred_save_path, row.names = FALSE)

      message("Wrote ST-GPR predictions into ", pred_save_path, " (original prediction file left unchanged)")

      # Prediction maps by year
      map_dir <- file.path(pred_save_dir, "predictionmaps")
      if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)

      # Try to load South Africa districts (GADM level 2)
      gadm_path <- NULL
      if (file.exists("data/spatial/gadm41_ZAF_2_sf.gpkg")) gadm_path <- "data/spatial/gadm41_ZAF_2_sf.gpkg"
      if (is.null(gadm_path) && file.exists("data/spatial/gadm41_ZAF_2_sf.rds")) gadm_path <- "data/spatial/gadm41_ZAF_2_sf.rds"

      if (is.null(gadm_path)) {
        dir.create("data/spatial", recursive = TRUE, showWarnings = FALSE)
        message("GADM not found locally; attempting download via geodata::gadm ...")
        try({
          gadm_v <- geodata::gadm("ZAF", level = 2, path = "data/spatial")
          sa_gadm2_tmp <- sf::st_as_sf(gadm_v)
          sf::st_write(sa_gadm2_tmp, "data/spatial/gadm41_ZAF_2_sf.gpkg", quiet = TRUE, delete_dsn = TRUE)
          saveRDS(sa_gadm2_tmp, "data/spatial/gadm41_ZAF_2_sf.rds")
          gadm_path <- "data/spatial/gadm41_ZAF_2_sf.gpkg"
        }, silent = TRUE)
      }

      if (!is.null(gadm_path)) {
        sa_gadm2 <- if (grepl("\\.rds$", gadm_path)) readRDS(gadm_path) else sf::st_read(gadm_path, quiet = TRUE)

        if ("CC_2" %in% names(sa_gadm2) && all(c("district_code", "year", "stgpr_pred") %in% names(pred_out))) {
          years <- sort(unique(pred_out$year))
          plot_list <- list()

          for (yr in years) {
            df_year <- dplyr::filter(pred_out, year == yr)
            if (nrow(df_year) == 0) next

            map_data <- sa_gadm2 %>%
              dplyr::left_join(df_year, by = c("CC_2" = "district_code"))

            p <- ggplot(map_data) +
              geom_sf(aes(fill = stgpr_pred), color = "grey50", size = 0.1) +
              scale_fill_gradientn(
                colors = c("#FFE52A", "#F79A19", "#CF0F0F"),
                limits = c(0, 1),
                na.value = "white",
                name = "Predicted"
              ) +
              labs(
                title = paste(str_to_title(gsub("_", " ", bug_name)), "-", yr),
                subtitle = "ST-GPR predicted resistance (district-level)"
              ) +
              theme_minimal(base_size = 12) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.box = "horizontal",
                panel.border = element_rect(color = "#fff", fill = NA, size = 0)
              )

            ggsave(
              filename = file.path(map_dir, paste0(bug_name, "_prediction_map_", yr, ".png")),
              plot = p,
              width = 10, height = 8, dpi = 300
            )
            plot_list[[length(plot_list) + 1]] <- p
          }

          if (length(plot_list) > 0) {
            combined_title <- paste(
              str_to_title(gsub("_", " ", bug_name)),
              "ST-GPR predicted resistance by year"
            )
            combined <- cowplot::plot_grid(
              cowplot::ggdraw() +
                cowplot::draw_label(combined_title, fontface = "bold", hjust = 0.5),
              cowplot::plot_grid(plotlist = plot_list, ncol = 1),
              ncol = 1,
              rel_heights = c(0.06, 0.94)
            )
            per_plot_height <- 6
            panel_height <- min(per_plot_height * length(plot_list), 45)
            ggsave(
              filename = file.path(map_dir, paste0(bug_name, "_prediction_map_panel.png")),
              plot = combined,
              width = 10, height = panel_height, dpi = 300
            )
          }
        } else {
          warning("Missing CC_2 in shapefile or district_code/year/stgpr_pred in prediction data for ", basename(train_path), "; skipping maps")
        }
      } else {
        warning("No GADM shapefile found (expected data/spatial/gadm41_ZAF_2_sf.gpkg or .rds); skipping maps for ", basename(train_path))
      }
    }
  }
}
