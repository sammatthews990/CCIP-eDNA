#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(posterior)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_operational_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
brms_dir <- file.path(project_root, "analysis", "brms", "output")
output_dir <- file.path(project_root, "analysis", "ensemble", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

arguments <- commandArgs(trailingOnly = TRUE)
folds_requested <- if (length(arguments)) as.integer(arguments) else 1:5
model_data <- readRDS(file.path(brms_dir, "brms_model_data.rds"))
fold_map <- read.csv(file.path(brms_dir, "brms_reef_fold_map.csv"))
model_data <- left_join(
  mutate(model_data, Reef = as.character(Reef)), fold_map, by = "Reef"
) |>
  mutate(Reef = factor(Reef, levels = levels(readRDS(
    file.path(brms_dir, "brms_model_data.rds")
  )$Reef)))
template <- readRDS(file.path(brms_dir, "brms_best_additive.rds"))
thresholds <- c(0.02, 0.04, 0.08)

for (fold_id in folds_requested) {
  message("BRMS out-of-fold fit ", fold_id, "/5")
  fit_path <- file.path(output_dir, paste0("brms_operational_fold_v2_", fold_id, ".rds"))
  prediction_path <- file.path(
    output_dir, paste0("brms_operational_oof_fold_", fold_id, ".csv")
  )
  train <- filter(model_data, fold != fold_id)
  test <- filter(model_data, fold == fold_id)

  if (file.exists(fit_path)) {
    fit <- readRDS(fit_path)
  } else {
    fit <- update(
      template,
      newdata = train,
      recompile = FALSE,
      chains = 4, cores = 4,
      iter = 1000, warmup = 500,
      seed = 20261000 + fold_id,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      save_pars = save_pars(all = TRUE),
      refresh = 200
    )
    saveRDS(fit, fit_path)
  }

  expected_conditional <- posterior_epred(
    fit, newdata = test, re_formula = NA, allow_new_levels = TRUE
  )
  draws <- posterior::as_draws_df(fit)
  reef_sd <- draws[["sd_Reef__Intercept"]]
  shape <- draws[["shape"]]
  stopifnot(nrow(expected_conditional) == length(reef_sd))

  expected_marginal <- sweep(
    expected_conditional, 1, exp(0.5 * reef_sd^2), "*"
  ) / rep(test$bottom_time, each = nrow(expected_conditional))
  expected_marginal <- matrix(
    expected_marginal, nrow = nrow(expected_conditional)
  )

  set.seed(20262000 + fold_id)
  reef_effect_draws <- matrix(
    rnorm(length(expected_conditional), 0, rep(reef_sd, ncol(expected_conditional))),
    nrow = nrow(expected_conditional)
  )
  future_mu <- expected_conditional * exp(reef_effect_draws)
  prediction_counts <- matrix(
    rnbinom(length(future_mu), mu = as.vector(future_mu), size = rep(shape, ncol(future_mu))),
    nrow = nrow(future_mu)
  )
  prediction_cpue <- sweep(prediction_counts, 2, test$bottom_time, "/")

  output <- test |>
    transmute(
      cull_id, Reef = as.character(Reef), fold,
      observed_cpue = cots_count / bottom_time,
      bottom_time,
      predicted_cpue = apply(expected_marginal, 2, median),
      expected_lower = apply(expected_marginal, 2, quantile, probs = 0.025),
      expected_upper = apply(expected_marginal, 2, quantile, probs = 0.975),
      prediction_lower = apply(prediction_cpue, 2, quantile, probs = 0.025),
      prediction_upper = apply(prediction_cpue, 2, quantile, probs = 0.975)
    )
  for (threshold in thresholds) {
    minimum_count <- ceiling(threshold * test$bottom_time)
    tail_probability <- matrix(1 - pnbinom(
      rep(minimum_count - 1, each = nrow(future_mu)),
      mu = as.vector(future_mu),
      size = rep(shape, ncol(future_mu))
    ), nrow = nrow(future_mu))
    probability <- colMeans(tail_probability)
    output[[paste0("prob_", gsub("\\.", "", sprintf("%.2f", threshold)))]] <- probability
  }
  write.csv(output, prediction_path, row.names = FALSE)
}

prediction_files <- file.path(
  output_dir, paste0("brms_operational_oof_fold_", 1:5, ".csv")
)
if (all(file.exists(prediction_files))) {
  all_predictions <- bind_rows(lapply(prediction_files, read.csv)) |>
    arrange(cull_id)
  write.csv(
    all_predictions,
    file.path(output_dir, "brms_operational_oof_predictions.csv"),
    row.names = FALSE
  )
}
