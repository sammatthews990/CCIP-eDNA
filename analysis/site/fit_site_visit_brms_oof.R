#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/site_visit_brms_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
template <- readRDS(file.path(output_dir, "site_visit_brms_final.rds"))
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds")) |>
  mutate(
    Reef = as.character(Reef), site_id = as.character(site_id),
    edna_campaign_id = as.character(edna_campaign_id),
    cpue = cots_count / bottom_time
  )
fold_map <- read.csv(file.path(output_dir, "site_visit_campaign_fold_map.csv"))
data <- left_join(data, select(fold_map, edna_campaign_id, fold), by = "edna_campaign_id")
stopifnot(!anyNA(data$fold), !anyDuplicated(data$site_visit_id))

arguments <- commandArgs(trailingOnly = TRUE)
folds_requested <- if (length(arguments)) as.integer(arguments) else 1:5
for (fold_id in folds_requested) {
  message("Site-visit BRMS validation fold ", fold_id, "/5")
  fit_path <- file.path(
    output_dir,
    if (fold_id == 5L) "site_visit_brms_fold_5_v2.rds" else {
      paste0("site_visit_brms_fold_", fold_id, ".rds")
    }
  )
  prediction_path <- file.path(output_dir, paste0("site_visit_brms_oof_fold_", fold_id, ".csv"))
  train <- filter(data, fold != fold_id) |>
    mutate(
      Reef = factor(Reef), site_id = factor(site_id),
      edna_campaign_id = factor(edna_campaign_id)
    )
  test <- filter(data, fold == fold_id)
  if (file.exists(fit_path)) {
    fit <- readRDS(fit_path)
  } else {
    fit <- update(
      template, newdata = train, recompile = FALSE,
      chains = 4, cores = 4,
      iter = if (fold_id == 5L) 1800 else 1200,
      warmup = if (fold_id == 5L) 900 else 600,
      seed = 20261100 + fold_id,
      control = list(
        adapt_delta = if (fold_id == 5L) 0.995 else 0.97,
        max_treedepth = 14
      ),
      save_pars = save_pars(all = TRUE), refresh = 200
    )
    saveRDS(fit, fit_path)
  }
  set.seed(20261200 + fold_id)
  expected_counts <- posterior_epred(
    fit, newdata = test, re_formula = NULL,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  )
  predicted_counts <- posterior_predict(
    fit, newdata = test, re_formula = NULL,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  )
  expected_draws_cpue <- sweep(expected_counts, 2, test$bottom_time, "/")
  predicted_draws_cpue <- sweep(predicted_counts, 2, test$bottom_time, "/")
  output <- test |>
    transmute(
      site_visit_id, Reef, site_id, edna_campaign_id, fold,
      observed_cpue = cpue, bottom_time,
      predicted_cpue = apply(expected_draws_cpue, 2, mean),
      expected_lower = apply(expected_draws_cpue, 2, quantile, probs = 0.025),
      expected_upper = apply(expected_draws_cpue, 2, quantile, probs = 0.975),
      prediction_lower = apply(predicted_draws_cpue, 2, quantile, probs = 0.025),
      prediction_upper = apply(predicted_draws_cpue, 2, quantile, probs = 0.975)
    )
  for (threshold in c(0.02, 0.04, 0.08)) {
    suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
    output[[paste0("prob_", suffix)]] <- colMeans(predicted_draws_cpue >= threshold)
  }
  write.csv(output, prediction_path, row.names = FALSE)
}

paths <- file.path(output_dir, paste0("site_visit_brms_oof_fold_", 1:5, ".csv"))
if (all(file.exists(paths))) {
  bind_rows(lapply(paths, read.csv)) |>
    arrange(site_visit_id) |>
    write.csv(
      file.path(output_dir, "site_visit_brms_oof.csv"), row.names = FALSE
    )
}
