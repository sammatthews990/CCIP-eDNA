#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "brms_scale_spec.rds"))
models <- list(
  no_concentration = readRDS(file.path(output_dir, "brms_best_additive.rds")),
  concentration_additive = readRDS(file.path(output_dir, "brms_full_additive.rds")),
  concentration_distance = readRDS(file.path(output_dir, "brms_selected_final.rds"))
)

effort <- median(model_data$bottom_time)
prediction_grid <- tibble(
  Reef = factor(
    rep(levels(model_data$Reef)[1], 101), levels = levels(model_data$Reef)
  ),
  bottom_time = effort,
  log_effort = log(effort),
  edna_pct = 0:100,
  edna_prop_z = ((edna_pct / 100) - scale_spec$edna_prop_model["center"]) /
    scale_spec$edna_prop_model["scale"],
  edna_conc_z = 0,
  distance_z = 0,
  lag_z = 0
)

find_crossing <- function(y, threshold) {
  if (!any(y >= threshold)) return(NA_real_)
  stats::approx(y, prediction_grid$edna_pct, xout = threshold)$y
}

curve_results <- bind_rows(lapply(names(models), function(model_name) {
  fit <- models[[model_name]]
  expected_counts <- posterior_epred(fit, newdata = prediction_grid, re_formula = NA)
  draws <- posterior::as_draws_df(fit)
  reef_sd <- draws[["sd_Reef__Intercept"]]
  stopifnot(nrow(expected_counts) == length(reef_sd))

  conditional_draws <- expected_counts / effort
  marginal_draws <- sweep(
    expected_counts, 1, exp(0.5 * reef_sd^2), "*"
  ) / effort

  bind_rows(lapply(c("conditional_typical_reef", "marginal_new_reef"), function(scale_name) {
    matrix_used <- if (scale_name == "conditional_typical_reef") {
      conditional_draws
    } else {
      marginal_draws
    }
    median_curve <- apply(matrix_used, 2, median)
    tibble(
      model = model_name,
      prediction_scale = scale_name,
      edna_pct = prediction_grid$edna_pct,
      estimate = median_curve,
      lower = apply(matrix_used, 2, quantile, probs = 0.025),
      upper = apply(matrix_used, 2, quantile, probs = 0.975),
      crossing_002 = find_crossing(median_curve, 0.02),
      crossing_004 = find_crossing(median_curve, 0.04),
      crossing_008 = find_crossing(median_curve, 0.08)
    )
  }))
}))

write.csv(
  curve_results,
  file.path(output_dir, "brms_curve_scale_diagnosis.csv"),
  row.names = FALSE
)
print(curve_results |>
  filter(edna_pct %in% c(0, 60, 80, 100)) |>
  select(model, prediction_scale, edna_pct, estimate, crossing_002, crossing_004))
