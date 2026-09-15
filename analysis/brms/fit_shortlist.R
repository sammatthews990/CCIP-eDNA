#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/inla/output/distance_time_additive_screen.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

screen <- readRDS(file.path(
  project_root, "analysis", "inla", "output", "distance_time_additive_screen.rds"
))
model_data <- screen$data |>
  transmute(
    cull_id, Reef = factor(Reef), cots_count = as.integer(cots_count),
    bottom_time, log_effort,
    edna_prop_z, edna_conc_z, distance_z, lag_z,
    edna_prop_positive, edna_conc_mean, edna_distance_m, edna_lag_days
  )
stopifnot(
  nrow(model_data) == dplyr::n_distinct(model_data$cull_id),
  all(model_data$bottom_time > 0),
  all(model_data$cots_count >= 0)
)
saveRDS(model_data, file.path(output_dir, "brms_model_data.rds"))
saveRDS(screen$scale_spec, file.path(output_dir, "brms_scale_spec.rds"))

model_formulas <- list(
  null = bf(cots_count ~ 1 + offset(log_effort) + (1 | Reef)),
  best_additive = bf(
    cots_count ~ edna_prop_z + distance_z + lag_z +
      offset(log_effort) + (1 | Reef)
  ),
  full_additive = bf(
    cots_count ~ edna_prop_z + edna_conc_z + distance_z + lag_z +
      offset(log_effort) + (1 | Reef)
  ),
  concentration_distance = bf(
    cots_count ~ edna_prop_z + edna_conc_z + distance_z + lag_z +
      edna_conc_z:distance_z + offset(log_effort) + (1 | Reef)
  ),
  dual_interaction = bf(
    cots_count ~ edna_prop_z + edna_conc_z + distance_z + lag_z +
      edna_conc_z:distance_z + edna_prop_z:edna_conc_z +
      offset(log_effort) + (1 | Reef)
  ),
  all_pairwise = bf(
    cots_count ~ (edna_prop_z + edna_conc_z + distance_z + lag_z)^2 +
      offset(log_effort) + (1 | Reef)
  )
)

common_priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(-3.5, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)

options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))
rstan::rstan_options(auto_write = TRUE)

fits <- vector("list", length(model_formulas))
names(fits) <- names(model_formulas)
for (i in seq_along(model_formulas)) {
  model_name <- names(model_formulas)[[i]]
  fit_path <- file.path(output_dir, paste0("brms_", model_name, ".rds"))
  if (file.exists(fit_path)) {
    message("Loading cached model: ", model_name)
    fits[[i]] <- readRDS(fit_path)
    next
  }
  message("Fitting model ", i, "/", length(model_formulas), ": ", model_name)
  fit <- brm(
    formula = model_formulas[[i]],
    data = model_data,
    family = negbinomial(link = "log"),
    prior = if (model_name == "null") common_priors[-1, ] else common_priors,
    backend = "rstan",
    chains = 4, cores = 4,
    iter = 2000, warmup = 1000,
    seed = 20260909 + i,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE),
    refresh = 200,
    silent = 0
  )
  saveRDS(fit, fit_path)
  fits[[i]] <- fit
}

saveRDS(fits, file.path(output_dir, "brms_shortlist_fits.rds"))
message("All shortlisted brms models completed.")
