#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/cull_first/output/cull_first_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "cull_first", "output")
data <- readRDS(file.path(output_dir, "cull_first_model_data.rds")) |>
  mutate(
    Reef = factor(Reef), site_id = factor(site_id), event_id = factor(event_id)
  )
fit_path <- file.path(output_dir, "cull_first_brms_final.rds")

formula <- bf(
  cots_count ~ edna_prop_z + distance_z + lag_z + offset(log_effort) +
    (1 | Reef) + (1 | site_id) + (1 | event_id)
)
priors <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(-3.5, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))
rstan::rstan_options(auto_write = TRUE)
if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
} else {
  fit <- brm(
    formula = formula, data = data,
    family = negbinomial(link = "log"), prior = priors,
    backend = "rstan", chains = 4, cores = 4,
    iter = 3000, warmup = 1000, seed = 20261021,
    control = list(adapt_delta = 0.97, max_treedepth = 13),
    save_pars = save_pars(all = TRUE), refresh = 200
  )
  saveRDS(fit, fit_path)
}
write.csv(
  tibble::rownames_to_column(as.data.frame(fixef(fit)), "term"),
  file.path(output_dir, "cull_first_brms_fixed_effects.csv"), row.names = FALSE
)
write.csv(
  tibble::rownames_to_column(as.data.frame(posterior_summary(fit)), "parameter"),
  file.path(output_dir, "cull_first_brms_posterior_summary.csv"), row.names = FALSE
)
print(summary(fit))
