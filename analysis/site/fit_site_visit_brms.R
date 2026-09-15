#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(loo)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/site_visit_glmm_screen.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
screen <- readRDS(file.path(output_dir, "site_visit_glmm_screen.rds"))
data <- screen$data |>
  mutate(
    Reef = factor(Reef),
    site_id = factor(site_id),
    edna_campaign_id = factor(edna_campaign_id)
  )

fit_path <- file.path(output_dir, "site_visit_brms_final.rds")
loo_path <- file.path(output_dir, "site_visit_brms_loo.rds")
formula <- bf(
  cots_count ~ edna_prop_z + distance_z + lag_z + offset(log_effort) +
    (1 | Reef) + (1 | site_id) + (1 | edna_campaign_id)
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
    formula = formula,
    data = data,
    family = negbinomial(link = "log"),
    prior = priors,
    backend = "rstan",
    chains = 4, cores = 4,
    iter = 3000, warmup = 1000,
    seed = 20261011,
    control = list(adapt_delta = 0.97, max_treedepth = 13),
    save_pars = save_pars(all = TRUE),
    refresh = 200
  )
  saveRDS(fit, fit_path)
}

if (file.exists(loo_path)) {
  fit_loo <- readRDS(loo_path)
} else {
  fit_loo <- brms::loo(fit, cores = 4)
  saveRDS(fit_loo, loo_path)
}

saveRDS(data, file.path(output_dir, "site_visit_model_data.rds"))
saveRDS(screen$scale_spec, file.path(output_dir, "site_visit_scale_spec.rds"))
write.csv(
  tibble::rownames_to_column(as.data.frame(fixef(fit)), "term"),
  file.path(output_dir, "site_visit_brms_fixed_effects.csv"), row.names = FALSE
)
draw_summary <- as.data.frame(posterior_summary(fit)) |>
  tibble::rownames_to_column("parameter")
write.csv(
  draw_summary,
  file.path(output_dir, "site_visit_brms_posterior_summary.csv"), row.names = FALSE
)
print(summary(fit))
print(fit_loo)
