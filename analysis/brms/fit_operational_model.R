#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_best_additive.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
fit_path <- file.path(output_dir, "brms_operational_final.rds")
loo_path <- file.path(output_dir, "brms_operational_final_loo.rds")

if (file.exists(fit_path)) {
  fit <- readRDS(fit_path)
} else {
  template <- readRDS(file.path(output_dir, "brms_best_additive.rds"))
  model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
  fit <- update(
    template,
    newdata = model_data,
    recompile = FALSE,
    chains = 4, cores = 4,
    iter = 4000, warmup = 1000,
    seed = 20260921,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE),
    refresh = 200
  )
  saveRDS(fit, fit_path)
}

if (!file.exists(loo_path)) {
  fit_loo <- brms::loo(fit, cores = 4, moment_match = TRUE)
  saveRDS(fit_loo, loo_path)
} else {
  fit_loo <- readRDS(loo_path)
}

print(summary(fit))
print(fit_loo)
