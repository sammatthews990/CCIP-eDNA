#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(loo)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_concentration_distance.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
fit_path <- file.path(output_dir, "brms_concentration_distance_raw_conc.rds")
loo_path <- file.path(output_dir, "brms_concentration_distance_raw_conc_loo.rds")

model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
raw_data <- model_data |>
  mutate(edna_conc_z = as.numeric(scale(edna_conc_mean)))
saveRDS(
  c(center = mean(model_data$edna_conc_mean), scale = sd(model_data$edna_conc_mean)),
  file.path(output_dir, "brms_raw_concentration_scale.rds")
)

if (file.exists(fit_path)) {
  raw_fit <- readRDS(fit_path)
} else {
  template_fit <- readRDS(file.path(output_dir, "brms_concentration_distance.rds"))
  raw_fit <- update(
    template_fit,
    newdata = raw_data,
    recompile = FALSE,
    chains = 4, cores = 4,
    iter = 2000, warmup = 1000,
    seed = 20260920,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 200
  )
  saveRDS(raw_fit, fit_path)
}

if (file.exists(loo_path)) {
  raw_loo <- readRDS(loo_path)
} else {
  raw_loo <- brms::loo(raw_fit, cores = 4, moment_match = TRUE)
  saveRDS(raw_loo, loo_path)
}

log_loo <- readRDS(file.path(output_dir, "brms_loo_final.rds"))$concentration_distance
comparison <- as.data.frame(loo::loo_compare(list(
  log1p_concentration = log_loo,
  raw_concentration = raw_loo
))) |>
  rownames_to_column("model")
write.csv(
  comparison,
  file.path(output_dir, "brms_concentration_transform_comparison.csv"),
  row.names = FALSE
)

distribution_summary <- tibble(
  scale = c("Raw", "log1p"),
  skewness = c(
    mean((model_data$edna_conc_mean - mean(model_data$edna_conc_mean))^3) /
      sd(model_data$edna_conc_mean)^3,
    mean((log1p(model_data$edna_conc_mean) -
      mean(log1p(model_data$edna_conc_mean)))^3) /
      sd(log1p(model_data$edna_conc_mean))^3
  ),
  maximum = c(max(model_data$edna_conc_mean), max(log1p(model_data$edna_conc_mean)))
)
write.csv(
  distribution_summary,
  file.path(output_dir, "brms_concentration_transform_summary.csv"),
  row.names = FALSE
)

print(comparison)
print(distribution_summary)
