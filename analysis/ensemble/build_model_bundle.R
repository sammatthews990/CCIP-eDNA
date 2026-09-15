#!/usr/bin/env Rscript

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_operational_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
brms_dir <- file.path(project_root, "analysis", "brms", "output")
brt_dir <- file.path(project_root, "analysis", "brt", "output")
output_dir <- file.path(project_root, "analysis", "ensemble", "output")

bundle <- list(
  brms_fit = readRDS(file.path(brms_dir, "brms_operational_final.rds")),
  brt_model = readRDS(file.path(brt_dir, "reefDNA_brt_model.rds")),
  ensemble_spec = readRDS(file.path(output_dir, "reefDNA_ensemble_spec.rds")),
  scale_spec = readRDS(file.path(brms_dir, "brms_scale_spec.rds")),
  metadata = list(
    version = "1.0.0",
    response = "COTS CPUE per minute",
    predictors = c("perc_pos", "distance_m", "lag_days"),
    thresholds = c(0.02, 0.04, 0.08),
    brms_formula = paste(
      "cots_count ~ edna_prop_z + distance_z + lag_z +",
      "offset(log_effort) + (1 | Reef)"
    ),
    brt_objective = "Poisson boosted trees with log-effort base margin",
    validation = "five-fold whole-reef cross-validation",
    note = paste(
      "Concentration excluded from operational models. Predictions for unseen",
      "reefs integrate reef-level heterogeneity in BRMS."
    ),
    created = Sys.time()
  )
)
saveRDS(
  bundle,
  file.path(output_dir, "reefDNA_operational_model_bundle.rds"),
  compress = "xz"
)
message("Operational model bundle saved.")
