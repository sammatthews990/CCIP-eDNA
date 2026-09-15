#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(loo)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_dual_interaction.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
fit <- readRDS(file.path(output_dir, "brms_dual_interaction.rds"))
base_loo <- readRDS(file.path(output_dir, "brms_loo_final.rds"))
target_loo_path <- file.path(output_dir, "brms_dual_interaction_loo_moment_match.rds")

if (file.exists(target_loo_path)) {
  target_loo <- readRDS(target_loo_path)
} else {
  target_loo <- brms::loo(fit, cores = 4, moment_match = TRUE)
  saveRDS(target_loo, target_loo_path)
}

combined_loo <- c(base_loo, list(dual_interaction = target_loo))
comparison <- as.data.frame(loo::loo_compare(combined_loo)) |>
  rownames_to_column("model")
write.csv(
  comparison,
  file.path(output_dir, "brms_loo_targeted_comparison.csv"),
  row.names = FALSE
)

target_vs_selected <- as.data.frame(loo::loo_compare(list(
  dual_interaction = target_loo,
  concentration_distance = base_loo$concentration_distance
))) |>
  rownames_to_column("model")
write.csv(
  target_vs_selected,
  file.path(output_dir, "brms_targeted_vs_selected.csv"),
  row.names = FALSE
)

target_effects <- as.data.frame(
  fixef(fit, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
) |>
  rownames_to_column("term")
write.csv(
  target_effects,
  file.path(output_dir, "brms_dual_interaction_effects.csv"),
  row.names = FALSE
)

fit_summary <- summary(fit)
diagnostic_rows <- bind_rows(
  as.data.frame(fit_summary$fixed),
  as.data.frame(fit_summary$spec_pars),
  as.data.frame(fit_summary$random$Reef)
)
nuts <- nuts_params(fit)
target_diagnostics <- tibble(
  max_rhat = max(diagnostic_rows$Rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostic_rows$Bulk_ESS, na.rm = TRUE),
  min_tail_ess = min(diagnostic_rows$Tail_ESS, na.rm = TRUE),
  divergences = sum(nuts$Value[nuts$Parameter == "divergent__"]),
  max_treedepth_hits = sum(nuts$Value[nuts$Parameter == "treedepth__"] >= 12),
  max_pareto_k = max(loo::pareto_k_values(target_loo))
)
write.csv(
  target_diagnostics,
  file.path(output_dir, "brms_dual_interaction_diagnostics.csv"),
  row.names = FALSE
)

print(target_vs_selected)
print(target_diagnostics)
