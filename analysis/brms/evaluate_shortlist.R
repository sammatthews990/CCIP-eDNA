#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(loo)
  library(posterior)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_null.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
model_names <- c(
  "null", "best_additive", "full_additive",
  "concentration_distance", "all_pairwise"
)
fits <- setNames(lapply(model_names, function(model_name) {
  readRDS(file.path(output_dir, paste0("brms_", model_name, ".rds")))
}), model_names)

diagnostics <- bind_rows(lapply(model_names, function(model_name) {
  fit <- fits[[model_name]]
  summary_all <- posterior::summarise_draws(
    posterior::as_draws_array(fit),
    posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  )
  fixed_rows <- grepl("^b_", summary_all$variable) |
    summary_all$variable %in% c("shape", "sd_Reef__Intercept")
  nuts <- brms::nuts_params(fit)
  tibble(
    model = model_name,
    max_rhat_all = max(summary_all$rhat, na.rm = TRUE),
    min_bulk_ess_all = min(summary_all$ess_bulk, na.rm = TRUE),
    min_tail_ess_all = min(summary_all$ess_tail, na.rm = TRUE),
    max_rhat_key = max(summary_all$rhat[fixed_rows], na.rm = TRUE),
    min_bulk_ess_key = min(summary_all$ess_bulk[fixed_rows], na.rm = TRUE),
    min_tail_ess_key = min(summary_all$ess_tail[fixed_rows], na.rm = TRUE),
    divergences = sum(nuts$Value[nuts$Parameter == "divergent__"]),
    max_treedepth_hits = sum(nuts$Value[nuts$Parameter == "treedepth__"] >= 12)
  )
}))
write.csv(diagnostics, file.path(output_dir, "brms_sampler_diagnostics.csv"), row.names = FALSE)

loo_objects <- vector("list", length(fits))
names(loo_objects) <- names(fits)
for (model_name in names(fits)) {
  message("Computing PSIS-LOO: ", model_name)
  loo_objects[[model_name]] <- brms::loo(fits[[model_name]], cores = 4)
}
saveRDS(loo_objects, file.path(output_dir, "brms_loo_objects.rds"))

comparison <- as.data.frame(loo::loo_compare(loo_objects)) |>
  rownames_to_column("model")
pareto <- bind_rows(lapply(names(loo_objects), function(model_name) {
  k <- loo::pareto_k_values(loo_objects[[model_name]])
  tibble(
    model = model_name,
    pareto_k_max = max(k),
    pareto_k_gt_0_7 = sum(k > 0.7),
    pareto_k_gt_1 = sum(k > 1)
  )
}))
comparison <- left_join(comparison, pareto, by = "model")
weights <- tibble(
  model = names(fits),
  stacking_weight = as.numeric(loo::loo_model_weights(loo_objects, method = "stacking")),
  pseudobma_weight = as.numeric(loo::loo_model_weights(loo_objects, method = "pseudobma"))
)
comparison <- left_join(comparison, weights, by = "model")
write.csv(comparison, file.path(output_dir, "brms_loo_comparison.csv"), row.names = FALSE)
write.csv(weights, file.path(output_dir, "brms_loo_weights.csv"), row.names = FALSE)

fixed_effects <- bind_rows(lapply(model_names, function(model_name) {
  out <- as.data.frame(brms::fixef(fits[[model_name]], probs = c(0.025, 0.1, 0.5, 0.9, 0.975))) |>
    rownames_to_column("term")
  mutate(out, model = model_name, .before = 1)
}))
write.csv(fixed_effects, file.path(output_dir, "brms_fixed_effects.csv"), row.names = FALSE)

print(diagnostics)
print(comparison)
print(weights)
