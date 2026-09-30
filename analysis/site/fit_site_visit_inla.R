#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(INLA)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/site_visit_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds")) |>
  mutate(
    sampling_design = factor(
      as.character(sampling_design),
      levels = c("3x12", "4x6", "other_or_mixed")
    ),
    reef_index = as.integer(factor(Reef)),
    site_index = as.integer(factor(site_id)),
    campaign_index = as.integer(factor(edna_campaign_id))
  )

required <- c(
  "cots_count", "log_effort", "edna_prop_z", "distance_z", "lag_z",
  "sampling_design", "reef_index", "site_index", "campaign_index"
)
if (!all(required %in% names(data))) {
  stop("Site-model data are missing: ", paste(setdiff(required, names(data)), collapse = ", "))
}
if (anyNA(data[, required])) stop("INLA site-model inputs contain missing values.")
if (anyDuplicated(data$site_visit_id)) stop("Site visits must be unique model responses.")

fixed_terms <- list(
  base = "edna_prop_z + distance_z + lag_z",
  sampling_design = "edna_prop_z + distance_z + lag_z + sampling_design",
  sampling_design_by_signal = paste(
    "edna_prop_z * sampling_design + distance_z + lag_z"
  )
)
random_terms <- paste(
  "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))",
  "+ f(site_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))",
  "+ f(campaign_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))"
)

fits <- lapply(fixed_terms, function(terms) {
  formula <- as.formula(paste(
    "cots_count ~", terms, "+ offset(log_effort) +", random_terms
  ))
  INLA::inla(
    formula,
    family = "nbinomial",
    data = data,
    control.predictor = list(compute = TRUE),
    control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE),
    control.family = list(
      hyper = list(theta = list(prior = "pc.mgamma", param = 7))
    ),
    num.threads = "2:1",
    verbose = FALSE
  )
})

ranking <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  valid_cpo <- is.finite(fit$cpo$cpo) & fit$cpo$cpo > 0
  tibble(
    model = model_name,
    n_site_visits = nrow(data),
    n_fixed = nrow(fit$summary.fixed),
    waic = fit$waic$waic,
    dic = fit$dic$dic,
    mean_log_cpo = mean(log(fit$cpo$cpo[valid_cpo])),
    cpo_failures = sum(!valid_cpo)
  )
})) |>
  arrange(waic) |>
  mutate(delta_waic = waic - min(waic), rank_waic = row_number())

fixed_effects <- bind_rows(lapply(names(fits), function(model_name) {
  rownames_to_column(as.data.frame(fits[[model_name]]$summary.fixed), "term") |>
    mutate(model = model_name, .before = 1)
}))

design_counts <- data |>
  count(sampling_design, name = "n_site_visits") |>
  mutate(proportion = n_site_visits / sum(n_site_visits))

write.csv(ranking, file.path(output_dir, "site_visit_inla_ranking.csv"), row.names = FALSE)
write.csv(fixed_effects, file.path(output_dir, "site_visit_inla_fixed_effects.csv"), row.names = FALSE)
write.csv(design_counts, file.path(output_dir, "site_visit_sampling_design_counts.csv"), row.names = FALSE)
saveRDS(
  list(
    fits = fits, ranking = ranking, fixed_effects = fixed_effects,
    design_counts = design_counts, n_site_visits = nrow(data),
    response_key = data$site_visit_id
  ),
  file.path(output_dir, "site_visit_inla_models.rds")
)

print(ranking, width = Inf)
print(design_counts, width = Inf)
