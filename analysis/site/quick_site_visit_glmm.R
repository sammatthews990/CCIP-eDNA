#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/edna_first_site_visit_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
data <- readRDS(file.path(output_dir, "edna_first_site_visit_data.rds")) |>
  mutate(
    Reef = factor(Reef),
    site_id = factor(paste(Reef, cull_site_name, sep = "__")),
    edna_campaign_id = factor(edna_campaign_id),
    sampling_design = factor(
      sampling_design,
      levels = c("3x12", "4x6", "other_or_mixed")
    ),
    edna_prop_model = edna_prop_positive,
    distance_model = log1p(distance_m / 200),
    lag_model = log1p(lag_days),
    log_effort = log(bottom_time)
  )

scale_spec <- lapply(
  data[c("edna_prop_model", "distance_model", "lag_model")],
  function(x) c(center = mean(x), scale = sd(x))
)
for (variable in names(scale_spec)) {
  short_name <- sub("_model$", "", variable)
  data[[paste0(short_name, "_z")]] <-
    (data[[variable]] - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}

# This is the shared, keyed input for the INLA, BRT and BRMS site models.
# Write it here so downstream models do not depend on a BRMS fit having run.
saveRDS(data, file.path(output_dir, "site_visit_model_data.rds"))
saveRDS(scale_spec, file.path(output_dir, "site_visit_scale_spec.rds"))

fixed <- cots_count ~ edna_prop_z + distance_z + lag_z + offset(log_effort)
formulas <- list(
  reef = update(fixed, . ~ . + (1 | Reef)),
  reef_site = update(fixed, . ~ . + (1 | Reef) + (1 | site_id)),
  reef_campaign = update(fixed, . ~ . + (1 | Reef) + (1 | edna_campaign_id)),
  reef_site_campaign = update(
    fixed, . ~ . + (1 | Reef) + (1 | site_id) + (1 | edna_campaign_id)
  )
)

fits <- lapply(formulas, function(formula) {
  glmmTMB(
    formula, data = data, family = nbinom2,
    control = glmmTMBControl(
      optimizer = optim,
      optArgs = list(method = "BFGS")
    )
  )
})

random_sd <- function(fit, group) {
  variance <- VarCorr(fit)$cond
  if (!group %in% names(variance)) return(NA_real_)
  unname(attr(variance[[group]], "stddev")[[1]])
}

performance <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  tibble(
    model = model_name,
    n = nobs(fit),
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    delta_AIC = NA_real_,
    reef_sd = random_sd(fit, "Reef"),
    site_sd = random_sd(fit, "site_id"),
    campaign_sd = random_sd(fit, "edna_campaign_id"),
    dispersion = sigma(fit),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = fit$sdr$pdHess
  )
})) |>
  mutate(delta_AIC = AIC - min(AIC)) |>
  arrange(AIC)

effects <- bind_rows(lapply(names(fits), function(model_name) {
  coefficients <- summary(fits[[model_name]])$coefficients$cond
  tibble(
    model = model_name,
    term = rownames(coefficients),
    estimate = coefficients[, "Estimate"],
    std_error = coefficients[, "Std. Error"],
    z_value = coefficients[, "z value"],
    p_value = coefficients[, "Pr(>|z|)"]
  )
}))

crossing <- function(x, y, threshold) {
  if (all(y < threshold) || all(y > threshold)) return(NA_real_)
  approx(y, x, xout = threshold)$y
}

edna_pct <- seq(0, 100, by = 0.25)
prediction_data <- data.frame(
  edna_prop_z = ((edna_pct / 100) - scale_spec$edna_prop_model["center"]) /
    scale_spec$edna_prop_model["scale"],
  distance_z = 0,
  lag_z = 0,
  log_effort = 0,
  Reef = data$Reef[[1]],
  site_id = data$site_id[[1]],
  edna_campaign_id = data$edna_campaign_id[[1]]
)

curves <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  conditional <- exp(predict(
    fit, newdata = prediction_data, type = "link", re.form = NA
  ))
  total_random_variance <- sum(
    c(
      random_sd(fit, "Reef"), random_sd(fit, "site_id"),
      random_sd(fit, "edna_campaign_id")
    )^2,
    na.rm = TRUE
  )
  tibble(
    model = model_name,
    edna_pct = edna_pct,
    conditional_cpue = conditional,
    marginal_new_group_cpue = conditional * exp(0.5 * total_random_variance)
  )
}))

crossings <- curves |>
  group_by(model) |>
  summarise(
    crossing_002 = crossing(edna_pct, marginal_new_group_cpue, 0.02),
    crossing_004 = crossing(edna_pct, marginal_new_group_cpue, 0.04),
    crossing_008 = crossing(edna_pct, marginal_new_group_cpue, 0.08),
    .groups = "drop"
  )

performance <- left_join(performance, crossings, by = "model")
saveRDS(
  list(fits = fits, data = data, scale_spec = scale_spec, performance = performance),
  file.path(output_dir, "site_visit_glmm_screen.rds")
)
write.csv(performance, file.path(output_dir, "site_visit_glmm_performance.csv"), row.names = FALSE)
write.csv(effects, file.path(output_dir, "site_visit_glmm_effects.csv"), row.names = FALSE)
write.csv(curves, file.path(output_dir, "site_visit_glmm_curves.csv"), row.names = FALSE)
print(performance, width = Inf)
print(effects |> filter(model == performance$model[[1]]), width = Inf)
