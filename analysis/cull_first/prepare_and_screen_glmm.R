#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
})

project_root <- normalizePath(
  if (file.exists("analysis/inla/output/distance_time_additive_screen.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "cull_first", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
screen <- readRDS(file.path(
  project_root, "analysis", "inla", "output", "distance_time_additive_screen.rds"
))
data <- screen$data |>
  transmute(
    cull_id, Reef = factor(Reef), cull_site_name = site_name,
    site_id = factor(paste(Reef, site_name, sep = "__")),
    event_id = factor(event_id), date_cull, edna_date, edna_site_name,
    cots_count = as.integer(cots_count), bottom_time, cpue,
    log_effort, edna_prop_z, distance_z, lag_z,
    edna_prop_positive, edna_distance_m, edna_lag_days
  )
stopifnot(nrow(data) == 2117L, !anyDuplicated(data$cull_id))

fixed <- cots_count ~ edna_prop_z + distance_z + lag_z + offset(log_effort)
formulas <- list(
  reef = update(fixed, . ~ . + (1 | Reef)),
  reef_site = update(fixed, . ~ . + (1 | Reef) + (1 | site_id)),
  reef_event = update(fixed, . ~ . + (1 | Reef) + (1 | event_id)),
  reef_site_event = update(
    fixed, . ~ . + (1 | Reef) + (1 | site_id) + (1 | event_id)
  )
)
fits <- lapply(formulas, function(formula) {
  glmmTMB(
    formula, data = data, family = nbinom2,
    control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS")
    )
  )
})

random_sd <- function(fit, group) {
  component <- VarCorr(fit)$cond
  if (!group %in% names(component)) return(NA_real_)
  unname(attr(component[[group]], "stddev")[[1]])
}
performance <- bind_rows(lapply(names(fits), function(model) {
  fit <- fits[[model]]
  tibble(
    model, n = nobs(fit), logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
    reef_sd = random_sd(fit, "Reef"), site_sd = random_sd(fit, "site_id"),
    event_sd = random_sd(fit, "event_id"), dispersion = sigma(fit),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = fit$sdr$pdHess
  )
})) |>
  mutate(delta_AIC = AIC - min(AIC)) |>
  arrange(AIC)
effects <- bind_rows(lapply(names(fits), function(model) {
  coefficients <- summary(fits[[model]])$coefficients$cond
  tibble(
    model, term = rownames(coefficients), estimate = coefficients[, "Estimate"],
    std_error = coefficients[, "Std. Error"], p_value = coefficients[, "Pr(>|z|)"]
  )
}))

saveRDS(
  list(data = data, scale_spec = screen$scale_spec, fits = fits, performance = performance),
  file.path(output_dir, "cull_first_glmm_screen.rds")
)
saveRDS(data, file.path(output_dir, "cull_first_model_data.rds"))
saveRDS(screen$scale_spec, file.path(output_dir, "cull_first_scale_spec.rds"))
write.csv(performance, file.path(output_dir, "cull_first_glmm_performance.csv"), row.names = FALSE)
write.csv(effects, file.path(output_dir, "cull_first_glmm_effects.csv"), row.names = FALSE)
print(performance, width = Inf)
print(effects |> filter(model %in% c("reef_site", "reef_site_event")), width = Inf)
