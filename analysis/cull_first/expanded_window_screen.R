#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readxl)
  library(sf)
  library(splines)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("data/eDNA data_ALL_20260528.xlsx")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
devtools::load_all(file.path(project_root, "reefDNA"), quiet = TRUE)
output_dir <- file.path(project_root, "analysis", "cull_first", "output")
plot_dir <- file.path(project_root, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

edna_raw <- read_excel(
  file.path(project_root, "data", "eDNA data_ALL_20260528.xlsx"),
  sheet = "eDNA_data_ALL"
)
cull_raw <- read_excel(
  file.path(project_root, "data", "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"),
  sheet = "Cull"
)
edna <- prepare_inla_edna(edna_raw, crs_projected = 3112)
culls <- prepare_inla_culls(
  cull_raw, crs_projected = 3112,
  min_date = min(edna$date_edna), max_date = max(edna$date_edna) + 365
)
links <- build_inla_space_time_links(
  edna, culls, max_distance_m = 5000,
  min_lag_days = 0, max_lag_days = 365, same_reef = TRUE
)
benchmark <- summarise_inla_event_benchmark(links, edna, culls)
baseline_scale <- readRDS(file.path(output_dir, "cull_first_scale_spec.rds"))
prepared <- prepare_inla_screen_data(benchmark, scale_spec = baseline_scale)
data <- prepared$data |>
  mutate(
    Reef = factor(Reef),
    site_id = factor(paste(Reef, site_name, sep = "__")),
    event_id = factor(event_id),
    cpue = cots_count / bottom_time
  )

coverage <- tibble(
  window = c("2,000 m / 183 days", "5,000 m / 365 days"),
  n_cull_dives = c(2117L, nrow(data)),
  n_reefs = c(36L, n_distinct(data$Reef)),
  n_cull_sites = c(332L, n_distinct(data$site_id)),
  n_edna_events = c(165L, n_distinct(data$event_id)),
  median_distance_m = c(
    median(readRDS(file.path(output_dir, "cull_first_model_data.rds"))$edna_distance_m),
    median(data$edna_distance_m)
  ),
  median_lag_days = c(
    median(readRDS(file.path(output_dir, "cull_first_model_data.rds"))$edna_lag_days),
    median(data$edna_lag_days)
  )
)

random_terms <- "(1 | Reef) + (1 | site_id) + (1 | event_id)"
formulas <- list(
  additive_log = as.formula(paste(
    "cots_count ~ edna_prop_z + distance_z + lag_z + offset(log_effort) +",
    random_terms
  )),
  linear_degradation = as.formula(paste(
    "cots_count ~ edna_prop_z * distance_z + edna_prop_z * lag_z +",
    "offset(log_effort) +", random_terms
  )),
  spline_additive = as.formula(paste(
    "cots_count ~ edna_prop_z + ns(distance_model, 3) + ns(lag_model, 3) +",
    "offset(log_effort) +", random_terms
  )),
  spline_degradation = as.formula(paste(
    "cots_count ~ edna_prop_z * ns(distance_model, 3) +",
    "edna_prop_z * ns(lag_model, 3) + offset(log_effort) +", random_terms
  ))
)
fits <- lapply(formulas, function(formula) {
  glmmTMB(
    formula, data = data, family = nbinom2,
    control = glmmTMBControl(
      optimizer = optim, optArgs = list(method = "BFGS")
    )
  )
})
performance <- bind_rows(lapply(names(fits), function(model) {
  fit <- fits[[model]]
  tibble(
    model, n = nobs(fit), logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = fit$sdr$pdHess
  )
})) |>
  mutate(delta_AIC = AIC - min(AIC)) |>
  arrange(AIC)
best_name <- performance$model[[1]]
best_fit <- fits[[best_name]]

random_variance <- sum(vapply(
  VarCorr(best_fit)$cond,
  function(component) unname(attr(component, "stddev")[[1]])^2,
  numeric(1)
))
grid <- crossing(
  edna_pct = seq(0, 100, by = 2),
  distance_m = c(200, 500, 1000, 2000, 3000, 5000),
  lag_days = c(30, 90, 183, 270, 365)
) |>
  mutate(
    edna_prop_model = edna_pct / 100,
    distance_model = log1p(distance_m / 200),
    lag_model = log1p(lag_days),
    edna_prop_z = (edna_prop_model - baseline_scale$edna_prop_model["center"]) /
      baseline_scale$edna_prop_model["scale"],
    distance_z = (distance_model - baseline_scale$distance_model["center"]) /
      baseline_scale$distance_model["scale"],
    lag_z = (lag_model - baseline_scale$lag_model["center"]) /
      baseline_scale$lag_model["scale"],
    log_effort = 0,
    Reef = data$Reef[[1]], site_id = data$site_id[[1]], event_id = data$event_id[[1]]
  )
grid_prediction <- predict(
  best_fit, newdata = grid, type = "link", re.form = NA, se.fit = TRUE
)
grid <- grid |>
  mutate(
    expected_cpue = exp(grid_prediction$fit + 0.5 * random_variance),
    lower = exp(grid_prediction$fit - 1.96 * grid_prediction$se.fit + 0.5 * random_variance),
    upper = exp(grid_prediction$fit + 1.96 * grid_prediction$se.fit + 0.5 * random_variance)
  )

degradation <- grid |>
  filter(edna_pct %in% c(0, 100)) |>
  select(edna_pct, distance_m, lag_days, expected_cpue) |>
  pivot_wider(names_from = edna_pct, values_from = expected_cpue, names_prefix = "pct_") |>
  mutate(
    absolute_uplift = pct_100 - pct_0,
    rate_ratio = pct_100 / pct_0
  )

write.csv(coverage, file.path(output_dir, "expanded_window_coverage.csv"), row.names = FALSE)
write.csv(performance, file.path(output_dir, "expanded_window_glmm_performance.csv"), row.names = FALSE)
write.csv(grid, file.path(output_dir, "expanded_window_prediction_grid.csv"), row.names = FALSE)
write.csv(degradation, file.path(output_dir, "expanded_window_degradation.csv"), row.names = FALSE)
saveRDS(
  list(
    fits = fits, best_model = best_name, data = data,
    scale_spec = baseline_scale, coverage = coverage
  ),
  file.path(output_dir, "expanded_window_glmm_screen.rds")
)

surface_plot <- ggplot(
  filter(grid, edna_pct %in% c(25, 50, 75, 100)),
  aes(distance_m, expected_cpue, colour = factor(edna_pct), group = edna_pct)
) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey65") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(edna_pct)),
    colour = NA, alpha = 0.08, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~lag_days, labeller = label_both) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  labs(
    title = "Expanded-window signal",
    subtitle = paste("Best GLMM:", best_name),
    x = "Distance from eDNA site (m)", y = "Expected CPUE (min^-1)",
    colour = "% positive"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(
  file.path(plot_dir, "cull_first_expanded_window.png"), surface_plot,
  width = 12, height = 7, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "cull_first_expanded_window.pdf"), surface_plot,
  width = 12, height = 7, units = "in", device = cairo_pdf
)
print(coverage)
print(performance)
print(degradation)
