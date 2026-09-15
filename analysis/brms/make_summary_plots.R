#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggdist)
  library(ggplot2)
  library(patchwork)
  library(posterior)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_selected_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
plot_dir <- file.path(project_root, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path(output_dir, "brms_selected_final.rds"))
model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "brms_scale_spec.rds"))
loo_path <- file.path(output_dir, "brms_loo_targeted_comparison.csv")
if (!file.exists(loo_path)) {
  loo_path <- file.path(output_dir, "brms_loo_final_comparison.csv")
}
loo_comparison <- read.csv(loo_path)
kfold_comparison <- read.csv(file.path(output_dir, "brms_reef_kfold_comparison.csv"))

blue <- "#2166AC"
light_blue <- "#92C5DE"
orange <- "#D95F02"
light_orange <- "#FDB863"
dark_grey <- "#303030"

paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8.5, colour = dark_grey),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5),
    plot.margin = margin(7, 7, 7, 7)
  )

model_labels <- c(
  null = "Null",
  best_additive = "% pos. + distance + time",
  full_additive = "All additive",
  concentration_distance = "Conc. x distance",
  dual_interaction = "Two interactions",
  all_pairwise = "All pairwise"
)
loo_plot_data <- loo_comparison |>
  mutate(
    selected = model == "concentration_distance",
    model_label = factor(
      model_labels[model],
      levels = model_labels[c(
        "null", "full_additive", "best_additive",
        "concentration_distance", "dual_interaction", "all_pairwise"
      )]
    )
  )

panel_a <- ggplot(
  loo_plot_data,
  aes(x = elpd_diff, y = model_label, colour = selected)
) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff),
    orientation = "y", width = 0.18, linewidth = 0.55
  ) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(`FALSE` = blue, `TRUE` = orange), guide = "none") +
  labs(
    title = "Model performance",
    subtitle = "PSIS-LOO difference +/- SE; orange = selected",
    x = "Difference in expected log predictive density",
    y = NULL
  ) +
  paper_theme

coef_lookup <- c(
  "b_edna_prop_z" = "% positive",
  "b_edna_conc_z" = "Mean concentration",
  "b_distance_z" = "Distance",
  "b_lag_z" = "Time since sample",
  "b_edna_conc_z:distance_z" = "Concentration x distance"
)
coef_draws <- posterior::as_draws_df(fit) |>
  as.data.frame() |>
  select(any_of(names(coef_lookup))) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "estimate") |>
  mutate(
    term = factor(
      unname(coef_lookup[parameter]),
      levels = rev(unname(coef_lookup))
    )
  )

panel_b <- ggplot(coef_draws, aes(x = estimate, y = term, fill = term)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.45) +
  ggdist::stat_halfeye(
    .width = c(0.8, 0.95), point_interval = median_qi,
    slab_alpha = 0.8, normalize = "panels", colour = dark_grey
  ) +
  scale_fill_manual(
    values = c(orange, light_orange, "#BDBDBD", light_blue, blue),
    guide = "none"
  ) +
  labs(
    title = "Effect sizes",
    subtitle = "Posterior median, 80% and 95% credible intervals",
    x = "Standardized log-rate effect",
    y = NULL
  ) +
  paper_theme

z_from_raw <- function(x, variable) {
  (x - unname(scale_spec[[variable]]["center"])) /
    unname(scale_spec[[variable]]["scale"])
}

newdata_template <- function(n) {
  tibble(
    Reef = factor(rep(levels(model_data$Reef)[1], n), levels = levels(model_data$Reef)),
    bottom_time = rep(stats::median(model_data$bottom_time), n),
    log_effort = log(bottom_time),
    edna_conc_z = 0,
    distance_z = 0,
    lag_z = 0
  )
}

summarise_prediction_matrix <- function(x, effort) {
  x <- sweep(x, 2, effort, "/")
  tibble(
    estimate = apply(x, 2, median),
    lower = apply(x, 2, quantile, probs = 0.025),
    upper = apply(x, 2, quantile, probs = 0.975)
  )
}

set.seed(20260909)
positive_grid <- seq(0, 100, by = 1)
marginal_data <- newdata_template(length(positive_grid)) |>
  mutate(
    edna_pct = positive_grid,
    edna_prop_z = z_from_raw(edna_pct / 100, "edna_prop_model")
  )

expected_draws <- posterior_epred(
  fit, newdata = marginal_data, re_formula = NA, ndraws = 2000
)
predictive_draws <- posterior_predict(
  fit, newdata = marginal_data, re_formula = NA, ndraws = 2000
)
marginal_summary <- marginal_data |>
  select(edna_pct, bottom_time) |>
  bind_cols(summarise_prediction_matrix(expected_draws, marginal_data$bottom_time)) |>
  rename(expected = estimate, expected_lower = lower, expected_upper = upper) |>
  bind_cols(summarise_prediction_matrix(predictive_draws, marginal_data$bottom_time)) |>
  rename(predicted = estimate, prediction_lower = lower, prediction_upper = upper)

write.csv(
  marginal_summary,
  file.path(output_dir, "selected_model_percent_positive_curve.csv"),
  row.names = FALSE
)

panel_c <- ggplot(marginal_summary, aes(x = edna_pct)) +
  geom_ribbon(
    aes(ymin = prediction_lower, ymax = prediction_upper, fill = "Prediction"),
    alpha = 0.25, colour = NA
  ) +
  geom_ribbon(
    aes(ymin = expected_lower, ymax = expected_upper, fill = "Expected"),
    alpha = 0.55, colour = NA
  ) +
  geom_line(aes(y = expected), colour = blue, linewidth = 0.9) +
  scale_fill_manual(
    values = c(Expected = light_blue, Prediction = light_orange),
    breaks = c("Expected", "Prediction"),
    labels = c("95% credible", "95% prediction")
  ) +
  labs(
    title = "% positive",
    subtitle = "Other predictors at their transformed means",
    x = "eDNA samples positive (%)",
    y = "Predicted COTS CPUE (min^-1)",
    fill = NULL
  ) +
  paper_theme +
  theme(legend.position = "top", legend.justification = "left")

expected_rate_draws <- sweep(
  expected_draws, 2, marginal_data$bottom_time, "/"
)
increment_draws <- sweep(
  expected_rate_draws, 1, expected_rate_draws[, 1], "-"
)
increment_summary <- tibble(edna_pct = positive_grid) |>
  bind_cols(summarise_prediction_matrix(
    increment_draws, rep(1, length(positive_grid))
  ))
write.csv(
  increment_summary,
  file.path(output_dir, "selected_model_percent_positive_uplift.csv"),
  row.names = FALSE
)

panel_c_operational <- ggplot(
  increment_summary, aes(x = edna_pct, y = estimate)
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper), fill = light_blue,
    alpha = 0.6, colour = NA
  ) +
  geom_line(colour = blue, linewidth = 0.9) +
  labs(
    title = "% positive uplift",
    subtitle = "Increase above the modelled 0% baseline",
    x = "eDNA samples positive (%)",
    y = "Additional COTS CPUE (min^-1)"
  ) +
  paper_theme

context_grid <- tidyr::crossing(
  edna_pct = seq(0, 100, by = 2),
  distance_m = c(200, 2000),
  lag_days = c(30, 183)
)
context_data <- newdata_template(nrow(context_grid)) |>
  bind_cols(context_grid) |>
  mutate(
    edna_prop_z = z_from_raw(edna_pct / 100, "edna_prop_model"),
    distance_z = z_from_raw(log1p(distance_m / 200), "distance_model"),
    lag_z = z_from_raw(log1p(lag_days), "lag_model"),
    distance_label = factor(
      distance_m, levels = c(200, 2000), labels = c("200 m", "2,000 m")
    ),
    lag_label = factor(
      lag_days, levels = c(30, 183), labels = c("30 days", "183 days")
    )
  )

context_draws <- posterior_epred(
  fit, newdata = context_data, re_formula = NA, ndraws = 2000
)
context_summary <- context_data |>
  select(edna_pct, distance_m, lag_days, distance_label, lag_label, bottom_time) |>
  bind_cols(summarise_prediction_matrix(context_draws, context_data$bottom_time))

write.csv(
  context_summary,
  file.path(output_dir, "selected_model_context_curves.csv"),
  row.names = FALSE
)

panel_d <- ggplot(context_summary, aes(x = edna_pct, y = estimate)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper), fill = light_blue,
    alpha = 0.55, colour = NA
  ) +
  geom_line(colour = blue, linewidth = 0.75) +
  facet_grid(lag_label ~ distance_label) +
  labs(
    title = "Context dependence",
    subtitle = "Rows: time since sample; columns: distance",
    x = "eDNA samples positive (%)",
    y = "Expected COTS CPUE (min^-1)"
  ) +
  paper_theme

individual_plots <- list(
  panel_model_performance = panel_a,
  panel_effect_halfeye = panel_b,
  panel_percent_positive = panel_c,
  panel_percent_positive_uplift = panel_c_operational,
  panel_context_curves = panel_d
)
for (plot_name in names(individual_plots)) {
  ggsave(
    file.path(output_dir, paste0(plot_name, ".png")),
    individual_plots[[plot_name]], width = 7, height = 5.2,
    units = "in", dpi = 300, bg = "white"
  )
}

four_panel <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(
  file.path(plot_dir, "brms_cpue_four_panel.png"), four_panel,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "brms_cpue_four_panel.pdf"), four_panel,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)

four_panel_operational <-
  (panel_a | panel_b) / (panel_c_operational | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "brms_cpue_four_panel_operational.png"),
  four_panel_operational,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "brms_cpue_four_panel_operational.pdf"),
  four_panel_operational,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)

fixed_summary <- as.data.frame(
  fixef(fit, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
) |>
  tibble::rownames_to_column("term") |>
  mutate(term = recode(
    term,
    Intercept = "Intercept",
    edna_prop_z = "% positive",
    edna_conc_z = "Mean concentration",
    distance_z = "Distance",
    lag_z = "Time since sample",
    `edna_conc_z:distance_z` = "Concentration x distance"
  ))
write.csv(
  fixed_summary,
  file.path(output_dir, "selected_model_fixed_effects.csv"),
  row.names = FALSE
)

fit_summary <- summary(fit)
diagnostic_rows <- bind_rows(
  as.data.frame(fit_summary$fixed) |>
    tibble::rownames_to_column("parameter") |>
    mutate(block = "fixed", .before = 1),
  as.data.frame(fit_summary$spec_pars) |>
    tibble::rownames_to_column("parameter") |>
    mutate(block = "distributional", .before = 1),
  as.data.frame(fit_summary$random$Reef) |>
    tibble::rownames_to_column("parameter") |>
    mutate(block = "reef", .before = 1)
)
nuts <- nuts_params(fit)
final_diagnostics <- tibble(
  max_rhat = max(diagnostic_rows$Rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostic_rows$Bulk_ESS, na.rm = TRUE),
  min_tail_ess = min(diagnostic_rows$Tail_ESS, na.rm = TRUE),
  divergences = sum(nuts$Value[nuts$Parameter == "divergent__"]),
  max_treedepth_hits = sum(nuts$Value[nuts$Parameter == "treedepth__"] >= 12)
)
write.csv(
  final_diagnostics,
  file.path(output_dir, "selected_model_diagnostics.csv"),
  row.names = FALSE
)

selection_summary <- tibble(
  selected_model = "concentration_distance",
  selection_rule = "whole-reef K-fold performance, then parsimony",
  loo_delta_vs_all_pairwise = loo_comparison$elpd_diff[
    loo_comparison$model == "concentration_distance"
  ],
  loo_delta_se = loo_comparison$se_diff[
    loo_comparison$model == "concentration_distance"
  ],
  reef_kfold_delta_vs_all_pairwise = kfold_comparison$elpd_diff[
    kfold_comparison$model == "concentration_distance"
  ],
  reef_kfold_delta_se = kfold_comparison$se_diff[
    kfold_comparison$model == "concentration_distance"
  ]
)
write.csv(
  selection_summary,
  file.path(output_dir, "selected_model_summary.csv"),
  row.names = FALSE
)

message("Saved the four-panel paper figure and all supporting summaries.")
