#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggdist)
  library(ggplot2)
  library(glmmTMB)
  library(patchwork)
  library(posterior)
  library(splines)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/cull_first/output/cull_first_brms_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "cull_first", "output")
site_dir <- file.path(project_root, "analysis", "site", "output")
plot_dir <- file.path(project_root, "plots")
fit <- readRDS(file.path(output_dir, "cull_first_brms_final.rds"))
data <- readRDS(file.path(output_dir, "cull_first_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "cull_first_scale_spec.rds"))

blue <- "#2166AC"
light_blue <- "#92C5DE"
orange <- "#D95F02"
green <- "#1B9E77"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(), legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5)
  )
z_value <- function(raw, variable, transform = identity) {
  transformed <- transform(raw)
  (transformed - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}
effort <- median(data$bottom_time)
template_data <- function(n) {
  tibble(
    Reef = rep(as.character(data$Reef[[1]]), n),
    site_id = rep(as.character(data$site_id[[1]]), n),
    event_id = rep(as.character(data$event_id[[1]]), n),
    bottom_time = effort, log_effort = log(effort),
    distance_z = 0, lag_z = 0
  )
}
summarise_draws <- function(values) {
  tibble(
    estimate = apply(values, 2, median),
    lower = apply(values, 2, quantile, probs = 0.025),
    upper = apply(values, 2, quantile, probs = 0.975)
  )
}
draws <- as_draws_df(fit)
total_variance <- draws$sd_Reef__Intercept^2 +
  draws$sd_site_id__Intercept^2 + draws$sd_event_id__Intercept^2
shape <- draws$shape

positive_data <- template_data(101) |>
  mutate(
    edna_pct = 0:100,
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model")
  )
conditional_counts <- posterior_epred(fit, newdata = positive_data, re_formula = NA)
marginal_cpue <- sweep(conditional_counts, 1, exp(0.5 * total_variance), "*") / effort
set.seed(20261022)
future_effect <- rnorm(nrow(draws), 0, sqrt(total_variance))
future_mu <- sweep(conditional_counts, 1, exp(future_effect), "*")
future_cpue <- matrix(
  rnbinom(
    length(future_mu), mu = as.vector(future_mu),
    size = rep(shape, ncol(future_mu))
  ),
  nrow = nrow(future_mu)
) / effort
cull_curve <- positive_data |>
  select(edna_pct) |>
  bind_cols(summarise_draws(marginal_cpue)) |>
  rename(expected = estimate, expected_lower = lower, expected_upper = upper) |>
  bind_cols(summarise_draws(future_cpue)) |>
  rename(predicted = estimate, prediction_lower = lower, prediction_upper = upper)
write.csv(
  cull_curve,
  file.path(output_dir, "cull_first_brms_percent_positive.csv"), row.names = FALSE
)

site_curve <- read.csv(file.path(site_dir, "site_visit_brms_percent_positive.csv"))
linkage_curves <- bind_rows(
  transmute(cull_curve, edna_pct, expected, lower = expected_lower, upper = expected_upper,
    linkage = "Cull-first dives"),
  transmute(site_curve, edna_pct, expected, lower = expected_lower, upper = expected_upper,
    linkage = "eDNA-first visits")
)
write.csv(
  linkage_curves,
  file.path(output_dir, "linkage_assumption_curves.csv"), row.names = FALSE
)
panel_a <- ggplot(linkage_curves, aes(edna_pct, expected, colour = linkage, fill = linkage)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.04, linetype = 2, colour = "grey45") +
  scale_colour_manual(values = c("Cull-first dives" = orange, "eDNA-first visits" = blue)) +
  scale_fill_manual(values = c("Cull-first dives" = orange, "eDNA-first visits" = blue)) +
  labs(
    title = "Linkage assumptions", subtitle = "New-group marginal predictions",
    x = "eDNA samples positive (%)", y = "Expected CPUE (min^-1)",
    colour = NULL, fill = NULL
  ) + paper_theme

effect_lookup <- c(
  b_edna_prop_z = "% positive", b_distance_z = "Distance",
  b_lag_z = "Time since sample"
)
effect_draws <- draws |>
  as.data.frame() |>
  select(all_of(names(effect_lookup))) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "estimate") |>
  mutate(term = factor(unname(effect_lookup[parameter]), levels = rev(unname(effect_lookup))))
panel_b <- ggplot(effect_draws, aes(estimate, term, fill = term)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
  stat_halfeye(
    .width = c(0.8, 0.95), point_interval = median_qi,
    slab_alpha = 0.8, normalize = "panels"
  ) +
  scale_fill_manual(values = c(orange, light_blue, blue), guide = "none") +
  labs(
    title = "Cull-first effects", subtitle = "Median, 80% and 95% intervals",
    x = "Standardized log-rate effect", y = NULL
  ) + paper_theme

panel_c <- ggplot(cull_curve, aes(edna_pct)) +
  geom_ribbon(aes(ymin = prediction_lower, ymax = prediction_upper, fill = "Prediction"), alpha = 0.22) +
  geom_ribbon(aes(ymin = expected_lower, ymax = expected_upper, fill = "Expected"), alpha = 0.55) +
  geom_line(aes(y = expected), colour = orange, linewidth = 0.9) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey50") +
  scale_fill_manual(values = c(Expected = light_blue, Prediction = "#FDB863")) +
  labs(
    title = "% positive", subtitle = "New reef, site and eDNA event",
    x = "eDNA samples positive (%)", y = "Predicted CPUE (min^-1)", fill = NULL
  ) + paper_theme

context_grid <- crossing(
  edna_pct = seq(0, 100, by = 2),
  distance_m = c(200, 2000), lag_days = c(30, 183)
)
context_data <- template_data(nrow(context_grid)) |>
  select(-distance_z, -lag_z) |>
  bind_cols(context_grid) |>
  mutate(
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model"),
    distance_z = z_value(distance_m, "distance_model", function(x) log1p(x / 200)),
    lag_z = z_value(lag_days, "lag_model", log1p),
    distance_label = factor(distance_m, c(200, 2000), c("200 m", "2,000 m")),
    lag_label = factor(lag_days, c(30, 183), c("30 days", "183 days"))
  )
context_conditional <- posterior_epred(fit, newdata = context_data, re_formula = NA)
context_marginal <- sweep(context_conditional, 1, exp(0.5 * total_variance), "*") / effort
context_summary <- context_data |>
  select(edna_pct, distance_m, lag_days, distance_label, lag_label) |>
  bind_cols(summarise_draws(context_marginal))
write.csv(context_summary, file.path(output_dir, "cull_first_brms_context.csv"), row.names = FALSE)
panel_d <- ggplot(context_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.55) +
  geom_line(colour = orange, linewidth = 0.75) +
  geom_hline(yintercept = 0.04, linetype = 2, colour = "#7B3294", linewidth = 0.4) +
  facet_grid(lag_label ~ distance_label) +
  labs(
    title = "Context", subtitle = "Rows: time; columns: distance",
    x = "eDNA samples positive (%)", y = "Expected CPUE (min^-1)"
  ) + paper_theme

figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "cull_first_brms_four_panel.png"), figure,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "cull_first_brms_four_panel.pdf"), figure,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)

# Add fixed-effect uncertainty to the already fitted expanded-window surface.
expanded <- readRDS(file.path(output_dir, "expanded_window_glmm_screen.rds"))
expanded_grid <- read.csv(file.path(output_dir, "expanded_window_prediction_grid.csv"))
expanded_fit <- expanded$fits[[expanded$best_model]]
expanded_prediction <- predict(
  expanded_fit, newdata = expanded_grid,
  type = "link", re.form = NA, se.fit = TRUE
)
expanded_random_variance <- sum(vapply(
  VarCorr(expanded_fit)$cond,
  function(component) unname(attr(component, "stddev")[[1]])^2,
  numeric(1)
))
expanded_grid <- expanded_grid |>
  mutate(
    expected_cpue = exp(expanded_prediction$fit + 0.5 * expanded_random_variance),
    lower = exp(expanded_prediction$fit - 1.96 * expanded_prediction$se.fit + 0.5 * expanded_random_variance),
    upper = exp(expanded_prediction$fit + 1.96 * expanded_prediction$se.fit + 0.5 * expanded_random_variance)
  )
write.csv(
  expanded_grid,
  file.path(output_dir, "expanded_window_prediction_grid.csv"), row.names = FALSE
)
expanded_plot <- ggplot(
  filter(expanded_grid, edna_pct %in% c(25, 50, 75, 100)),
  aes(distance_m, expected_cpue, colour = factor(edna_pct), fill = factor(edna_pct), group = edna_pct)
) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey65") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.07, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~lag_days, labeller = label_both) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  scale_fill_viridis_d(option = "C", end = 0.9) +
  labs(
    title = "Expanded-window signal", subtitle = "Exploratory 5 km / 365-day GLMM",
    x = "Distance from eDNA site (m)", y = "Expected CPUE (min^-1)", colour = "% positive"
  ) + paper_theme
ggsave(
  file.path(plot_dir, "cull_first_expanded_window.png"), expanded_plot,
  width = 12, height = 7, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "cull_first_expanded_window.pdf"), expanded_plot,
  width = 12, height = 7, units = "in", device = cairo_pdf
)

expanded_support <- expanded$data |>
  mutate(
    distance_band = cut(
      edna_distance_m,
      breaks = c(0, 200, 500, 1000, 2000, 5000),
      include.lowest = TRUE,
      labels = c("0-200", "200-500", "500-1,000", "1,000-2,000", "2,000-5,000")
    ),
    lag_band = cut(
      edna_lag_days,
      breaks = c(-1, 30, 90, 183, 270, 365),
      labels = c("0-30", "31-90", "91-183", "184-270", "271-365")
    )
  ) |>
  count(lag_band, distance_band, name = "n_dives") |>
  complete(lag_band, distance_band, fill = list(n_dives = 0))
write.csv(
  expanded_support,
  file.path(output_dir, "expanded_window_support.csv"), row.names = FALSE
)
support_plot <- ggplot(
  expanded_support, aes(distance_band, lag_band, fill = n_dives)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_dives), size = 3) +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Observation support", x = "Distance band (m)",
    y = "Elapsed-time band (days)", fill = "Cull dives"
  ) +
  paper_theme
expanded_assessment <- expanded_plot / support_plot +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "cull_first_expanded_window_assessment.png"),
  expanded_assessment,
  width = 12, height = 11, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "cull_first_expanded_window_assessment.pdf"),
  expanded_assessment,
  width = 12, height = 11, units = "in", device = cairo_pdf
)

cat("Cull-first 0.04 crossing:", approx(cull_curve$expected, cull_curve$edna_pct, xout = 0.04)$y, "\n")
cat("eDNA-first 0.04 crossing:", approx(site_curve$expected, site_curve$edna_pct, xout = 0.04)$y, "\n")
