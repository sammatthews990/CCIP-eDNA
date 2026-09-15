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
  if (file.exists("analysis/site/output/site_visit_brms_final.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
plot_dir <- file.path(project_root, "plots")
fit <- readRDS(file.path(output_dir, "site_visit_brms_final.rds"))
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "site_visit_scale_spec.rds"))
glmm_curves <- read.csv(file.path(output_dir, "site_visit_glmm_curves.csv")) |>
  mutate(model = recode(
    model,
    reef = "Reef", reef_site = "Reef + site",
    reef_campaign = "Reef + campaign",
    reef_site_campaign = "Reef + site + campaign"
  ))

blue <- "#2166AC"
light_blue <- "#92C5DE"
orange <- "#D95F02"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(), legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5)
  )

panel_a <- ggplot(
  glmm_curves,
  aes(edna_pct, marginal_new_group_cpue, colour = model)
) +
  geom_hline(yintercept = 0.04, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.8) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title = "Grouping structure", subtitle = "Fast GLMM screen",
    x = "eDNA samples positive (%)", y = "Expected CPUE (min^-1)", colour = NULL
  ) + paper_theme

effect_lookup <- c(
  b_edna_prop_z = "% positive", b_distance_z = "Distance",
  b_lag_z = "Time since sample"
)
effect_draws <- as_draws_df(fit) |>
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
    title = "Effect sizes", subtitle = "Median, 80% and 95% credible intervals",
    x = "Standardized log-rate effect", y = NULL
  ) + paper_theme

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
    edna_campaign_id = rep(as.character(data$edna_campaign_id[[1]]), n),
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
  draws$sd_site_id__Intercept^2 + draws$sd_edna_campaign_id__Intercept^2
shape <- draws$shape

positive_data <- template_data(101) |>
  mutate(
    edna_pct = 0:100,
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model")
  )
conditional_counts <- posterior_epred(fit, newdata = positive_data, re_formula = NA)
marginal_cpue <- sweep(conditional_counts, 1, exp(0.5 * total_variance), "*") / effort
set.seed(20261013)
new_group_effect <- rnorm(nrow(draws), 0, sqrt(total_variance))
future_mu <- sweep(conditional_counts, 1, exp(new_group_effect), "*")
future_cpue <- matrix(
  rnbinom(
    length(future_mu), mu = as.vector(future_mu),
    size = rep(shape, ncol(future_mu))
  ),
  nrow = nrow(future_mu)
) / effort
positive_summary <- positive_data |>
  select(edna_pct) |>
  bind_cols(summarise_draws(marginal_cpue)) |>
  rename(expected = estimate, expected_lower = lower, expected_upper = upper) |>
  bind_cols(summarise_draws(future_cpue)) |>
  rename(predicted = estimate, prediction_lower = lower, prediction_upper = upper)
write.csv(
  positive_summary,
  file.path(output_dir, "site_visit_brms_percent_positive.csv"), row.names = FALSE
)
panel_c <- ggplot(positive_summary, aes(edna_pct)) +
  geom_ribbon(
    aes(ymin = prediction_lower, ymax = prediction_upper, fill = "Prediction"),
    alpha = 0.22
  ) +
  geom_ribbon(
    aes(ymin = expected_lower, ymax = expected_upper, fill = "Expected"),
    alpha = 0.55
  ) +
  geom_line(aes(y = expected), colour = blue, linewidth = 0.9) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey50") +
  scale_fill_manual(values = c(Expected = light_blue, Prediction = "#FDB863")) +
  labs(
    title = "% positive", subtitle = "New reef, site and campaign",
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
context_marginal <- sweep(
  context_conditional, 1, exp(0.5 * total_variance), "*"
) / effort
context_summary <- context_data |>
  select(edna_pct, distance_m, lag_days, distance_label, lag_label) |>
  bind_cols(summarise_draws(context_marginal))
write.csv(
  context_summary,
  file.path(output_dir, "site_visit_brms_context.csv"), row.names = FALSE
)
panel_d <- ggplot(context_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.55) +
  geom_line(colour = blue, linewidth = 0.75) +
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
  file.path(plot_dir, "site_visit_brms_four_panel.png"), figure,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "site_visit_brms_four_panel.pdf"), figure,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)
cat("0.04 crossing:", approx(positive_summary$expected, positive_summary$edna_pct, xout = 0.04)$y, "\n")
