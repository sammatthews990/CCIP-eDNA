#!/usr/bin/env Rscript

# Script to generate 4-panel composite figure comparing BRMS (top row) and BRT (bottom row) site-visit models
# Top-Left (A): BRMS Effect sizes half-eye plot (80% & 95% CI)
# Top-Right (B): BRMS % positive expected CPUE (60 days vs 270 days, 80% CI)
# Bottom-Left (C): BRT Variable importance bar plot
# Bottom-Right (D): BRT % positive partial-dependence curve (95% bootstrap CI)

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

fit_brms <- readRDS(file.path(output_dir, "site_visit_brms_final.rds"))
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "site_visit_scale_spec.rds"))

z_value <- function(raw, variable, transform = identity) {
  transformed <- transform(raw)
  (transformed - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}
effort <- median(data$bottom_time)

paper_theme <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

# ==========================================
# PANEL A (Top-Left): BRMS Effect Sizes
# ==========================================
effect_lookup <- c(
  b_edna_prop_z = "% positive",
  b_distance_z = "Distance",
  b_lag_z = "Time since sample"
)
effect_draws <- as_draws_df(fit_brms) |>
  as.data.frame() |>
  select(all_of(names(effect_lookup))) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "estimate") |>
  mutate(term = factor(unname(effect_lookup[parameter]), levels = rev(unname(effect_lookup))))

panel_a <- ggplot(effect_draws, aes(estimate, term, fill = term)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
  stat_halfeye(
    .width = c(0.8, 0.95), point_interval = median_qi,
    slab_alpha = 0.8, normalize = "panels"
  ) +
  scale_fill_manual(values = c("#D95F02", "#92C5DE", "#2166AC"), guide = "none") +
  labs(
    title = "A. BRMS GLMM: Effect Sizes",
    subtitle = "Median, 80% and 95% credible intervals",
    x = "Standardized log-rate effect", y = NULL
  ) + paper_theme

# ==========================================
# PANEL B (Top-Right): BRMS % Positive (60 vs 270 Days, 80% CI)
# ==========================================
draws <- as_draws_df(fit_brms)
total_variance <- draws$sd_Reef__Intercept^2 +
  draws$sd_site_id__Intercept^2 + draws$sd_edna_campaign_id__Intercept^2

get_brms_curve <- function(lag_days, label) {
  lag_z_val <- z_value(lag_days, "lag_model", log1p)
  df_grid <- tibble(
    Reef = rep(as.character(data$Reef[[1]]), 101),
    site_id = rep(as.character(data$site_id[[1]]), 101),
    edna_campaign_id = rep(as.character(data$edna_campaign_id[[1]]), 101),
    bottom_time = effort, log_effort = log(effort),
    distance_z = 0, lag_z = lag_z_val,
    edna_pct = 0:100,
    edna_prop_z = z_value((0:100) / 100, "edna_prop_model")
  )
  cond_counts <- posterior_epred(fit_brms, newdata = df_grid, re_formula = NA)
  marg_cpue <- sweep(cond_counts, 1, exp(0.5 * total_variance), "*") / effort
  
  tibble(
    edna_pct = 0:100,
    lag_group = label,
    expected = apply(marg_cpue, 2, median),
    lower80 = apply(marg_cpue, 2, quantile, probs = 0.10),
    upper80 = apply(marg_cpue, 2, quantile, probs = 0.90)
  )
}

brms_curves <- bind_rows(
  get_brms_curve(60, "60 days (Recent)"),
  get_brms_curve(270, "270 days (Aged)")
) |> mutate(lag_group = factor(lag_group, levels = c("60 days (Recent)", "270 days (Aged)")))

panel_b <- ggplot(brms_curves, aes(x = edna_pct, y = expected, color = lag_group, fill = lag_group)) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.22, color = NA) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey50") +
  scale_color_manual(values = c("60 days (Recent)" = "#1f78b4", "270 days (Aged)" = "#e31a1c"), name = "Time since sample") +
  scale_fill_manual(values = c("60 days (Recent)" = "#1f78b4", "270 days (Aged)" = "#e31a1c"), name = "Time since sample") +
  labs(
    title = "B. BRMS GLMM: % Positive Response by Sample Age",
    subtitle = "Marginal expected CPUE with 80% credible intervals",
    x = "eDNA samples positive (%)", y = "Predicted CPUE (min^-1)"
  ) + paper_theme +
  theme(
    legend.position = c(0.28, 0.78),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80", linewidth = 0.3),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 8)
  )

# ==========================================
# PANEL C (Bottom-Left): BRT Variable Importance
# ==========================================
brt_importance <- read.csv(file.path(output_dir, "site_visit_brt_importance.csv")) |>
  mutate(variable = factor(variable, levels = rev(c("% positive", "Distance", "Time since sample"))))

panel_c <- ggplot(brt_importance, aes(importance, variable, fill = variable)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c("#D95F02", "#92C5DE", "#2166AC"), guide = "none") +
  labs(
    title = "C. BRT: Variable Importance",
    subtitle = "Relative influence (%)",
    x = "Relative influence (%)", y = NULL
  ) + paper_theme

# ==========================================
# PANEL D (Bottom-Right): BRT % Positive Response
# ==========================================
brt_pos_summary <- read.csv(file.path(output_dir, "site_visit_brt_percent_positive.csv"))

panel_d <- ggplot(brt_pos_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#92C5DE", alpha = 0.5) +
  geom_line(colour = "#2166AC", linewidth = 0.9) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey50") +
  labs(
    title = "D. BRT: % Positive Response",
    subtitle = "Whole-reef bootstrap 95% confidence interval",
    x = "eDNA samples positive (%)", y = "Partial-dependence CPUE (min^-1)"
  ) + paper_theme

# ==========================================
# COMBINE COMPOSITE FIGURE (2x2)
# ==========================================
composite_fig <- (panel_a | panel_b) / (panel_c | panel_d)

ggsave(
  file.path(plot_dir, "site_visit_composite_brms_brt.png"), composite_fig,
  width = 13, height = 9.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "site_visit_composite_brms_brt.pdf"), composite_fig,
  width = 13, height = 9.5, units = "in", device = cairo_pdf
)

cat("Successfully saved site_visit_composite_brms_brt.png and .pdf!\n")
