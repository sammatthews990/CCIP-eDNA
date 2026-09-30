# Script to generate 4-panel figure comparing F1 vs F2 optimization surfaces and confusion matrices
# Using 5-fold CV (R=100) and Isotonic Monotonic Trajectories (isoreg) to ensure operational monotonicity
# Panels A & B: F1 and F2 Optimization Surfaces with Monotonic CV Threshold Trajectories
# Panels C & D: F1-optimised vs F2-optimised Confusion Matrices at CPUE >= 0.04
# Also exports plots/cv_performance_summary_table.csv aligned with 5-fold CV monotonic thresholds

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(fuzzyjoin)
  library(broom)
  library(MASS)
  library(glmmTMB)
  library(tidyr)
  library(patchwork)
  library(rsample)
  library(purrr)
  library(tibble)
  library(viridis)
  library(scales)
  library(readr)
})

# Load reefDNA package functions safely
pkg_path <- if (dir.exists("reefDNA")) "reefDNA" else if (dir.exists("../reefDNA")) "../reefDNA" else "."
if (dir.exists(pkg_path)) {
  suppressMessages(devtools::load_all(pkg_path, quiet = TRUE))
}

# Set clean global plot theme
theme_set(theme_bw(base_family = "Helvetica") + theme(panel.grid.minor = element_blank()))

# CPUE target for the confusion matrices (panels C & D); default 0.04. Usage: Rscript <script> 0.02
cli_args <- commandArgs(trailingOnly = TRUE)
cpue_target <- if (length(cli_args) > 0) as.numeric(cli_args[1]) else 0.04
target_suffix <- if (cpue_target == 0.04) "" else sprintf("_cpue%03d", round(cpue_target * 100))

# 1. Load Data
cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

# Clean Cull Data
cull <- cull.dat %>%
  rename(Reef = ReefName) %>%
  mutate(date_cull = as.Date(SurveyDate))

# Aggregate eDNA data
edna_agg <- edna.dat %>%
  filter(!is.na(Year)) %>%
  rename(Collection.organisation = `Collection organisation`) %>%
  mutate(
    Reef = ReefName,
    Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
    date_edna = as.Date(Date),
    Conc_mean = as.numeric(Conc_mean)
  ) %>%
  arrange(Reef, Year, date_edna) %>%
  group_by(Reef, Collection.org, Year) %>%
  mutate(
    grp = cumsum(
      if_else(
        is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > 7,
        1L, 0L
      )
    )
  ) %>%
  ungroup() %>%
  group_by(Reef, Collection.org, Year, grp) %>%
  summarise(
    date_edna = min(date_edna),
    conc_mean = mean(Conc_mean, na.rm = TRUE),
    perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
    n_samples = n(),
    .groups   = "drop"
  )

# Extract multi-horizon cohorts
reef_prior_3m  <- get_prior_cohort(cull, edna_agg, 91,  "0-3 Months", min_days = 0)
reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "3-6 Months", min_days = 92)
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "6-12 Months", min_days = 184)

dat_glmm <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
  mutate(
    conc_t    = log1p(conc_mean_reef),
    obs_cpue  = total_cots / total_bottom,
    horizon   = factor(horizon, levels = c("0-3 Months", "3-6 Months", "6-12 Months"))
  ) %>%
  filter(
    !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
    !is.na(perc_pos_reef), !is.na(conc_t), !is.na(Reef),
    conc_mean_reef < 5000
  ) %>%
  mutate(Reef = factor(Reef)) %>%
  droplevels()

# Define event dataset (using prior horizons <= 6 months aggregated by Reef, Collection.org, Year)
dat_evt <- dat_glmm %>%
  filter(horizon %in% c("0-3 Months", "3-6 Months")) %>%
  group_by(Reef, Collection.org, Year) %>%
  summarise(
    counts = sum(total_cots, na.rm = TRUE),
    total_bottom = sum(total_bottom, na.rm = TRUE),
    perc_pos = mean(perc_pos_reef, na.rm = TRUE),
    conc_mean_reef = mean(conc_mean_reef, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    reef = Reef,
    cpue = counts / total_bottom
  )

# ----------------------------------------------------
# 2D RASTER GRID: F1 & F2 Surfaces
# ----------------------------------------------------
perc_grid <- seq(0, 100, by = 2)
cpue_grid <- seq(0, quantile(dat_evt$cpue, 0.99, na.rm = TRUE), length.out = 40)

f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (beta^2 * prec + rec == 0)) return(0)
  (1 + beta^2) * (prec * rec) / ((beta^2 * prec) + rec)
}

grid <- tidyr::expand_grid(perc_thresh = perc_grid, cpue_thresh = cpue_grid) %>%
  mutate(
    out = purrr::pmap(., ~ metrics_for(..1, ..2, dat_evt)),
    F1 = vapply(out, `[[`, numeric(1), "F1"),
    Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
    Precision = vapply(out, `[[`, numeric(1), "Precision"),
    Recall = vapply(out, `[[`, numeric(1), "Recall"),
    F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2))
  ) %>%
  dplyr::select(-out)

# ----------------------------------------------------
# REPEATED 5-FOLD CROSS VALIDATION & MONOTONIC SMOOTHING
# ----------------------------------------------------
cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
k <- 5
R <- 100
set.seed(42)

run_5fold_tuning <- function(dat_evt, cpue_thr, perc_grid, k = 5, R = 100) {
  dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_long <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    dplyr::select(id, id2, assess) %>%
    tidyr::expand_grid(perc_thresh = perc_grid) %>%
    mutate(
      out = pmap(list(perc_thresh, assess), \(p, d) metrics_for(p, cpue_thr, d)),
      F1 = vapply(out, `[[`, numeric(1), "F1"),
      Precision = vapply(out, `[[`, numeric(1), "Precision"),
      Recall = vapply(out, `[[`, numeric(1), "Recall")
    ) %>%
    mutate(F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2))) %>%
    dplyr::select(-out, -assess)

  cv_sum <- cv_long %>%
    group_by(perc_thresh) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F2_mean = mean(F2, na.rm = TRUE),
      Rec_mean = mean(Recall, na.rm = TRUE),
      .groups = "drop"
    )

  best_f1 <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
  best_f2 <- cv_sum %>% slice_max(F2_mean, n = 1, with_ties = FALSE)

  tibble(
    cpue_thr = cpue_thr,
    perc_star_f1_raw = best_f1$perc_thresh,
    f1_val = best_f1$F1_mean,
    perc_star_f2_raw = best_f2$perc_thresh,
    f2_val = best_f2$F2_mean
  )
}

tuning_5fold <- map_dfr(cpue_levels, ~ run_5fold_tuning(dat_evt, .x, perc_grid, k = k, R = R)) %>%
  mutate(
    # Apply Isotonic Monotonic Regression to guarantee operational consistency
    perc_star_f1 = isoreg(cpue_thr, perc_star_f1_raw)$yf,
    perc_star_f2 = isoreg(cpue_thr, perc_star_f2_raw)$yf
  )

# ----------------------------------------------------
# GENERATE & EXPORT CV PERFORMANCE SUMMARY TABLE
# ----------------------------------------------------
eval_targets <- c(0.02, 0.04, 0.08)

get_cv_stats <- function(c_thr, strategy, p_thresh) {
  dat2 <- dat_evt %>% mutate(actual = if (c_thr == 0) (cpue > 0) else (cpue >= c_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_eval <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    mutate(
      out = map(assess, ~ metrics_for(p_thresh, c_thr, .x)),
      F1 = vapply(out, `[[`, numeric(1), "F1"),
      Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
      Precision = vapply(out, `[[`, numeric(1), "Precision"),
      Recall = vapply(out, `[[`, numeric(1), "Recall")
    ) %>%
    mutate(F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2)))

  se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  fmt <- function(m, s) sprintf("%.3f ± %.3f", m, s)

  tibble(
    `CPUE Target` = c_thr,
    `Optimization Strategy` = strategy,
    `% eDNA Positive Threshold` = sprintf("%d%%", round(p_thresh)),
    `F1 Score` = fmt(mean(cv_eval$F1, na.rm = TRUE), se(cv_eval$F1)),
    `F2 Score` = fmt(mean(cv_eval$F2, na.rm = TRUE), se(cv_eval$F2)),
    `Accuracy` = fmt(mean(cv_eval$Accuracy, na.rm = TRUE), se(cv_eval$Accuracy)),
    `Precision` = fmt(mean(cv_eval$Precision, na.rm = TRUE), se(cv_eval$Precision)),
    `Recall` = fmt(mean(cv_eval$Recall, na.rm = TRUE), se(cv_eval$Recall))
  )
}

cv_table_summary <- map_dfr(eval_targets, function(ct) {
  row_tune <- tuning_5fold %>% filter(cpue_thr == ct)
  p1 <- row_tune$perc_star_f1
  p2 <- row_tune$perc_star_f2

  bind_rows(
    get_cv_stats(ct, "F1-Optimized (Balanced)", p1),
    get_cv_stats(ct, "F2-Optimized (Recall-Weighted)", p2)
  )
})

dir.create("plots", showWarnings = FALSE)
write_csv(cv_table_summary, "plots/cv_performance_summary_table.csv")

cat("\n=======================================================\n")
cat(" 5-FOLD CV MONOTONIC THRESHOLD PERFORMANCE TABLE\n")
cat("=======================================================\n")
print(cv_table_summary, n = 20)
cat("=======================================================\n\n")

# ----------------------------------------------------
# PANEL A: F1 Optimization Surface
# ----------------------------------------------------
panel_a <- ggplot(grid, aes(x = perc_thresh, y = cpue_thresh, fill = F1)) +
  geom_raster() +
  geom_contour(aes(z = F1), breaks = seq(0.2, 1, by = 0.2), colour = "white", linewidth = 0.35, alpha = 0.8) +
  geom_path(data = tuning_5fold, aes(x = perc_star_f1, y = cpue_thr), color = "#e31a1c", linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  geom_point(data = tuning_5fold, aes(x = perc_star_f1, y = cpue_thr), color = "#e31a1c", size = 2.8, inherit.aes = FALSE) +
  scale_fill_viridis_c(limits = c(0, 1), option = "viridis", name = "F1 Score") +
  labs(
    title = "A. F1-optimised threshold surface (Equal weight)",
    x = "eDNA % positive threshold",
    y = "CPUE threshold"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 9.5),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.position = "right"
  )

# ----------------------------------------------------
# PANEL B: F2 Optimization Surface (Recall-Weighted)
# ----------------------------------------------------
panel_b <- ggplot(grid, aes(x = perc_thresh, y = cpue_thresh, fill = F2)) +
  geom_raster() +
  geom_contour(aes(z = F2), breaks = seq(0.2, 1, by = 0.2), colour = "white", linewidth = 0.35, alpha = 0.8) +
  geom_path(data = tuning_5fold, aes(x = perc_star_f2, y = cpue_thr), color = "#ff7f00", linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  geom_point(data = tuning_5fold, aes(x = perc_star_f2, y = cpue_thr), color = "#ff7f00", size = 2.8, inherit.aes = FALSE) +
  scale_fill_viridis_c(limits = c(0, 1), option = "plasma", name = "F2 Score") +
  labs(
    title = "B. F2-optimised threshold surface (Recall-weighted)",
    x = "eDNA % positive threshold",
    y = "CPUE threshold"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 9.5),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.position = "right"
  )

# ----------------------------------------------------
# PANELS C & D: F1 vs F2 Confusion Matrices (CPUE >= 0.04)
# ----------------------------------------------------
perc_grid_fine <- seq(0, 100, by = 1)

grid_eval <- tibble(perc_thresh = perc_grid_fine) %>%
  mutate(
    out = map(perc_thresh, ~ metrics_for(.x, cpue_target, dat_evt)),
    F1 = vapply(out, `[[`, numeric(1), "F1"),
    Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
    Precision = vapply(out, `[[`, numeric(1), "Precision"),
    Recall = vapply(out, `[[`, numeric(1), "Recall"),
    F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2))
  )

# Extract Monotonic CV thresholds for CPUE >= 0.04
p_star_f1 <- round((tuning_5fold %>% filter(cpue_thr == cpue_target))$perc_star_f1)
p_star_f2 <- round((tuning_5fold %>% filter(cpue_thr == cpue_target))$perc_star_f2)

cm_f1_obj <- make_confusion(dat_evt, perc_thresh = p_star_f1, cpue_thresh = cpue_target)
cm_f2_obj <- make_confusion(dat_evt, perc_thresh = p_star_f2, cpue_thresh = cpue_target)

f1_metrics_row <- grid_eval %>% filter(perc_thresh == p_star_f1) %>% slice(1)
f2_metrics_row <- grid_eval %>% filter(perc_thresh == p_star_f2) %>% slice(1)

df_cm_f1 <- cm_f1_obj$cm %>%
  mutate(panel = sprintf("C. F1-Optimised Matrix (CPUE ≥ %.2f)\n(%%pos* = %d%% | F1 = %.3f | Recall = %.1f%%)", cpue_target, p_star_f1, f1_metrics_row$F1, f1_metrics_row$Recall * 100))

df_cm_f2 <- cm_f2_obj$cm %>%
  mutate(panel = sprintf("D. F2-Optimised Matrix (CPUE ≥ %.2f)\n(%%pos* = %d%% | F2 = %.3f | Recall = %.1f%%)", cpue_target, p_star_f2, f2_metrics_row$F2, f2_metrics_row$Recall * 100))

panel_c <- ggplot(df_cm_f1, aes(x = pred, y = actual, fill = row_prop)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4.0, fontface = "bold", color = "black") +
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Prop.") +
  labs(
    title = sprintf("C. F1-Optimised Matrix (CPUE ≥ %.2f, %%pos* = %d%%)", cpue_target, p_star_f1),
    x = "Prediction from eDNA (% pos)",
    y = "Reference Culling CPUE"
  ) +
  coord_fixed() +
  theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.title = element_text(face = "bold", size = 9),
    legend.position = "none"
  )

panel_d <- ggplot(df_cm_f2, aes(x = pred, y = actual, fill = row_prop)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4.0, fontface = "bold", color = "black") +
  scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Prop.") +
  labs(
    title = sprintf("D. F2-Optimised Matrix (CPUE ≥ %.2f, %%pos* = %d%%)", cpue_target, p_star_f2),
    x = "Prediction from eDNA (% pos)",
    y = "Reference Culling CPUE"
  ) +
  coord_fixed() +
  theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.title = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

# ----------------------------------------------------
# ASSEMBLE 4-PANEL COMPOSITE GRAPH
# ----------------------------------------------------
four_panel_plot <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(heights = c(1.1, 1))

out_file <- function(ext) sprintf("plots/multi_horizon_CPUE%s.%s", target_suffix, ext)
ggsave(out_file("png"), four_panel_plot, width = 11, height = 9.5, dpi = 300, bg = "white")
ggsave(out_file("pdf"), four_panel_plot, width = 11, height = 9.5, dpi = 300, device = cairo_pdf)
ggsave(out_file("eps"), four_panel_plot, width = 11, height = 9.5, dpi = 300)

cat(sprintf("Successfully generated %s and exported CSV table to plots/cv_performance_summary_table.csv!\n", out_file("png")))
