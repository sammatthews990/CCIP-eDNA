# Script to generate 4-panel figure using SITE-VISIT data within 1000m and 6 months (183 days)
# Uses the same eDNA-first site-visit linkage as the BRMS/BRT models (analysis/site/audit_edna_first_linkage.R):
# one row per first subsequent cull-site visit, linked to the nearest eDNA site in its sampling campaign,
# restricted here to visits within 1000m of that eDNA site
# Mirrors plots/generate_two_panel_multi_horizon.R (reef level), but CV folds hold out whole eDNA campaigns
# because several visits share one campaign's eDNA value
# Using 5-fold CV (R=100) and Isotonic Monotonic Trajectories (isoreg) to ensure operational monotonicity
# Panels A & B: F1 and F2 Optimization Surfaces with Monotonic CV Threshold Trajectories
# Panels C & D: F1-optimised vs F2-optimised Confusion Matrices at CPUE >= 0.04
# Also exports plots/cv_performance_summary_table_site_1000m.csv aligned with 5-fold CV monotonic thresholds

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
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

out_suffix <- "site_1000m"
level_label <- "Site visit: ≤ 1000m, ≤ 6 months"

# 1. Load site-visit model data (shared input for the INLA, BRT and BRMS site models)
dat_evt <- readRDS("analysis/site/output/site_visit_model_data.rds") %>%
    filter(distance_m <= 1000) %>%
    transmute(
        Reef = as.character(Reef),
        reef = Reef,
        edna_campaign_id = as.character(edna_campaign_id),
        site_visit_id,
        perc_pos = edna_prop_positive * 100,
        counts = cots_count,
        total_bottom = bottom_time,
        cpue = counts / total_bottom
    )

cat(sprintf(
    "Site-visit dataset: N = %d site visits within 1000m (%d campaigns, %d reefs)\n",
    nrow(dat_evt), n_distinct(dat_evt$edna_campaign_id), n_distinct(dat_evt$Reef)
))

# ----------------------------------------------------
# 2D RASTER GRID: F1 & F2 Surfaces
# ----------------------------------------------------
perc_grid <- seq(0, 100, by = 2)
# Match the reef-level figure's y-axis: 99th percentile of reef-level CPUE
# (quantile(dat_evt$cpue, 0.99) in plots/generate_two_panel_multi_horizon.R)
cpue_axis_max <- 0.1636491
cpue_grid <- seq(0, cpue_axis_max, length.out = 40)

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

make_folds <- function(dat2) {
    # Several visits share one campaign's eDNA value, so hold out whole campaigns
    group_vfold_cv(dat2, group = edna_campaign_id, v = k, repeats = R)
}
set.seed(42)

run_5fold_tuning <- function(dat_evt, cpue_thr, perc_grid, k = 5, R = 100) {
    dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
    folds <- make_folds(dat2)

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
    folds <- make_folds(dat2)

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
write_csv(cv_table_summary, sprintf("plots/cv_performance_summary_table_%s.csv", out_suffix))

cat("\n=======================================================\n")
cat(" SITE-LEVEL (1000m) 5-FOLD CV MONOTONIC THRESHOLD TABLE\n")
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
        title = sprintf("A. F1-optimised threshold surface (%s)", level_label),
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

# Extract Monotonic CV thresholds for CPUE >= 0.04
p_star_f1 <- round((tuning_5fold %>% filter(cpue_thr == cpue_target))$perc_star_f1)
p_star_f2 <- round((tuning_5fold %>% filter(cpue_thr == cpue_target))$perc_star_f2)

cm_f1_obj <- make_confusion(dat_evt, perc_thresh = p_star_f1, cpue_thresh = cpue_target)
cm_f2_obj <- make_confusion(dat_evt, perc_thresh = p_star_f2, cpue_thresh = cpue_target)

panel_c <- ggplot(cm_f1_obj$cm, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.0, fontface = "bold", color = "black") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Prop.") +
    labs(
        title = sprintf("C. F1-Optimised Matrix (CPUE ≥ %.2f, %%pos* = %d%%, N = %d)", cpue_target, p_star_f1, nrow(dat_evt)),
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

panel_d <- ggplot(cm_f2_obj$cm, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.0, fontface = "bold", color = "black") +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Prop.") +
    labs(
        title = sprintf("D. F2-Optimised Matrix (CPUE ≥ %.2f, %%pos* = %d%%, N = %d)", cpue_target, p_star_f2, nrow(dat_evt)),
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

out_file <- function(ext) sprintf("plots/multi_horizon_CPUE_%s%s.%s", out_suffix, target_suffix, ext)
ggsave(out_file("png"), four_panel_plot, width = 11, height = 9.5, dpi = 300, bg = "white")
ggsave(out_file("pdf"), four_panel_plot, width = 11, height = 9.5, dpi = 300, device = cairo_pdf)
ggsave(out_file("eps"), four_panel_plot, width = 11, height = 9.5, dpi = 300)

cat(sprintf("Successfully generated %s and plots/cv_performance_summary_table_%s.csv!\n", out_file("png"), out_suffix))
