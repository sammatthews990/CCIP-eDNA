library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(purrr)
library(tidyr)
devtools::load_all("reefDNA")

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

cull <- cull.dat %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

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

reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "6 Months", min_days = 0)

dat_glmm <- reef_prior_6m %>%
    mutate(
        conc_t    = log1p(conc_mean_reef),
        obs_cpue  = total_cots / total_bottom
    ) %>%
    filter(
        !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
        !is.na(perc_pos_reef), !is.na(conc_t), !is.na(Reef),
        conc_mean_reef < 5000
    ) %>%
    mutate(Reef = factor(Reef)) %>%
    droplevels()

dat_evt <- dat_glmm %>%
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
        cpue = counts / total_bottom,
        conc_t = log1p(conc_mean_reef)
    )

cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
conc_t_max <- quantile(dat_evt$conc_t, 0.99, na.rm = TRUE)
conc_grid_t <- seq(min(dat_evt$conc_t, na.rm = TRUE), conc_t_max, length.out = 30)
cpue_star <- 0.04

# 1. Grid search for heatmaps at CPUE = 0.04
perc_grid <- seq(0, 100, by = 2)
grid_both_and <- tidyr::expand_grid(perc_thresh = perc_grid, conc_thresh_t = conc_grid_t) %>%
    mutate(
        out = pmap(list(perc_thresh, conc_thresh_t), \(p, c) metrics_for_both(p, c, cpue_star, rule = "and", dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
        Precision = vapply(out, `[[`, numeric(1), "Precision"),
        Recall = vapply(out, `[[`, numeric(1), "Recall")
    )

grid_both_or <- tidyr::expand_grid(perc_thresh = perc_grid, conc_thresh_t = conc_grid_t) %>%
    mutate(
        out = pmap(list(perc_thresh, conc_thresh_t), \(p, c) metrics_for_both(p, c, cpue_star, rule = "or", dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
        Precision = vapply(out, `[[`, numeric(1), "Precision"),
        Recall = vapply(out, `[[`, numeric(1), "Recall")
    )

p_2dhm_and <- ggplot(grid_both_and, aes(x = perc_thresh, y = expm1(conc_thresh_t), fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 0.9, by = 0.2), colour = "black", linewidth = 0.3) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(x = "% eDNA Positive Threshold", y = "Concentration Threshold (copies/µL)", title = sprintf("2D Search: AND Rule F1 (CPUE ≥ %.3f)", cpue_star)) +
    theme_bw()

p_2dhm_or <- ggplot(grid_both_or, aes(x = perc_thresh, y = expm1(conc_thresh_t), fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 0.9, by = 0.2), colour = "black", linewidth = 0.3) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(x = "% eDNA Positive Threshold", y = "Concentration Threshold (copies/µL)", title = sprintf("2D Search: OR Rule F1 (CPUE ≥ %.3f)", cpue_star)) +
    theme_bw()

# 2. Run 2D CV for AND rule
results_2d_and <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "and"))
summary_table_2d_and <- map_dfr(compact(results_2d_and), "perf") %>%
    mutate(
        rule = "AND",
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    )

# 3. Run 2D CV for OR rule
results_2d_or <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "or"))
summary_table_2d_or <- map_dfr(compact(results_2d_or), "perf") %>%
    mutate(
        rule = "OR",
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    )

# 4. Plot 2D Multi-CM for AND
thr_df_and <- map_dfr(compact(results_2d_and), "thr") %>% mutate(cpue_thr = sapply(compact(results_2d_and), function(x) x$perf$cpue_thr))
cm_all_2d_and <- pmap_dfr(list(thr_df_and$cpue_thr, thr_df_and$perc_star, thr_df_and$conc_t_star), function(cp, p, c) {
    obj <- make_confusion_2d(dat_evt, p, c, cp, rule = "and")
    obj$cm %>% mutate(panel = paste0("CPUE ≥ ", cp, "\n%pos*=", round(p, 1), "  conc*≈", round(pmax(expm1(c), 0), 1)))
})

p_2dcm_and <- ggplot(cm_all_2d_and, aes(pred, actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.5, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Rate") +
    labs(x = "Prediction (AND Rule)", y = "Ground Truth (CPUE)", title = "2D AND Rule Confusion Matrices (CV-Tuned Thresholds)") +
    theme_bw() +
    coord_fixed()

# 5. Plot 2D Multi-CM for OR
thr_df_or <- map_dfr(compact(results_2d_or), "thr") %>% mutate(cpue_thr = sapply(compact(results_2d_or), function(x) x$perf$cpue_thr))
cm_all_2d_or <- pmap_dfr(list(thr_df_or$cpue_thr, thr_df_or$perc_star, thr_df_or$conc_t_star), function(cp, p, c) {
    obj <- make_confusion_2d(dat_evt, p, c, cp, rule = "or")
    obj$cm %>% mutate(panel = paste0("CPUE ≥ ", cp, "\n%pos*=", round(p, 1), "  conc*≈", round(pmax(expm1(c), 0), 1)))
})

p_2dcm_or <- ggplot(cm_all_2d_or, aes(pred, actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.5, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Rate") +
    labs(x = "Prediction (OR Rule)", y = "Ground Truth (CPUE)", title = "2D OR Rule Confusion Matrices (CV-Tuned Thresholds)") +
    theme_bw() +
    coord_fixed()

cat("Section test script completed successfully!\n")
