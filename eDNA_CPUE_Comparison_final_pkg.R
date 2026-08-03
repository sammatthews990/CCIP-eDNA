library(dplyr)
pkg_path <- if(dir.exists("reefDNA")) "reefDNA" else if(dir.exists("../reefDNA")) "../reefDNA" else "."
devtools::load_all(pkg_path)
library(ggplot2)
library(readxl)
library(fuzzyjoin)
library(broom)
library(MASS)
library(glmmTMB)
library(mgcv)
library(tidyr)
library(patchwork)
library(rsample)
library(purrr)
library(tibble)
library(viridis)
library(scales)

# Setup a clean plotting theme
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
# Formatting labels for linear models

# Confusion matrices and Evaluation Metrics


# Load Cull data
cull.dat <- read_excel("data/260201_COTS-Cull-Data-Ewels.xlsx", sheet = "Cull")

# Load eDNA data
edna.dat <- read_excel("data/eDNA data_ALL_20260225.xlsx", sheet = "eDNA_data_ALL")

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
        date_edna = as.Date(Date), # Excel supplies native dates
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    arrange(Reef, Year, date_edna) %>%
    group_by(Reef, Collection.org, Year) %>%
    # cluster rows into groups where consecutive samples are ≤ 7 days apart
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
win_days <- 183

same6_any <- fuzzy_inner_join(
    cull, edna_agg,
    by = c("Reef" = "Reef", "date_cull" = "date_edna"),
    match_fun = list(
        `==`,
        function(cull_date, edna_date) {
            abs(as.numeric(cull_date - edna_date)) <= win_days
        }
    )
) %>%
    mutate(diff_days = as.numeric(date_cull - date_edna), Reef = Reef.x)

reef_any <- same6_any %>%
    group_by(Reef, Collection.org, Year) %>%
    summarise(
        total_cots      = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
        total_bottom    = sum(Bottomtime, na.rm = TRUE),
        cpue_reef       = total_cots / total_bottom,
        perc_pos_reef   = mean(perc_pos, na.rm = TRUE),
        conc_mean_reef  = mean(conc_mean, na.rm = TRUE),
        n_matches       = n(),
        .groups         = "drop"
    )

reef_prior_3m <- get_prior_cohort(cull, edna_agg, 91, "3 Months")
reef_prior_6m <- get_prior_cohort(cull, edna_agg, 183, "6 Months")
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "12 Months")

reef_before <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(horizon = factor(horizon, levels = c("3 Months", "6 Months", "12 Months")))
p1 <- ggplot(reef_before, aes(x = perc_pos_reef, y = cpue_reef, color = horizon)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    labs(x = "% eDNA positive", y = "CPUE", title = "Prior eDNA (% Pos) vs CPUE") +
    facet_wrap(~horizon)

p2 <- ggplot(reef_before, aes(x = log1p(conc_mean_reef), y = cpue_reef, color = horizon)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    labs(x = "log1p(Mean Concentration)", y = "CPUE", title = "Prior eDNA (Conc) vs CPUE") +
    facet_wrap(~horizon)

p1 / p2
# Filter outliers systematically at modeling onset
dat_glmm <- reef_before %>%
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

effort_ref <- median(dat_glmm$total_bottom, na.rm = TRUE)

horizons <- levels(dat_glmm$horizon)
nd_all <- map_dfr(horizons, ~ fit_horizon_glmm(.x, dat_glmm)) %>%
    mutate(horizon = factor(horizon, levels = horizons))

ggplot(dat_glmm, aes(x = perc_pos_reef, y = obs_cpue)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_all, aes(x = perc_pos_reef, ymin = lcl_cpue, ymax = ucl_cpue, fill = horizon), alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = nd_all, aes(x = perc_pos_reef, y = fit_cpue, color = horizon), linewidth = 1, inherit.aes = FALSE) +
    facet_wrap(~horizon) +
    labs(x = "% eDNA positive", y = "CPUE", title = "NegBin GLM Offset Modeling: % Pos Degradation") +
    theme(legend.position = "none")
m_nb_base <- glmmTMB(
    total_cots ~ perc_pos_reef + offset(log(total_bottom)) + (1 | Reef),
    family = nbinom2(),
    data = dat_glmm
)

m_nb_both <- glmmTMB(
    total_cots ~ perc_pos_reef * conc_t + offset(log(total_bottom)) + (1 | Reef),
    family = nbinom2(),
    data = dat_glmm
)

# Interaction Anova Validation
anova(m_nb_base, m_nb_both)
g_perc <- gam(
    total_cots ~ s(perc_pos_reef, k = 3) + s(Reef, bs = "re") + offset(log(total_bottom)),
    family = nb(), method = "REML",
    data = dat_glmm
)

nd_gam <- tibble(
    perc_pos_reef = seq(min(dat_glmm$perc_pos_reef), max(dat_glmm$perc_pos_reef), length.out = 200),
    total_bottom = effort_ref,
    Reef = dat_glmm$Reef[1]
)

pp_perc <- predict(g_perc, newdata = nd_gam, type = "link", se.fit = TRUE, exclude = "s(Reef)")
nd_gam <- nd_gam %>%
    mutate(
        fit_cpue  = exp(pp_perc$fit) / total_bottom,
        lcl_cpue  = exp(pp_perc$fit - 1.96 * pp_perc$se.fit) / total_bottom,
        ucl_cpue  = exp(pp_perc$fit + 1.96 * pp_perc$se.fit) / total_bottom
    )

ggplot(dat_glmm, aes(x = perc_pos_reef, y = obs_cpue)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_gam, aes(x = perc_pos_reef, ymin = lcl_cpue, ymax = ucl_cpue), alpha = 0.3, fill = "skyblue", inherit.aes = FALSE) +
    geom_line(data = nd_gam, aes(x = perc_pos_reef, y = fit_cpue), linewidth = 1, color = "black", inherit.aes = FALSE) +
    labs(title = "GAM Spline NegBin Estimate: % Pos")
# Set CPUE Management Threshold
cpue_thresh <- 0.04

# Baseline uses 6 Months max window; Aggregating by Reef, Collection.org, and Year
dat_evt <- dat_glmm %>%
    filter(horizon == "6 Months") %>%
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
    ) %>%
    dplyr::select(reef, Collection.org, Year, perc_pos, conc_t, conc_mean_reef, cpue)

res_all <- make_confusion(dat_evt, perc_thresh = 48, cpue_thresh = cpue_thresh)
plot_confusion(res_all, title = sprintf("Full Data Matrix (eDNA ≥ 48%%, CPUE ≥ %.3f)", cpue_thresh))
# 1) Evaluate on a grid of %pos AND cpue thresholds (Heatmap)
perc_grid <- seq(0, 100, by = 2)
cpue_grid <- seq(0, quantile(dat_evt$cpue, 0.99, na.rm = TRUE), length.out = 40)

grid <- tidyr::expand_grid(perc_thresh = perc_grid, cpue_thresh = cpue_grid) %>%
    mutate(
        out = purrr::pmap(., ~ metrics_for(..1, ..2, dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy")
    ) %>%
    dplyr::select(-out)

# F1 Heatmap
p_hm <- ggplot(grid, aes(x = perc_thresh, y = cpue_thresh, fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 1, by = 0.2), colour = "black", linewidth = 0.3, alpha = 0.7) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = "eDNA % positive threshold", y = "CPUE threshold", fill = "F1", title = "F1 across %pos thresholds with iso-F1 contours")

# 2) Cross validate across multiple CPUE levels
cpue_levels <- c(0.01, 0.02, 0.04, 0.06, 0.08, 0.10)
k <- 5
R <- 100 # Reduced from 200 for rendering speed
set.seed(1)


results <- map(cpue_levels, ~ run_cv_for_cpue(dat_evt, .x, perc_grid, k = k, R = R))

summary_table <- map_dfr(results, "perf_sum") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(cpue_thr, perc_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

# Show the optimal threshold summary
summary_table

# Plot multi-panel confusion matrices
cm_all <- map_dfr(results, function(r) {
    r$cm_full$cm %>% mutate(cpue_thr = r$cpue_thr, perc_star = r$perc_star)
}) %>% mutate(panel = paste0("CPUE ≥ ", cpue_thr, "\n%pos* = ", perc_star))

p_cm <- ggplot(cm_all, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
    labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = "Confusion matrices with CV-tuned thresholds") +
    coord_fixed()

p_hm
p_cm
cpue_target <- 0.02

# Ensure we aggregate by Reef, Collection.org, and Year to preserve annual sampling events
dat_evt_hor <- dat_glmm %>%
    group_by(Reef, Collection.org, Year, horizon) %>%
    summarise(
        counts = sum(total_cots, na.rm = TRUE),
        total_bottom = sum(total_bottom, na.rm = TRUE),
        perc_pos = mean(perc_pos_reef, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        reef = Reef,
        cpue = counts / total_bottom
    ) %>%
    dplyr::select(reef, Collection.org, Year, perc_pos, cpue, horizon)


horizons <- levels(dat_evt_hor$horizon)
results_hor <- map(horizons, ~ run_cv_for_horizon(dat_evt_hor, .x, cpue_target, perc_grid, k = k, R = R))

summary_hor <- map_dfr(results_hor, "perf_sum") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(horizon, cpue_thr, perc_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

summary_hor

cm_all_hor <- map_dfr(results_hor, function(r) {
    r$cm_full$cm %>% mutate(cpue_thr = r$cpue_thr, perc_star = r$perc_star)
}) %>% mutate(panel = factor(paste0(horizon, "\n%pos* = ", perc_star), levels = paste0(levels(dat_evt_hor$horizon), "\n%pos* = ", map_dbl(results_hor, "perc_star"))))

p_cm_hor <- ggplot(cm_all_hor, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
    labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = sprintf("Confusion Matrices across Horizons (CPUE ≥ %.2f)", cpue_target)) +
    coord_fixed()

p_cm_hor

# Tune on transformed concentration to avoid skewness clustering 0s all the way across
conc_t_max <- quantile(dat_evt$conc_t, 0.99, na.rm = TRUE)
conc_grid_t <- seq(min(dat_evt$conc_t, na.rm = TRUE), conc_t_max, length.out = 30)

# Evaluate OR and AND grids over 2D thresholds for a baseline cpue
cpue_star <- 0.04

grid_both_and <- tidyr::expand_grid(perc_thresh = seq(0, 100, by = 2), conc_thresh_t = conc_grid_t) %>%
    mutate(
        out = purrr::pmap(., ~ metrics_for_both(..1, ..2, cpue_star, "and", dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1")
    )

p_2dhm <- ggplot(grid_both_and, aes(x = perc_thresh, y = expm1(conc_thresh_t), fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 0.9, by = 0.2), colour = "black", linewidth = 0.3) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(x = "% eDNA Positive Threshold", y = "Concentration Threshold (pseudo-log)", title = sprintf("2D Boolean Search: AND Rule F1 (CPUE ≥ %.3f)", cpue_star))


# Full 2D Loop over CPUE levels

results_2d <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "and")) # %>% discard(is.null)

summary_table_2d <- map_dfr(results_2d, "perf") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(cpue_thr, perc_star, conc_mean_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

# Show operational 2D thresholds
summary_table_2d

# Plot 2D Multi-CM

thr_df <- map_dfr(results_2d, "thr") %>% mutate(cpue_thr = sapply(results_2d, function(x) x$perf$cpue_thr))
cm_all_2d <- pmap_dfr(list(thr_df$cpue_thr, thr_df$perc_star, thr_df$conc_t_star), function(cp, p, c) {
    obj <- make_confusion_2d(dat_evt, p, c, cp)
    obj$cm %>% mutate(panel = paste0("CPUE ≥ ", cp, "\n%pos*=", round(p, 1), "  conc*≈", round(pmax(expm1(c), 0), 1)))
})

p_2dcm <- ggplot(cm_all_2d, aes(pred, actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Rate") +
    labs(x = "Prediction (2D rule)", y = "Ground truth (CPUE)", title = "2D AND Rule Confusion Matrices") +
    coord_fixed()

p_2dhm
p_2dcm
library(sf)

# 1. Standardize eDNA sample cohorts by their physical 'Site' level rather than 'Reef' layer
edna_site <- edna.dat %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Year)) %>%
    mutate(date_edna = as.Date(Date)) %>%
    group_by(ReefName, Site_name, Year) %>%
    summarise(
        Lat = mean(Lat),
        Long = mean(Long),
        date_edna = min(date_edna),
        conc_mean = mean(as.numeric(Conc_mean), na.rm = TRUE),
        perc_pos = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        .groups = "drop"
    ) %>%
    rename(Reef = ReefName)

# 2. Tag Sampling Regimes based on local density
# (avg_samples per site usually bounds around ~6 or ~12)
regime_df <- edna_site %>%
    group_by(Reef, Year) %>%
    summarise(
        avg_samples = mean(n_samples),
        .groups = "drop"
    ) %>%
    mutate(
        regime = case_when(
            avg_samples >= 10 ~ "3x12 Density",
            TRUE ~ "4x6 Density"
        )
    )

edna_site <- edna_site %>% left_join(regime_df, by = c("Reef", "Year"))

# 3. Project coordinate matrices into native EPSG:3112 (Lambert Conformal Conic for AUS)
edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)
buffer_radii <- c(200, 500, 1000, 2000)
spatial_results <- list()

for (dist_m in buffer_radii) {
    # Radial expansion around each eDNA point
    edna_buf <- st_buffer(edna_sf, dist = dist_m)

    # Isolate Cull coordinates bounding inside the circle
    intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>%
        filter(!is.na(date_cull))

    # Constrain to 6 month temporal overlaps
    valid_encounters <- intersect_cull %>%
        mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
        filter(diff_days >= 0 & diff_days <= 183)

    # Collapse cull metadata back onto the respective eDNA Site ID
    site_cpue <- valid_encounters %>%
        st_drop_geometry() %>%
        rename(any_of(c(Reef = "Reef.x", Year = "Year.x"))) %>%
        group_by(Reef, Site_name, Year, regime, perc_pos) %>%
        summarise(
            counts = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
            total_bottom = sum(Bottomtime, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            cpue = counts / total_bottom,
            radius = paste0(dist_m, "m Radius")
        )

    spatial_results[[as.character(dist_m)]] <- site_cpue
}

df_spatial <- bind_rows(spatial_results) %>%
    mutate(radius = factor(radius, levels = paste0(buffer_radii, "m Radius")))

# Run CV Threshold tuning evaluating spatial density bounds 
# We evaluate both 0.04 and 0.02 CPUE Targets, using unified global thresholds per radius.

res_04 <- evaluate_unified_spatial(df_spatial, 0.04, seq(0, 100, by = 5))
if (!is.null(res_04)) {
    p_spatial_04 <- res_04$p
    print(p_spatial_04)
}

res_02 <- evaluate_unified_spatial(df_spatial, 0.02, seq(0, 100, by = 5))
if (!is.null(res_02)) {
    p_spatial_02 <- res_02$p
    print(p_spatial_02)
}
metrics_list <- list()
if (exists("res_04") && !is.null(res_04)) metrics_list[["0.04"]] <- res_04$cm
if (exists("res_02") && !is.null(res_02)) metrics_list[["0.02"]] <- res_02$cm

if (length(metrics_list) > 0) {
    df_metrics <- bind_rows(metrics_list, .id = "CPUE") %>%
        mutate(class = case_when(
            pred == "Pred +" & actual == "Actual +" ~ "TP",
            pred == "Pred +" & actual == "Actual -" ~ "FP",
            pred == "Pred -" & actual == "Actual +" ~ "FN",
            pred == "Pred -" & actual == "Actual -" ~ "TN"
        )) %>%
        dplyr::select(CPUE, radius, regime, perc_star, class, n) %>%
        pivot_wider(names_from = class, values_from = n, values_fill = list(n = 0)) %>%
        mutate(
            TPR = TP / (TP + FN + 1e-9),
            Sensitivity = TPR,
            FPR = FP / (FP + TN + 1e-9),
            Precision = TP / (TP + FP + 1e-9),
            F1 = 2 * (Precision * TPR) / (Precision + TPR + 1e-9)
        ) %>%
        arrange(CPUE, radius, regime) %>%
        dplyr::select(
            CPUE, Radius = radius, Regime = regime, `% Pos Target` = perc_star, 
            TPR, FPR, Sensitivity, Precision, F1
        )
        
    # Expose the formatted statistics natively inside the rendered workflow
    df_metrics %>%
        mutate(across(c(TPR, FPR, Sensitivity, Precision, F1), ~ round(.x, 3))) %>%
        knitr::kable(caption = "Site-Level Performance Statistics bounded by Spatial Extents & Regimes")
}


# Safely load ggplot and export defined variables
library(ggplot2)
dir.create("plots", showWarnings=FALSE)
if(exists("p_hm")) ggsave("plots/temporal_degradation_GLM.png", plot = p_hm, width=10, height=8, dpi=300)
if(exists("p_cm")) ggsave("plots/multi_horizon_CPUE.png", plot = p_cm, width=10, height=8, dpi=300)
if(exists("p_cm_hor")) ggsave("plots/cv_horizons_02_CPUE.png", plot = p_cm_hor, width=10, height=4, dpi=300)
if(exists("p_2dhm")) ggsave("plots/2D_AND_F1_Heatmap.png", plot = p_2dhm, width=10, height=8, dpi=300)
if(exists("p_2dcm")) ggsave("plots/2D_AND_Multi_CM.png", plot = p_2dcm, width=10, height=8, dpi=300)
if(exists("p_spatial_04")) ggsave("plots/spatial_buffer_CM_04.png", plot = p_spatial_04, width=12, height=8, dpi=300)
if(exists("p_spatial_02")) ggsave("plots/spatial_buffer_CM_02.png", plot = p_spatial_02, width=12, height=8, dpi=300)
if(exists("df_metrics")) write.csv(df_metrics, "plots/spatial_metrics.csv", row.names=FALSE)
