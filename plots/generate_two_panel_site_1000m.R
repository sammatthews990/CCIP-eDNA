# Script to generate 2-panel figure using SITE-LEVEL data matched within 1000m and 6 months (183 days)
# Panel A: F1 Heatmap across %pos threshold grid & CV Trajectory
# Panel B: Side-by-side Confusion Matrices (F1-Optimised vs Recall-Optimised)

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
library(sf)

# Load reefDNA package functions safely
pkg_path <- if (dir.exists("reefDNA")) "reefDNA" else if (dir.exists("../reefDNA")) "../reefDNA" else "."
devtools::load_all(pkg_path)

# Set clean global plot theme
theme_set(theme_bw(base_family = "Helvetica") + theme(panel.grid.minor = element_blank()))

# 1. Load Data
cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

# Clean Cull Data
cull <- cull.dat %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

# Aggregation of eDNA at Site Level (within Reef, Site_name, and Year)
edna_site <- edna.dat %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Year)) %>%
    rename(Collection.organisation = `Collection organisation`) %>%
    mutate(
        Reef = ReefName,
        Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
        date_edna = as.Date(Date),
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    group_by(Reef, Site_name, Collection.org, Year) %>%
    summarise(
        Lat = mean(Lat, na.rm = TRUE),
        Long = mean(Long, na.rm = TRUE),
        date_edna = min(date_edna),
        conc_mean = mean(Conc_mean, na.rm = TRUE),
        perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        .groups   = "drop"
    )

# 2. Spatial (1000m) & Temporal (6 months / 183 days) Matching
edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

# Buffer site eDNA locations by 1000m
edna_buf <- st_buffer(edna_sf, dist = 1000)
intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>%
    filter(!is.na(date_cull))

# Filter culls occurring within 0 to 183 days AFTER eDNA sampling (prior predictive window)
valid_encounters <- intersect_cull %>%
    mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
    filter(diff_days >= 0 & diff_days <= 183)

dat_site_1000m <- valid_encounters %>%
    st_drop_geometry() %>%
    rename(any_of(c(Reef = "Reef.x", Year = "Year.x"))) %>%
    group_by(Reef, Site_name, Collection.org, Year, perc_pos, conc_mean) %>%
    summarise(
        counts = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
        total_bottom = sum(Bottomtime, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(total_bottom > 0) %>%
    mutate(
        cpue = counts / total_bottom,
        reef = Reef
    )

cat(sprintf("Site-level dataset constructed: N = %d site-year encounters (1000m radius, <= 6 months prior)\n", nrow(dat_site_1000m)))

# ----------------------------------------------------
# PANEL A: F1 Heatmap across Threshold Grid & CV Trajectory (Site Level)
# ----------------------------------------------------
perc_grid <- seq(0, 100, by = 2)
cpue_grid <- seq(0, quantile(dat_site_1000m$cpue, 0.99, na.rm = TRUE), length.out = 40)

grid <- tidyr::expand_grid(perc_thresh = perc_grid, cpue_thresh = cpue_grid) %>%
    mutate(
        out = purrr::pmap(., ~ metrics_for(..1, ..2, dat_site_1000m)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy")
    ) %>%
    dplyr::select(-out)

# Cross-validate across multiple CPUE levels
cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
k <- 5
R <- 50
set.seed(42)

results <- map(cpue_levels, ~ run_cv_for_cpue(dat_site_1000m, .x, perc_grid, k = k, R = R))
perf_summary <- map_dfr(results, "perf_sum")

panel_a <- ggplot(grid, aes(x = perc_thresh, y = cpue_thresh, fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 1, by = 0.2), colour = "white", linewidth = 0.35, alpha = 0.8) +
    geom_path(data = perf_summary, aes(x = perc_star, y = cpue_thr), color = "red", linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
    geom_point(data = perf_summary, aes(x = perc_star, y = cpue_thr), color = "red", size = 2.8, inherit.aes = FALSE) +
    scale_fill_viridis_c(limits = c(0, 1), name = "F1 Score") +
    labs(
        title = "A. F1 optimised % positive thresholds (Site Level: 1000m & ≤ 6 Months)",
        x = "eDNA % positive threshold",
        y = "CPUE threshold"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 9),
        legend.position = "right"
    )

# ----------------------------------------------------
# PANEL B: Two Confusion Matrices (F1 vs Recall Optimised at CPUE >= 0.04)
# ----------------------------------------------------
cpue_target <- 0.04
perc_grid_1 <- seq(0, 100, by = 1)

grid_eval <- tibble(perc_thresh = perc_grid_1) %>%
    mutate(
        out = map(perc_thresh, ~ metrics_for(.x, cpue_target, dat_site_1000m)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
        Precision = vapply(out, `[[`, numeric(1), "Precision"),
        Recall = vapply(out, `[[`, numeric(1), "Recall")
    )

# F1-Optimised Threshold
best_f1 <- grid_eval %>% slice_max(F1, n = 1, with_ties = FALSE)
p_star_f1 <- best_f1$perc_thresh

# Recall-Optimised Threshold (Minimises False Negatives)
max_rec <- max(grid_eval$Recall)
best_rec <- grid_eval %>% filter(Recall == max_rec) %>% slice_max(perc_thresh, n = 1, with_ties = FALSE)
p_star_rec <- best_rec$perc_thresh

# Generate confusion matrices
cm_f1_obj  <- make_confusion(dat_site_1000m, perc_thresh = p_star_f1, cpue_thresh = cpue_target)
cm_rec_obj <- make_confusion(dat_site_1000m, perc_thresh = p_star_rec, cpue_thresh = cpue_target)

df_cm_f1 <- cm_f1_obj$cm %>%
    mutate(panel = sprintf("F1-Optimised Threshold (Site Level)\n(%%pos* = %d%% | F1 = %.3f)", p_star_f1, best_f1$F1))

df_cm_rec <- cm_rec_obj$cm %>%
    mutate(panel = sprintf("Recall-Optimised Threshold (Minimises FN)\n(%%pos* = %d%% | Recall = %.1f%% | FN = 0)", p_star_rec, best_rec$Recall * 100))

cm_combined <- bind_rows(df_cm_f1, df_cm_rec) %>%
    mutate(panel = factor(panel, levels = c(
        sprintf("F1-Optimised Threshold (Site Level)\n(%%pos* = %d%% | F1 = %.3f)", p_star_f1, best_f1$F1),
        sprintf("Recall-Optimised Threshold (Minimises FN)\n(%%pos* = %d%% | Recall = %.1f%% | FN = 0)", p_star_rec, best_rec$Recall * 100)
    )))

panel_b <- ggplot(cm_combined, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.2, fontface = "bold", color = "black") +
    facet_wrap(~panel, ncol = 2) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Proportion") +
    labs(
        title = sprintf("B. Outbreak Classification Performance (CPUE ≥ %.2f, Site Level N = %d): F1 vs Recall Optimised", cpue_target, nrow(dat_site_1000m)),
        x = "Prediction from eDNA (% pos)",
        y = "Ground Truth from COTS CPUE"
    ) +
    coord_fixed() +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey92", color = "grey70"),
        strip.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 9),
        legend.position = "right"
    )

# ----------------------------------------------------
# STACK VERTICALLY: Combined 2-Panel Graph
# ----------------------------------------------------
combined_plot <- panel_a / panel_b +
    plot_layout(heights = c(1, 1.05))

# Save to plots directory
dir.create("plots", showWarnings = FALSE)
ggsave("plots/multi_horizon_CPUE_site_1000m.png", combined_plot, width = 10, height = 9.5, dpi = 300)
ggsave("plots/multi_horizon_CPUE_site_1000m.pdf", combined_plot, width = 10, height = 9.5, dpi = 300)
ggsave("plots/multi_horizon_CPUE_site_1000m.eps", combined_plot, width = 10, height = 9.5, dpi = 300)

# Also update the primary multi_horizon_CPUE.png if desired
ggsave("plots/multi_horizon_CPUE.png", combined_plot, width = 10, height = 9.5, dpi = 300)

cat("Successfully generated site-level 2-panel figure (1000m & <=6 months)!\n")
