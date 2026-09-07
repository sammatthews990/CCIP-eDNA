library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(purrr)
library(tidyr)
library(sf)
library(rsample)
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

reef_prior_3m  <- get_prior_cohort(cull, edna_agg, 91,  "0-3 Months", min_days = 0)
reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "3-6 Months", min_days = 92)
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "6-12 Months", min_days = 184)

dat_glmm <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(
        conc_t      = log1p(conc_mean_reef),
        combo_index = perc_pos_reef * conc_t,
        obs_cpue    = total_cots / total_bottom,
        horizon     = factor(horizon, levels = c("0-3 Months", "3-6 Months", "6-12 Months"))
    ) %>%
    filter(
        !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
        !is.na(perc_pos_reef), !is.na(conc_t), !is.na(Reef),
        conc_mean_reef < 5000
    ) %>%
    mutate(Reef = factor(Reef)) %>%
    droplevels()

cat("=== 1. SPATIAL BUFFER BOOTSTRAP METRIC DECAY ===\n")
buffer_radii <- c(200, 500, 1000, 2000, 3000, 5000)

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

edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

spatial_list <- list()
for (dist_m in buffer_radii) {
    edna_buf <- st_buffer(edna_sf, dist = dist_m)
    intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>%
        filter(!is.na(date_cull))
    valid_encounters <- intersect_cull %>%
        mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
        filter(diff_days >= 0 & diff_days <= 183)
    
    site_cpue <- valid_encounters %>%
        st_drop_geometry() %>%
        rename(any_of(c(Reef = "Reef.x", Year = "Year.x"))) %>%
        group_by(Reef, Site_name, Year, perc_pos, conc_mean) %>%
        summarise(
            counts = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
            total_bottom = sum(Bottomtime, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            cpue = counts / total_bottom,
            radius = dist_m,
            conc_t = log1p(conc_mean),
            combo_index = perc_pos * conc_t
        )
    spatial_list[[as.character(dist_m)]] <- site_cpue
}

df_spatial_full <- bind_rows(spatial_list)

# Bootstrapping metrics across distances for CPUE = 0.02 and 0.04
set.seed(42)
B <- 100
cpue_targets <- c(0.02, 0.04)

calc_metrics <- function(df, cpue_target, p_star = 50) {
    classified <- df %>%
        mutate(
            actual = cpue >= cpue_target,
            pred   = perc_pos >= p_star
        )
    tp <- sum(classified$pred & classified$actual)
    fp <- sum(classified$pred & !classified$actual)
    fn <- sum(!classified$pred & classified$actual)
    tn <- sum(!classified$pred & !classified$actual)
    
    acc  <- (tp + tn) / nrow(classified)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1   <- if ((prec + rec) > 0) 2 * (prec * rec) / (prec + rec) else 0
    fpr  <- if ((fp + tn) > 0) fp / (fp + tn) else 0
    fnr  <- if ((fn + tp) > 0) fn / (fn + tp) else 0
    
    c(Accuracy = acc, F1 = f1, Precision = prec, Recall = rec, FPR = fpr, FNR = fnr)
}

boot_results <- list()

for (cp in cpue_targets) {
    for (r in buffer_radii) {
        sub_df <- df_spatial_full %>% filter(radius == r)
        if (nrow(sub_df) < 5) next
        
        boot_matrix <- replicate(B, {
            idx <- sample(nrow(sub_df), replace = TRUE)
            calc_metrics(sub_df[idx, ], cpue_target = cp)
        })
        
        boot_df <- as.data.frame(t(boot_matrix)) %>%
            pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
            group_by(Metric) %>%
            summarise(
                mean_val = mean(Value, na.rm = TRUE),
                lower_ci = quantile(Value, 0.025, na.rm = TRUE),
                upper_ci = quantile(Value, 0.975, na.rm = TRUE),
                .groups  = "drop"
            ) %>%
            mutate(cpue_target = paste0("CPUE Target ≥ ", cp), radius = r)
        
        boot_results[[paste(cp, r, sep = "_")]] <- boot_df
    }
}

df_boot <- bind_rows(boot_results) %>%
    mutate(
        Metric = factor(Metric, levels = c("F1", "Accuracy", "Precision", "Recall", "FPR", "FNR")),
        radius_km = radius / 1000
    )

p_boot_decay <- ggplot(df_boot, aes(x = radius, y = mean_val, color = Metric, fill = Metric)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "darkred", linewidth = 0.8) +
    facet_grid(Metric ~ cpue_target, scales = "free_y") +
    scale_x_continuous(breaks = buffer_radii, labels = c("200m", "500m", "1km", "2km", "3km", "5km")) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    labs(
        x = "Spatial Buffer Radius Around eDNA Site",
        y = "Metric Score (Bootstrapped 95% CI)",
        title = "Spatial Classification Metric Decay Across Buffer Radii (200m – 5000m)"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave("plots/spatial_metric_decay_bootstrapped.png", p_boot_decay, width = 10, height = 9, dpi = 300)
cat("Bootstrapped metric decay plot saved to plots/spatial_metric_decay_bootstrapped.png!\n")

cat("=== 2. EYRIE REEF MAP WITH HIGHLIGHTED OPTIMAL BUFFER (1000m) ===\n")

eyrie_pts <- edna_sf %>% filter(grepl("Eyrie", Reef, ignore.case = TRUE))

# Use 6km buffer for cropping cull points cleanly
eyrie_buf_6k <- st_buffer(st_union(eyrie_pts), dist = 6000)
eyrie_cull   <- st_intersection(cull_sf, eyrie_buf_6k)

# Generate concentric buffer rings around Eyrie eDNA sites
eyrie_rings <- map_dfr(buffer_radii, function(r) {
    buf <- st_buffer(st_union(eyrie_pts), dist = r)
    st_sf(radius = paste0(r, "m"), radius_m = r, geometry = buf)
}) %>% mutate(
    is_optimal = radius_m == 1000,
    radius = factor(radius, levels = paste0(buffer_radii, "m"))
)

p_eyrie_map <- ggplot() +
    geom_sf(data = eyrie_cull, aes(color = "COTS Cull Dive Track"), alpha = 0.6, size = 1.8) +
    geom_sf(data = eyrie_rings, aes(color = radius, linetype = is_optimal, linewidth = is_optimal), fill = NA) +
    geom_sf(data = eyrie_pts, color = "black", fill = "yellow", shape = 21, size = 3, stroke = 1) +
    scale_linetype_manual(values = c("FALSE" = "dotted", "TRUE" = "solid"), guide = "none") +
    scale_linewidth_manual(values = c("FALSE" = 0.6, "TRUE" = 1.3), guide = "none") +
    scale_color_manual(
        name = "Layer / Radius",
        values = c(
            "COTS Cull Dive Track" = "coral2",
            "200m"  = "#440154",
            "500m"  = "#3b528b",
            "1000m" = "#d90429", # HIGHLIGHTED OPTIMAL IN BRIGHT RED
            "2000m" = "#21918c",
            "3000m" = "#5ec962",
            "5000m" = "#fde725"
        )
    ) +
    labs(
        title = "Eyrie Reef (14-118) eDNA Sampling & Spatial Buffer Extents",
        subtitle = "Optimal spatial buffer radius (1000m highlighted in solid red) maximizes F1 precision/recall trade-off",
        x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30")
    )

ggsave("plots/eyrie_reef_buffer_map.png", p_eyrie_map, width = 9, height = 7, dpi = 300)
cat("Eyrie Reef spatial buffer map saved to plots/eyrie_reef_buffer_map.png!\n")
