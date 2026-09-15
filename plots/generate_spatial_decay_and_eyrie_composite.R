# Script to generate 2-panel composite figure:
# Panel A: 4-quadrant Quantified Spatial Signal Degradation (Max F1, AUC, Beta, Precision across 200m - 5000m)
# Panel B: Eyrie Reef (14-118) High-Res Satellite Map (with tight legend & zero dead whitespace)

library(dplyr)
library(ggplot2)
library(readxl)
library(purrr)
library(tidyr)
library(sf)
library(patchwork)
library(maptiles)
library(tidyterra)
library(glmmTMB)
library(pROC)

devtools::load_all("reefDNA")

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"
cull_gpkg <- "data/Eotr_CotsCullSites_2025_11_19_1_58_PM.gpkg"

cull <- read_excel(cull_file, sheet = "Cull") %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

edna_site <- read_excel(edna_file, sheet = "eDNA_data_ALL") %>%
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

edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

# ==========================================
# PANEL A: Quantified Degradation Across 200m - 5000m
# ==========================================
buffer_radii <- c(200, 500, 1000, 1500, 2000, 3000, 5000)
cpue_target <- 0.04
set.seed(42)

degradation_results <- map_dfr(buffer_radii, function(r) {
    edna_buf <- st_buffer(edna_sf, dist = r)
    intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>% filter(!is.na(date_cull))
    
    valid_encounters <- intersect_cull %>%
        mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
        filter(diff_days >= 0 & diff_days <= 183)
    
    site_cpue <- valid_encounters %>%
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
            actual = cpue >= cpue_target
        )
    
    # Tune optimal F1 for this radius
    p_grid <- seq(0, 100, by = 1)
    grid_eval <- map_dfr(p_grid, function(p) {
        m <- metrics_for(p, cpue_target, site_cpue)
        tibble(
            p_star = p, 
            F1 = unname(m["F1"]), 
            Precision = unname(m["Precision"]), 
            Recall = unname(m["Recall"]), 
            Accuracy = unname(m["Accuracy"]), 
            FPR = unname(m["FPR"])
        )
    })
    
    best_row <- grid_eval %>% slice_max(F1, n = 1, with_ties = FALSE)
    
    # Calculate AUC-ROC
    roc_obj <- tryCatch(pROC::roc(site_cpue$actual, site_cpue$perc_pos, quiet = TRUE), error = function(e) NULL)
    auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA
    
    # GLMM fit effect size
    glmm_fit <- tryCatch(
        glmmTMB(counts ~ perc_pos + offset(log(total_bottom)) + (1|Reef), data = site_cpue, family = nbinom2),
        error = function(e) NULL
    )
    beta_perc <- if (!is.null(glmm_fit)) fixef(glmm_fit)$cond["perc_pos"] else NA
    
    tibble(
        radius_m = r,
        Max_F1 = best_row$F1,
        Precision = best_row$Precision,
        Recall = best_row$Recall,
        Accuracy = best_row$Accuracy,
        FPR = best_row$FPR,
        AUC = auc_val,
        Beta_perc_pos = beta_perc
    )
})

# Construct 4 degradation subplots
p1 <- ggplot(degradation_results, aes(x = radius_m, y = Max_F1)) +
    geom_line(color = "#2b5c8f", linewidth = 1.1) +
    geom_point(color = "#2b5c8f", size = 2.5) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = 1100, y = min(degradation_results$Max_F1, na.rm=TRUE) + 0.02, label = "Optimal Threshold (1000m)", color = "red", hjust = 0, fontface = "bold", size = 3) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "A1. Max F1 Score vs Buffer Radius", x = NULL, y = "Max F1 Score") +
    theme_bw(base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 10.5), plot.margin = margin(3, 3, 3, 3),
        axis.text.x = element_blank())

p2 <- ggplot(degradation_results, aes(x = radius_m, y = AUC)) +
    geom_line(color = "#2a9d8f", linewidth = 1.1) +
    geom_point(color = "#2a9d8f", size = 2.5) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "A2. Area Under ROC Curve (AUC)", x = NULL, y = "AUC-ROC") +
    theme_bw(base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 10.5), plot.margin = margin(3, 3, 3, 3),
        axis.text.x = element_blank())

p3 <- ggplot(degradation_results, aes(x = radius_m, y = Beta_perc_pos)) +
    geom_line(color = "#e76f51", linewidth = 1.1) +
    geom_point(color = "#e76f51", size = 2.5) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "A3. GLMM Effect Size (Beta)", x = "Spatial Buffer Radius", y = "Beta (%pos coefficient)") +
    theme_bw(base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 10.5), plot.margin = margin(3, 3, 3, 3),
        axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(degradation_results, aes(x = radius_m, y = Precision)) +
    geom_line(color = "#e63946", linewidth = 1.1) +
    geom_point(color = "#e63946", size = 2.5) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "A4. Classification Precision", x = "Spatial Buffer Radius", y = "Precision") +
    theme_bw(base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 10.5), plot.margin = margin(3, 3, 3, 3),
        axis.text.x = element_text(angle = 45, hjust = 1))

panel_a_quantified <- (p1 + p2) / (p3 + p4) + 
    plot_annotation(
        title = "A. Quantified eDNA Spatial Signal Degradation Across Buffer Radii (200m – 5000m)",
        theme = theme(plot.title = element_text(face = "bold", size = 12, family = "Helvetica"))
    )

# ==========================================
# PANEL B: Eyrie Reef High-Res Satellite Map (Inset Legend)
# ==========================================
eyrie_site_agg_projected <- edna_sf %>% filter(grepl("Eyrie", Reef, ignore.case = TRUE))
eyrie_site_agg_projected <- st_transform(eyrie_site_agg_projected, 3112)

cull_sites_gpkg <- st_read(cull_gpkg, quiet = TRUE) %>% st_transform(3112)
eyrie_buffer_radii <- c(200, 500, 1000, 2000)
eyrie_buffer_labels <- c("0–200m", "200–500m", "500–1000m", "1000–2000m")

eyrie_cumulative_buffers <- lapply(eyrie_buffer_radii, function(r) {
    st_union(st_buffer(eyrie_site_agg_projected, dist = r))
})

eyrie_buffer_bands <- Map(
    function(current, previous) {
        if (is.null(previous)) current else st_difference(current, previous)
    },
    eyrie_cumulative_buffers,
    c(list(NULL), eyrie_cumulative_buffers[-length(eyrie_cumulative_buffers)])
)

eyrie_buffer_bands_projected <- st_sf(
    buffer_band = factor(eyrie_buffer_labels, levels = eyrie_buffer_labels),
    geometry = do.call(c, eyrie_buffer_bands)
)

eyrie_site_agg <- st_transform(eyrie_site_agg_projected, 4326)
eyrie_buffer_bands <- st_transform(eyrie_buffer_bands_projected, 4326)

# eyrie_map_extent_3112 <- st_as_sfc(st_bbox(st_buffer(eyrie_cumulative_buffers[[4]], 250)))
base_bbox <- st_bbox(eyrie_cumulative_buffers[[4]])
# Add a massive horizontal buffer (e.g., 2000m left/right) and smaller vertical buffer (e.g., 250m top/bottom)
eyrie_map_extent_3112 <- st_as_sfc(st_bbox(c(
    base_bbox["xmin"] - 0,  # Expand left
    base_bbox["xmax"] + 1800,  # Expand right
    base_bbox["ymin"] - 250,   # Keep top/bottom relatively tight
    base_bbox["ymax"] + 250
), crs = 3112))
eyrie_map_extent <- st_transform(eyrie_map_extent_3112, 4326)
eyrie_map_bbox <- st_bbox(eyrie_map_extent)

cull_eyrie_boxes_3112 <- suppressWarnings(st_intersection(cull_sites_gpkg, eyrie_map_extent_3112))
cull_eyrie_boxes <- st_transform(cull_eyrie_boxes_3112, 4326)

eyrie_satellite <- get_tiles(
    eyrie_map_extent,
    provider = "Esri.WorldImagery",
    crop = TRUE,
    zoom = 15
)

panel_b <- ggplot() +
    geom_spatraster_rgb(data = eyrie_satellite) +
    geom_sf(
        data = eyrie_buffer_bands,
        aes(fill = buffer_band),
        color = "white",
        linewidth = 0.3,
        alpha = 0.20
    ) +
    geom_sf(data = cull_eyrie_boxes, fill = alpha("orange", 0.35), color = "darkorange", linewidth = 0.4, show.legend = FALSE) +
    geom_sf(data = eyrie_site_agg, color = "black", size = 4.5, shape = 16, show.legend = FALSE) +
    geom_sf(data = eyrie_site_agg, aes(color = perc_pos), size = 3.2, shape = 16) +
    scale_fill_manual(
        name = "Dissolved buffer band",
        values = c(
            "0–200m" = "#ffffcc",
            "200–500m" = "#c2e699",
            "500–1000m" = "#78c679",
            "1000–2000m" = "#238443"
        )
    ) +
    scale_color_gradient2(
        name = "eDNA samples\n% positive",
        low = "#2166ac",
        mid = "#f7f7f7",
        high = "#b2182b",
        midpoint = 50,
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100)
    ) +
    coord_sf(
        crs = st_crs(4326),
        default_crs = st_crs(4326),
        xlim = unname(c(eyrie_map_bbox["xmin"], eyrie_map_bbox["xmax"])),
        ylim = unname(c(eyrie_map_bbox["ymin"], eyrie_map_bbox["ymax"])),
        expand = FALSE
    ) +
    labs(
        title = "B. Eyrie Reef (14-118): eDNA Signal & Spatial Buffer Extents (200m – 2000m)",
        subtitle = "Site % positive (blue–red), COTS cull polygons (orange), and dissolved buffer bands",
        x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        legend.position = c(0.85, 0.25),
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "grey60", linewidth = 0.4),
        legend.box.margin = margin(4, 4, 4, 4),
        legend.title = element_text(face = "bold", size = 8.5),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey20"),
        plot.margin = margin(3, 3, 3, 3)
    )

# Save individual Eyrie satellite map with inset legend
ggsave("plots/eyrie_reef_satellite_2000m.png", panel_b, width = 9, height = 7.5, dpi = 300)

# Combine Panel A and Panel B vertically
composite_2panel <- wrap_elements(panel_a_quantified) / panel_b + plot_layout(heights = c(1.1, 1.15))

ggsave("plots/spatial_degradation_and_eyrie_composite.png", composite_2panel, width = 8.5, height = 13.5, dpi = 300)
ggsave("plots/spatial_degradation_and_eyrie_composite.pdf", composite_2panel, width = 9, height = 13.5)

cat("Successfully generated updated composite figure with inset legend and minimal white space!\n")
