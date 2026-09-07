library(sf)
library(ggplot2)
library(dplyr)
library(readxl)
library(lubridate)
library(rnaturalearth)
library(ggspatial)
library(ggrepel)
library(cowplot)

sf_use_s2(FALSE)

font_family <- "Helvetica"

cat("1. Loading datasets...\n")
cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"
gpkg_path <- "data/rrap_canonical_2025-03-20-T15-18-17.gpkg"
shp_path <- "data/SectorShapefile"

gbr_reefs <- st_read(gpkg_path, quiet = TRUE) %>% st_transform(4326)
sector_sf <- st_read(shp_path, quiet = TRUE) %>% st_transform(4326)
aus_land <- ne_countries(country = "australia", scale = "medium", returnclass = "sf")

cull.dat <- read_excel(cull_file, sheet = "Cull")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

cat("2. Processing eDNA sampling locations and temporal alignment...\n")
edna_reef <- edna.dat %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Date)) %>%
    mutate(
        Reef = ReefName,
        date_edna = as.Date(Date),
        Year = year(date_edna)
    ) %>%
    group_by(Reef) %>%
    summarise(
        Lat = mean(Lat, na.rm = TRUE),
        Long = mean(Long, na.rm = TRUE),
        n_years = n_distinct(Year),
        n_samples = n(),
        .groups = "drop"
    )

edna_dates <- edna.dat %>%
    filter(!is.na(Date)) %>%
    transmute(Reef = ReefName, date_edna = as.Date(Date)) %>%
    distinct()

all_surveys <- bind_rows(
    cull.dat %>% filter(!is.na(SurveyDate), !is.na(ReefName)) %>% transmute(Reef = ReefName, date_survey = as.Date(SurveyDate))
) %>% distinct()

# Minimum absolute day difference for each eDNA date to nearest program survey on same reef
pairs <- edna_dates %>%
    inner_join(all_surveys, by = "Reef", relationship = "many-to-many") %>%
    mutate(diff_days = abs(as.numeric(date_survey - date_edna)))

abs_match <- pairs %>%
    group_by(Reef, date_edna) %>%
    summarise(min_abs_diff = min(diff_days), .groups = "drop") %>%
    group_by(Reef) %>%
    summarise(reef_min_diff = min(min_abs_diff), .groups = "drop")

edna_map_df <- edna_reef %>%
    left_join(abs_match, by = "Reef") %>%
    mutate(
        window = case_when(
            is.na(reef_min_diff) ~ "> 12 Months",
            reef_min_diff <= 91 ~ "Within 3 Months",
            reef_min_diff <= 183 ~ "Within 6 Months",
            reef_min_diff <= 365 ~ "Within 12 Months",
            TRUE ~ "> 12 Months"
        ),
        window = factor(window, levels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months"))
    )

edna_sf <- st_as_sf(edna_map_df, coords = c("Long", "Lat"), crs = 4326)

# Coastal towns & key locations
towns <- tibble::tribble(
    ~name, ~lat, ~long, ~nudge_x, ~nudge_y,
    "Thursday Island", -10.58, 142.22, -0.70, 0.05,
    "Mer Island", -9.91, 144.05, 0.50, 0.05,
    "Lockhart River", -12.78, 143.34, -0.70, 0.00,
    "Cooktown", -15.47, 145.25, -0.70, 0.05,
    "Cairns", -16.92, 145.77, -0.70, 0.00,
    "Townsville", -19.26, 146.81, -0.75, -0.10,
    "Mackay", -21.14, 149.19, -0.70, -0.10,
    "Gladstone", -23.84, 151.26, -0.75, -0.10
)
towns_sf <- st_as_sf(towns, coords = c("long", "lat"), crs = 4326)

# 3. Rotate spatial layers by 30 degrees (counter-clockwise)
cat("3. Rotating spatial layers (30 degrees counter-clockwise)...\n")
center_pt <- c(147.5, -17.5)
rot_angle_deg <- 30

rotate_sf <- function(sf_obj, angle_deg, center = center_pt) {
    angle_rad <- angle_deg * pi / 180
    rot <- matrix(c(cos(angle_rad), sin(angle_rad), -sin(angle_rad), cos(angle_rad)), 2, 2)
    geom <- st_geometry(sf_obj)
    crs_orig <- st_crs(sf_obj)
    geom_rot <- (geom - center) * rot + center
    st_geometry(sf_obj) <- geom_rot
    st_crs(sf_obj) <- crs_orig
    return(sf_obj)
}

gbr_reefs_rot <- rotate_sf(gbr_reefs, rot_angle_deg)
sector_sf_rot <- rotate_sf(sector_sf, rot_angle_deg)
aus_land_rot <- rotate_sf(aus_land, rot_angle_deg)
edna_sf_rot <- rotate_sf(edna_sf, rot_angle_deg)
towns_sf_rot <- rotate_sf(towns_sf, rot_angle_deg)

towns_rot_coords <- st_coordinates(towns_sf_rot)
towns_df_rot <- towns %>%
    mutate(
        rot_long = towns_rot_coords[, 1],
        rot_lat  = towns_rot_coords[, 2]
    )

bbox_rot <- st_bbox(edna_sf_rot)

# Okabe-Ito Color-Blind Safe Palette
window_colors <- c(
    "Within 3 Months"  = "#009E73", # Okabe-Ito Bluish Green
    "Within 6 Months"  = "#E69F00", # Okabe-Ito Amber Orange
    "Within 12 Months" = "#D55E00", # Okabe-Ito Vermilion Red
    "> 12 Months"      = "#555555" # Slate Dark Grey
)

cat("4. Building Main GBR Map...\n")
p_main <- ggplot() +
    # Sector shapefile layer (white background with black sector borders)
    geom_sf(data = sector_sf_rot, fill = "white", color = "black", linewidth = 0.55) +
    # GBR Reef outlines background
    geom_sf(data = gbr_reefs_rot, fill = "#d8d8d8", color = "#c4c4c4", linewidth = 0.08, alpha = 0.75) +
    # Australian land mass
    # geom_sf(data = aus_land_rot, fill = "#f0f0ec", color = "#a0a0a0", linewidth = 0.35) +
    # eDNA Reef Points with outline
    geom_sf(
        data = edna_sf_rot,
        aes(fill = window, size = n_years),
        shape = 21, color = "black", stroke = 0.35, alpha = 0.90
    ) +
    # Coastal Towns & Key Locations
    geom_sf(data = towns_sf_rot, size = 2.2, color = "black", shape = 16) +
    geom_text_repel(
        data = towns_df_rot,
        aes(x = rot_long, y = rot_lat, label = name),
        size = 3.6, fontface = "plain", family = font_family, color = "grey10",
        nudge_x = towns$nudge_x, nudge_y = towns$nudge_y,
        segment.color = "grey40"
    ) +
    # Rotated North Arrow
    annotation_north_arrow(
        location = "br", which_north = "true",
        pad_x = unit(0.5, "in"), pad_y = unit(0.75, "in"),
        rotation = -rot_angle_deg,
        style = north_arrow_fancy_orienteering(text_size = 10)
    ) +
    # Scale bar
    annotation_scale(
        location = "br", width_hint = 0.22,
        pad_x = unit(0.5, "in"), pad_y = unit(0.4, "in"),
        text_cex = 1
    ) +
    # Color & Size scales (Okabe-Ito)
    scale_fill_manual(
        name = "eDNA & Survey Window",
        values = window_colors,
        labels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months / Unmatched")
    ) +
    scale_size_continuous(
        name = "Sampling Effort",
        range = c(2.2, 6.8),
        breaks = c(1, 2, 4, 7),
        labels = c("1 Year", "2 Years", "4 Years", "7 Years")
    ) +
    guides(
        fill = guide_legend(order = 1, override.aes = list(size = 4, shape = 21, stroke = 0.4)),
        size = guide_legend(order = 2, override.aes = list(fill = "grey40", shape = 21, stroke = 0.4))
    ) +
    coord_sf(
        xlim = c(bbox_rot$xmin - 2, bbox_rot$xmax + 3.5),
        ylim = c(bbox_rot$ymin - 1.5, bbox_rot$ymax + 0.5),
        expand = FALSE
    ) +
    theme_bw(base_family = font_family) +
    theme(
        text = element_text(family = font_family),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.78, 0.42),
        legend.background = element_rect(fill = alpha("white", 0.93), color = "grey50", linewidth = 0.5),
        legend.key = element_blank(),
        legend.title = element_text(family = font_family, face = "bold", size = 9),
        legend.text = element_text(family = font_family, size = 8.5),
        legend.spacing.y = unit(2, "pt"),
        legend.margin = margin(6, 8, 6, 8)
    )

p_main

cat("5. Building Inset Map with Rotated & Narrow Rectangle (-30 degrees)...\n")
rect_center <- c(147.5, -17.5)
rect_half_w <- 1.8 # narrowed width (3.6 deg total width)
rect_half_h <- 8.0 # length along GBR

box_pts <- matrix(c(
    -rect_half_w, -rect_half_h,
    rect_half_w, -rect_half_h,
    rect_half_w, rect_half_h,
    -rect_half_w, rect_half_h,
    -rect_half_w, -rect_half_h
), ncol = 2, byrow = TRUE)

box_angle_rad <- (-rot_angle_deg) * pi / 180
box_rot_matrix <- matrix(c(cos(box_angle_rad), sin(box_angle_rad), -sin(box_angle_rad), cos(box_angle_rad)), 2, 2)
rot_box_pts <- t(t(box_pts %*% box_rot_matrix) + rect_center)

inset_rect_sf <- st_sf(
    geometry = st_sfc(st_polygon(list(rot_box_pts)), crs = 4326)
)

p_inset <- ggplot() +
    geom_sf(data = aus_land, fill = "#e5e5e0", color = "#888888", linewidth = 0.2) +
    geom_sf(data = inset_rect_sf, fill = NA, color = "#D55E00", linewidth = 0.8) +
    coord_sf(xlim = c(110, 155), ylim = c(-44, -9), expand = FALSE) +
    theme_void(base_family = font_family) +
    theme(
        text = element_text(family = font_family),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        plot.margin = margin(0, 0, 0, 0)
    )

cat("6. Combining and rendering final figure...\n")
final_map <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset, x = 0.49, y = 0.65, width = 0.29, height = 0.26)

final_map

dir.create("plots", showWarnings = FALSE)

cat("Saving PNG...\n")
ggsave("plots/gbr_edna_sample_map.png", final_map, width = 8.5, height = 11, dpi = 300)

cat("Saving PDF...\n")
ggsave("plots/gbr_edna_sample_map.pdf", final_map, width = 8.5, height = 11, dpi = 300)

cat("Saving EPS...\n")
ggsave("plots/gbr_edna_sample_map.eps", final_map, width = 8.5, height = 11, dpi = 300)

cat("Saving TIFF...\n")
ggsave("plots/gbr_edna_sample_map.tiff", final_map, width = 8.5, height = 11, dpi = 300, compression = "lzw")

cat("All map files saved successfully!\n")
