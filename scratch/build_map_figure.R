library(sf)
library(ggplot2)
library(dplyr)
library(readxl)
library(lubridate)
library(rnaturalearth)
library(ggspatial)
library(ggrepel)
library(cowplot)

# 1. Load Data
cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"
gpkg_path <- "data/Eotr_CotsCullSites_2025_11_19_1_58_PM.gpkg"

cat("Loading reef boundaries...\n")
gbr_reefs <- st_read(gpkg_path, quiet = TRUE)

cat("Loading Australia land...\n")
aus_land <- ne_countries(country = "australia", scale = "medium", returnclass = "sf")

cat("Loading survey & eDNA datasets...\n")
cull.dat <- read_excel(cull_file, sheet = "Cull")
manta.dat <- read_excel(cull_file, sheet = "Manta")
rhis.dat <- read_excel(cull_file, sheet = "RHIS")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

# 2. Process eDNA Reef Locations & Year counts
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

# Combine all survey dates (Cull, Manta, RHIS)
all_surveys <- bind_rows(
    cull.dat %>% filter(!is.na(SurveyDate), !is.na(ReefName)) %>% transmute(Reef = ReefName, date_survey = as.Date(SurveyDate)),
    manta.dat %>% filter(!is.na(SurveyTime), !is.na(ReefName)) %>% transmute(Reef = ReefName, date_survey = as.Date(SurveyTime)),
    rhis.dat %>% filter(!is.na(SurveyTime), !is.na(ReefName)) %>% transmute(Reef = ReefName, date_survey = as.Date(SurveyTime))
) %>% distinct()

# Calculate min absolute day difference for each eDNA date to nearest survey date on same reef
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
            reef_min_diff <= 91 ~ "3 Months",
            reef_min_diff <= 183 ~ "6 Months",
            reef_min_diff <= 365 ~ "12 Months",
            TRUE ~ "> 12 Months"
        ),
        window = factor(window, levels = c("3 Months", "6 Months", "12 Months", "> 12 Months")),
        # Binned years for discrete sizing in legend
        years_cat = case_when(
            n_years == 1 ~ "1 Year",
            n_years == 2 ~ "2 Years",
            n_years == 3 ~ "3 Years",
            n_years >= 4 ~ "4+ Years"
        ),
        years_cat = factor(years_cat, levels = c("1 Year", "2 Years", "3 Years", "4+ Years"))
    )

edna_sf <- st_as_sf(edna_map_df, coords = c("Long", "Lat"), crs = 4326)

# Coastal towns
towns <- tibble::tribble(
    ~name, ~lat, ~long, ~nudge_x, ~nudge_y,
    "Cooktown", -15.47, 145.25, -0.6, 0.0,
    "Cairns", -16.92, 145.77, -0.6, 0.0,
    "Townsville", -19.26, 146.81, -0.7, -0.1,
    "Mackay", -21.14, 149.19, -0.6, -0.1,
    "Gladstone", -23.84, 151.26, -0.7, -0.1
)
towns_sf <- st_as_sf(towns, coords = c("long", "lat"), crs = 4326)

# Main GBR Map
cat("Building main plot...\n")

# Color palette: Vibrant, publication-quality
window_colors <- c(
    "3 Months"   = "#1b9e77", # Dark green/teal
    "6 Months"   = "#d95f02", # Bright orange
    "12 Months"  = "#7570b3", # Purple/indigo
    "> 12 Months"= "#e7298a"  # Magenta/pink or dark gray
)

p_main <- ggplot() +
    # GBR Reef boundaries background
    geom_sf(data = gbr_reefs, fill = "#e0e0e0", color = "#d0d0d0", size = 0.05, alpha = 0.6) +
    # Australian land mass
    geom_sf(data = aus_land, fill = "#f5f5f3", color = "#b0b0b0", size = 0.3) +
    # eDNA Reef Points
    geom_sf(
        data = edna_sf,
        aes(color = window, size = n_years),
        alpha = 0.85
    ) +
    # Coastal Towns points & text
    geom_sf(data = towns_sf, size = 2, color = "black") +
    geom_text_repel(
        data = towns,
        aes(x = long, y = lat, label = name),
        size = 3.5, fontface = "bold", color = "black",
        nudge_x = towns$nudge_x, nudge_y = towns$nudge_y,
        segment.color = "grey50", segment.size = 0.3
    ) +
    # Scale bar & North arrow
    annotation_scale(
        location = "bl", width_hint = 0.25,
        pad_x = unit(0.6, "in"), pad_y = unit(0.4, "in"),
        text_size = 8
    ) +
    annotation_north_arrow(
        location = "bl", which_north = "true",
        pad_x = unit(0.7, "in"), pad_y = unit(0.8, "in"),
        style = north_arrow_fancy_orienteering(text_size = 7)
    ) +
    # Scales
    scale_color_manual(
        name = "eDNA & Survey Window",
        values = c("3 Months" = "#2ca02c", "6 Months" = "#ff7f0e", "12 Months" = "#d62728", "> 12 Months" = "#7f7f7f"),
        labels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months / Unmatched")
    ) +
    scale_size_continuous(
        name = "Sampling Effort (Years)",
        range = c(2, 6.5),
        breaks = c(1, 2, 4, 7),
        labels = c("1 Year", "2 Years", "4 Years", "7 Years")
    ) +
    coord_sf(xlim = c(142.5, 153.5), ylim = c(-24.5, -10.5), expand = FALSE) +
    theme_bw() +
    theme(
        panel.background = element_rect(fill = "#edf4f9"), # Light ocean blue tint
        panel.grid.major = element_line(color = "#d9e4ec", linetype = "dashed", size = 0.3),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        legend.position = c(0.82, 0.35),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey60", size = 0.5),
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 8.5),
        legend.margin = margin(6, 8, 6, 8)
    )

# Inset map of Australia
cat("Building inset map...\n")
p_inset <- ggplot() +
    geom_sf(data = aus_land, fill = "#e8e8e8", color = "#999999", size = 0.2) +
    # Bounding box for GBR region
    geom_rect(
        aes(xmin = 142.5, xmax = 153.5, ymin = -24.5, ymax = -10.5),
        fill = NA, color = "red", size = 0.6
    ) +
    coord_sf(xlim = c(110, 155), ylim = c(-44, -9), expand = FALSE) +
    theme_void() +
    theme(
        panel.background = element_rect(fill = "white", color = "black", size = 0.5),
        plot.margin = margin(0, 0, 0, 0)
    )

# Combine main plot and inset using cowplot
cat("Combining main plot and inset...\n")
final_map <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset, x = 0.70, y = 0.72, width = 0.26, height = 0.24)

# Save test PNG
dir.create("plots", showWarnings = FALSE)
cat("Saving test plot...\n")
ggsave("plots/test_gbr_edna_map.png", final_map, width = 8.5, height = 11, dpi = 300)
cat("Saved plots/test_gbr_edna_map.png successfully!\n")
