library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(purrr)
library(tidyr)
library(sf)
library(rsample)
library(tibble)

devtools::load_all("reefDNA")

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull <- read_excel(cull_file, sheet = "Cull") %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

edna_site <- read_excel(edna_file, sheet = "eDNA_data_ALL") %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Year)) %>%
    mutate(
        date_edna = as.Date(Date),
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    group_by(ReefName, Site_name, Year) %>%
    summarise(
        Lat = mean(Lat, na.rm = TRUE),
        Long = mean(Long, na.rm = TRUE),
        date_edna = min(date_edna),
        conc_mean = mean(Conc_mean, na.rm = TRUE),
        perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        .groups   = "drop"
    ) %>%
    rename(Reef = ReefName)

edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

edna_buf <- st_buffer(edna_sf, dist = 1000)
intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>%
    filter(!is.na(date_cull))

valid_encounters <- intersect_cull %>%
    mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
    filter(diff_days >= 0 & diff_days <= 183)

dat_site_1000m <- valid_encounters %>%
    st_drop_geometry() %>%
    rename(any_of(c(Reef = "Reef.x", Year = "Year.x"))) %>%
    group_by(Reef, Site_name, Year, perc_pos, conc_mean) %>%
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

cpue_target <- 0.04
perc_grid <- seq(0, 100, by = 1)

grid_eval <- tibble(perc_thresh = perc_grid) %>%
    mutate(
        out = map(perc_thresh, ~ metrics_for(.x, cpue_target, dat_site_1000m)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
        Precision = vapply(out, `[[`, numeric(1), "Precision"),
        Recall = vapply(out, `[[`, numeric(1), "Recall")
    )

print(grid_eval %>% filter(Recall >= 0.90), n = 30)
