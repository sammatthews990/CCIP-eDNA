#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(sf)
})

project_root <- normalizePath(
  if (file.exists("analysis/inla/output/space_time_links.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
devtools::load_all(file.path(project_root, "reefDNA"), quiet = TRUE)
output_dir <- file.path(project_root, "analysis", "site", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

edna_raw <- read_excel(
  file.path(project_root, "data", "eDNA data_ALL_20260528.xlsx"),
  sheet = "eDNA_data_ALL"
)
cull_raw <- read_excel(
  file.path(project_root, "data", "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"),
  sheet = "Cull"
)
edna <- prepare_inla_edna(edna_raw, crs_projected = 3112)
voyage_design <- summarise_edna_voyage_design(edna_raw)
culls <- prepare_inla_culls(
  cull_raw, crs_projected = 3112,
  min_date = min(edna$date_edna), max_date = max(edna$date_edna) + 183
)
links <- readRDS(file.path(project_root, "analysis", "inla", "output", "space_time_links.rds"))

event_metadata <- edna |>
  st_drop_geometry() |>
  group_by(Reef, site_name, event_id, date_edna) |>
  summarise(
    Voyage = dplyr::first(as.character(Voyage)),
    n_voyages = n_distinct(as.character(Voyage)),
    edna_n_replicates = n_distinct(edna_id),
    edna_prop_positive = mean(detection, na.rm = TRUE),
    edna_conc_mean = mean(concentration, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(voyage_design, by = "Voyage") |>
  arrange(Reef, date_edna, event_id) |>
  group_by(Reef) |>
  mutate(
    new_campaign = is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > 7,
    campaign_number = cumsum(new_campaign),
    edna_campaign_id = paste(Reef, campaign_number, sep = "__campaign_")
  ) |>
  ungroup() |>
  select(-new_campaign, -campaign_number)

if (any(event_metadata$n_voyages != 1L)) {
  stop("At least one eDNA event contains samples from multiple voyages.")
}
campaign_voyages <- event_metadata |>
  count(edna_campaign_id, Voyage, name = "n_events") |>
  count(edna_campaign_id, name = "n_voyages")
if (any(campaign_voyages$n_voyages != 1L)) {
  stop("At least one seven-day eDNA campaign spans multiple voyages.")
}
if (anyNA(event_metadata$sampling_design)) {
  stop("Sampling design could not be assigned to every eDNA event.")
}

cull_metadata <- culls |>
  st_drop_geometry() |>
  transmute(
    cull_id, Reef, cull_site_name = site_name, date_cull,
    cots_count, bottom_time
  )

# Collapse replicate-level links first so every cull response contributes once.
event_cull <- links |>
  group_by(event_id, cull_id) |>
  summarise(
    lag_days = first(lag_days),
    distance_m = median(distance_m),
    .groups = "drop"
  ) |>
  inner_join(cull_metadata, by = "cull_id")

# A visit is all cull dives at one named site on one date.
event_site_visits <- event_cull |>
  group_by(
    event_id, Reef, cull_site_name, date_cull, lag_days
  ) |>
  summarise(
    distance_m = median(distance_m),
    cots_count = sum(cots_count),
    bottom_time = sum(bottom_time),
    n_cull_dives = n_distinct(cull_id),
    .groups = "drop"
  ) |>
  mutate(
    cpue = cots_count / bottom_time,
    site_visit_id = paste(Reef, cull_site_name, date_cull, sep = "__")
  ) |>
  inner_join(event_metadata, by = c("event_id", "Reef"))

# Literal event-first version: first later visit at each cull site for each
# individual eDNA site/date event. Outcomes can be reused by several events.
event_first_by_site <- event_site_visits |>
  arrange(event_id, cull_site_name, lag_days, distance_m, site_visit_id) |>
  group_by(event_id, cull_site_name) |>
  slice_head(n = 1L) |>
  ungroup()

event_first_overall <- event_site_visits |>
  arrange(event_id, lag_days, distance_m, site_visit_id) |>
  group_by(event_id) |>
  slice_head(n = 1L) |>
  ungroup()

# Operational campaign version:
# 1. treat eDNA sites sampled within a seven-day reef campaign together;
# 2. assign each candidate visit to the nearest eDNA site in that campaign;
# 3. retain the first subsequent visit to each cull site; and
# 4. if campaigns compete for the same outcome, retain the most recent one.
campaign_visit <- event_site_visits |>
  arrange(edna_campaign_id, site_visit_id, distance_m, lag_days, event_id) |>
  group_by(edna_campaign_id, site_visit_id) |>
  slice_head(n = 1L) |>
  ungroup()

campaign_first_overall <- campaign_visit |>
  arrange(edna_campaign_id, lag_days, distance_m, site_visit_id) |>
  group_by(edna_campaign_id) |>
  slice_head(n = 1L) |>
  ungroup()

campaign_first_by_site <- campaign_visit |>
  arrange(edna_campaign_id, cull_site_name, lag_days, distance_m, site_visit_id) |>
  group_by(edna_campaign_id, cull_site_name) |>
  slice_head(n = 1L) |>
  ungroup()

unique_campaign_site_visits <- campaign_first_by_site |>
  arrange(site_visit_id, lag_days, distance_m, edna_campaign_id) |>
  group_by(site_visit_id) |>
  slice_head(n = 1L) |>
  ungroup()

if (anyDuplicated(unique_campaign_site_visits$site_visit_id)) {
  stop("The final site-model data contain duplicated site visits.")
}

summarise_design <- function(data, design) {
  tibble(
    design = design,
    n_rows = nrow(data),
    n_reefs = n_distinct(data$Reef),
    n_edna_events = n_distinct(data$event_id),
    n_edna_campaigns = n_distinct(data$edna_campaign_id),
    n_cull_sites = n_distinct(paste(data$Reef, data$cull_site_name)),
    n_site_visits = n_distinct(data$site_visit_id),
    duplicated_outcome_rows = nrow(data) - n_distinct(data$site_visit_id),
    median_lag_days = median(data$lag_days),
    median_distance_m = median(data$distance_m)
  )
}

audit <- bind_rows(
  summarise_design(event_first_overall, "event_first_single_visit"),
  summarise_design(event_first_by_site, "event_first_each_site"),
  summarise_design(campaign_first_overall, "campaign_first_single_visit"),
  summarise_design(campaign_first_by_site, "campaign_first_each_site_before_deduplication"),
  summarise_design(unique_campaign_site_visits, "unique_campaign_to_next_site_visit")
)

write.csv(audit, file.path(output_dir, "edna_first_linkage_audit.csv"), row.names = FALSE)
saveRDS(
  unique_campaign_site_visits,
  file.path(output_dir, "edna_first_site_visit_data.rds")
)
write.csv(
  st_drop_geometry(unique_campaign_site_visits),
  file.path(output_dir, "edna_first_site_visit_data.csv"), row.names = FALSE
)
print(audit, width = Inf)
