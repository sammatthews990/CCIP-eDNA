suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

edna_path <- "data/eDNA data_ALL_20260528.xlsx"
output_dir <- "analysis/manuscript"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
devtools::load_all("reefDNA", quiet = TRUE)

edna <- read_excel(edna_path, sheet = "eDNA_data_ALL") %>%
  mutate(
    Date = as.Date(Date),
    Voyage = as.character(Voyage),
    ReefName = as.character(ReefName),
    Site_name = as.character(Site_name)
  )

voyages <- read_excel(edna_path, sheet = "eDNAVoyages") %>%
  transmute(
    Vessel = as.character(Vessel),
    Voyage = as.character(Voyage),
    start_date = as.Date(`Departure Date`),
    end_date = as.Date(`Return Date`)
  )

stopifnot(n_distinct(voyages$Voyage) == nrow(voyages))
stopifnot(setequal(unique(edna$Voyage), voyages$Voyage))

# Use the exact deduplicated response rows consumed by the site model.
voyage_design <- summarise_edna_voyage_design(edna) %>%
  mutate(
    sampling_method = recode(
      sampling_design,
      `3x12` = "3 sites x 12 reps",
      `4x6` = "4 sites x 6 reps",
      other_or_mixed = "Other/mixed design"
    )
  )

model_data_path <- "analysis/site/output/edna_first_site_visit_data.rds"
if (!file.exists(model_data_path)) {
  stop("Run analysis/site/audit_edna_first_linkage.R before building the voyage table.")
}
model_rows <- readRDS(model_data_path)
required_model_columns <- c(
  "site_visit_id", "Voyage", "sampling_design", "n_cull_dives"
)
if (!all(required_model_columns %in% names(model_rows))) {
  stop("The site-model data predate voyage/design linkage; rebuild them first.")
}
if (anyDuplicated(model_rows$site_visit_id)) {
  stop("The site-model response contains duplicated site_visit_id values.")
}

match_summary <- model_rows %>%
  group_by(Voyage) %>%
  summarise(
    culling_matches = n(),
    unique_cull_dives = sum(n_cull_dives),
    .groups = "drop"
  )

sample_summary <- edna %>%
  group_by(Voyage) %>%
  summarise(
    reefs_visited = n_distinct(ReefName, na.rm = TRUE),
    sites = n_distinct(Site_name, na.rm = TRUE),
    funding_program = paste(sort(unique(na.omit(as.character(Project)))), collapse = "; "),
    total_samples = n(),
    sample_first_date = min(Date, na.rm = TRUE),
    sample_last_date = max(Date, na.rm = TRUE),
    .groups = "drop"
  )

voyage_summary <- voyages %>%
  left_join(sample_summary, by = "Voyage") %>%
  left_join(voyage_design, by = "Voyage") %>%
  left_join(match_summary, by = "Voyage") %>%
  mutate(
    culling_matches = replace_na(culling_matches, 0L),
    unique_cull_dives = replace_na(unique_cull_dives, 0L)
  ) %>%
  arrange(start_date, Voyage)

stopifnot(nrow(voyage_summary) == 36L)
stopifnot(sum(voyage_summary$total_samples) == nrow(edna))
stopifnot(sum(voyage_summary$culling_matches) == nrow(model_rows))
stopifnot(sum(voyage_summary$unique_cull_dives) == sum(model_rows$n_cull_dives))

manuscript_table <- voyage_summary %>%
  select(
    Vessel,
    Voyage,
    start_date,
    end_date,
    reefs_visited,
    sites,
    funding_program,
    sampling_method,
    total_samples,
    culling_matches
  )

audit_table <- voyage_summary %>%
  select(
    Voyage,
    observed_designs,
    sample_first_date,
    sample_last_date,
    unique_cull_dives,
    culling_matches
  )

write_csv(manuscript_table, file.path(output_dir, "voyage_summary_table.csv"), na = "")
write_csv(audit_table, file.path(output_dir, "voyage_summary_audit.csv"), na = "")

cat(sprintf(
  "Created %d-voyage summary: %d samples, %d site-model matches, %d underlying cull dives.\n",
  nrow(manuscript_table), sum(manuscript_table$total_samples),
  sum(manuscript_table$culling_matches), sum(model_rows$n_cull_dives)
))
