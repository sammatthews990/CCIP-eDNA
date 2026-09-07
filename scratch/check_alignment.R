library(dplyr)
library(readxl)
library(fuzzyjoin)
library(lubridate)

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
manta.dat <- read_excel(cull_file, sheet = "Manta")
rhis.dat <- read_excel(cull_file, sheet = "RHIS")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

# Clean eDNA
edna_reefs <- edna.dat %>%
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
        years_list = paste(sort(unique(Year)), collapse = ", "),
        n_samples = n(),
        min_date_edna = min(date_edna),
        max_date_edna = max(date_edna),
        .groups = "drop"
    )

# Prepare all survey dates (Cull, Manta, RHIS)
cull_dates <- cull.dat %>%
    filter(!is.na(SurveyDate), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyDate), source = "Cull")

manta_dates <- manta.dat %>%
    filter(!is.na(SurveyTime), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyTime), source = "Manta")

rhis_dates <- rhis.dat %>%
    filter(!is.na(SurveyTime), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyTime), source = "RHIS")

all_surveys <- bind_rows(cull_dates, manta_dates, rhis_dates)

# Also prepare cull dates alone
cull_only <- cull_dates

cat(sprintf("Total eDNA reefs with coords: %d\n", nrow(edna_reefs)))
cat("eDNA Reefs by n_years:\n")
print(table(edna_reefs$n_years))

# For each eDNA reef, find the minimum difference to any Cull/Manta/RHIS survey
# We can check both directional (cull after/with eDNA, or eDNA prior to cull) and bidirectional (closest survey within X months)

# Let's check for each eDNA sampling event (aggregated by Reef & date_edna)
edna_events <- edna.dat %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Date)) %>%
    mutate(
        Reef = ReefName,
        date_edna = as.Date(Date)
    ) %>%
    group_by(Reef, date_edna) %>%
    summarise(
        Lat = mean(Lat, na.rm = TRUE),
        Long = mean(Long, na.rm = TRUE),
        .groups = "drop"
    )

# Method A: Matching Cull only (as in current .qmd script)
# Find min absolute diff in days for each edna event to cull survey on same reef
match_cull <- edna_events %>%
    left_join(cull_only, by = "Reef", relationship = "many-to-many") %>%
    mutate(diff_days = as.numeric(date_survey - date_edna)) %>%
    # let's look at absolute difference or positive difference
    group_by(Reef, date_edna, Lat, Long) %>%
    summarise(
        min_abs_diff = if (all(is.na(diff_days))) NA_real_ else min(abs(diff_days), na.rm = TRUE),
        min_pos_diff = if (all(is.na(diff_days[diff_days >= 0]))) NA_real_ else min(diff_days[diff_days >= 0], na.rm = TRUE),
        .groups = "drop"
    )

# At Reef level, determine best window match across its eDNA sampling events
reef_window_cull_abs <- match_cull %>%
    group_by(Reef, Lat, Long) %>%
    summarise(
        min_diff = min(min_abs_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(edna_reefs %>% select(Reef, n_years), by = "Reef") %>%
    mutate(
        window = case_when(
            is.na(min_diff) ~ "> 12 Months / No match",
            min_diff <= 91 ~ "Within 3 Months",
            min_diff <= 183 ~ "Within 6 Months",
            min_diff <= 365 ~ "Within 12 Months",
            TRUE ~ "> 12 Months / No match"
        ),
        window = factor(window, levels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months / No match"))
    )

cat("\n=== Reef level matching with Cull data (Absolute Days) ===\n")
print(table(reef_window_cull_abs$window))

# Method B: Matching All Surveys (Cull, Manta, RHIS)
match_all <- edna_events %>%
    left_join(all_surveys, by = "Reef", relationship = "many-to-many") %>%
    mutate(diff_days = as.numeric(date_survey - date_edna)) %>%
    group_by(Reef, date_edna, Lat, Long) %>%
    summarise(
        min_abs_diff = if (all(is.na(diff_days))) NA_real_ else min(abs(diff_days), na.rm = TRUE),
        min_pos_diff = if (all(is.na(diff_days[diff_days >= 0]))) NA_real_ else min(diff_days[diff_days >= 0], na.rm = TRUE),
        .groups = "drop"
    )

reef_window_all_abs <- match_all %>%
    group_by(Reef, Lat, Long) %>%
    summarise(
        min_diff = min(min_abs_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(edna_reefs %>% select(Reef, n_years), by = "Reef") %>%
    mutate(
        window = case_when(
            is.na(min_diff) ~ "> 12 Months / No match",
            min_diff <= 91 ~ "Within 3 Months",
            min_diff <= 61 ~ "Within 3 Months",
            min_diff <= 91 ~ "Within 3 Months",
            min_diff <= 183 ~ "Within 6 Months",
            min_diff <= 365 ~ "Within 12 Months",
            TRUE ~ "> 12 Months / No match"
        ),
        window = factor(window, levels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months / No match"))
    )

cat("\n=== Reef level matching with ALL Surveys (Cull + Manta + RHIS) (Absolute Days) ===\n")
print(table(reef_window_all_abs$window))

# Directional matching (eDNA prior to survey, 0 to N days)
reef_window_all_dir <- match_all %>%
    group_by(Reef, Lat, Long) %>%
    summarise(
        min_diff = if (all(is.na(min_pos_diff))) NA_real_ else min(min_pos_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(edna_reefs %>% select(Reef, n_years), by = "Reef") %>%
    mutate(
        window = case_when(
            is.na(min_diff) ~ "> 12 Months / No match",
            min_diff <= 91 ~ "Within 3 Months",
            min_diff <= 183 ~ "Within 6 Months",
            min_diff <= 365 ~ "Within 12 Months",
            TRUE ~ "> 12 Months / No match"
        ),
        window = factor(window, levels = c("Within 3 Months", "Within 6 Months", "Within 12 Months", "> 12 Months / No match"))
    )

cat("\n=== Reef level matching with ALL Surveys (Directional eDNA -> Survey) ===\n")
print(table(reef_window_all_dir$window))
