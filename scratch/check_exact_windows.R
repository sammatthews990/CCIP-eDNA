library(dplyr)
library(readxl)
library(lubridate)

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
manta.dat <- read_excel(cull_file, sheet = "Manta")
rhis.dat <- read_excel(cull_file, sheet = "RHIS")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

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
        years = paste(sort(unique(Year)), collapse = ", "),
        n_samples = n(),
        .groups = "drop"
    )

edna_dates <- edna.dat %>%
    filter(!is.na(Date)) %>%
    transmute(Reef = ReefName, date_edna = as.Date(Date)) %>%
    distinct()

# 1. Cull dates
cull_dates <- cull.dat %>%
    filter(!is.na(SurveyDate), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyDate), type = "Cull") %>%
    distinct()

# 2. Manta dates
manta_dates <- manta.dat %>%
    filter(!is.na(SurveyTime), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyTime), type = "Manta") %>%
    distinct()

# 3. RHIS dates
rhis_dates <- rhis.dat %>%
    filter(!is.na(SurveyTime), !is.na(ReefName)) %>%
    transmute(Reef = ReefName, date_survey = as.Date(SurveyTime), type = "RHIS") %>%
    distinct()

all_surveys <- bind_rows(cull_dates, manta_dates, rhis_dates) %>% distinct()

# Calculate min difference for each eDNA date to any survey date on the same reef
get_reef_window <- function(survey_df, label_suffix = "") {
    pairs <- edna_dates %>%
        inner_join(survey_df, by = "Reef", relationship = "many-to-many") %>%
        mutate(diff_days = as.numeric(date_survey - date_edna))
    
    # We test both:
    # 1) Prior/concurrent eDNA: diff_days >= -30 & diff_days <= 365 (eDNA before or up to 30 days after survey)
    # 2) Absolute days: abs(diff_days)
    
    # Absolute days nearest match per eDNA sampling date
    abs_match <- pairs %>%
        group_by(Reef, date_edna) %>%
        summarise(min_abs_diff = min(abs(diff_days)), .groups = "drop") %>%
        group_by(Reef) %>%
        summarise(reef_min_diff = min(min_abs_diff), .groups = "drop")
    
    # Prior match (eDNA prior to cull within 0..365 days)
    prior_match <- pairs %>%
        filter(diff_days >= 0) %>%
        group_by(Reef, date_edna) %>%
        summarise(min_pos_diff = min(diff_days), .groups = "drop") %>%
        group_by(Reef) %>%
        summarise(reef_min_pos = min(min_pos_diff), .groups = "drop")
    
    res <- edna_reef %>%
        left_join(abs_match, by = "Reef") %>%
        left_join(prior_match, by = "Reef") %>%
        mutate(
            window_abs = case_when(
                is.na(reef_min_diff) ~ "> 12 Months / No match",
                reef_min_diff <= 91 ~ "Within 3 Months",
                reef_min_diff <= 183 ~ "Within 6 Months",
                reef_min_diff <= 365 ~ "Within 12 Months",
                TRUE ~ "> 12 Months / No match"
            ),
            window_prior = case_when(
                is.na(reef_min_pos) ~ "> 12 Months / No match",
                reef_min_pos <= 91 ~ "Within 3 Months",
                reef_min_pos <= 183 ~ "Within 6 Months",
                reef_min_pos <= 365 ~ "Within 12 Months",
                TRUE ~ "> 12 Months / No match"
            )
        )
    return(res)
}

cat("--- Cull Only Alignment ---\n")
res_cull <- get_reef_window(cull_dates)
print("Absolute diff window:")
print(table(res_cull$window_abs))
print("Prior eDNA window:")
print(table(res_cull$window_prior))

cat("\n--- All Surveys (Cull + Manta + RHIS) Alignment ---\n")
res_all <- get_reef_window(all_surveys)
print("Absolute diff window:")
print(table(res_all$window_abs))
print("Prior eDNA window:")
print(table(res_all$window_prior))
