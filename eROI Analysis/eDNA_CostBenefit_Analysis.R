# eDNA Cost-Benefit Analysis for COTS Culling (eROI Analysis)
# This script processes historical COTS Cull and Manta Tow data to identify 
# "wasted", "borderline", and "needed" diver hours, and context from eDNA & Manta Tow.
#
# TEMPORAL CORRECTION (2026-04-23):
# All cull metrics are now computed using ONLY dives that occurred ON or AFTER
# the earliest eDNA sample date for that reef-year. This ensures:
#   1. We only count hours that could actually have been avoided via eDNA decisions
#   2. Ghost-cull exoneration (>0.04 CPUE) is not contaminated by pre-eDNA culling
# For reefs WITHOUT eDNA, all dives are retained (they are excluded from the
# confusion matrix analysis anyway).

library(dplyr)
library(readxl)
library(lubridate)
library(readr)
library(tidyr)
library(stringr)

# Set base path (assuming this script is run from the project root or eROI folder)
root_dir <- if(dir.exists("data")) "." else ".."

cull_file <- file.path(root_dir, "data/260201_COTS-Cull-Data-Ewels.xlsx")
manta_file <- file.path(root_dir, "data/COTS Program  Manta Tow Data-2026-02-04.csv")
edna_file <- file.path(root_dir, "data/eDNA data_ALL_20260225.xlsx")

# =============================================================================
# 1. Load and prepare Cull Data (all dives, before temporal filtering)
# =============================================================================
cull_raw <- read_excel(cull_file, sheet = "Cull")

cull_data <- cull_raw %>%
  mutate(
    Year = year(SurveyDate),
    Cohort1 = replace_na(as.numeric(Cohort1), 0),
    Cohort2 = replace_na(as.numeric(Cohort2), 0),
    Cohort3 = replace_na(as.numeric(Cohort3), 0),
    Cohort4 = replace_na(as.numeric(Cohort4), 0),
    Total_COTS = Cohort1 + Cohort2 + Cohort3 + Cohort4,
    Bottomtime = as.numeric(Bottomtime),
    date_cull = as.Date(SurveyDate)
  ) %>%
  filter(!is.na(Bottomtime) & Bottomtime > 0, Year >= 2020) %>%
  mutate(
    CPUE = Total_COTS / Bottomtime,
    Wasted_Dive = CPUE < 0.02,
    Borderline_Dive = CPUE >= 0.02 & CPUE <= 0.04,
    Needed_Dive = CPUE > 0.04,
    NeededLow_Dive = CPUE > 0.04 & CPUE <= 0.08,
    NeededHigh_Dive = CPUE > 0.08
  )

# =============================================================================
# 2. Load eDNA Data FIRST to get temporal anchor dates
# =============================================================================
edna_raw <- read_excel(edna_file)

# Extract earliest eDNA sample date per reef-year
edna_dates <- edna_raw %>%
  mutate(
    ReefLabel = str_extract(ReefName, "[0-9]+-[0-9a-zA-Z]+"),
    date_edna = as.Date(Date)
  ) %>%
  filter(!is.na(ReefLabel), !is.na(date_edna)) %>%
  group_by(Year, ReefLabel) %>%
  summarise(date_edna = min(date_edna), .groups = "drop")

message("eDNA temporal anchors computed: ", nrow(edna_dates), " reef-year combinations")

# =============================================================================
# 3. TEMPORAL FILTER: Only keep cull dives ON or AFTER the eDNA sample date
#    Reefs without eDNA (date_edna = NA) retain all their dives
# =============================================================================
cull_data_temporal <- cull_data %>%
  left_join(edna_dates, by = c("Year", "ReefLabel")) %>%
  filter(is.na(date_edna) | date_cull >= date_edna)

n_removed <- nrow(cull_data) - nrow(cull_data_temporal)
message("Temporal filter: removed ", n_removed, " pre-eDNA cull dives (",
        round(n_removed / nrow(cull_data) * 100, 1), "% of total)")

# =============================================================================
# 4. Aggregate Cull Data per Reef and Year (POST-eDNA dives only)
# =============================================================================
cull_summary <- cull_data_temporal %>%
  group_by(Year, ReefName, ReefLabel) %>%
  summarise(
    Total_Dives = n(),
    Total_Diver_Hours = sum(Bottomtime) / 60,
    Wasted_Diver_Hours = sum(Bottomtime[Wasted_Dive], na.rm = TRUE) / 60,
    Borderline_Diver_Hours = sum(Bottomtime[Borderline_Dive], na.rm = TRUE) / 60,
    Needed_Diver_Hours = sum(Bottomtime[Needed_Dive], na.rm = TRUE) / 60,
    NeededLow_Diver_Hours = sum(Bottomtime[NeededLow_Dive], na.rm = TRUE) / 60,
    NeededHigh_Diver_Hours = sum(Bottomtime[NeededHigh_Dive], na.rm = TRUE) / 60,
    Total_COTS_Culled = sum(Total_COTS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Proportion_Wasted = ifelse(Total_Diver_Hours > 0, Wasted_Diver_Hours / Total_Diver_Hours, 0),
    Proportion_Borderline = ifelse(Total_Diver_Hours > 0, Borderline_Diver_Hours / Total_Diver_Hours, 0),
    Proportion_Needed = ifelse(Total_Diver_Hours > 0, Needed_Diver_Hours / Total_Diver_Hours, 0),
    Proportion_NeededLow = ifelse(Total_Diver_Hours > 0, NeededLow_Diver_Hours / Total_Diver_Hours, 0),
    Proportion_NeededHigh = ifelse(Total_Diver_Hours > 0, NeededHigh_Diver_Hours / Total_Diver_Hours, 0)
  )

# =============================================================================
# 5. Scout Culling Data: First Day Metrics (POST-eDNA only)
#    The "scout" is the first cull day AFTER the eDNA sample was taken
# =============================================================================
first_day_summary <- cull_data_temporal %>%
  group_by(Year, ReefLabel) %>%
  arrange(SurveyDate) %>%
  filter(SurveyDate == min(SurveyDate, na.rm = TRUE)) %>%
  summarise(
    Scout_Diver_Hours = sum(Bottomtime, na.rm = TRUE) / 60,
    Scout_Total_COTS = sum(Total_COTS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Initial_CPUE = Scout_Total_COTS / (Scout_Diver_Hours * 60)
  ) %>%
  select(Year, ReefLabel, Scout_Diver_Hours, Initial_CPUE)

# =============================================================================
# 6. Load and prepare Manta Tow Data (unchanged - not temporally filtered)
# =============================================================================
manta_raw <- read_csv(manta_file, show_col_types = FALSE)

manta_data <- manta_raw %>%
  mutate(
    Year = year(SurveyTime),
    ScarsCount = as.numeric(ScarsCount),
    CrownOfThornsStarfishCount = as.numeric(CrownOfThornsStarfishCount),
    ScarsPresent = (replace_na(ScarsCount, 0) > 0) | 
                   (tolower(FeedingScarCountRangeCode) %in% c("p", "c"))
  ) %>%
  filter(Year >= 2020)

manta_summary <- manta_data %>%
  group_by(Year, ReefName, ReefLabel) %>%
  summarise(
    Total_Tows = n(),
    Mean_COTS_per_Tow = mean(CrownOfThornsStarfishCount, na.rm = TRUE),
    Tows_With_Scars = sum(ScarsPresent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Pct_Tows_With_Scars = ifelse(Total_Tows > 0, (Tows_With_Scars / Total_Tows) * 100, 0)
  )

# =============================================================================
# 7. Prepare eDNA Summary (now includes eDNA_Date for transparency)
# =============================================================================
edna_summary <- edna_raw %>%
  mutate(
    ReefLabel = str_extract(ReefName, "[0-9]+-[0-9a-zA-Z]+")
  ) %>%
  filter(!is.na(ReefLabel)) %>%
  group_by(Year, ReefLabel) %>%
  summarise(
    Total_eDNA_Samples = n(),
    eDNA_Positives = sum(LOD_sample_positive == 1 | LOD_sample_positive == "1", na.rm = TRUE),
    eDNA_Mean_Conc = mean(as.numeric(Conc_mean), na.rm = TRUE),
    eDNA_Date = min(as.Date(Date), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eDNA_Pct_Positive = (eDNA_Positives / Total_eDNA_Samples) * 100
  )

# =============================================================================
# 8. Join ALL data together
# =============================================================================
combined_data <- cull_summary %>%
  left_join(first_day_summary, by = c("Year", "ReefLabel")) %>%
  left_join(manta_summary, by = c("Year", "ReefLabel", "ReefName")) %>%
  left_join(edna_summary, by = c("Year", "ReefLabel"))

# 9. Extract Full Matched Dataset (eDNA + Cull both exist)
full_matched_data <- combined_data %>%
  filter(!is.na(eDNA_Pct_Positive))

message("Full matched dataset: ", nrow(full_matched_data), " reef-year records with both eDNA and post-eDNA culling")

# Keep the top 20 logic for the original deep dive section
top_20_annual <- combined_data %>%
  group_by(Year) %>%
  slice_max(order_by = Wasted_Diver_Hours, n = 20, with_ties = FALSE) %>%
  ungroup()

# =============================================================================
# 10. Save datasets
# =============================================================================
write_csv(combined_data, file.path(root_dir, "eROI Analysis", "combined_reef_analysis_v2.csv"))
write_csv(top_20_annual, file.path(root_dir, "eROI Analysis", "top_20_annual_wasted_v2.csv"))
write_csv(full_matched_data, file.path(root_dir, "eROI Analysis", "full_matched_data_v2.csv"))

message("Done. All CSVs written with temporal correction applied.")
