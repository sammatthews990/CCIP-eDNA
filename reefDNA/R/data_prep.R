#' Prepare Cull Data
#'
#' @param cull_file Path to the cull data excel file.
#' @return Cleaned cull dataset with CPUE and dive class labels.
#' @export
#' @import dplyr
#' @importFrom readxl read_excel
#' @importFrom lubridate year
#' @importFrom tidyr replace_na
prepare_cull_data <- function(cull_file) {
  read_excel(cull_file, sheet = "Cull") %>%
    mutate(
      Year = year(SurveyDate),
      Cohort1 = replace_na(as.numeric(Cohort1), 0),
      Cohort2 = replace_na(as.numeric(Cohort2), 0),
      Cohort3 = replace_na(as.numeric(Cohort3), 0),
      Cohort4 = replace_na(as.numeric(Cohort4), 0),
      Total_COTS = Cohort1 + Cohort2 + Cohort3 + Cohort4,
      Bottomtime = as.numeric(Bottomtime)
    ) %>%
    filter(!is.na(Bottomtime) & Bottomtime > 0, Year >= 2020) %>%
    mutate(
      CPUE = Total_COTS / Bottomtime,
      Wasted_Dive = CPUE < 0.02,
      Borderline_Dive = CPUE >= 0.02 & CPUE <= 0.04,
      Needed_Dive = CPUE > 0.04,
      NeededLow_Dive = CPUE > 0.04 & CPUE <= 0.08,
      NeededHigh_Dive = CPUE > 0.08,
      date_cull = as.Date(SurveyDate),
      Reef = ReefName
    )
}

#' Aggregate Cull Data
#'
#' @param cull_data Prepared cull data.
#' @return Summary dataset by Year and ReefLabel.
#' @export
#' @import dplyr
aggregate_cull_data <- function(cull_data) {
  cull_data %>%
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
}

#' Get Scout Metrics
#'
#' @param cull_data Prepared cull data.
#' @return Early scout metrics for the first day on a reef.
#' @export
#' @import dplyr
get_scout_metrics <- function(cull_data) {
  cull_data %>%
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
}

#' Prepare Manta Tow Data
#'
#' @param manta_file Path to manta tow csv file.
#' @return Cleaned manta tow dataset.
#' @export
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom lubridate year
#' @importFrom tidyr replace_na
prepare_manta_data <- function(manta_file) {
  read_csv(manta_file, show_col_types = FALSE) %>%
    mutate(
      Year = year(SurveyTime),
      ScarsCount = as.numeric(ScarsCount),
      CrownOfThornsStarfishCount = as.numeric(CrownOfThornsStarfishCount),
      ScarsPresent = (replace_na(ScarsCount, 0) > 0) | 
                     (tolower(FeedingScarCountRangeCode) %in% c("p", "c"))
    ) %>%
    filter(Year >= 2020)
}

#' Summarize Manta Tow Data
#'
#' @param manta_data Cleaned manta tow data
#' @return Summary by Reef and Year.
#' @export
#' @import dplyr
summarize_manta_data <- function(manta_data) {
  manta_data %>%
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
}

#' Summarize eDNA Data
#'
#' @param edna_file Path to eDNA excel file
#' @param sheet Sheet name, defaults to "eDNA_data_ALL"
#' @return Summarised eDNA data per reef/year, including earliest eDNA_Date.
#' @export
#' @import dplyr
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
summarize_edna_data <- function(edna_file, sheet = "eDNA_data_ALL") {
  read_excel(edna_file, sheet = sheet) %>%
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
}

#' Get Earliest eDNA Dates per Reef-Year
#'
#' Extracts the earliest eDNA sample date for each reef-year combination.
#' Used as a temporal anchor to filter cull data so that only post-eDNA
#' dives are counted in the cost-benefit analysis.
#'
#' @param edna_file Path to eDNA excel file.
#' @param sheet Sheet name, defaults to "eDNA_data_ALL".
#' @return Data frame with columns: Year, ReefLabel, date_edna.
#' @export
#' @import dplyr
#' @importFrom readxl read_excel
#' @importFrom stringr str_extract
get_edna_dates <- function(edna_file, sheet = "eDNA_data_ALL") {
  read_excel(edna_file, sheet = sheet) %>%
    mutate(
      ReefLabel = str_extract(ReefName, "[0-9]+-[0-9a-zA-Z]+"),
      date_edna = as.Date(Date)
    ) %>%
    filter(!is.na(ReefLabel), !is.na(date_edna)) %>%
    group_by(Year, ReefLabel) %>%
    summarise(date_edna = min(date_edna), .groups = "drop")
}

#' Filter Cull Data to Post-eDNA Dives Only
#'
#' Joins eDNA temporal anchors to cull data and removes any cull dives that
#' occurred before the earliest eDNA sample date for that reef-year.
#' Reefs without eDNA data (date_edna NA after join) retain all dives.
#'
#' @param cull_data Prepared cull data from [prepare_cull_data()].
#' @param edna_dates eDNA date anchors from [get_edna_dates()].
#' @return Filtered cull data containing only post-eDNA dives.
#' @export
#' @import dplyr
filter_cull_post_edna <- function(cull_data, edna_dates) {
  cull_data %>%
    left_join(edna_dates, by = c("Year", "ReefLabel")) %>%
    filter(is.na(date_edna) | date_cull >= date_edna) %>%
    select(-date_edna)
}

#' Process eDNA Survey Events
#'
#' @param edna_data Raw eDNA data.
#' @param gap_days Maximum days between samples to be considered the same event.
#' @return Aggregated event-level dataset.
#' @export
#' @import dplyr
process_edna_events <- function(edna_data, gap_days = 7) {
  edna_data %>%
    mutate(
      Reef = ReefName,
      date_edna = as.Date(Date),
      Conc_mean = as.numeric(Conc_mean)
    ) %>%
    arrange(Reef, Year, date_edna) %>%
    group_by(Reef, Year) %>%
    mutate(
      grp = cumsum(
        ifelse(is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > gap_days,
               1L, 0L
        )
      )
    ) %>%
    ungroup() %>%
    group_by(Reef, Year, grp) %>%
    summarise(
      date_edna = min(date_edna),
      conc_mean = mean(Conc_mean, na.rm = TRUE),
      perc_pos  = mean(LOD_sample_positive == 1 | LOD_sample_positive == "1", na.rm = TRUE) * 100,
      n_samples = n(),
      lat_edna  = mean(Lat, na.rm = TRUE),
      lon_edna  = mean(Long, na.rm = TRUE),
      .groups   = "drop"
    )
}

#' Combine All Reef Data
#'
#' @param cull_summary Checked and evaluated summarizations for cull.
#'   For temporally correct analysis, cull data should be pre-filtered
#'   via [filter_cull_post_edna()] before aggregation.
#' @param first_day_summary Scout dive metrics.
#' @param manta_summary Evaluated Manta metrics.
#' @param edna_summary Summary of eDNA traces.
#' @return Joined wide dataset combining all field campaigns.
#' @export
#' @import dplyr
combine_reef_data <- function(cull_summary, first_day_summary, manta_summary, edna_summary) {
  cull_summary %>%
    left_join(first_day_summary, by = c("Year", "ReefLabel")) %>%
    left_join(manta_summary, by = c("Year", "ReefLabel", "ReefName")) %>%
    left_join(edna_summary, by = c("Year", "ReefLabel"))
}
