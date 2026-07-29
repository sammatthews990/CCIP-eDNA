#' Analyze Cull Site Closure Episodes
#'
#' Takes prepared cull data and builds site-year episodes, tracking the
#' chronological CPUE journey for each cull site within a year. A "closure"
#' is defined as a site that had at least one dive with CPUE >= 0.04 and
#' whose final dive(s) dropped below 0.04.
#'
#' @param cull_data Prepared cull data from \code{\link{prepare_cull_data}}.
#'   Must contain columns: CullSiteName, Year, SurveyDate, Bottomtime,
#'   CPUE, Total_COTS, ReefName, ReefLabel.
#' @param cpue_threshold CPUE threshold for "active" outbreak status
#'   (default 0.04).
#' @param exclude_years Integer vector of years to exclude (default 2026
#'   for incomplete data).
#'
#' @return A tibble with one row per site-year episode containing:
#' \describe{
#'   \item{Year, CullSiteName, ReefName, ReefLabel}{Identifiers}
#'   \item{n_dives}{Total number of dives at the site in that year}
#'   \item{total_hrs}{Total diver hours (Bottomtime / 60)}
#'   \item{initial_cpue}{CPUE of the first dive (chronologically)}
#'   \item{max_cpue}{Maximum CPUE observed across all dives}
#'   \item{last_cpue}{CPUE of the final dive}
#'   \item{first_date, last_date}{Date range of activity}
#'   \item{ever_above_threshold}{Logical — did any dive exceed the threshold?}
#'   \item{closed}{Logical — site was above threshold and last dive fell below}
#'   \item{initial_density_class}{Factor — density class based on initial CPUE}
#'   \item{max_density_class}{Factor — density class based on max CPUE}
#' }
#'
#' @export
#' @import dplyr
analyze_site_closures <- function(cull_data,
                                   cpue_threshold = 0.04,
                                   exclude_years = 2026) {

  density_levels <- c("<0.02", "0.02-0.04", "0.04-0.08", ">0.08")

  cull_data %>%
    filter(
      !is.na(CullSiteName),
      !(Year %in% exclude_years)
    ) %>%
    arrange(CullSiteName, Year, SurveyDate) %>%
    group_by(Year, CullSiteName, ReefName, ReefLabel) %>%
    summarise(
      n_dives = n(),
      total_hrs = sum(Bottomtime, na.rm = TRUE) / 60,
      initial_cpue = first(CPUE),
      max_cpue = max(CPUE, na.rm = TRUE),
      last_cpue = last(CPUE),
      first_date = min(as.Date(SurveyDate), na.rm = TRUE),
      last_date = max(as.Date(SurveyDate), na.rm = TRUE),
      ever_above_threshold = any(CPUE >= cpue_threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      closed = ever_above_threshold & (last_cpue < cpue_threshold),
      initial_density_class = factor(
        case_when(
          initial_cpue >= 0.08 ~ ">0.08",
          initial_cpue >= 0.04 ~ "0.04-0.08",
          initial_cpue >= 0.02 ~ "0.02-0.04",
          TRUE ~ "<0.02"
        ),
        levels = density_levels
      ),
      max_density_class = factor(
        case_when(
          max_cpue >= 0.08 ~ ">0.08",
          max_cpue >= 0.04 ~ "0.04-0.08",
          max_cpue >= 0.02 ~ "0.02-0.04",
          TRUE ~ "<0.02"
        ),
        levels = density_levels
      )
    )
}


#' Summarize Cull Site Effort
#'
#' Takes the output of \code{\link{analyze_site_closures}} and produces
#' summary tables for closure effort by density class and annual diver
#' hour allocation percentiles.
#'
#' @param site_episodes Output from \code{\link{analyze_site_closures}}.
#' @param ha_per_site Hectares managed per cull site (default 10).
#'
#' @return A list with components:
#' \describe{
#'   \item{closure_by_max_density}{Tibble — hours to close by max CPUE
#'     density class (n, mean, median, IQR, mean dives).}
#'   \item{closure_by_initial_density}{Tibble — same, stratified by
#'     initial CPUE density class.}
#'   \item{closure_overall}{Tibble — single-row summary across all
#'     closures (median, IQR, percentiles).}
#'   \item{annual_percentiles}{Tibble — hours-per-site-year at
#'     p10/p25/p50/p75/p90 across all site-years.}
#'   \item{annual_by_year}{Tibble — same percentiles broken down by year.}
#'   \item{ha_per_site}{Numeric — hectares per site constant.}
#' }
#'
#' @export
#' @import dplyr
summarize_site_effort <- function(site_episodes, ha_per_site = 10) {

  closed_sites <- site_episodes %>% filter(closed)

  # ── Closures by max density class ──
  closure_by_max <- closed_sites %>%
    group_by(max_density_class) %>%
    summarise(
      n_closures = n(),
      mean_hrs = round(mean(total_hrs, na.rm = TRUE), 1),
      median_hrs = round(median(total_hrs, na.rm = TRUE), 1),
      p25_hrs = round(quantile(total_hrs, 0.25, na.rm = TRUE), 1),
      p75_hrs = round(quantile(total_hrs, 0.75, na.rm = TRUE), 1),
      mean_dives = round(mean(n_dives, na.rm = TRUE), 1),
      .groups = "drop"
    )

  # ── Closures by initial density class ──
  closure_by_initial <- closed_sites %>%
    group_by(initial_density_class) %>%
    summarise(
      n_closures = n(),
      mean_hrs = round(mean(total_hrs, na.rm = TRUE), 1),
      median_hrs = round(median(total_hrs, na.rm = TRUE), 1),
      p25_hrs = round(quantile(total_hrs, 0.25, na.rm = TRUE), 1),
      p75_hrs = round(quantile(total_hrs, 0.75, na.rm = TRUE), 1),
      mean_dives = round(mean(n_dives, na.rm = TRUE), 1),
      .groups = "drop"
    )

  # ── Overall closure summary ──
  closure_overall <- tibble::tibble(
    n_closures = nrow(closed_sites),
    mean_hrs = round(mean(closed_sites$total_hrs, na.rm = TRUE), 1),
    median_hrs = round(median(closed_sites$total_hrs, na.rm = TRUE), 1),
    p10_hrs = round(quantile(closed_sites$total_hrs, 0.10, na.rm = TRUE), 1),
    p25_hrs = round(quantile(closed_sites$total_hrs, 0.25, na.rm = TRUE), 1),
    p75_hrs = round(quantile(closed_sites$total_hrs, 0.75, na.rm = TRUE), 1),
    p90_hrs = round(quantile(closed_sites$total_hrs, 0.90, na.rm = TRUE), 1)
  )

  # ── Annual percentiles: hours per site-year (all sites, not just closures) ──
  annual_percentiles <- tibble::tibble(
    p10_hrs = round(quantile(site_episodes$total_hrs, 0.10, na.rm = TRUE), 1),
    p25_hrs = round(quantile(site_episodes$total_hrs, 0.25, na.rm = TRUE), 1),
    p50_hrs = round(quantile(site_episodes$total_hrs, 0.50, na.rm = TRUE), 1),
    p75_hrs = round(quantile(site_episodes$total_hrs, 0.75, na.rm = TRUE), 1),
    p90_hrs = round(quantile(site_episodes$total_hrs, 0.90, na.rm = TRUE), 1),
    n_site_years = nrow(site_episodes)
  )

  # ── Annual percentiles broken down by year ──
  annual_by_year <- site_episodes %>%
    group_by(Year) %>%
    summarise(
      n_sites = n(),
      n_closed = sum(closed, na.rm = TRUE),
      total_hrs = round(sum(total_hrs, na.rm = TRUE), 0),
      p10_hrs = round(quantile(total_hrs, 0.10, na.rm = TRUE), 1),
      p25_hrs = round(quantile(total_hrs, 0.25, na.rm = TRUE), 1),
      p50_hrs = round(quantile(total_hrs, 0.50, na.rm = TRUE), 1),
      p75_hrs = round(quantile(total_hrs, 0.75, na.rm = TRUE), 1),
      p90_hrs = round(quantile(total_hrs, 0.90, na.rm = TRUE), 1),
      .groups = "drop"
    )

  list(
    closure_by_max_density = closure_by_max,
    closure_by_initial_density = closure_by_initial,
    closure_overall = closure_overall,
    annual_percentiles = annual_percentiles,
    annual_by_year = annual_by_year,
    ha_per_site = ha_per_site
  )
}
