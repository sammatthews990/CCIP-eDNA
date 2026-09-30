# Helpers for assigning the voyage-level eDNA sampling design.

#' Summarise eDNA sampling design by voyage
#'
#' Classifies voyages into the two principal protocols used in the analysis.
#' Minor shortfalls are accepted as incomplete replicates (5-6 for the 4 x 6
#' protocol and 11-12 for the 3 x 12 protocol). Voyages containing both
#' protocols or older designs are retained as `other_or_mixed`.
#'
#' @param edna Raw eDNA data containing Voyage, ReefName and Site_name.
#' @return A tibble with one row per voyage and an auditable design summary.
#' @export
summarise_edna_voyage_design <- function(edna) {
  assert_inla_columns(edna, c("Voyage", "ReefName", "Site_name"), "edna")

  site_reps <- edna |>
    dplyr::transmute(
      Voyage = as.character(.data$Voyage),
      ReefName = as.character(.data$ReefName),
      Site_name = as.character(.data$Site_name)
    ) |>
    dplyr::filter(
      !is.na(.data$Voyage), !is.na(.data$ReefName), !is.na(.data$Site_name)
    ) |>
    dplyr::count(.data$Voyage, .data$ReefName, .data$Site_name, name = "n_reps")

  reef_design <- site_reps |>
    dplyr::group_by(.data$Voyage, .data$ReefName) |>
    dplyr::summarise(
      n_sites = dplyr::n_distinct(.data$Site_name),
      rep_min = min(.data$n_reps),
      rep_max = max(.data$n_reps),
      rep_median = stats::median(.data$n_reps),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      reef_design = dplyr::case_when(
        .data$n_sites == 4L & .data$rep_min >= 5L & .data$rep_max <= 6L ~ "4x6",
        .data$n_sites <= 3L & .data$rep_min >= 11L & .data$rep_max <= 12L ~ "3x12",
        TRUE ~ "other"
      ),
      exact_pattern = dplyr::if_else(
        .data$rep_min == .data$rep_max,
        paste0(.data$n_sites, "x", .data$rep_min),
        paste0(
          .data$n_sites, "x", round(.data$rep_median),
          " [", .data$rep_min, "-", .data$rep_max, "]"
        )
      )
    )

  pattern_summary <- reef_design |>
    dplyr::count(.data$Voyage, .data$exact_pattern, name = "n_reefs") |>
    dplyr::arrange(.data$Voyage, .data$exact_pattern) |>
    dplyr::summarise(
      observed_designs = paste0(
        .data$exact_pattern, " (", .data$n_reefs, ")", collapse = "; "
      ),
      .by = .data$Voyage
    )

  reef_design |>
    dplyr::group_by(.data$Voyage) |>
    dplyr::summarise(
      sampling_design = dplyr::case_when(
        all(.data$reef_design == "3x12") ~ "3x12",
        all(.data$reef_design == "4x6") ~ "4x6",
        TRUE ~ "other_or_mixed"
      ),
      .groups = "drop"
    ) |>
    dplyr::left_join(pattern_summary, by = "Voyage")
}
