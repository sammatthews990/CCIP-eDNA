# Helpers for the INLA data-linkage workflow.

assert_inla_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

#' Normalize a cull-site identifier
#'
#' @param x Character vector containing cull-site names.
#' @return Upper-case alphanumeric identifiers.
#' @export
normalize_cull_site_id <- function(x) {
  x <- trimws(as.character(x))
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

extract_reef_key <- function(x) {
  x <- as.character(x)
  hit <- stringr::str_extract(x, "[0-9]{1,2}-[0-9]{2,3}[A-Za-z]?")
  normalize_cull_site_id(hit)
}

#' Prepare raw eDNA observations for INLA
#'
#' Keeps one row per raw replicate and creates stable identifiers, detection and
#' concentration fields, time indices, and projected point geometry.
#'
#' @param edna Raw eDNA data frame.
#' @param crs_projected Projected CRS used for distance and mesh calculations.
#' @return An `sf` point object with one row per valid raw eDNA observation.
#' @export
prepare_inla_edna <- function(edna, crs_projected = 3112) {
  assert_inla_columns(
    edna,
    c("Date", "ReefName", "Site_name", "Lat", "Long", "LOD_sample_positive", "Conc_mean"),
    "edna"
  )

  source_row <- seq_len(nrow(edna))
  detection_chr <- toupper(trimws(as.character(edna$LOD_sample_positive)))

  out <- dplyr::mutate(
    edna,
    source_row = source_row,
    edna_id = sprintf("edna_%07d", source_row),
    date_edna = as.Date(.data$Date),
    Reef = as.character(.data$ReefName),
    site_name = as.character(.data$Site_name),
    detection = dplyr::case_when(
      detection_chr %in% c("1", "TRUE", "T", "YES", "Y") ~ 1L,
      detection_chr %in% c("0", "FALSE", "F", "NO", "N") ~ 0L,
      TRUE ~ NA_integer_
    ),
    concentration = suppressWarnings(as.numeric(.data$Conc_mean)),
    year = as.integer(format(.data$date_edna, "%Y")),
    quarter = paste0(.data$year, "-Q", ((as.integer(format(.data$date_edna, "%m")) - 1L) %/% 3L) + 1L),
    event_id = paste(.data$Reef, .data$site_name, .data$date_edna, sep = "__")
  )

  out <- dplyr::filter(
    out,
    !is.na(.data$date_edna),
    !is.na(.data$Lat), !is.na(.data$Long),
    dplyr::between(.data$Lat, -90, 90),
    dplyr::between(.data$Long, -180, 180),
    !is.na(.data$detection)
  )

  sf::st_as_sf(out, coords = c("Long", "Lat"), crs = 4326, remove = FALSE) |>
    sf::st_transform(crs_projected)
}

#' Prepare cull-dive observations for INLA
#'
#' @param culls Raw Cull-sheet data frame.
#' @param crs_projected Projected CRS used for distance and mesh calculations.
#' @param min_date Optional earliest retained survey date.
#' @param max_date Optional latest retained survey date.
#' @return An `sf` point object with one row per valid cull dive.
#' @export
prepare_inla_culls <- function(culls, crs_projected = 3112, min_date = NULL, max_date = NULL) {
  required <- c(
    "SurveyDate", "ReefName", "ReefLabel", "CullSiteName", "Latitude", "Longitude",
    "Bottomtime", "Depth", "Cohort1", "Cohort2", "Cohort3", "Cohort4"
  )
  assert_inla_columns(culls, required, "culls")

  source_row <- seq_len(nrow(culls))
  cohorts <- lapply(culls[c("Cohort1", "Cohort2", "Cohort3", "Cohort4")], function(x) {
    suppressWarnings(as.numeric(x))
  })
  counts <- rowSums(as.data.frame(cohorts), na.rm = TRUE)
  raw_id <- if ("CrownOfThornsStarfishCullDiveId" %in% names(culls)) {
    as.character(culls$CrownOfThornsStarfishCullDiveId)
  } else {
    rep(NA_character_, nrow(culls))
  }
  raw_id[is.na(raw_id) | raw_id == ""] <- sprintf("row_%07d", source_row[is.na(raw_id) | raw_id == ""])

  out <- dplyr::mutate(
    culls,
    source_row = source_row,
    cull_id = make.unique(raw_id, sep = "__dup"),
    date_cull = as.Date(.data$SurveyDate),
    Reef = as.character(.data$ReefName),
    site_name = as.character(.data$CullSiteName),
    site_key = normalize_cull_site_id(.data$site_name),
    reef_key = dplyr::coalesce(extract_reef_key(.data$ReefLabel), extract_reef_key(.data$site_name)),
    bottom_time = suppressWarnings(as.numeric(.data$Bottomtime)),
    cots_count = counts,
    cpue = .data$cots_count / .data$bottom_time,
    depth_m = suppressWarnings(as.numeric(.data$Depth)),
    year = as.integer(format(.data$date_cull, "%Y")),
    quarter = paste0(.data$year, "-Q", ((as.integer(format(.data$date_cull, "%m")) - 1L) %/% 3L) + 1L)
  )

  out <- dplyr::filter(
    out,
    !is.na(.data$date_cull),
    !is.na(.data$Latitude), !is.na(.data$Longitude),
    dplyr::between(.data$Latitude, -90, 90),
    dplyr::between(.data$Longitude, -180, 180),
    !is.na(.data$bottom_time), .data$bottom_time > 0,
    !is.na(.data$cots_count), .data$cots_count >= 0
  )
  if (!is.null(min_date)) out <- dplyr::filter(out, .data$date_cull >= as.Date(min_date))
  if (!is.null(max_date)) out <- dplyr::filter(out, .data$date_cull <= as.Date(max_date))

  sf::st_as_sf(out, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) |>
    sf::st_transform(crs_projected)
}

#' Prepare and dissolve static cull-site polygons
#'
#' @param sites An `sf` object containing cull-site polygons.
#' @param crs_projected Projected CRS used for distance calculations.
#' @return An `sf` object with one valid dissolved geometry per normalized site.
#' @export
prepare_inla_cull_sites <- function(sites, crs_projected = 3112) {
  if (!inherits(sites, "sf")) stop("sites must be an sf object", call. = FALSE)
  assert_inla_columns(sites, "Name", "sites")

  sites |>
    dplyr::mutate(
      site_name_polygon = as.character(.data$Name),
      site_key = normalize_cull_site_id(.data$site_name_polygon),
      reef_key = extract_reef_key(.data$site_name_polygon)
    ) |>
    dplyr::filter(!is.na(.data$site_key), .data$site_key != "") |>
    sf::st_make_valid() |>
    sf::st_transform(crs_projected) |>
    dplyr::group_by(.data$site_key, .data$reef_key) |>
    dplyr::summarise(
      site_name_polygon = dplyr::first(.data$site_name_polygon),
      .groups = "drop"
    )
}

#' Link cull dives to static cull-site polygons
#'
#' Uses normalized site name first and then a same-reef nearest-polygon fallback.
#'
#' @param culls Prepared cull points from [prepare_inla_culls()].
#' @param sites Prepared polygons from [prepare_inla_cull_sites()].
#' @param max_fallback_m Maximum accepted nearest-polygon fallback distance.
#' @return The cull `sf` object with link method, polygon name, and distance.
#' @export
link_inla_culls_to_sites <- function(culls, sites, max_fallback_m = 2000) {
  if (!inherits(culls, "sf") || !inherits(sites, "sf")) {
    stop("culls and sites must both be sf objects", call. = FALSE)
  }
  if (sf::st_crs(culls) != sf::st_crs(sites)) stop("culls and sites must share a CRS", call. = FALSE)

  out <- culls
  out$polygon_site_key <- NA_character_
  out$polygon_site_name <- NA_character_
  out$site_link_method <- "unmatched"
  out$distance_to_site_m <- NA_real_

  direct <- match(out$site_key, sites$site_key)
  direct_rows <- which(!is.na(direct))
  if (length(direct_rows) > 0L) {
    direct_distance <- sf::st_distance(
      out[direct_rows, ], sites[direct[direct_rows], ], by_element = TRUE
    )
    out$polygon_site_key[direct_rows] <- sites$site_key[direct[direct_rows]]
    out$polygon_site_name[direct_rows] <- sites$site_name_polygon[direct[direct_rows]]
    out$site_link_method[direct_rows] <- "normalized_name"
    out$distance_to_site_m[direct_rows] <- as.numeric(direct_distance)
  }

  unmatched <- which(is.na(direct))
  reef_groups <- split(unmatched, out$reef_key[unmatched], drop = TRUE)
  for (reef_key in names(reef_groups)) {
    rows <- reef_groups[[reef_key]]
    candidates <- which(!is.na(sites$reef_key) & sites$reef_key == reef_key)
    if (length(rows) == 0L || length(candidates) == 0L) next

    local_nearest <- sf::st_nearest_feature(out[rows, ], sites[candidates, ])
    polygon_rows <- candidates[local_nearest]
    distance <- as.numeric(sf::st_distance(out[rows, ], sites[polygon_rows, ], by_element = TRUE))
    accepted <- is.finite(distance) & distance <= max_fallback_m
    accepted_rows <- rows[accepted]
    accepted_polygons <- polygon_rows[accepted]
    if (length(accepted_rows) == 0L) next

    out$polygon_site_key[accepted_rows] <- sites$site_key[accepted_polygons]
    out$polygon_site_name[accepted_rows] <- sites$site_name_polygon[accepted_polygons]
    out$site_link_method[accepted_rows] <- "nearest_same_reef"
    out$distance_to_site_m[accepted_rows] <- distance[accepted]
  }

  out
}

#' Build raw eDNA-to-cull space-time links
#'
#' @param edna Prepared raw eDNA points.
#' @param culls Prepared cull-dive points.
#' @param max_distance_m Maximum spatial distance in metres.
#' @param min_lag_days Minimum signed lag, where positive means eDNA precedes cull.
#' @param max_lag_days Maximum signed lag.
#' @param same_reef Require exact normalized reef-name equality.
#' @return A tibble with one row per admissible raw eDNA/cull pair.
#' @export
build_inla_space_time_links <- function(
    edna, culls, max_distance_m = 2000, min_lag_days = 0,
    max_lag_days = 183, same_reef = TRUE) {
  if (!inherits(edna, "sf") || !inherits(culls, "sf")) {
    stop("edna and culls must both be sf objects", call. = FALSE)
  }
  if (sf::st_crs(edna) != sf::st_crs(culls)) stop("edna and culls must share a CRS", call. = FALSE)
  if (sf::st_is_longlat(edna)) stop("A projected CRS is required for metric distances", call. = FALSE)

  candidates <- sf::st_is_within_distance(culls, edna, dist = max_distance_m)
  n_candidates <- lengths(candidates)
  if (sum(n_candidates) == 0L) return(tibble::tibble())

  links <- tibble::tibble(
    cull_index = rep.int(seq_along(candidates), n_candidates),
    edna_index = unlist(candidates, use.names = FALSE)
  )
  edna_data <- sf::st_drop_geometry(edna)
  cull_data <- sf::st_drop_geometry(culls)
  links$lag_days <- as.numeric(
    cull_data$date_cull[links$cull_index] - edna_data$date_edna[links$edna_index]
  )
  links$same_reef <- cull_data$Reef[links$cull_index] == edna_data$Reef[links$edna_index]
  links <- dplyr::filter(
    links,
    .data$lag_days >= min_lag_days,
    .data$lag_days <= max_lag_days
  )
  if (isTRUE(same_reef)) links <- dplyr::filter(links, .data$same_reef)
  if (nrow(links) == 0L) return(links)

  links$distance_m <- as.numeric(sf::st_distance(
    culls[links$cull_index, ], edna[links$edna_index, ], by_element = TRUE
  ))
  links$cull_id <- cull_data$cull_id[links$cull_index]
  links$edna_id <- edna_data$edna_id[links$edna_index]
  links$event_id <- edna_data$event_id[links$edna_index]
  dplyr::filter(links, .data$distance_m <= max_distance_m)
}

#' Build one benchmark row per cull dive
#'
#' Aggregates linked raw eDNA observations before joining them to the cull
#' response, preventing multiplication of cull counts by replicate count.
#'
#' @param links Output from [build_inla_space_time_links()].
#' @param edna Prepared raw eDNA points.
#' @param culls Prepared cull-dive points.
#' @return An `sf` object with one unique row per linked cull dive.
#' @export
summarise_inla_benchmark <- function(links, edna, culls) {
  if (nrow(links) == 0L) stop("No admissible eDNA/cull links were supplied", call. = FALSE)
  edna_data <- sf::st_drop_geometry(edna)

  linked <- dplyr::left_join(
    links,
    dplyr::transmute(
      edna_data,
      edna_id = .data$edna_id,
      detection = .data$detection,
      concentration = .data$concentration,
      edna_date = .data$date_edna
    ),
    by = "edna_id"
  )

  summaries <- linked |>
    dplyr::group_by(.data$cull_index, .data$cull_id) |>
    dplyr::summarise(
      n_edna = dplyr::n(),
      n_edna_events = dplyr::n_distinct(.data$event_id),
      edna_prop_positive = mean(.data$detection, na.rm = TRUE),
      edna_conc_mean = mean(.data$concentration, na.rm = TRUE),
      edna_conc_median = stats::median(.data$concentration, na.rm = TRUE),
      edna_conc_positive_mean = mean(.data$concentration[.data$detection == 1L], na.rm = TRUE),
      edna_distance_min_m = min(.data$distance_m, na.rm = TRUE),
      edna_distance_mean_m = mean(.data$distance_m, na.rm = TRUE),
      edna_lag_min_days = min(.data$lag_days, na.rm = TRUE),
      edna_lag_median_days = stats::median(.data$lag_days, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("edna_conc_mean", "edna_conc_median", "edna_conc_positive_mean")),
        ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
      )
    )

  benchmark <- culls[summaries$cull_index, ]
  benchmark$cull_index <- summaries$cull_index
  benchmark <- dplyr::left_join(benchmark, summaries, by = c("cull_index", "cull_id"))
  if (anyDuplicated(benchmark$cull_id)) stop("Benchmark construction duplicated cull responses", call. = FALSE)
  benchmark
}

#' Summarize INLA linkage coverage
#'
#' @param edna Prepared raw eDNA points.
#' @param culls Prepared cull-dive points, optionally polygon-linked.
#' @param links Raw space-time link table.
#' @param benchmark One-row-per-cull benchmark table.
#' @return A named list of audit tibbles.
#' @export
audit_inla_linkage <- function(edna, culls, links, benchmark) {
  cull_data <- sf::st_drop_geometry(culls)
  benchmark_data <- sf::st_drop_geometry(benchmark)
  link_method <- if ("site_link_method" %in% names(cull_data)) {
    dplyr::count(cull_data, .data$site_link_method, name = "n_culls")
  } else {
    tibble::tibble(site_link_method = "not_run", n_culls = nrow(cull_data))
  }

  list(
    totals = tibble::tibble(
      n_edna_rows = nrow(edna),
      n_cull_rows = nrow(culls),
      n_raw_links = nrow(links),
      n_linked_culls = nrow(benchmark),
      n_linked_reefs = dplyr::n_distinct(benchmark_data$Reef),
      first_cull_date = min(benchmark_data$date_cull),
      last_cull_date = max(benchmark_data$date_cull)
    ),
    by_year = benchmark_data |>
      dplyr::count(.data$year, name = "n_linked_culls") |>
      dplyr::arrange(.data$year),
    polygon_link_method = link_method,
    link_distance = tibble::tibble(
      minimum_m = min(links$distance_m),
      median_m = stats::median(links$distance_m),
      mean_m = mean(links$distance_m),
      maximum_m = max(links$distance_m)
    ),
    replicate_check = tibble::tibble(
      duplicated_cull_ids = sum(duplicated(benchmark_data$cull_id)),
      max_rows_per_cull_id = max(table(benchmark_data$cull_id))
    )
  )
}
