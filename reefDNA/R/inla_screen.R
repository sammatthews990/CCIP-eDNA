# Lightweight INLA formulation screening for eDNA-informed COTS CPUE.

#' Build a one-event-per-cull INLA screening table
#'
#' Raw eDNA replicates are first summarized within sampling events. Each cull
#' dive is then linked to its most recent eligible preceding event; ties are
#' resolved by the shortest distance. This retains one response row per cull.
#'
#' @param links Output from [build_inla_space_time_links()].
#' @param edna Prepared raw eDNA observations from [prepare_inla_edna()].
#' @param culls Prepared culls, or a linked-cull subset containing unique
#'   `cull_id` values.
#' @return An `sf` object with one selected eDNA event per cull response.
#' @export
summarise_inla_event_benchmark <- function(links, edna, culls) {
  if (nrow(links) == 0L) stop("No admissible eDNA/cull links were supplied", call. = FALSE)
  if (!inherits(edna, "sf") || !inherits(culls, "sf")) {
    stop("edna and culls must both be sf objects", call. = FALSE)
  }
  if (anyDuplicated(culls$cull_id)) stop("culls must contain unique cull_id values", call. = FALSE)

  edna_data <- sf::st_drop_geometry(edna)
  linked <- dplyr::left_join(
    links,
    dplyr::transmute(
      edna_data,
      edna_id = .data$edna_id,
      event_id = .data$event_id,
      edna_date = .data$date_edna,
      edna_site_name = .data$site_name,
      detection = .data$detection,
      concentration = .data$concentration
    ),
    by = c("edna_id", "event_id")
  )

  event_links <- linked |>
    dplyr::group_by(.data$cull_id, .data$event_id) |>
    dplyr::summarise(
      edna_date = dplyr::first(.data$edna_date),
      edna_site_name = dplyr::first(.data$edna_site_name),
      edna_n_replicates = dplyr::n_distinct(.data$edna_id),
      edna_prop_positive = mean(.data$detection, na.rm = TRUE),
      edna_conc_mean = mean(.data$concentration, na.rm = TRUE),
      edna_conc_median = stats::median(.data$concentration, na.rm = TRUE),
      edna_distance_m = stats::median(.data$distance_m, na.rm = TRUE),
      edna_lag_days = stats::median(.data$lag_days, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("edna_conc_mean", "edna_conc_median")),
        ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
      )
    )

  selected <- event_links |>
    dplyr::arrange(.data$cull_id, .data$edna_lag_days, .data$edna_distance_m, .data$event_id) |>
    dplyr::group_by(.data$cull_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  cull_rows <- match(selected$cull_id, culls$cull_id)
  if (anyNA(cull_rows)) stop("Some selected cull IDs are absent from culls", call. = FALSE)
  benchmark <- culls[cull_rows, ]
  replacement_columns <- setdiff(names(selected), "cull_id")
  benchmark <- benchmark[, setdiff(names(benchmark), replacement_columns), drop = FALSE]
  benchmark <- dplyr::left_join(benchmark, selected, by = "cull_id")
  if (anyDuplicated(benchmark$cull_id)) stop("Event selection duplicated cull responses", call. = FALSE)
  benchmark
}

inla_screen_predictors <- function() {
  c("edna_prop_z", "edna_conc_z", "distance_z", "lag_z")
}

#' Generate all additive eDNA CPUE candidate formulations
#'
#' @param predictors Model-ready predictor names.
#' @param include_null Include the random-intercept-only formulation.
#' @return A named list of character vectors containing formula terms.
#' @export
build_inla_additive_candidates <- function(
    predictors = inla_screen_predictors(), include_null = TRUE) {
  labels <- c(
    edna_prop_z = "prop", edna_conc_z = "conc",
    distance_z = "distance", lag_z = "lag"
  )
  candidates <- list()
  if (include_null) candidates$null <- character()
  for (size in seq_along(predictors)) {
    combinations <- utils::combn(predictors, size, simplify = FALSE)
    for (terms in combinations) {
      name <- paste(unname(labels[terms]), collapse = "+")
      candidates[[name]] <- terms
    }
  }
  candidates
}

prepare_inla_screen_data <- function(benchmark, scale_spec = NULL) {
  data <- sf::st_drop_geometry(benchmark)
  required <- c(
    "cots_count", "bottom_time", "Reef", "edna_prop_positive",
    "edna_conc_mean", "edna_distance_m", "edna_lag_days"
  )
  assert_inla_columns(data, required, "benchmark")
  complete <- stats::complete.cases(data[, required]) & data$bottom_time > 0
  data <- data[complete, , drop = FALSE]
  if (nrow(data) < 30L) stop("At least 30 complete cull observations are required", call. = FALSE)

  data$edna_prop_model <- data$edna_prop_positive
  data$edna_conc_model <- log1p(pmax(data$edna_conc_mean, 0))
  data$distance_model <- log1p(pmax(data$edna_distance_m, 0) / 200)
  data$lag_model <- log1p(pmax(data$edna_lag_days, 0))
  variables <- c("edna_prop_model", "edna_conc_model", "distance_model", "lag_model")
  if (is.null(scale_spec)) scale_spec <- inla_scale_spec(data, variables)
  data <- apply_inla_scale(data, scale_spec)
  names(data)[match(paste0(variables, "_z"), names(data))] <- inla_screen_predictors()
  data$reef_index <- as.integer(factor(data$Reef))
  data$log_effort <- log(data$bottom_time)
  if ("VesselName_ID" %in% names(data)) {
    data$vessel_index <- as.integer(factor(data$VesselName_ID))
  }
  list(data = data, scale_spec = scale_spec, complete = complete)
}

#' Generate hierarchical pairwise-interaction candidates
#'
#' Each candidate retains all four main effects. Six models add one pairwise
#' interaction at a time and a final model adds all pairwise interactions.
#'
#' @return A named list of character vectors containing formula terms.
#' @export
build_inla_interaction_candidates <- function() {
  main <- inla_screen_predictors()
  pairs <- utils::combn(main, 2, simplify = FALSE)
  labels <- c(
    edna_prop_z = "prop", edna_conc_z = "conc",
    distance_z = "distance", lag_z = "lag"
  )
  candidates <- lapply(pairs, function(pair) c(main, paste(pair, collapse = ":")))
  names(candidates) <- vapply(
    pairs,
    function(pair) paste0("full+", paste(unname(labels[pair]), collapse = ":")),
    character(1)
  )
  candidates$full_all_pairwise <- c(
    main,
    vapply(pairs, paste, collapse = ":", FUN.VALUE = character(1))
  )
  candidates
}

#' Fit and rank lightweight negative-binomial INLA formulations
#'
#' All candidate models use the same complete rows, log effort offset, and reef
#' random intercept. An optional vessel random intercept is held constant across
#' candidates. No SPDE or calendar-year field is included: distance and lag
#' provide the requested spatial and temporal context directly.
#'
#' @param benchmark One-event-per-cull data from
#'   [summarise_inla_event_benchmark()].
#' @param candidates Named list of formula-term character vectors.
#' @param include_vessel Include a vessel iid intercept when multiple vessels
#'   are present.
#' @param num_threads INLA thread specification.
#' @return An object of class `reefdna_inla_screen` containing rankings and fits.
#' @export
fit_inla_cpue_candidates <- function(
    benchmark,
    candidates = build_inla_additive_candidates(),
    include_vessel = TRUE,
    num_threads = "1:1") {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required to fit candidate models", call. = FALSE)
  }
  if (is.null(names(candidates)) || any(names(candidates) == "")) {
    stop("candidates must be a fully named list", call. = FALSE)
  }
  allowed <- c(
    inla_screen_predictors(),
    vapply(utils::combn(inla_screen_predictors(), 2, simplify = FALSE),
      paste, collapse = ":", FUN.VALUE = character(1)
    )
  )
  supplied <- unique(unlist(candidates, use.names = FALSE))
  invalid <- setdiff(supplied, allowed)
  if (length(invalid) > 0L) {
    stop("Unsupported candidate terms: ", paste(invalid, collapse = ", "), call. = FALSE)
  }

  prepared <- prepare_inla_screen_data(benchmark)
  data <- prepared$data
  use_vessel <- include_vessel && "vessel_index" %in% names(data) &&
    dplyr::n_distinct(data$vessel_index) > 1L
  random_terms <- paste0(
    "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))",
    if (use_vessel) {
      " + f(vessel_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))"
    } else ""
  )

  fits <- vector("list", length(candidates))
  names(fits) <- names(candidates)
  ranking_rows <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    terms <- candidates[[i]]
    fixed <- if (length(terms) > 0L) paste(terms, collapse = " + ") else "1"
    formula <- stats::as.formula(paste(
      "cots_count ~", fixed, "+ offset(log_effort) +", random_terms
    ))
    fit <- INLA::inla(
      formula,
      family = "nbinomial",
      data = data,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE),
      control.family = list(
        hyper = list(theta = list(prior = "pc.mgamma", param = 7))
      ),
      num.threads = num_threads,
      verbose = FALSE
    )
    fits[[i]] <- fit
    cpo <- fit$cpo$cpo
    valid_cpo <- is.finite(cpo) & cpo > 0
    ranking_rows[[i]] <- tibble::tibble(
      model = names(candidates)[[i]],
      terms = if (length(terms) > 0L) paste(terms, collapse = " + ") else "intercept only",
      n = nrow(data),
      n_fixed = nrow(fit$summary.fixed),
      waic = fit$waic$waic,
      effective_parameters_waic = fit$waic$p.eff,
      dic = fit$dic$dic,
      mean_log_cpo = mean(log(cpo[valid_cpo])),
      cpo_failures = sum(!valid_cpo)
    )
  }
  ranking <- dplyr::bind_rows(ranking_rows) |>
    dplyr::arrange(.data$waic) |>
    dplyr::mutate(
      delta_waic = .data$waic - min(.data$waic),
      waic_weight = exp(-0.5 * .data$delta_waic) /
        sum(exp(-0.5 * .data$delta_waic)),
      rank_waic = dplyr::row_number()
    )

  structure(
    list(
      ranking = ranking,
      fits = fits,
      candidates = candidates,
      data = data,
      scale_spec = prepared$scale_spec,
      complete = prepared$complete,
      include_vessel = use_vessel,
      reef_levels = levels(factor(data$Reef)),
      vessel_levels = if (use_vessel) levels(factor(data$VesselName_ID)) else NULL
    ),
    class = "reefdna_inla_screen"
  )
}

#' Extract fixed-effect summaries from an INLA formulation screen
#'
#' @param screen A fitted `reefdna_inla_screen` object.
#' @return A tibble containing fixed-effect summaries for every candidate.
#' @export
extract_inla_screen_coefficients <- function(screen) {
  dplyr::bind_rows(lapply(names(screen$fits), function(model_name) {
    out <- tibble::rownames_to_column(
      as.data.frame(screen$fits[[model_name]]$summary.fixed), "term"
    )
    dplyr::mutate(out, model = model_name, .before = 1)
  }))
}

#' Build distance-by-time scenarios from the latest eDNA site events
#'
#' @param edna Prepared raw eDNA observations from [prepare_inla_edna()].
#' @param distances_m Distance-from-sample scenarios in metres.
#' @param lag_days Days-since-sample scenarios.
#' @return An `sf` object with one row per site, distance, and lag combination.
#' @export
build_inla_distance_time_scenarios <- function(
    edna, distances_m = c(0, 200, 500, 1000, 2000),
    lag_days = c(0, 30, 90, 183)) {
  if (!inherits(edna, "sf")) stop("edna must be an sf object", call. = FALSE)
  events <- edna |>
    dplyr::group_by(.data$Reef, .data$site_name, .data$event_id, .data$date_edna) |>
    dplyr::summarise(
      edna_prop_positive = mean(.data$detection, na.rm = TRUE),
      edna_conc_mean = mean(.data$concentration, na.rm = TRUE),
      edna_n_replicates = dplyr::n(),
      geometry = sf::st_centroid(sf::st_union(.data$geometry)),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Reef, .data$site_name, dplyr::desc(.data$date_edna)) |>
    dplyr::group_by(.data$Reef, .data$site_name) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  scenario_grid <- expand.grid(
    site_row = seq_len(nrow(events)),
    edna_distance_m = distances_m,
    edna_lag_days = lag_days,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenarios <- events[scenario_grid$site_row, ]
  scenarios$edna_distance_m <- scenario_grid$edna_distance_m
  scenarios$edna_lag_days <- scenario_grid$edna_lag_days
  scenarios$bottom_time <- 1
  scenarios$cots_count <- NA_real_
  scenarios$scenario_id <- paste(
    scenarios$Reef, scenarios$site_name,
    paste0(scenarios$edna_distance_m, "m"),
    paste0(scenarios$edna_lag_days, "d"),
    sep = "__"
  )
  scenarios
}

#' Refit a screened formulation with prediction scenarios
#'
#' @param screen A fitted `reefdna_inla_screen` object.
#' @param model_name Candidate name in `screen$candidates`.
#' @param newdata Scenario data containing raw eDNA, distance, lag, reef, and
#'   effort fields.
#' @param num_threads INLA thread specification.
#' @return A list containing scenario posterior summaries and the prediction fit.
#' @export
predict_inla_cpue_screen <- function(
    screen, model_name, newdata, num_threads = "1:1") {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required for prediction", call. = FALSE)
  }
  if (!model_name %in% names(screen$candidates)) {
    stop("Unknown model_name: ", model_name, call. = FALSE)
  }
  scenario <- sf::st_drop_geometry(newdata)
  required <- c(
    "bottom_time", "Reef", "edna_prop_positive", "edna_conc_mean",
    "edna_distance_m", "edna_lag_days"
  )
  assert_inla_columns(scenario, required, "newdata")
  if (any(!is.finite(scenario$bottom_time) | scenario$bottom_time <= 0)) {
    stop("newdata bottom_time must be positive", call. = FALSE)
  }

  scenario$edna_prop_model <- scenario$edna_prop_positive
  scenario$edna_conc_model <- log1p(pmax(scenario$edna_conc_mean, 0))
  scenario$distance_model <- log1p(pmax(scenario$edna_distance_m, 0) / 200)
  scenario$lag_model <- log1p(pmax(scenario$edna_lag_days, 0))
  scenario <- apply_inla_scale(scenario, screen$scale_spec)
  raw_names <- paste0(names(screen$scale_spec), "_z")
  names(scenario)[match(raw_names, names(scenario))] <- inla_screen_predictors()
  scenario$reef_index <- match(as.character(scenario$Reef), screen$reef_levels)
  if (anyNA(scenario$reef_index)) stop("newdata contain reefs absent from the screen", call. = FALSE)
  scenario$log_effort <- log(scenario$bottom_time)

  if (screen$include_vessel) {
    assert_inla_columns(scenario, "VesselName_ID", "newdata")
    scenario$vessel_index <- match(as.character(scenario$VesselName_ID), screen$vessel_levels)
    if (anyNA(scenario$vessel_index)) stop("newdata contain unknown vessels", call. = FALSE)
  }
  scenario$cots_count <- NA_real_
  prediction_rows <- nrow(screen$data) + seq_len(nrow(scenario))
  combined <- dplyr::bind_rows(screen$data, scenario)

  terms <- screen$candidates[[model_name]]
  fixed <- if (length(terms) > 0L) paste(terms, collapse = " + ") else "1"
  random_terms <- paste0(
    "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))",
    if (screen$include_vessel) {
      " + f(vessel_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))"
    } else ""
  )
  formula <- stats::as.formula(paste(
    "cots_count ~", fixed, "+ offset(log_effort) +", random_terms
  ))
  fit <- INLA::inla(
    formula,
    family = "nbinomial",
    data = combined,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(waic = TRUE, dic = TRUE),
    control.family = list(
      hyper = list(theta = list(prior = "pc.mgamma", param = 7))
    ),
    num.threads = num_threads,
    verbose = FALSE
  )

  fitted <- fit$summary.fitted.values[prediction_rows, , drop = FALSE]
  linear <- fit$summary.linear.predictor[prediction_rows, , drop = FALSE]
  effort <- scenario$bottom_time
  predictions <- newdata
  predictions$cpue_mean <- fitted$mean / effort
  predictions$cpue_q025 <- fitted$`0.025quant` / effort
  predictions$cpue_q10 <- exp(linear$mean + stats::qnorm(0.10) * linear$sd) / effort
  predictions$cpue_q50 <- fitted$`0.5quant` / effort
  predictions$cpue_q90 <- exp(linear$mean + stats::qnorm(0.90) * linear$sd) / effort
  predictions$cpue_q975 <- fitted$`0.975quant` / effort
  predictions$prob_cpue_gt_002 <- stats::pnorm(
    log(0.02 * effort), linear$mean, linear$sd, lower.tail = FALSE
  )
  predictions$prob_cpue_gt_004 <- stats::pnorm(
    log(0.04 * effort), linear$mean, linear$sd, lower.tail = FALSE
  )
  list(predictions = predictions, fit = fit, model_name = model_name)
}
