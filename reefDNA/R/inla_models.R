inla_scale_spec <- function(data, variables) {
  stats <- lapply(variables, function(variable) {
    value <- data[[variable]]
    center <- mean(value, na.rm = TRUE)
    scale <- stats::sd(value, na.rm = TRUE)
    if (!is.finite(scale) || scale == 0) scale <- 1
    c(center = center, scale = scale)
  })
  names(stats) <- variables
  stats
}

apply_inla_scale <- function(data, spec) {
  for (variable in names(spec)) {
    data[[paste0(variable, "_z")]] <-
      (data[[variable]] - spec[[variable]][["center"]]) / spec[[variable]][["scale"]]
  }
  data
}

prepare_inla_fixed_effects <- function(data, scale_spec, reef_levels, time_levels, time_variable) {
  data$edna_prop_positive_model <- data$edna_prop_positive
  data$edna_log_concentration <- log1p(pmax(data$edna_conc_mean, 0))
  data$edna_log_distance <- log1p(pmax(data$edna_distance_mean_m, 0) / 200)
  data$edna_log_lag <- log1p(pmax(data$edna_lag_median_days, 0))
  data <- apply_inla_scale(data, scale_spec)
  data$reef_index <- match(as.character(data$Reef), reef_levels)
  data$time_index <- match(as.character(data[[time_variable]]), time_levels)
  if (anyNA(data$reef_index)) stop("Prediction data contain an unknown reef level", call. = FALSE)
  if (anyNA(data$time_index)) stop("Prediction data contain an unknown time level", call. = FALSE)

  data.frame(
    intercept = 1,
    edna_prop_z = data$edna_prop_positive_model_z,
    edna_conc_z = data$edna_log_concentration_z,
    distance_z = data$edna_log_distance_z,
    lag_z = data$edna_log_lag_z,
    reef_index = data$reef_index,
    log_effort = log(data$bottom_time)
  )
}

#' Build site-level distance scenarios for INLA prediction
#'
#' @param edna Prepared raw eDNA points from [prepare_inla_edna()].
#' @param distances_m Distance scenarios in metres.
#' @param prediction_date Date represented by the prediction scenarios.
#' @param lag_days Assumed eDNA-to-cull lag for the scenario.
#' @return An `sf` object with one row per site and distance scenario.
#' @export
build_inla_distance_scenarios <- function(
    edna, distances_m = c(0, 200, 500, 1000, 2000),
    prediction_date = max(edna$date_edna), lag_days = 0) {
  if (!inherits(edna, "sf")) stop("edna must be an sf object", call. = FALSE)
  prediction_date <- as.Date(prediction_date)

  sites <- edna |>
    dplyr::group_by(.data$Reef, .data$site_name) |>
    dplyr::summarise(
      edna_prop_positive = mean(.data$detection, na.rm = TRUE),
      edna_conc_mean = mean(.data$concentration, na.rm = TRUE),
      n_edna = dplyr::n(),
      geometry = sf::st_centroid(sf::st_union(.data$geometry)),
      .groups = "drop"
    )
  site_rows <- rep(seq_len(nrow(sites)), each = length(distances_m))
  scenarios <- sites[site_rows, ]
  scenarios$distance_scenario_m <- rep(distances_m, times = nrow(sites))
  scenarios$edna_distance_mean_m <- scenarios$distance_scenario_m
  scenarios$edna_lag_median_days <- lag_days
  scenarios$bottom_time <- 1
  scenarios$date_cull <- prediction_date
  scenarios$year <- as.integer(format(prediction_date, "%Y"))
  scenarios$quarter <- paste0(
    scenarios$year, "-Q",
    ((as.integer(format(prediction_date, "%m")) - 1L) %/% 3L) + 1L
  )
  scenarios$scenario_id <- paste(
    scenarios$Reef, scenarios$site_name, scenarios$distance_scenario_m,
    sep = "__"
  )
  scenarios
}

#' Fit the Stage A spatio-temporal INLA CPUE benchmark
#'
#' Fits negative-binomial cull counts with log bottom time as an offset, eDNA
#' summary covariates, a reef iid effect, and a Matérn SPDE field grouped through
#' time with an AR(1) dependence model.
#'
#' @param benchmark One-row-per-cull `sf` data from [summarise_inla_benchmark()].
#' @param mesh An INLA mesh, normally from [build_inla_spde_mesh()].
#' @param prediction Optional `sf` prediction scenarios.
#' @param time_resolution Either `year` or `quarter`.
#' @param prior_range_m Spatial PC-prior range threshold in metres.
#' @param prior_range_prob Probability that range is below `prior_range_m`.
#' @param prior_sigma Spatial PC-prior marginal SD threshold.
#' @param prior_sigma_prob Probability that SD exceeds `prior_sigma`.
#' @param num_threads INLA thread specification.
#' @return An object of class `reefdna_inla_benchmark`.
#' @export
fit_inla_cpue_benchmark <- function(
    benchmark, mesh, prediction = NULL, time_resolution = c("year", "quarter"),
    prior_range_m = 25000, prior_range_prob = 0.05,
    prior_sigma = 1, prior_sigma_prob = 0.05, num_threads = "1:1") {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required to fit this model", call. = FALSE)
  }
  if (!inherits(benchmark, "sf")) stop("benchmark must be an sf object", call. = FALSE)
  time_resolution <- match.arg(time_resolution)
  time_variable <- time_resolution

  model_data <- sf::st_drop_geometry(benchmark)
  required <- c(
    "cots_count", "bottom_time", "Reef", time_variable,
    "edna_prop_positive", "edna_conc_mean", "edna_distance_mean_m",
    "edna_lag_median_days"
  )
  assert_inla_columns(model_data, required, "benchmark")
  complete <- stats::complete.cases(model_data[, required])
  benchmark <- benchmark[complete, ]
  model_data <- model_data[complete, ]
  if (nrow(model_data) < 30L) stop("At least 30 complete cull observations are required", call. = FALSE)

  transform_variables <- c(
    "edna_prop_positive_model", "edna_log_concentration",
    "edna_log_distance", "edna_log_lag"
  )
  model_data$edna_prop_positive_model <- model_data$edna_prop_positive
  model_data$edna_log_concentration <- log1p(pmax(model_data$edna_conc_mean, 0))
  model_data$edna_log_distance <- log1p(pmax(model_data$edna_distance_mean_m, 0) / 200)
  model_data$edna_log_lag <- log1p(pmax(model_data$edna_lag_median_days, 0))
  scale_spec <- inla_scale_spec(model_data, transform_variables)

  prediction_data <- if (is.null(prediction)) NULL else sf::st_drop_geometry(prediction)
  reef_levels <- sort(unique(c(as.character(model_data$Reef), as.character(prediction_data$Reef))))
  time_levels <- sort(unique(c(as.character(model_data[[time_variable]]), as.character(prediction_data[[time_variable]]))))
  fixed_est <- prepare_inla_fixed_effects(model_data, scale_spec, reef_levels, time_levels, time_variable)

  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,
    prior.range = c(prior_range_m, prior_range_prob),
    prior.sigma = c(prior_sigma, prior_sigma_prob)
  )
  spatial_index <- INLA::inla.spde.make.index(
    "spatial", n.spde = spde$n.spde, n.group = length(time_levels)
  )
  A_est <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = sf::st_coordinates(benchmark),
    group = match(as.character(model_data[[time_variable]]), time_levels),
    n.group = length(time_levels)
  )
  stack_est <- INLA::inla.stack(
    data = list(y = model_data$cots_count),
    A = list(A_est, 1),
    effects = list(spatial_index, fixed_est),
    tag = "est"
  )

  stack <- stack_est
  if (!is.null(prediction)) {
    assert_inla_columns(prediction_data, required[-1], "prediction")
    fixed_pred <- prepare_inla_fixed_effects(
      prediction_data, scale_spec, reef_levels, time_levels, time_variable
    )
    A_pred <- INLA::inla.spde.make.A(
      mesh = mesh,
      loc = sf::st_coordinates(prediction),
      group = match(as.character(prediction_data[[time_variable]]), time_levels),
      n.group = length(time_levels)
    )
    stack_pred <- INLA::inla.stack(
      data = list(y = rep(NA_real_, nrow(prediction_data))),
      A = list(A_pred, 1),
      effects = list(spatial_index, fixed_pred),
      tag = "pred"
    )
    stack <- INLA::inla.stack(stack_est, stack_pred)
  }

  formula <- y ~ 0 + intercept + edna_prop_z + edna_conc_z + distance_z + lag_z +
    offset(log_effort) +
    f(
      reef_index,
      model = "iid",
      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
    ) +
    f(
      spatial,
      model = spde,
      group = spatial.group,
      control.group = list(model = "ar1")
    )

  fit <- INLA::inla(
    formula,
    family = "nbinomial",
    data = INLA::inla.stack.data(stack),
    control.predictor = list(
      A = INLA::inla.stack.A(stack),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE),
    control.family = list(
      hyper = list(theta = list(prior = "pc.mgamma", param = 7))
    ),
    num.threads = num_threads,
    verbose = FALSE
  )

  structure(
    list(
      fit = fit,
      stack = stack,
      mesh = mesh,
      spde = spde,
      benchmark = benchmark,
      prediction = prediction,
      scale_spec = scale_spec,
      reef_levels = reef_levels,
      time_levels = time_levels,
      time_resolution = time_resolution,
      formula = formula,
      config = list(
        prior_range_m = prior_range_m,
        prior_range_prob = prior_range_prob,
        prior_sigma = prior_sigma,
        prior_sigma_prob = prior_sigma_prob
      )
    ),
    class = "reefdna_inla_benchmark"
  )
}

#' Extract fitted or scenario CPUE posterior summaries
#'
#' @param model A fitted `reefdna_inla_benchmark` object.
#' @param tag Either `est` for fitted observations or `pred` for scenarios.
#' @return An `sf` object containing posterior response summaries and CPUE.
#' @export
extract_inla_cpue_predictions <- function(model, tag = c("pred", "est")) {
  tag <- match.arg(tag)
  source <- if (tag == "pred") model$prediction else model$benchmark
  if (is.null(source)) stop("The requested stack tag is not present", call. = FALSE)
  index <- INLA::inla.stack.index(model$stack, tag = tag)$data
  summary <- model$fit$summary.fitted.values[index, , drop = FALSE]
  source$expected_count_mean <- summary$mean
  source$expected_count_sd <- summary$sd
  source$expected_count_q025 <- summary$`0.025quant`
  source$expected_count_q50 <- summary$`0.5quant`
  source$expected_count_q975 <- summary$`0.975quant`
  effort <- source$bottom_time
  source$cpue_mean <- source$expected_count_mean / effort
  source$cpue_sd <- source$expected_count_sd / effort
  source$cpue_q025 <- source$expected_count_q025 / effort
  source$cpue_q50 <- source$expected_count_q50 / effort
  source$cpue_q975 <- source$expected_count_q975 / effort
  marginals <- model$fit$marginals.fitted.values[index]
  if (length(marginals) == length(index)) {
    source$cpue_q10 <- vapply(
      seq_along(marginals),
      function(i) INLA::inla.qmarginal(0.10, marginals[[i]]) / effort[[i]],
      numeric(1)
    )
    source$cpue_q90 <- vapply(
      seq_along(marginals),
      function(i) INLA::inla.qmarginal(0.90, marginals[[i]]) / effort[[i]],
      numeric(1)
    )
    source$prob_cpue_gt_002 <- vapply(
      seq_along(marginals),
      function(i) 1 - INLA::inla.pmarginal(0.02 * effort[[i]], marginals[[i]]),
      numeric(1)
    )
    source$prob_cpue_gt_004 <- vapply(
      seq_along(marginals),
      function(i) 1 - INLA::inla.pmarginal(0.04 * effort[[i]], marginals[[i]]),
      numeric(1)
    )
    attr(source, "threshold_probability_method") <- "response marginal"
  } else {
    linear <- model$fit$summary.linear.predictor[index, , drop = FALSE]
    source$cpue_q10 <- exp(linear$mean + stats::qnorm(0.10) * linear$sd) / effort
    source$cpue_q90 <- exp(linear$mean + stats::qnorm(0.90) * linear$sd) / effort
    source$prob_cpue_gt_002 <- stats::pnorm(
      log(0.02 * effort), mean = linear$mean, sd = linear$sd,
      lower.tail = FALSE
    )
    source$prob_cpue_gt_004 <- stats::pnorm(
      log(0.04 * effort), mean = linear$mean, sd = linear$sd,
      lower.tail = FALSE
    )
    attr(source, "threshold_probability_method") <-
      "Gaussian approximation to the posterior linear predictor"
  }
  source
}

#' Summarize diagnostics for an INLA CPUE benchmark
#'
#' @param model A fitted `reefdna_inla_benchmark` object.
#' @return A named list containing information criteria and CPO/PIT summaries.
#' @export
summarise_inla_diagnostics <- function(model) {
  fit <- model$fit
  estimation_index <- INLA::inla.stack.index(model$stack, tag = "est")$data
  cpo <- fit$cpo$cpo[estimation_index]
  pit <- fit$cpo$pit[estimation_index]
  list(
    information_criteria = tibble::tibble(
      dic = fit$dic$dic,
      effective_parameters_dic = fit$dic$p.eff,
      waic = fit$waic$waic,
      effective_parameters_waic = fit$waic$p.eff
    ),
    predictive = tibble::tibble(
      mean_log_cpo = mean(log(cpo[is.finite(cpo) & cpo > 0])),
      cpo_failures = sum(!is.finite(cpo) | cpo <= 0),
      pit_mean = mean(pit, na.rm = TRUE),
      pit_sd = stats::sd(pit, na.rm = TRUE)
    )
  )
}
