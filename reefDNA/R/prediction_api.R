.resolve_prediction_column <- function(data, candidates, label, default = NULL) {
  available <- intersect(candidates, names(data))
  if (length(available)) return(data[[available[[1]]]])
  if (!is.null(default)) return(rep(default, nrow(data)))
  stop(
    "Missing ", label, ". Supply one of: ", paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

.prepare_operational_predictors <- function(newdata, scale_spec, default_effort = 216) {
  data <- as.data.frame(newdata)
  if (!nrow(data)) stop("newdata must contain at least one row.", call. = FALSE)
  percent_positive <- .resolve_prediction_column(
    data, c("perc_pos", "edna_percent_positive", "edna_pct"),
    "eDNA percent positive"
  )
  distance_m <- .resolve_prediction_column(
    data, c("distance_m", "edna_distance_m"), "distance in metres"
  )
  lag_days <- .resolve_prediction_column(
    data, c("lag_days", "edna_lag_days"), "time since eDNA sampling in days"
  )
  bottom_time <- .resolve_prediction_column(
    data, c("bottom_time", "survey_minutes"), "survey effort", default_effort
  )
  if (any(!is.finite(percent_positive) | percent_positive < 0 | percent_positive > 100)) {
    stop("eDNA percent positive must be finite and between 0 and 100.", call. = FALSE)
  }
  if (any(!is.finite(distance_m) | distance_m < 0)) {
    stop("Distance must be finite and non-negative.", call. = FALSE)
  }
  if (any(!is.finite(lag_days) | lag_days < 0)) {
    stop("Time since sampling must be finite and non-negative.", call. = FALSE)
  }
  if (any(!is.finite(bottom_time) | bottom_time <= 0)) {
    stop("Survey effort must be finite and greater than zero.", call. = FALSE)
  }

  standardize <- function(value, variable) {
    (value - unname(scale_spec[[variable]]["center"])) /
      unname(scale_spec[[variable]]["scale"])
  }
  data$edna_prop_z <- standardize(percent_positive / 100, "edna_prop_model")
  data$distance_z <- standardize(log1p(distance_m / 200), "distance_model")
  data$lag_z <- standardize(log1p(lag_days), "lag_model")
  design_column <- intersect(c("sampling_design", "sampling_method"), names(data))
  if (length(design_column)) {
    design_key <- tolower(gsub("[^a-z0-9]", "", as.character(data[[design_column[[1]]]])))
    data$sampling_design <- dplyr::case_when(
      design_key %in% c("3x12", "3sitesx12reps") ~ "3x12",
      design_key %in% c("4x6", "4sitesx6reps") ~ "4x6",
      design_key %in% c("other", "mixed", "othermixed", "otherormixed", "othermixeddesign") ~
        "other_or_mixed",
      TRUE ~ NA_character_
    )
    if (anyNA(data$sampling_design)) {
      stop("Sampling design must be 3x12, 4x6, or other_or_mixed.", call. = FALSE)
    }
    data$design_4x6 <- as.integer(data$sampling_design == "4x6")
    data$design_other_or_mixed <- as.integer(data$sampling_design == "other_or_mixed")
  }
  data$bottom_time <- bottom_time
  data$log_effort <- log(bottom_time)
  data$Reef <- if ("Reef" %in% names(data)) {
    as.character(data$Reef)
  } else {
    paste0("new_reef_", seq_len(nrow(data)))
  }
  data$site_id <- if ("site_id" %in% names(data)) {
    as.character(data$site_id)
  } else {
    site_column <- intersect(c("cull_site_name", "site_name"), names(data))
    if (length(site_column)) {
      paste(data$Reef, as.character(data[[site_column[[1]]]]), sep = "__")
    } else {
      paste0("new_site_", seq_len(nrow(data)))
    }
  }
  data$edna_campaign_id <- if ("edna_campaign_id" %in% names(data)) {
    as.character(data$edna_campaign_id)
  } else {
    paste0("new_campaign_", seq_len(nrow(data)))
  }
  data
}

#' Predict COTS CPUE with the operational BRMS model
#'
#' Generates new-reef or known-reef posterior predictions, including expected
#' CPUE intervals, observation-level prediction intervals, and probabilities
#' that CPUE exceeds operational thresholds.
#'
#' @param fit A fitted `brmsfit` object.
#' @param newdata Data frame with percent positive (`perc_pos`), distance in
#'   metres (`distance_m`), elapsed days (`lag_days`), and optionally
#'   `bottom_time` and `Reef`.
#' @param scale_spec Predictor scaling specification stored with the model.
#' @param thresholds CPUE thresholds for exceedance probabilities.
#' @param ndraws Number of posterior draws, or `NULL` for all draws.
#' @param seed Random seed.
#' @param default_effort Survey effort used when it is absent from newdata.
#' @return A data frame of CPUE summaries and exceedance probabilities.
#' @export
predict_brms_cpue <- function(fit, newdata, scale_spec,
                              thresholds = c(0.02, 0.04, 0.08),
                              ndraws = 2000, seed = 1,
                              default_effort = 216) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for BRMS predictions.", call. = FALSE)
  }
  prepared <- .prepare_operational_predictors(
    newdata, scale_spec, default_effort = default_effort
  )
  set.seed(seed)
  expected_counts <- brms::posterior_epred(
    fit, newdata = prepared, re_formula = NULL,
    allow_new_levels = TRUE, sample_new_levels = "gaussian", ndraws = ndraws
  )
  predicted_counts <- brms::posterior_predict(
    fit, newdata = prepared, re_formula = NULL,
    allow_new_levels = TRUE, sample_new_levels = "gaussian", ndraws = ndraws
  )
  expected_cpue <- sweep(expected_counts, 2, prepared$bottom_time, "/")
  predicted_cpue <- sweep(predicted_counts, 2, prepared$bottom_time, "/")
  output <- data.frame(
    row_id = seq_len(nrow(prepared)),
    cpue_estimate = apply(expected_cpue, 2, mean),
    cpue_lower = apply(expected_cpue, 2, stats::quantile, probs = 0.025),
    cpue_upper = apply(expected_cpue, 2, stats::quantile, probs = 0.975),
    prediction_lower = apply(predicted_cpue, 2, stats::quantile, probs = 0.025),
    prediction_upper = apply(predicted_cpue, 2, stats::quantile, probs = 0.975)
  )
  for (threshold in thresholds) {
    suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
    output[[paste0("prob_exceeds_", suffix)]] <- colMeans(predicted_cpue >= threshold)
  }
  output
}

#' Predict COTS CPUE with the operational boosted-tree model
#'
#' @param model A `reefDNA_brt_model` object produced by the analysis pipeline.
#' @inheritParams predict_brms_cpue
#' @return A data frame of CPUE summaries and threshold probabilities.
#' @export
predict_brt_cpue <- function(model, newdata,
                             thresholds = c(0.02, 0.04, 0.08)) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for BRT predictions.", call. = FALSE)
  }
  design_predictors <- c("design_4x6", "design_other_or_mixed")
  if (any(design_predictors %in% model$predictors) &&
      !any(c("sampling_design", "sampling_method") %in% names(newdata))) {
    stop(
      "This BRT requires sampling_design: 3x12, 4x6, or other_or_mixed.",
      call. = FALSE
    )
  }
  default_effort <- if (is.null(model$default_effort)) 216 else model$default_effort
  prepared <- .prepare_operational_predictors(
    newdata, model$scale_spec, default_effort = default_effort
  )
  matrix <- xgboost::xgb.DMatrix(
    as.matrix(prepared[model$predictors]), base_margin = prepared$log_effort
  )
  cpue <- pmax(
    predict(model$regression, matrix) / prepared$bottom_time, 0
  )
  residual_interval <- model$oof_residual_quantiles
  output <- data.frame(
    row_id = seq_len(nrow(prepared)),
    cpue_estimate = cpue,
    cpue_lower = cpue,
    cpue_upper = cpue,
    prediction_lower = pmax(cpue + residual_interval[[1]], 0),
    prediction_upper = pmax(cpue + residual_interval[[2]], 0)
  )
  predictor_matrix <- xgboost::xgb.DMatrix(as.matrix(prepared[model$predictors]))
  for (threshold in thresholds) {
    suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
    classifier_name <- paste0("above_", suffix)
    output[[paste0("prob_exceeds_", suffix)]] <- predict(
      model$classifiers[[classifier_name]], predictor_matrix
    )
  }
  output
}

#' Predict COTS CPUE with the reefDNA operational ensemble
#'
#' The ensemble combines the concentration-free BRMS and boosted-tree models
#' using weights estimated from whole-reef cross-validation. Returned alerts
#' use threshold-specific cutoffs selected by cross-validated F1 score.
#'
#' @param bundle A model bundle containing `brms_fit`, `brt_model`,
#'   `ensemble_spec`, and `scale_spec`.
#' @param method One of `"ensemble"`, `"brms"`, or `"brt"`.
#' @inheritParams predict_brms_cpue
#' @return A data frame with expected CPUE, intervals, threshold probabilities,
#'   and operational alert flags.
#' @export
predict_reef_cpue <- function(bundle, newdata,
                              method = c("ensemble", "brms", "brt"),
                              thresholds = c(0.02, 0.04, 0.08),
                              ndraws = 2000, seed = 1) {
  method <- match.arg(method)
  brms_prediction <- NULL
  brt_prediction <- NULL
  if (method %in% c("ensemble", "brms")) {
    default_effort <- if (is.null(bundle$metadata$default_effort)) {
      216
    } else {
      bundle$metadata$default_effort
    }
    brms_prediction <- predict_brms_cpue(
      bundle$brms_fit, newdata, bundle$scale_spec,
      thresholds = thresholds, ndraws = ndraws, seed = seed,
      default_effort = default_effort
    )
  }
  if (method %in% c("ensemble", "brt")) {
    brt_prediction <- predict_brt_cpue(
      bundle$brt_model, newdata, thresholds = thresholds
    )
  }
  if (method == "brms") return(brms_prediction)
  if (method == "brt") return(brt_prediction)

  brms_weight <- bundle$ensemble_spec$cpue_brms_weight
  brt_weight <- bundle$ensemble_spec$cpue_brt_weight
  output <- data.frame(
    row_id = brms_prediction$row_id,
    cpue_estimate = brms_weight * brms_prediction$cpue_estimate +
      brt_weight * brt_prediction$cpue_estimate,
    cpue_lower = brms_weight * brms_prediction$cpue_lower +
      brt_weight * brt_prediction$cpue_estimate,
    cpue_upper = brms_weight * brms_prediction$cpue_upper +
      brt_weight * brt_prediction$cpue_estimate,
    prediction_lower = brms_weight * brms_prediction$prediction_lower +
      brt_weight * brt_prediction$prediction_lower,
    prediction_upper = brms_weight * brms_prediction$prediction_upper +
      brt_weight * brt_prediction$prediction_upper
  )
  for (threshold in thresholds) {
    suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
    weight_row <- bundle$ensemble_spec$threshold_weights[
      abs(bundle$ensemble_spec$threshold_weights$threshold - threshold) < 1e-8,
      , drop = FALSE
    ]
    if (!nrow(weight_row)) {
      stop("No ensemble weight is stored for threshold ", threshold, ".", call. = FALSE)
    }
    probability_name <- paste0("prob_exceeds_", suffix)
    output[[probability_name]] <-
      weight_row$brms_weight * brms_prediction[[probability_name]] +
      weight_row$brt_weight * brt_prediction[[probability_name]]
    cutoff_row <- bundle$ensemble_spec$threshold_cutoffs[
      abs(bundle$ensemble_spec$threshold_cutoffs$threshold - threshold) < 1e-8,
      , drop = FALSE
    ]
    output[[paste0("alert_", suffix)]] <-
      output[[probability_name]] >= cutoff_row$cutoff
  }
  attr(output, "interval_note") <- paste(
    "Ensemble expected intervals combine BRMS parameter uncertainty with the",
    "BRT point prediction; prediction intervals combine BRMS posterior",
    "prediction bounds with empirical BRT out-of-fold residual bounds."
  )
  output
}
