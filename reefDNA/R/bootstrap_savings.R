#' Bootstrap Confidence Intervals for Reef Savings
#'
#' Performs case-resampling bootstrap on reef-year observations to produce
#' confidence intervals around per-reef and program-level savings estimates.
#' Each bootstrap replicate resamples reef-years with replacement, applies
#' the supplied classification function, and computes savings metrics.
#'
#' @param data Data frame with one row per reef-year.  Must contain at minimum
#'   the columns that \code{classify_fn} requires plus \code{Total_Diver_Hours}.
#' @param classify_fn A function that accepts a data frame (same schema as
#'   \code{data}) and returns it with three additional logical/numeric columns:
#'   \describe{
#'     \item{is_true_negative}{Logical — reef was correctly closed.}
#'     \item{net_saved_hrs}{Numeric — diver hours that would have been saved
#'       for this reef-year (0 for non-TN observations).}
#'     \item{is_false_negative}{Logical — reef was a missed outbreak.}
#'   }
#' @param n_boot Number of bootstrap replicates (default 10 000).
#' @param ci_level Confidence level for the percentile interval (default 0.80).
#' @param cost_per_diver_hour Dollar cost per diver hour (default 1000).
#' @param edna_cost_per_reef eDNA testing cost per reef (default 2500).
#' @param n_operational_reefs Number of reefs in the operational program to
#'   scale estimates to (default 40).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{per_reef}{Tibble — bootstrap summary statistics for per-reef metrics
#'     (hours saved, gross dollar savings, net dollar savings, false negative
#'     rate).}
#'   \item{scaled}{Tibble — same metrics scaled to \code{n_operational_reefs}.}
#'   \item{boot_dist}{Tibble — full bootstrap distribution with one row per
#'     replicate, useful for plotting.}
#' }
#'
#' @export
#' @import dplyr
bootstrap_savings <- function(data,
                              classify_fn,
                              n_boot = 10000,
                              ci_level = 0.80,
                              cost_per_diver_hour = 1000,
                              edna_cost_per_reef = 2500,
                              n_operational_reefs = 40,
                              seed = 42) {

  set.seed(seed)

  n <- nrow(data)
  alpha <- (1 - ci_level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)

  # Pre-compute: classify the original data once for structure validation
  test_out <- classify_fn(data[1, , drop = FALSE])
  required_cols <- c("is_true_negative", "net_saved_hrs", "is_false_negative")
  missing <- setdiff(required_cols, names(test_out))
  if (length(missing) > 0) {
    stop("classify_fn must return columns: ", paste(missing, collapse = ", "))
  }

  # Bootstrap replicates
  boot_results <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_data <- data[idx, , drop = FALSE]
    classified <- classify_fn(boot_data)

    tn_count <- sum(classified$is_true_negative, na.rm = TRUE)
    fn_count <- sum(classified$is_false_negative, na.rm = TRUE)
    total_saved_hrs <- sum(classified$net_saved_hrs, na.rm = TRUE)

    # Per-reef averages (averaged across all n resampled reefs)
    avg_hrs_saved <- total_saved_hrs / n
    avg_fn_rate <- fn_count / n

    boot_results[[b]] <- tibble::tibble(
      replicate = b,
      avg_hrs_saved_per_reef = avg_hrs_saved,
      avg_fn_rate_per_reef = avg_fn_rate,
      tn_count = tn_count,
      fn_count = fn_count,
      total_saved_hrs = total_saved_hrs
    )
  }

  boot_dist <- dplyr::bind_rows(boot_results)

  # Compute per-reef dollar metrics
  boot_dist <- boot_dist %>%
    mutate(
      gross_dollars_per_reef = avg_hrs_saved_per_reef * cost_per_diver_hour,
      net_dollars_per_reef = gross_dollars_per_reef - edna_cost_per_reef,
      scaled_hrs_saved = avg_hrs_saved_per_reef * n_operational_reefs,
      scaled_gross_savings = gross_dollars_per_reef * n_operational_reefs,
      scaled_edna_cost = edna_cost_per_reef * n_operational_reefs,
      scaled_net_savings = scaled_gross_savings - scaled_edna_cost,
      scaled_fn_count = avg_fn_rate_per_reef * n_operational_reefs
    )

  # Summarise per-reef metrics
  summarise_metric <- function(x, label) {
    qs <- quantile(x, probs = probs, na.rm = TRUE)
    tibble::tibble(
      metric = label,
      mean = mean(x, na.rm = TRUE),
      median = qs[2],
      ci_lower = qs[1],
      ci_upper = qs[3],
      sd = sd(x, na.rm = TRUE)
    )
  }

  per_reef <- dplyr::bind_rows(
    summarise_metric(boot_dist$avg_hrs_saved_per_reef, "Hours saved per reef"),
    summarise_metric(boot_dist$gross_dollars_per_reef, "Gross dollar savings per reef"),
    summarise_metric(boot_dist$net_dollars_per_reef, "Net dollar savings per reef (after eDNA cost)"),
    summarise_metric(boot_dist$avg_fn_rate_per_reef, "Missed outbreaks per reef")
  )

  scaled <- dplyr::bind_rows(
    summarise_metric(boot_dist$scaled_hrs_saved, "Annual hours saved"),
    summarise_metric(boot_dist$scaled_gross_savings, "Annual gross dollar savings"),
    summarise_metric(boot_dist$scaled_net_savings, "Annual net savings (after eDNA costs)"),
    summarise_metric(boot_dist$scaled_fn_count, "Missed outbreaks per year")
  )

  # Break-even probability
  breakeven_prob <- mean(boot_dist$scaled_net_savings > 0, na.rm = TRUE)

  list(
    per_reef = per_reef,
    scaled = scaled,
    boot_dist = boot_dist,
    breakeven_prob = breakeven_prob,
    ci_level = ci_level,
    n_boot = n_boot,
    n_obs = n,
    params = list(
      cost_per_diver_hour = cost_per_diver_hour,
      edna_cost_per_reef = edna_cost_per_reef,
      n_operational_reefs = n_operational_reefs
    )
  )
}


#' Sensitivity Analysis for Cost Parameters
#'
#' Given a bootstrap result, recalculates the median savings under different
#' cost assumptions.  This is a deterministic recalculation — the bootstrap
#' distribution of hours-saved is held fixed and only the dollar conversion
#' parameters vary.
#'
#' @param boot_result Output from \code{\link{bootstrap_savings}}.
#' @param cost_per_hour_range Numeric vector of cost-per-diver-hour values to
#'   evaluate (default: seq(600, 1400, 200)).
#' @param edna_program_cost_range Numeric vector of total eDNA program costs
#'   (default: seq(200000, 600000, 100000)).
#' @param n_reefs_range Numeric vector of operational reef counts
#'   (default: seq(30, 50, 5)).
#'
#' @return A tibble with one row per parameter combination, containing the
#'   median net annual savings for that scenario.
#'
#' @export
#' @import dplyr
sensitivity_analysis <- function(boot_result,
                                 cost_per_hour_range = seq(600, 1400, 200),
                                 edna_program_cost_range = seq(200000, 600000, 100000),
                                 n_reefs_range = seq(30, 50, 5)) {

  # Extract the median hours saved per reef from bootstrap
  median_hrs_per_reef <- median(boot_result$boot_dist$avg_hrs_saved_per_reef, na.rm = TRUE)
  base_params <- boot_result$params

  results <- list()

  # Vary cost per hour (hold others at base)
  for (cph in cost_per_hour_range) {
    gross_per_reef <- median_hrs_per_reef * cph
    edna_per_reef <- base_params$edna_cost_per_reef
    n_reefs <- base_params$n_operational_reefs
    net_annual <- (gross_per_reef - edna_per_reef) * n_reefs
    results <- c(results, list(tibble::tibble(
      parameter = "Cost per diver hour",
      value = cph,
      value_label = paste0("$", format(cph, big.mark = ",")),
      net_annual_savings = net_annual,
      is_base = (cph == base_params$cost_per_diver_hour)
    )))
  }

  # Vary eDNA program cost (hold others at base)
  for (epc in edna_program_cost_range) {
    gross_per_reef <- median_hrs_per_reef * base_params$cost_per_diver_hour
    edna_per_reef <- epc / base_params$n_operational_reefs
    n_reefs <- base_params$n_operational_reefs
    net_annual <- (gross_per_reef * n_reefs) - epc
    results <- c(results, list(tibble::tibble(
      parameter = "eDNA program cost",
      value = epc,
      value_label = paste0("$", format(epc, big.mark = ",")),
      net_annual_savings = net_annual,
      is_base = (epc == (base_params$edna_cost_per_reef * base_params$n_operational_reefs))
    )))
  }

  # Vary number of reefs (hold others at base)
  for (nr in n_reefs_range) {
    gross_per_reef <- median_hrs_per_reef * base_params$cost_per_diver_hour
    edna_per_reef <- base_params$edna_cost_per_reef
    net_annual <- (gross_per_reef - edna_per_reef) * nr
    results <- c(results, list(tibble::tibble(
      parameter = "Operational reefs",
      value = nr,
      value_label = as.character(nr),
      net_annual_savings = net_annual,
      is_base = (nr == base_params$n_operational_reefs)
    )))
  }

  dplyr::bind_rows(results)
}
