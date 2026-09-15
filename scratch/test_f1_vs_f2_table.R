suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rsample)
  library(glmmTMB)
  library(pROC)
})

# Load dat_glmm and helper functions directly
source("eDNA_CPUE_Comparison_final.R", local = FALSE)

cpue_levels <- c(0.02, 0.04, 0.08)
perc_grid <- seq(0, 100, by = 1)
k <- 5
R <- 100

f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (prec + rec == 0)) return(0)
  (1 + beta^2) * (prec * rec) / ((beta^2 * prec) + rec)
}

run_cv_comparison <- function(dat_evt, cpue_thr, perc_grid, k = 5, R = 100) {
  dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  set.seed(42)
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_long <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    dplyr::select(id, id2, assess) %>%
    tidyr::expand_grid(perc_thresh = perc_grid) %>%
    mutate(
      out = pmap(list(perc_thresh, assess), \(p, d) metrics_for(p, cpue_thr, d)),
      F1 = vapply(out, `[[`, numeric(1), "F1"),
      Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
      Precision = vapply(out, `[[`, numeric(1), "Precision"),
      Recall = vapply(out, `[[`, numeric(1), "Recall")
    ) %>%
    mutate(
      F2 = ifelse(is.na(Precision) | is.na(Recall) | (4 * Precision + Recall == 0), 0, 5 * (Precision * Recall) / (4 * Precision + Recall))
    ) %>%
    dplyr::select(-out, -assess)

  cv_sum <- cv_long %>%
    group_by(perc_thresh) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F2_mean = mean(F2, na.rm = TRUE),
      Rec_mean = mean(Recall, na.rm = TRUE),
      Prec_mean = mean(Precision, na.rm = TRUE),
      .groups = "drop"
    )

  # 1. Best F1 threshold
  best_f1 <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
  p_f1 <- best_f1$perc_thresh

  # 2. Best F2 threshold (Recall-weighted)
  best_f2 <- cv_sum %>% slice_max(F2_mean, n = 1, with_ties = FALSE)
  p_f2 <- best_f2$perc_thresh

  # Evaluate performance at both thresholds across all folds
  eval_thresh <- function(p_star, opt_label) {
    perf_split <- folds %>%
      mutate(
        assess = map(splits, assessment),
        m = map(assess, ~ {
          out <- metrics_for(p_star, cpue_thr, .x)
          prec <- out["Precision"]
          rec <- out["Recall"]
          f2_val <- if (is.na(prec) || is.na(rec) || (4 * prec + rec == 0)) 0 else (5 * prec * rec / (4 * prec + rec))
          tibble(
            F1 = out["F1"],
            Accuracy = out["Accuracy"],
            Precision = prec,
            Recall = rec,
            F2 = f2_val
          )
        })
      ) %>%
      dplyr::select(m) %>%
      unnest(m)

    perf_split %>%
      summarise(
        cpue_thr = cpue_thr,
        opt_criterion = opt_label,
        perc_star = p_star,
        F1_mean = mean(F1, na.rm = TRUE), F1_se = sd(F1, na.rm = TRUE) / sqrt(n()),
        F2_mean = mean(F2, na.rm = TRUE), F2_se = sd(F2, na.rm = TRUE) / sqrt(n()),
        Acc_mean = mean(Accuracy, na.rm = TRUE), Acc_se = sd(Accuracy, na.rm = TRUE) / sqrt(n()),
        Prec_mean = mean(Precision, na.rm = TRUE), Prec_se = sd(Precision, na.rm = TRUE) / sqrt(n()),
        Rec_mean = mean(Recall, na.rm = TRUE), Rec_se = sd(Recall, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
  }

  res_f1 <- eval_thresh(p_f1, "F1-Optimized (Equal Weight)")
  res_f2 <- eval_thresh(p_f2, "F2-Optimized (Recall-Weighted)")

  bind_rows(res_f1, res_f2)
}

results_cmp <- map_dfr(cpue_levels, ~ run_cv_comparison(dat_evt, .x, perc_grid, k = k, R = R))
print(results_cmp, width = Inf)

write.csv(results_cmp, "scratch/results_cmp_f1_f2.csv", row.names = FALSE)
