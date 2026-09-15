suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rsample)
  library(readxl)
  library(readr)
})

source("eDNA_CPUE_Comparison_final.R", local = FALSE)

f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (beta^2 * prec + rec == 0)) return(0)
  (1 + beta^2) * (prec * rec) / ((beta^2 * prec) + rec)
}

cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
perc_grid <- seq(0, 100, by = 2)
k <- 5
R <- 100
set.seed(42)

# Step 1: Tune raw & isotonic monotonic thresholds
tune_res <- map_dfr(cpue_levels, function(c_thr) {
  dat2 <- dat_evt %>% mutate(actual = if (c_thr == 0) (cpue > 0) else (cpue >= c_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_long <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    dplyr::select(id, id2, assess) %>%
    tidyr::expand_grid(perc_thresh = perc_grid) %>%
    mutate(
      out = pmap(list(perc_thresh, assess), \(p, d) metrics_for(p, c_thr, d)),
      F1 = vapply(out, `[[`, numeric(1), "F1"),
      Precision = vapply(out, `[[`, numeric(1), "Precision"),
      Recall = vapply(out, `[[`, numeric(1), "Recall")
    ) %>%
    mutate(F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2))) %>%
    dplyr::select(-out, -assess)

  cv_sum <- cv_long %>%
    group_by(perc_thresh) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F2_mean = mean(F2, na.rm = TRUE),
      .groups = "drop"
    )

  best_f1 <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
  best_f2 <- cv_sum %>% slice_max(F2_mean, n = 1, with_ties = FALSE)

  tibble(
    cpue_thr = c_thr,
    p_f1_raw = best_f1$perc_thresh,
    p_f2_raw = best_f2$perc_thresh
  )
}) %>%
  mutate(
    p_f1_iso = isoreg(cpue_thr, p_f1_raw)$yf,
    p_f2_iso = isoreg(cpue_thr, p_f2_raw)$yf
  )

# Step 2: Evaluate CV performance metrics at the tuned monotonic thresholds
eval_targets <- c(0.02, 0.04, 0.08)

get_cv_stats <- function(c_thr, strategy, p_thresh) {
  dat2 <- dat_evt %>% mutate(actual = if (c_thr == 0) (cpue > 0) else (cpue >= c_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_eval <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    mutate(
      out = map(assess, ~ metrics_for(p_thresh, c_thr, .x)),
      F1 = vapply(out, `[[`, numeric(1), "F1"),
      Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
      Precision = vapply(out, `[[`, numeric(1), "Precision"),
      Recall = vapply(out, `[[`, numeric(1), "Recall")
    ) %>%
    mutate(F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2)))

  se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

  fmt <- function(m, s) sprintf("%.3f ± %.3f", m, s)

  tibble(
    `CPUE Target` = c_thr,
    `Optimization Strategy` = strategy,
    `% eDNA Positive Threshold` = sprintf("%d%%", round(p_thresh)),
    `F1 Score` = fmt(mean(cv_eval$F1, na.rm = TRUE), se(cv_eval$F1)),
    `F2 Score` = fmt(mean(cv_eval$F2, na.rm = TRUE), se(cv_eval$F2)),
    `Accuracy` = fmt(mean(cv_eval$Accuracy, na.rm = TRUE), se(cv_eval$Accuracy)),
    `Precision` = fmt(mean(cv_eval$Precision, na.rm = TRUE), se(cv_eval$Precision)),
    `Recall` = fmt(mean(cv_eval$Recall, na.rm = TRUE), se(cv_eval$Recall))
  )
}

tbl_rows <- map_dfr(eval_targets, function(ct) {
  row_f1 <- tune_res %>% filter(cpue_thr == ct)
  p1 <- row_f1$p_f1_iso
  p2 <- row_f1$p_f2_iso

  bind_rows(
    get_cv_stats(ct, "F1-Optimized (Balanced)", p1),
    get_cv_stats(ct, "F2-Optimized (Recall-Weighted)", p2)
  )
})

print(tbl_rows)
write_csv(tbl_rows, "plots/cv_performance_summary_table.csv")
