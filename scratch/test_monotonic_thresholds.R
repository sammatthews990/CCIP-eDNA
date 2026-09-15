suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rsample)
  library(readxl)
  library(glmmTMB)
})

source("eDNA_CPUE_Comparison_final.R", local = FALSE)

cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
perc_grid <- seq(0, 100, by = 1)

# Compare 5-fold vs 10-fold stability across R=100 repeats
f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (beta^2 * prec + rec == 0)) return(0)
  (1 + beta^2) * (prec * rec) / ((beta^2 * prec) + rec)
}

run_cv_grid <- function(dat_evt, k_folds = 5, R_repeats = 100) {
  set.seed(42)
  res_list <- map(cpue_levels, function(c_thr) {
    dat2 <- dat_evt %>% mutate(actual = if (c_thr == 0) (cpue > 0) else (cpue >= c_thr))
    strat_ok <- length(unique(dat2$actual)) > 1
    folds <- if (strat_ok) vfold_cv(dat2, v = k_folds, repeats = R_repeats, strata = actual) else vfold_cv(dat2, v = k_folds, repeats = R_repeats)

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

    p_f1_raw <- (cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE))$perc_thresh
    p_f2_raw <- (cv_sum %>% slice_max(F2_mean, n = 1, with_ties = FALSE))$perc_thresh

    tibble(cpue_thr = c_thr, p_f1_raw = p_f1_raw, p_f2_raw = p_f2_raw)
  })

  df_raw <- bind_rows(res_list)

  # Apply Monotonic Smoothing (Isotonic Regression)
  df_raw %>%
    mutate(
      p_f1_iso = isoreg(cpue_thr, p_f1_raw)$yf,
      p_f2_iso = isoreg(cpue_thr, p_f2_raw)$yf
    )
}

cat("--- 5-FOLD CV RESULT ---\n")
res_5f <- run_cv_grid(dat_evt, k_folds = 5, R_repeats = 100)
print(res_5f, n = 20)

cat("\n--- 10-FOLD CV RESULT ---\n")
res_10f <- run_cv_grid(dat_evt, k_folds = 10, R_repeats = 100)
print(res_10f, n = 20)
