# Test 4-panel script with 5-fold CV (R=100) and Isotonic Monotonic Trajectory
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(fuzzyjoin)
  library(broom)
  library(MASS)
  library(glmmTMB)
  library(tidyr)
  library(patchwork)
  library(rsample)
  library(purrr)
  library(tibble)
  library(viridis)
  library(scales)
})

source("eDNA_CPUE_Comparison_final.R", local = FALSE)

perc_grid <- seq(0, 100, by = 2)
cpue_grid <- seq(0, quantile(dat_evt$cpue, 0.99, na.rm = TRUE), length.out = 40)

f_beta <- function(prec, rec, beta = 1) {
  if (is.na(prec) || is.na(rec) || (beta^2 * prec + rec == 0)) return(0)
  (1 + beta^2) * (prec * rec) / ((beta^2 * prec) + rec)
}

grid <- tidyr::expand_grid(perc_thresh = perc_grid, cpue_thresh = cpue_grid) %>%
  mutate(
    out = purrr::pmap(., ~ metrics_for(..1, ..2, dat_evt)),
    F1 = vapply(out, `[[`, numeric(1), "F1"),
    Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
    Precision = vapply(out, `[[`, numeric(1), "Precision"),
    Recall = vapply(out, `[[`, numeric(1), "Recall"),
    F2 = map2_dbl(Precision, Recall, ~ f_beta(.x, .y, beta = 2))
  ) %>%
  dplyr::select(-out)

cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
k <- 5
R <- 100
set.seed(42)

run_5fold_tuning <- function(dat_evt, cpue_thr, perc_grid, k = 5, R = 100) {
  dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
  strat_ok <- length(unique(dat2$actual)) > 1
  folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

  cv_long <- folds %>%
    mutate(assess = map(splits, assessment)) %>%
    dplyr::select(id, id2, assess) %>%
    tidyr::expand_grid(perc_thresh = perc_grid) %>%
    mutate(
      out = pmap(list(perc_thresh, assess), \(p, d) metrics_for(p, cpue_thr, d)),
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
      Rec_mean = mean(Recall, na.rm = TRUE),
      .groups = "drop"
    )

  best_f1 <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
  best_f2 <- cv_sum %>% slice_max(F2_mean, n = 1, with_ties = FALSE)

  tibble(
    cpue_thr = cpue_thr,
    perc_star_f1 = best_f1$perc_thresh,
    perc_star_f2 = best_f2$perc_thresh
  )
}

tuning_5fold <- map_dfr(cpue_levels, ~ run_5fold_tuning(dat_evt, .x, perc_grid, k = k, R = R)) %>%
  mutate(
    perc_star_f1_iso = isoreg(cpue_thr, perc_star_f1)$yf,
    perc_star_f2_iso = isoreg(cpue_thr, perc_star_f2)$yf
  )

print(tuning_5fold)
