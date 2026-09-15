#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/site_visit_brms_oof.csv")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
plot_dir <- file.path(project_root, "plots")
brms <- read.csv(file.path(output_dir, "site_visit_brms_oof.csv")) |>
  transmute(
    site_visit_id, observed_cpue, brms_cpue = predicted_cpue,
    brms_002 = prob_002, brms_004 = prob_004, brms_008 = prob_008
  )
brt <- read.csv(file.path(output_dir, "site_visit_brt_oof.csv")) |>
  transmute(site_visit_id, brt_cpue = predicted_cpue)
brt_thresholds <- read.csv(file.path(output_dir, "site_visit_brt_threshold_oof.csv")) |>
  mutate(suffix = gsub("\\.", "", sprintf("%.2f", threshold))) |>
  select(site_visit_id, suffix, probability) |>
  pivot_wider(names_from = suffix, values_from = probability, names_prefix = "brt_")
oof <- brms |>
  inner_join(brt, by = "site_visit_id") |>
  inner_join(brt_thresholds, by = "site_visit_id")
stopifnot(nrow(oof) == 519L, !anyDuplicated(oof$site_visit_id))

weight_search <- tibble(brms_weight = seq(0, 1, by = 0.01)) |>
  mutate(
    brt_weight = 1 - brms_weight,
    predicted = lapply(
      brms_weight,
      function(weight) weight * oof$brms_cpue + (1 - weight) * oof$brt_cpue
    ),
    RMSE = vapply(predicted, function(x) sqrt(mean((oof$observed_cpue - x)^2)), numeric(1)),
    MAE = vapply(predicted, function(x) mean(abs(oof$observed_cpue - x)), numeric(1))
  ) |>
  select(-predicted)
best_cpue_weight <- slice_min(weight_search, RMSE, n = 1, with_ties = FALSE)
oof$ensemble_cpue <- best_cpue_weight$brms_weight * oof$brms_cpue +
  best_cpue_weight$brt_weight * oof$brt_cpue

regression_metrics <- bind_rows(lapply(
  c(BRMS = "brms_cpue", BRT = "brt_cpue", Ensemble = "ensemble_cpue"),
  function(column) {
    estimate <- oof[[column]]
    tibble(
      RMSE = sqrt(mean((oof$observed_cpue - estimate)^2)),
      MAE = mean(abs(oof$observed_cpue - estimate)),
      correlation = suppressWarnings(cor(oof$observed_cpue, estimate)),
      bias = mean(estimate - oof$observed_cpue)
    )
  }
), .id = "model") |>
  arrange(RMSE)

auc <- function(observed, probability) {
  positives <- sum(observed == 1)
  negatives <- sum(observed == 0)
  if (!positives || !negatives) return(NA_real_)
  (sum(rank(probability)[observed == 1]) - positives * (positives + 1) / 2) /
    (positives * negatives)
}
classification_metrics <- function(observed, probability, cutoff) {
  predicted <- probability >= cutoff
  tp <- sum(predicted & observed == 1)
  fp <- sum(predicted & observed == 0)
  fn <- sum(!predicted & observed == 1)
  tn <- sum(!predicted & observed == 0)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  tibble(
    cutoff = cutoff,
    F1 = if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0,
    recall = recall, precision = precision,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else 0,
    false_negatives = fn
  )
}

threshold_weights <- list()
threshold_results <- list()
calibration <- list()
for (threshold in c(0.02, 0.04, 0.08)) {
  suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
  observed <- as.integer(oof$observed_cpue >= threshold)
  brms_probability <- oof[[paste0("brms_", suffix)]]
  brt_probability <- oof[[paste0("brt_", suffix)]]
  probability_weights <- tibble(brms_weight = seq(0, 1, by = 0.01)) |>
    mutate(
      brt_weight = 1 - brms_weight,
      Brier = vapply(
        brms_weight,
        function(weight) mean((observed - (weight * brms_probability + (1 - weight) * brt_probability))^2),
        numeric(1)
      )
    )
  best_probability_weight <- slice_min(
    probability_weights, Brier, n = 1, with_ties = FALSE
  ) |>
    mutate(threshold = threshold, .before = 1)
  threshold_weights[[suffix]] <- best_probability_weight
  ensemble_probability <- best_probability_weight$brms_weight * brms_probability +
    best_probability_weight$brt_weight * brt_probability
  probabilities <- list(
    BRMS = brms_probability, BRT = brt_probability,
    Ensemble = ensemble_probability
  )
  threshold_results[[suffix]] <- bind_rows(lapply(names(probabilities), function(model) {
    probability <- probabilities[[model]]
    cutoff_search <- bind_rows(lapply(
      seq(0, 1, by = 0.01),
      function(cutoff) classification_metrics(observed, probability, cutoff)
    ))
    best_cutoff <- slice_max(cutoff_search, F1, n = 1, with_ties = FALSE)
    best_cutoff |>
      mutate(
        threshold = threshold, model = model,
        AUC = auc(observed, probability),
        Brier = mean((observed - probability)^2),
        .before = 1
      )
  }))
  calibration[[suffix]] <- tibble(
    threshold = factor(threshold), observed = observed,
    probability = ensemble_probability
  ) |>
    mutate(bin = ntile(probability, 8)) |>
    group_by(threshold, bin) |>
    summarise(
      predicted = mean(probability), observed = mean(observed), n = n(),
      .groups = "drop"
    )
  oof[[paste0("ensemble_", suffix)]] <- ensemble_probability
}
threshold_weights <- bind_rows(threshold_weights)
threshold_results <- bind_rows(threshold_results)
calibration <- bind_rows(calibration)
write.csv(weight_search, file.path(output_dir, "site_visit_cpue_weight_search.csv"), row.names = FALSE)
write.csv(regression_metrics, file.path(output_dir, "site_visit_model_performance.csv"), row.names = FALSE)
write.csv(threshold_weights, file.path(output_dir, "site_visit_threshold_weights.csv"), row.names = FALSE)
write.csv(threshold_results, file.path(output_dir, "site_visit_threshold_performance.csv"), row.names = FALSE)
write.csv(oof, file.path(output_dir, "site_visit_ensemble_oof.csv"), row.names = FALSE)

ensemble_cutoffs <- threshold_results |>
  filter(model == "Ensemble") |>
  select(threshold, cutoff)
ensemble_spec <- list(
  cpue_brms_weight = best_cpue_weight$brms_weight,
  cpue_brt_weight = best_cpue_weight$brt_weight,
  threshold_weights = threshold_weights |>
    select(threshold, brms_weight, brt_weight, brier = Brier),
  threshold_cutoffs = ensemble_cutoffs,
  validation = "five-fold complete eDNA-campaign holdout"
)
saveRDS(ensemble_spec, file.path(output_dir, "site_visit_ensemble_spec.rds"))

blue <- "#2166AC"
orange <- "#D95F02"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(), legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5)
  )
limits <- range(c(oof$observed_cpue, oof$ensemble_cpue))
panel_a <- ggplot(oof, aes(observed_cpue, ensemble_cpue)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(colour = blue, alpha = 0.4, size = 1.2) +
  coord_equal(xlim = limits, ylim = limits) +
  labs(
    title = "Ensemble fit", subtitle = "Complete campaigns held out",
    x = "Observed CPUE (min^-1)", y = "Predicted CPUE (min^-1)"
  ) + paper_theme
performance_long <- regression_metrics |>
  select(model, RMSE, MAE) |>
  pivot_longer(c(RMSE, MAE), names_to = "metric", values_to = "value")
panel_b <- ggplot(performance_long, aes(model, value, fill = model)) +
  geom_col(width = 0.65) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(BRMS = blue, BRT = orange, Ensemble = "#1B9E77"), guide = "none") +
  labs(title = "Prediction error", x = NULL, y = NULL) + paper_theme
panel_c <- ggplot(
  threshold_results,
  aes(factor(threshold), AUC, colour = model, group = model)
) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c(BRMS = blue, BRT = orange, Ensemble = "#1B9E77")) +
  labs(title = "Threshold discrimination", x = "CPUE threshold", y = "AUC", colour = NULL) +
  paper_theme
panel_d <- ggplot(calibration, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_line(colour = blue, linewidth = 0.7) +
  geom_point(aes(size = n), colour = blue) +
  facet_wrap(~threshold, labeller = label_both) +
  scale_size(range = c(1.5, 4), guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Ensemble calibration", x = "Predicted probability", y = "Observed frequency") +
  paper_theme
figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "site_visit_ensemble_assessment.png"), figure,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "site_visit_ensemble_assessment.pdf"), figure,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)

bundle <- list(
  brms_fit = readRDS(file.path(output_dir, "site_visit_brms_final.rds")),
  brt_model = readRDS(file.path(output_dir, "site_visit_brt_model.rds")),
  ensemble_spec = ensemble_spec,
  scale_spec = readRDS(file.path(output_dir, "site_visit_scale_spec.rds")),
  metadata = list(
    version = "2.0.0-site-visit",
    response = "site-visit COTS CPUE per minute",
    predictors = c("perc_pos", "distance_m", "lag_days"),
    linkage = paste(
      "eDNA campaigns within seven days; nearest eDNA site; first subsequent",
      "visit at each cull site; unique response rows"
    ),
    grouping = c("Reef", "site_id", "edna_campaign_id"),
    validation = ensemble_spec$validation,
    default_effort = median(readRDS(file.path(output_dir, "site_visit_model_data.rds"))$bottom_time),
    created = Sys.time()
  )
)
saveRDS(
  bundle,
  file.path(output_dir, "reefDNA_site_visit_model_bundle.rds"), compress = "xz"
)
print(regression_metrics)
print(threshold_weights)
print(threshold_results)
