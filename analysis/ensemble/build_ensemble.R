#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(patchwork)
  library(pROC)
  library(tidyr)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
brms_dir <- file.path(project_root, "analysis", "brms", "output")
brt_dir <- file.path(project_root, "analysis", "brt", "output")
output_dir <- file.path(project_root, "analysis", "ensemble", "output")
plot_dir <- file.path(project_root, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- readRDS(file.path(brms_dir, "brms_model_data.rds")) |>
  mutate(Reef = as.character(Reef), observed_cpue = cots_count / bottom_time)
fold_map <- read.csv(file.path(brms_dir, "brms_reef_fold_map.csv"))
data <- left_join(data, fold_map, by = "Reef")
thresholds <- c(0.02, 0.04, 0.08)

message("Generating GLMM whole-reef out-of-fold predictions...")
glmm_oof <- bind_rows(lapply(sort(unique(data$fold)), function(fold_id) {
  train <- filter(data, fold != fold_id) |>
    mutate(Reef = factor(Reef))
  test <- filter(data, fold == fold_id) |>
    mutate(Reef = factor(Reef, levels = levels(train$Reef)))
  fit <- glmmTMB(
    cots_count ~ edna_prop_z + distance_z + lag_z +
      offset(log_effort) + (1 | Reef),
    family = nbinom2(), data = train,
    control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
  )
  test$Reef <- factor(levels(train$Reef)[1], levels = levels(train$Reef))
  conditional_link <- predict(fit, newdata = test, type = "link", re.form = NA)
  reef_sd <- sqrt(as.numeric(VarCorr(fit)$cond$Reef[1, 1]))
  shape <- sigma(fit)
  predicted_cpue <- exp(conditional_link + 0.5 * reef_sd^2) / test$bottom_time
  set.seed(20263000 + fold_id)
  probability_matrix <- sapply(seq_along(thresholds), function(j) {
    threshold <- thresholds[j]
    vapply(seq_len(nrow(test)), function(i) {
      random_effect <- rnorm(2000, 0, reef_sd)
      mu <- exp(conditional_link[i] + random_effect)
      minimum_count <- ceiling(threshold * test$bottom_time[i])
      mean(1 - pnbinom(minimum_count - 1, mu = mu, size = shape))
    }, numeric(1))
  })
  colnames(probability_matrix) <- paste0(
    "prob_", gsub("\\.", "", sprintf("%.2f", thresholds))
  )
  bind_cols(
    transmute(
      test, cull_id, Reef = as.character(Reef), fold,
      observed_cpue, predicted_cpue
    ),
    as.data.frame(probability_matrix)
  )
}))
write.csv(
  glmm_oof,
  file.path(output_dir, "glmm_operational_oof_predictions.csv"),
  row.names = FALSE
)

brms_oof <- read.csv(file.path(output_dir, "brms_operational_oof_predictions.csv")) |>
  mutate(cull_id = as.character(cull_id)) |>
  select(cull_id, observed_cpue, predicted_cpue, starts_with("prob_")) |>
  rename_with(~ paste0("brms_", .x), -c(cull_id, observed_cpue))
brt_oof <- read.csv(file.path(brt_dir, "brt_oof_predictions.csv")) |>
  mutate(cull_id = as.character(cull_id)) |>
  select(cull_id, predicted_cpue) |>
  rename(brt_predicted_cpue = predicted_cpue)
brt_class <- read.csv(file.path(brt_dir, "brt_threshold_oof_predictions.csv")) |>
  mutate(cull_id = as.character(cull_id)) |>
  mutate(probability_name = paste0(
    "brt_prob_", gsub("\\.", "", sprintf("%.2f", threshold))
  )) |>
  select(cull_id, probability_name, probability) |>
  pivot_wider(names_from = probability_name, values_from = probability)
glmm_join <- glmm_oof |>
  mutate(cull_id = as.character(cull_id)) |>
  select(cull_id, predicted_cpue, starts_with("prob_")) |>
  rename_with(~ paste0("glmm_", .x), -cull_id)

predictions <- brms_oof |>
  left_join(brt_oof, by = "cull_id") |>
  left_join(brt_class, by = "cull_id") |>
  left_join(glmm_join, by = "cull_id")
stopifnot(nrow(predictions) == nrow(data), !anyDuplicated(predictions$cull_id))

weight_grid <- seq(0, 1, by = 0.01)
cpue_weight_scores <- tibble(
  brms_weight = weight_grid,
  rmse = vapply(weight_grid, function(weight) {
    ensemble <- weight * predictions$brms_predicted_cpue +
      (1 - weight) * predictions$brt_predicted_cpue
    sqrt(mean((predictions$observed_cpue - ensemble)^2))
  }, numeric(1))
)
best_cpue_weight <- cpue_weight_scores$brms_weight[which.min(cpue_weight_scores$rmse)]
predictions$ensemble_predicted_cpue <-
  best_cpue_weight * predictions$brms_predicted_cpue +
  (1 - best_cpue_weight) * predictions$brt_predicted_cpue

cpue_long <- predictions |>
  select(
    cull_id, observed_cpue,
    BRMS = brms_predicted_cpue,
    BRT = brt_predicted_cpue,
    GLMM = glmm_predicted_cpue,
    Ensemble = ensemble_predicted_cpue
  ) |>
  pivot_longer(BRMS:Ensemble, names_to = "model", values_to = "prediction")
cpue_metrics <- cpue_long |>
  group_by(model) |>
  summarise(
    RMSE = sqrt(mean((observed_cpue - prediction)^2)),
    MAE = mean(abs(observed_cpue - prediction)),
    correlation = cor(observed_cpue, prediction),
    bias = mean(prediction - observed_cpue),
    .groups = "drop"
  ) |>
  arrange(RMSE)
write.csv(cpue_metrics, file.path(output_dir, "cpue_model_performance.csv"), row.names = FALSE)
write.csv(cpue_weight_scores, file.path(output_dir, "cpue_ensemble_weight_search.csv"), row.names = FALSE)

binary_metrics <- function(observed, probability, cutoff) {
  predicted <- probability >= cutoff
  tp <- sum(predicted & observed == 1)
  fp <- sum(predicted & observed == 0)
  fn <- sum(!predicted & observed == 1)
  tn <- sum(!predicted & observed == 0)
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  tibble(
    cutoff, F1 = f1, recall, precision,
    specificity = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp),
    false_negatives = fn
  )
}

threshold_metrics <- list()
probability_weights <- list()
for (threshold in thresholds) {
  suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
  observed <- as.integer(predictions$observed_cpue >= threshold)
  brms_probability <- predictions[[paste0("brms_prob_", suffix)]]
  brt_probability <- predictions[[paste0("brt_prob_", suffix)]]
  glmm_probability <- predictions[[paste0("glmm_prob_", suffix)]]

  weight_scores <- bind_rows(lapply(weight_grid, function(weight) {
    probability <- weight * brms_probability + (1 - weight) * brt_probability
    tibble(
      brms_weight = weight,
      brier = mean((observed - probability)^2),
      log_loss = -mean(
        observed * log(pmax(probability, 1e-8)) +
          (1 - observed) * log(pmax(1 - probability, 1e-8))
      )
    )
  }))
  best_probability_weight <- weight_scores$brms_weight[which.min(weight_scores$brier)]
  probability_weights[[suffix]] <- tibble(
    threshold, brms_weight = best_probability_weight,
    brt_weight = 1 - best_probability_weight,
    brier = min(weight_scores$brier)
  )
  predictions[[paste0("ensemble_prob_", suffix)]] <-
    best_probability_weight * brms_probability +
    (1 - best_probability_weight) * brt_probability

  probabilities <- list(
    BRMS = brms_probability,
    BRT = brt_probability,
    GLMM = glmm_probability,
    Ensemble = predictions[[paste0("ensemble_prob_", suffix)]]
  )
  threshold_metrics[[suffix]] <- bind_rows(lapply(names(probabilities), function(model_name) {
    probability <- probabilities[[model_name]]
    cutoff_results <- bind_rows(lapply(seq(0.01, 0.99, by = 0.01), function(cutoff) {
      binary_metrics(observed, probability, cutoff)
    }))
    best_cutoff <- cutoff_results |>
      arrange(desc(F1), desc(recall), cutoff) |>
      slice_head(n = 1)
    tibble(
      threshold, model = model_name,
      AUC = as.numeric(pROC::auc(observed, probability, quiet = TRUE)),
      Brier = mean((observed - probability)^2),
      cutoff = best_cutoff$cutoff,
      F1 = best_cutoff$F1,
      recall = best_cutoff$recall,
      precision = best_cutoff$precision,
      specificity = best_cutoff$specificity,
      false_negatives = best_cutoff$false_negatives
    )
  }))
}
threshold_metrics <- bind_rows(threshold_metrics)
probability_weights <- bind_rows(probability_weights)
write.csv(
  threshold_metrics,
  file.path(output_dir, "threshold_model_performance.csv"),
  row.names = FALSE
)
write.csv(
  probability_weights,
  file.path(output_dir, "threshold_ensemble_weights.csv"),
  row.names = FALSE
)
write.csv(
  predictions,
  file.path(output_dir, "ensemble_oof_predictions.csv"),
  row.names = FALSE
)

ensemble_spec <- list(
  cpue_brms_weight = best_cpue_weight,
  cpue_brt_weight = 1 - best_cpue_weight,
  threshold_weights = probability_weights,
  threshold_cutoffs = threshold_metrics |>
    filter(model == "Ensemble") |>
    select(threshold, cutoff),
  thresholds = thresholds,
  validation = "five-fold whole-reef cross-validation",
  created = Sys.time()
)
saveRDS(ensemble_spec, file.path(output_dir, "reefDNA_ensemble_spec.rds"))

blue <- "#2166AC"
orange <- "#D95F02"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin = margin(7, 7, 7, 7)
  )
panel_a <- ggplot(cpue_metrics, aes(RMSE, reorder(model, -RMSE), fill = model)) +
  geom_col(width = 0.65) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "CPUE error", subtitle = "Whole-reef validation", x = "RMSE", y = NULL) +
  paper_theme

plot_limit <- quantile(predictions$observed_cpue, 0.99)
panel_b <- ggplot(
  predictions, aes(observed_cpue, ensemble_predicted_cpue)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.25, size = 1, colour = blue) +
  coord_equal(xlim = c(0, plot_limit), ylim = c(0, plot_limit)) +
  labs(
    title = "Ensemble fit", subtitle = "Axes limited to observed 99th percentile",
    x = "Observed CPUE (min^-1)", y = "Predicted CPUE (min^-1)"
  ) +
  paper_theme

panel_c <- ggplot(
  threshold_metrics,
  aes(factor(threshold), F1, colour = model, group = model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Threshold F1", x = "CPUE threshold", y = "Tuned F1", colour = NULL) +
  paper_theme

panel_d <- ggplot(
  threshold_metrics,
  aes(factor(threshold), recall, colour = model, group = model)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Threshold recall", x = "CPUE threshold", y = "Tuned recall", colour = NULL) +
  paper_theme

assessment_plot <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "cpue_ensemble_assessment.png"), assessment_plot,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "cpue_ensemble_assessment.pdf"), assessment_plot,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)

print(cpue_metrics)
print(probability_weights)
print(threshold_metrics)
