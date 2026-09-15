#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(xgboost)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
brms_dir <- file.path(project_root, "analysis", "brms", "output")
output_dir <- file.path(project_root, "analysis", "brt", "output")
plot_dir <- file.path(project_root, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_data <- readRDS(file.path(brms_dir, "brms_model_data.rds")) |>
  mutate(Reef = as.character(Reef), cpue = cots_count / bottom_time)
scale_spec <- readRDS(file.path(brms_dir, "brms_scale_spec.rds"))
fold_map <- read.csv(file.path(brms_dir, "brms_reef_fold_map.csv"))
model_data <- left_join(model_data, fold_map, by = "Reef")
stopifnot(!anyNA(model_data$fold))

predictors <- c("edna_prop_z", "distance_z", "lag_z")
tuning_grid <- tidyr::crossing(
  interaction_depth = 1:3,
  shrinkage = c(0.02, 0.05),
  n_trees = 1200L
)

fit_regression <- function(data, depth, shrinkage, n_trees) {
  training_matrix <- xgb.DMatrix(
    as.matrix(data[predictors]),
    label = data$cots_count,
    base_margin = data$log_effort
  )
  xgb.train(
    params = list(
      objective = "count:poisson",
      eval_metric = "poisson-nloglik",
      eta = shrinkage,
      max_depth = depth,
      min_child_weight = 15,
      subsample = 0.7,
      colsample_bytree = 1,
      max_delta_step = 0.7,
      nthread = 4
    ),
    data = training_matrix,
    nrounds = n_trees,
    verbose = 0
  )
}

predict_regression <- function(fit, data) {
  matrix <- xgb.DMatrix(
    as.matrix(data[predictors]), base_margin = data$log_effort
  )
  pmax(predict(fit, matrix) / exp(data$log_effort), 0)
}

message("Tuning BRT over whole-reef folds...")
tuning_results <- bind_rows(lapply(seq_len(nrow(tuning_grid)), function(i) {
  setting <- tuning_grid[i, ]
  fold_metrics <- bind_rows(lapply(sort(unique(model_data$fold)), function(fold_id) {
    train <- filter(model_data, fold != fold_id)
    test <- filter(model_data, fold == fold_id)
    fit <- fit_regression(
      train, setting$interaction_depth, setting$shrinkage, setting$n_trees
    )
    prediction <- predict_regression(fit, test)
    tibble(
      fold = fold_id,
      rmse = sqrt(mean((test$cpue - prediction)^2)),
      mae = mean(abs(test$cpue - prediction)),
      correlation = suppressWarnings(cor(test$cpue, prediction))
    )
  }))
  bind_cols(setting, summarise(
    fold_metrics,
    rmse_se = sd(rmse) / sqrt(n()),
    rmse = mean(rmse),
    mae = mean(mae),
    correlation = mean(correlation, na.rm = TRUE)
  ))
})) |>
  arrange(rmse, mae)
write.csv(tuning_results, file.path(output_dir, "brt_tuning.csv"), row.names = FALSE)

best <- slice_head(tuning_results, n = 1)
message(
  "Selected depth=", best$interaction_depth,
  ", shrinkage=", best$shrinkage,
  ", trees=", best$n_trees
)

oof_predictions <- bind_rows(lapply(sort(unique(model_data$fold)), function(fold_id) {
  train <- filter(model_data, fold != fold_id)
  test <- filter(model_data, fold == fold_id)
  fit <- fit_regression(
    train, best$interaction_depth, best$shrinkage, best$n_trees
  )
  test |>
    transmute(
      cull_id, Reef, fold, observed_cpue = cpue,
      predicted_cpue = predict_regression(fit, test)
    )
}))
write.csv(
  oof_predictions,
  file.path(output_dir, "brt_oof_predictions.csv"),
  row.names = FALSE
)

final_regression <- fit_regression(
  model_data, best$interaction_depth, best$shrinkage, best$n_trees
)
importance <- xgb.importance(
  feature_names = predictors, model = final_regression
) |>
  as_tibble() |>
  transmute(
    variable = recode(
      Feature,
      edna_prop_z = "% positive",
      distance_z = "Distance",
      lag_z = "Time since sample"
    ),
    importance = 100 * Gain
  )
write.csv(importance, file.path(output_dir, "brt_importance.csv"), row.names = FALSE)

thresholds <- c(0.02, 0.04, 0.08)
classification_oof <- list()
final_classifiers <- list()
fit_classifier <- function(data) {
  matrix <- xgb.DMatrix(
    as.matrix(data[predictors]), label = data$above_threshold
  )
  xgb.train(
    params = list(
      objective = "binary:logistic", eval_metric = "logloss",
      eta = best$shrinkage, max_depth = best$interaction_depth,
      min_child_weight = 15, subsample = 0.7,
      colsample_bytree = 1, nthread = 4
    ),
    data = matrix, nrounds = best$n_trees, verbose = 0
  )
}
predict_classifier <- function(fit, data) {
  predict(fit, xgb.DMatrix(as.matrix(data[predictors])))
}
for (threshold in thresholds) {
  response_name <- paste0("above_", gsub("\\.", "", sprintf("%.2f", threshold)))
  data_threshold <- model_data |>
    mutate(above_threshold = as.integer(cpue >= threshold))
  classification_oof[[response_name]] <- bind_rows(lapply(
    sort(unique(data_threshold$fold)), function(fold_id) {
      train <- filter(data_threshold, fold != fold_id)
      test <- filter(data_threshold, fold == fold_id)
      classifier <- fit_classifier(train)
      tibble(
        cull_id = test$cull_id,
        threshold = threshold,
        observed = test$above_threshold,
        probability = predict_classifier(classifier, test)
      )
    }
  ))
  final_classifiers[[response_name]] <- fit_classifier(data_threshold)
}
classification_oof <- bind_rows(classification_oof)
write.csv(
  classification_oof,
  file.path(output_dir, "brt_threshold_oof_predictions.csv"),
  row.names = FALSE
)

z_value <- function(raw, variable, transform = identity) {
  transformed <- transform(raw)
  (transformed - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}
positive_grid <- tibble(
  edna_pct = 0:100,
  edna_prop_z = z_value((0:100) / 100, "edna_prop_model"),
  distance_z = 0,
  lag_z = 0,
  log_effort = log(median(model_data$bottom_time))
)
context_grid <- crossing(
  edna_pct = seq(0, 100, by = 2),
  distance_m = c(200, 2000),
  lag_days = c(30, 183)
) |>
  mutate(
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model"),
    distance_z = z_value(
      distance_m, "distance_model", function(x) log1p(x / 200)
    ),
    lag_z = z_value(lag_days, "lag_model", log1p),
    log_effort = log(median(model_data$bottom_time)),
    distance_label = factor(
      distance_m, c(200, 2000), c("200 m", "2,000 m")
    ),
    lag_label = factor(lag_days, c(30, 183), c("30 days", "183 days"))
  )

bootstrap_replicates <- 100L
set.seed(20260922)
reef_levels <- unique(model_data$Reef)
positive_boot <- matrix(NA_real_, bootstrap_replicates, nrow(positive_grid))
context_boot <- matrix(NA_real_, bootstrap_replicates, nrow(context_grid))
message("Running ", bootstrap_replicates, " whole-reef bootstrap fits...")
for (b in seq_len(bootstrap_replicates)) {
  sampled_reefs <- sample(reef_levels, length(reef_levels), replace = TRUE)
  bootstrap_data <- bind_rows(lapply(sampled_reefs, function(reef_name) {
    filter(model_data, Reef == reef_name)
  }))
  bootstrap_fit <- fit_regression(
    bootstrap_data,
    best$interaction_depth,
    best$shrinkage,
    best$n_trees
  )
  positive_boot[b, ] <- predict_regression(bootstrap_fit, positive_grid)
  context_boot[b, ] <- predict_regression(bootstrap_fit, context_grid)
  if (b %% 10 == 0) message("Bootstrap ", b, "/", bootstrap_replicates)
}
saveRDS(
  list(positive = positive_boot, context = context_boot),
  file.path(output_dir, "brt_bootstrap_partial_dependence.rds")
)

summarise_bootstrap <- function(matrix_used) {
  tibble(
    estimate = apply(matrix_used, 2, median),
    lower = apply(matrix_used, 2, quantile, probs = 0.025),
    upper = apply(matrix_used, 2, quantile, probs = 0.975)
  )
}
positive_summary <- bind_cols(positive_grid, summarise_bootstrap(positive_boot))
context_summary <- bind_cols(context_grid, summarise_bootstrap(context_boot))
write.csv(
  positive_summary,
  file.path(output_dir, "brt_percent_positive_partial_dependence.csv"),
  row.names = FALSE
)
write.csv(
  context_summary,
  file.path(output_dir, "brt_context_partial_dependence.csv"),
  row.names = FALSE
)

brt_object <- list(
  regression = final_regression,
  classifiers = final_classifiers,
  predictors = predictors,
  thresholds = thresholds,
  scale_spec = scale_spec,
  tuning = best,
  training_ranges = lapply(model_data[predictors], range),
  oof_residual_quantiles = quantile(
    oof_predictions$observed_cpue - oof_predictions$predicted_cpue,
    c(0.025, 0.975)
  ),
  created = Sys.time()
)
saveRDS(brt_object, file.path(output_dir, "reefDNA_brt_model.rds"))

blue <- "#2166AC"
light_blue <- "#92C5DE"
orange <- "#D95F02"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5),
    plot.margin = margin(7, 7, 7, 7)
  )

panel_a <- ggplot(oof_predictions, aes(observed_cpue, predicted_cpue)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_point(alpha = 0.25, size = 1, colour = blue) +
  coord_equal() +
  labs(
    title = "Reef-held-out fit",
    subtitle = paste0("RMSE = ", round(sqrt(mean(
      (oof_predictions$observed_cpue - oof_predictions$predicted_cpue)^2
    )), 3)),
    x = "Observed CPUE (min^-1)", y = "Predicted CPUE (min^-1)"
  ) +
  paper_theme

panel_b <- ggplot(
  importance,
  aes(x = importance, y = reorder(variable, importance), fill = variable)
) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(blue, light_blue, orange), guide = "none") +
  labs(title = "Variable importance", x = "Relative influence (%)", y = NULL) +
  paper_theme

panel_c <- ggplot(positive_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.6) +
  geom_line(colour = blue, linewidth = 0.9) +
  geom_hline(
    yintercept = c(0.02, 0.04, 0.08), linetype = 2,
    colour = c("grey55", "#7B3294", "grey55"), linewidth = 0.4
  ) +
  labs(
    title = "% positive",
    subtitle = "Whole-reef bootstrap 95% interval",
    x = "eDNA samples positive (%)", y = "Partial-dependence CPUE (min^-1)"
  ) +
  paper_theme

panel_d <- ggplot(context_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.6) +
  geom_line(colour = blue, linewidth = 0.75) +
  facet_grid(lag_label ~ distance_label) +
  labs(
    title = "Context dependence",
    subtitle = "Rows: time; columns: distance",
    x = "eDNA samples positive (%)", y = "Partial-dependence CPUE (min^-1)"
  ) +
  paper_theme

four_panel <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "brt_cpue_four_panel.png"), four_panel,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "brt_cpue_four_panel.pdf"), four_panel,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)
message("BRT model, bootstrap summaries, and paper figure saved.")
