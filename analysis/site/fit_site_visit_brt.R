#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(xgboost)
})

project_root <- normalizePath(
  if (file.exists("analysis/site/output/site_visit_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "site", "output")
plot_dir <- file.path(project_root, "plots")
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds")) |>
  mutate(
    Reef = as.character(Reef), site_id = as.character(site_id),
    edna_campaign_id = as.character(edna_campaign_id),
    sampling_design = factor(
      as.character(sampling_design),
      levels = c("3x12", "4x6", "other_or_mixed")
    ),
    design_4x6 = as.integer(sampling_design == "4x6"),
    design_other_or_mixed = as.integer(sampling_design == "other_or_mixed"),
    cpue = cots_count / bottom_time
  )
scale_spec <- readRDS(file.path(output_dir, "site_visit_scale_spec.rds"))
set.seed(20261012)

# Balance complete eDNA campaigns across five validation folds. No campaign is
# allowed to contribute to both training and assessment data.
campaign_fold_map <- data |>
  count(edna_campaign_id, sort = TRUE, name = "n_rows") |>
  mutate(fold = rep(1:5, length.out = n())) |>
  select(edna_campaign_id, fold, n_rows)
data <- left_join(data, campaign_fold_map, by = "edna_campaign_id")
write.csv(
  campaign_fold_map,
  file.path(output_dir, "site_visit_campaign_fold_map.csv"), row.names = FALSE
)

predictors <- c(
  "edna_prop_z", "distance_z", "lag_z",
  "design_4x6", "design_other_or_mixed"
)
fit_regression <- function(training, depth, eta, n_trees, predictor_names = predictors) {
  matrix <- xgb.DMatrix(
    as.matrix(training[predictor_names]), label = training$cots_count,
    base_margin = training$log_effort
  )
  xgb.train(
    params = list(
      objective = "count:poisson", eval_metric = "poisson-nloglik",
      eta = eta, max_depth = depth, min_child_weight = 8,
      subsample = 0.8, colsample_bytree = 1,
      max_delta_step = 0.7, nthread = 4
    ),
    data = matrix, nrounds = n_trees, verbose = 0
  )
}
predict_regression <- function(fit, newdata, predictor_names = predictors) {
  matrix <- xgb.DMatrix(
    as.matrix(newdata[predictor_names]), base_margin = newdata$log_effort
  )
  pmax(predict(fit, matrix) / newdata$bottom_time, 0)
}

tuning_grid <- crossing(
  interaction_depth = 1:3,
  shrinkage = c(0.02, 0.05),
  n_trees = 800L
)
tuning <- bind_rows(lapply(seq_len(nrow(tuning_grid)), function(i) {
  setting <- tuning_grid[i, ]
  fold_results <- bind_rows(lapply(1:5, function(fold_id) {
    train <- filter(data, fold != fold_id)
    test <- filter(data, fold == fold_id)
    fit <- fit_regression(
      train, setting$interaction_depth, setting$shrinkage, setting$n_trees
    )
    estimate <- predict_regression(fit, test)
    tibble(
      fold = fold_id,
      rmse = sqrt(mean((test$cpue - estimate)^2)),
      mae = mean(abs(test$cpue - estimate)),
      correlation = suppressWarnings(cor(test$cpue, estimate))
    )
  }))
  bind_cols(
    setting,
    summarise(
      fold_results,
      rmse_se = sd(rmse) / sqrt(n()), rmse = mean(rmse),
      mae = mean(mae), correlation = mean(correlation, na.rm = TRUE)
    )
  )
})) |>
  arrange(rmse, mae)
write.csv(tuning, file.path(output_dir, "site_visit_brt_tuning.csv"), row.names = FALSE)
best <- slice_head(tuning, n = 1)

oof <- bind_rows(lapply(1:5, function(fold_id) {
  train <- filter(data, fold != fold_id)
  test <- filter(data, fold == fold_id)
  fit <- fit_regression(train, best$interaction_depth, best$shrinkage, best$n_trees)
  test |>
    transmute(
      site_visit_id, Reef, site_id, edna_campaign_id, fold,
      sampling_design,
      observed_cpue = cpue,
      predicted_cpue = predict_regression(fit, test)
    )
}))
write.csv(oof, file.path(output_dir, "site_visit_brt_oof.csv"), row.names = FALSE)

ablation_predictors <- list(
  base = c("edna_prop_z", "distance_z", "lag_z"),
  sampling_design = predictors
)
design_ablation <- bind_rows(lapply(names(ablation_predictors), function(model_name) {
  predictor_names <- ablation_predictors[[model_name]]
  fold_results <- bind_rows(lapply(1:5, function(fold_id) {
    train <- filter(data, fold != fold_id)
    test <- filter(data, fold == fold_id)
    fit <- fit_regression(
      train, best$interaction_depth, best$shrinkage, best$n_trees,
      predictor_names = predictor_names
    )
    estimate <- predict_regression(fit, test, predictor_names = predictor_names)
    tibble(
      fold = fold_id,
      rmse = sqrt(mean((test$cpue - estimate)^2)),
      mae = mean(abs(test$cpue - estimate)),
      correlation = suppressWarnings(cor(test$cpue, estimate))
    )
  }))
  summarise(
    fold_results,
    model = model_name,
    rmse_se = sd(rmse) / sqrt(n()), rmse = mean(rmse),
    mae = mean(mae), correlation = mean(correlation, na.rm = TRUE)
  )
})) |>
  arrange(rmse)
write.csv(
  design_ablation,
  file.path(output_dir, "site_visit_brt_design_ablation.csv"), row.names = FALSE
)

final_regression <- fit_regression(
  data, best$interaction_depth, best$shrinkage, best$n_trees
)
importance <- tibble(Feature = predictors) |>
  left_join(
    xgb.importance(feature_names = predictors, model = final_regression) |>
      as_tibble() |>
      select(Feature, Gain),
    by = "Feature"
  ) |>
  mutate(Gain = replace_na(Gain, 0)) |>
  transmute(
    variable = recode(
      Feature,
      edna_prop_z = "% positive", distance_z = "Distance",
      lag_z = "Time since sample", design_4x6 = "Sampling design: 4x6",
      design_other_or_mixed = "Sampling design: other/mixed"
    ),
    importance = 100 * Gain
  )
write.csv(importance, file.path(output_dir, "site_visit_brt_importance.csv"), row.names = FALSE)

thresholds <- c(0.02, 0.04, 0.08)
fit_classifier <- function(training) {
  matrix <- xgb.DMatrix(
    as.matrix(training[predictors]), label = training$above_threshold
  )
  xgb.train(
    params = list(
      objective = "binary:logistic", eval_metric = "logloss",
      eta = best$shrinkage, max_depth = best$interaction_depth,
      min_child_weight = 8, subsample = 0.8,
      colsample_bytree = 1, nthread = 4
    ),
    data = matrix, nrounds = best$n_trees, verbose = 0
  )
}
classification_oof <- list()
classifiers <- list()
for (threshold in thresholds) {
  suffix <- gsub("\\.", "", sprintf("%.2f", threshold))
  threshold_data <- mutate(data, above_threshold = as.integer(cpue >= threshold))
  classification_oof[[suffix]] <- bind_rows(lapply(1:5, function(fold_id) {
    train <- filter(threshold_data, fold != fold_id)
    test <- filter(threshold_data, fold == fold_id)
    classifier <- fit_classifier(train)
    tibble(
      site_visit_id = test$site_visit_id,
      threshold = threshold,
      observed = test$above_threshold,
      probability = predict(
        classifier, xgb.DMatrix(as.matrix(test[predictors]))
      )
    )
  }))
  classifiers[[paste0("above_", suffix)]] <- fit_classifier(threshold_data)
}
classification_oof <- bind_rows(classification_oof)
write.csv(
  classification_oof,
  file.path(output_dir, "site_visit_brt_threshold_oof.csv"), row.names = FALSE
)

z_value <- function(raw, variable, transform = identity) {
  transformed <- transform(raw)
  (transformed - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}
effort <- median(data$bottom_time)
positive_grid <- tibble(
  edna_pct = 0:100,
  edna_prop_z = z_value((0:100) / 100, "edna_prop_model"),
  distance_z = 0, lag_z = 0,
  design_4x6 = 0L, design_other_or_mixed = 0L,
  bottom_time = effort, log_effort = log(effort),
  sampling_design = factor("3x12", levels = levels(data$sampling_design))
)
context_grid <- crossing(
  edna_pct = seq(0, 100, by = 2),
  distance_m = c(200, 2000), lag_days = c(30, 183)
) |>
  mutate(
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model"),
    distance_z = z_value(distance_m, "distance_model", function(x) log1p(x / 200)),
    lag_z = z_value(lag_days, "lag_model", log1p),
    design_4x6 = 0L, design_other_or_mixed = 0L,
    bottom_time = effort, log_effort = log(effort),
    distance_label = factor(distance_m, c(200, 2000), c("200 m", "2,000 m")),
    lag_label = factor(lag_days, c(30, 183), c("30 days", "183 days"))
  )

design_grid <- crossing(
  edna_pct = seq(0, 100, by = 2),
  sampling_design = factor(levels(data$sampling_design), levels = levels(data$sampling_design))
) |>
  mutate(
    edna_prop_z = z_value(edna_pct / 100, "edna_prop_model"),
    distance_z = 0, lag_z = 0,
    design_4x6 = as.integer(sampling_design == "4x6"),
    design_other_or_mixed = as.integer(sampling_design == "other_or_mixed"),
    bottom_time = effort, log_effort = log(effort)
  )

# Resample complete reefs so bootstrap intervals retain both site and campaign
# clustering within each reef.
bootstrap_replicates <- 100L
set.seed(20261012)
reef_levels <- unique(data$Reef)
positive_boot <- matrix(NA_real_, bootstrap_replicates, nrow(positive_grid))
context_boot <- matrix(NA_real_, bootstrap_replicates, nrow(context_grid))
design_boot <- matrix(NA_real_, bootstrap_replicates, nrow(design_grid))
for (b in seq_len(bootstrap_replicates)) {
  sampled <- sample(reef_levels, length(reef_levels), replace = TRUE)
  bootstrap_data <- bind_rows(lapply(seq_along(sampled), function(i) {
    filter(data, Reef == sampled[[i]])
  }))
  fit <- fit_regression(
    bootstrap_data, best$interaction_depth, best$shrinkage, best$n_trees
  )
  positive_boot[b, ] <- predict_regression(fit, positive_grid)
  context_boot[b, ] <- predict_regression(fit, context_grid)
  design_boot[b, ] <- predict_regression(fit, design_grid)
}
summarise_bootstrap <- function(values) {
  tibble(
    estimate = apply(values, 2, median),
    lower = apply(values, 2, quantile, probs = 0.025),
    upper = apply(values, 2, quantile, probs = 0.975)
  )
}
positive_summary <- bind_cols(positive_grid, summarise_bootstrap(positive_boot))
context_summary <- bind_cols(context_grid, summarise_bootstrap(context_boot))
design_summary <- bind_cols(design_grid, summarise_bootstrap(design_boot))
write.csv(
  positive_summary,
  file.path(output_dir, "site_visit_brt_percent_positive.csv"), row.names = FALSE
)
write.csv(
  context_summary,
  file.path(output_dir, "site_visit_brt_context.csv"), row.names = FALSE
)
write.csv(
  design_summary,
  file.path(output_dir, "site_visit_brt_sampling_design.csv"), row.names = FALSE
)
saveRDS(
  list(positive = positive_boot, context = context_boot, design = design_boot),
  file.path(output_dir, "site_visit_brt_bootstrap.rds")
)

model <- list(
  regression = final_regression, classifiers = classifiers,
  predictors = predictors, thresholds = thresholds,
  scale_spec = scale_spec, tuning = best,
  sampling_design_levels = levels(data$sampling_design),
  training_ranges = lapply(data[predictors], range),
  default_effort = effort,
  oof_residual_quantiles = quantile(
    oof$observed_cpue - oof$predicted_cpue, c(0.025, 0.975)
  ),
  linkage = "unique eDNA campaign to first subsequent visit at each cull site",
  created = Sys.time()
)
saveRDS(model, file.path(output_dir, "site_visit_brt_model.rds"))

blue <- "#2166AC"
light_blue <- "#92C5DE"
orange <- "#D95F02"
paper_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    panel.grid.minor = element_blank(), legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 8.5)
  )
limits <- range(c(oof$observed_cpue, oof$predicted_cpue))
panel_a <- ggplot(oof, aes(observed_cpue, predicted_cpue)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_point(alpha = 0.4, size = 1.2, colour = blue) +
  coord_equal(xlim = limits, ylim = limits) +
  labs(
    title = "Campaign-held-out fit",
    subtitle = paste0("RMSE = ", round(sqrt(mean((oof$observed_cpue - oof$predicted_cpue)^2)), 3)),
    x = "Observed CPUE (min^-1)", y = "Predicted CPUE (min^-1)"
  ) + paper_theme
panel_b <- ggplot(importance, aes(importance, reorder(variable, importance), fill = variable)) +
  geom_col(width = 0.65) +
  scale_fill_manual(
    values = c(blue, light_blue, orange, "#7B3294", "#008837"), guide = "none"
  ) +
  labs(title = "Variable importance", x = "Relative influence (%)", y = NULL) +
  paper_theme
panel_c <- ggplot(positive_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.6) +
  geom_line(colour = blue, linewidth = 0.9) +
  geom_hline(yintercept = c(0.02, 0.04, 0.08), linetype = 2, colour = "grey50") +
  labs(
    title = "% positive", subtitle = "Whole-reef bootstrap 95% interval",
    x = "eDNA samples positive (%)", y = "Partial-dependence CPUE (min^-1)"
  ) + paper_theme
panel_d <- ggplot(context_summary, aes(edna_pct, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = light_blue, alpha = 0.6) +
  geom_line(colour = blue, linewidth = 0.75) +
  facet_grid(lag_label ~ distance_label) +
  labs(
    title = "Context", subtitle = "Rows: time; columns: distance",
    x = "eDNA samples positive (%)", y = "Partial-dependence CPUE (min^-1)"
  ) + paper_theme
figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(
  file.path(plot_dir, "site_visit_brt_four_panel.png"), figure,
  width = 14, height = 10.5, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(plot_dir, "site_visit_brt_four_panel.pdf"), figure,
  width = 14, height = 10.5, units = "in", device = cairo_pdf
)
print(best)
print(design_ablation)
print(importance)
