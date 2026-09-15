#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
  library(dplyr)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_all_pairwise.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
set.seed(20260909)
folds <- loo::kfold_split_grouped(K = 5, x = model_data$Reef)
fold_map <- tibble(Reef = model_data$Reef, fold = folds) |>
  distinct() |>
  arrange(fold, Reef)
write.csv(fold_map, file.path(output_dir, "brms_reef_fold_map.csv"), row.names = FALSE)

model_names <- c("all_pairwise", "concentration_distance")
kfold_objects <- vector("list", length(model_names))
names(kfold_objects) <- model_names
for (i in seq_along(model_names)) {
  model_name <- model_names[[i]]
  cache_path <- file.path(output_dir, paste0("brms_", model_name, "_reef_kfold.rds"))
  if (file.exists(cache_path)) {
    message("Loading cached grouped K-fold: ", model_name)
    kfold_objects[[i]] <- readRDS(cache_path)
    next
  }
  message("Running grouped K-fold: ", model_name)
  fit <- readRDS(file.path(output_dir, paste0("brms_", model_name, ".rds")))
  validation <- brms::kfold(
    fit,
    folds = folds,
    group = "Reef",
    joint = "group",
    save_fits = FALSE,
    recompile = FALSE,
    chains = 2,
    cores = 2,
    iter = 1200,
    warmup = 600,
    seed = 20260920 + i,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  saveRDS(validation, cache_path)
  kfold_objects[[i]] <- validation
}

saveRDS(kfold_objects, file.path(output_dir, "brms_reef_kfold_objects.rds"))
comparison_matrix <- loo::loo_compare(kfold_objects)
comparison <- data.frame(unclass(comparison_matrix))
comparison$model <- rownames(comparison)
comparison <- comparison[, c("model", setdiff(names(comparison), "model"))]
write.csv(comparison, file.path(output_dir, "brms_reef_kfold_comparison.csv"), row.names = FALSE)
print(comparison_matrix)
