#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(fuzzyjoin)
  library(glmmTMB)
  library(readxl)
  library(tibble)
})

project_root <- normalizePath(
  if (file.exists("analysis/brms/output/brms_model_data.rds")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
output_dir <- file.path(project_root, "analysis", "brms", "output")
model_data <- readRDS(file.path(output_dir, "brms_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "brms_scale_spec.rds"))

formulas <- list(
  prop_only = cots_count ~ edna_prop_z + offset(log_effort) + (1 | Reef),
  prop_time = cots_count ~ edna_prop_z + lag_z + offset(log_effort) + (1 | Reef),
  prop_distance_time = cots_count ~ edna_prop_z + distance_z + lag_z +
    offset(log_effort) + (1 | Reef),
  concentration_additive = cots_count ~ edna_prop_z + edna_conc_z +
    distance_z + lag_z + offset(log_effort) + (1 | Reef),
  concentration_distance = cots_count ~ edna_prop_z + edna_conc_z +
    distance_z + lag_z + edna_conc_z:distance_z +
    offset(log_effort) + (1 | Reef)
)

fits <- lapply(formulas, function(formula_used) {
  glmmTMB(
    formula_used, family = nbinom2(), data = model_data,
    control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
  )
})
saveRDS(fits, file.path(output_dir, "quick_glmm_ablation_models.rds"))

effort <- median(model_data$bottom_time)
prediction_grid <- tibble(
  Reef = factor(
    rep(levels(model_data$Reef)[1], 101), levels = levels(model_data$Reef)
  ),
  bottom_time = effort,
  log_effort = log(effort),
  edna_pct = 0:100,
  edna_prop_z = ((edna_pct / 100) - scale_spec$edna_prop_model["center"]) /
    scale_spec$edna_prop_model["scale"],
  edna_conc_z = 0,
  distance_z = 0,
  lag_z = 0
)

find_crossing <- function(y, threshold) {
  if (!any(y >= threshold)) return(NA_real_)
  stats::approx(y, prediction_grid$edna_pct, xout = threshold)$y
}

new_glmm_results <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  link <- predict(fit, newdata = prediction_grid, type = "link", re.form = NA)
  conditional <- exp(link) / effort
  reef_sd <- sqrt(as.numeric(VarCorr(fit)$cond$Reef[1, 1]))
  marginal <- conditional * exp(0.5 * reef_sd^2)
  bind_rows(lapply(c("conditional_typical_reef", "marginal_new_reef"), function(scale_name) {
    curve <- if (scale_name == "conditional_typical_reef") conditional else marginal
    tibble(
      model = model_name,
      prediction_scale = scale_name,
      edna_pct = prediction_grid$edna_pct,
      estimate = curve,
      AIC = AIC(fit),
      reef_sd = reef_sd,
      crossing_002 = find_crossing(curve, 0.02),
      crossing_004 = find_crossing(curve, 0.04),
      crossing_008 = find_crossing(curve, 0.08)
    )
  }))
}))
write.csv(
  new_glmm_results,
  file.path(output_dir, "quick_glmm_ablation_curves.csv"),
  row.names = FALSE
)

# Recreate the old reef/event-aggregated GLMM data for an exact scale audit.
cull_raw <- read_excel(
  file.path(project_root, "data", "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"),
  sheet = "Cull"
)
edna_raw <- read_excel(
  file.path(project_root, "data", "eDNA data_ALL_20260528.xlsx"),
  sheet = "eDNA_data_ALL"
)
cull <- cull_raw |>
  rename(Reef = ReefName) |>
  mutate(date_cull = as.Date(SurveyDate))
edna_agg <- edna_raw |>
  filter(!is.na(Year)) |>
  rename(Collection.organisation = `Collection organisation`) |>
  mutate(
    Reef = ReefName,
    Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
    date_edna = as.Date(Date),
    Conc_mean = as.numeric(Conc_mean)
  ) |>
  arrange(Reef, Year, date_edna) |>
  group_by(Reef, Collection.org, Year) |>
  mutate(grp = cumsum(if_else(
    is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > 7,
    1L, 0L
  ))) |>
  ungroup() |>
  group_by(Reef, Collection.org, Year, grp) |>
  summarise(
    date_edna = min(date_edna),
    conc_mean = mean(Conc_mean, na.rm = TRUE),
    perc_pos = mean(LOD_sample_positive, na.rm = TRUE) * 100,
    .groups = "drop"
  )

get_old_cohort <- function(max_days, min_days, label) {
  fuzzy_inner_join(
    cull, edna_agg,
    by = c("Reef" = "Reef", "date_cull" = "date_edna"),
    match_fun = list(`==`, function(cull_date, edna_date) {
      difference <- as.numeric(cull_date - edna_date)
      difference >= min_days & difference <= max_days
    })
  ) |>
    mutate(Reef = Reef.x) |>
    group_by(Reef, Collection.org, Year, grp) |>
    summarise(
      total_cots = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
      total_bottom = sum(Bottomtime, na.rm = TRUE),
      perc_pos_reef = mean(perc_pos, na.rm = TRUE),
      conc_mean_reef = mean(conc_mean, na.rm = TRUE),
      horizon = label,
      .groups = "drop"
    )
}

old_data <- bind_rows(
  get_old_cohort(91, 0, "0-3 Months"),
  get_old_cohort(183, 92, "3-6 Months"),
  get_old_cohort(365, 184, "6-12 Months")
) |>
  filter(
    total_bottom > 0, is.finite(perc_pos_reef),
    is.finite(conc_mean_reef), conc_mean_reef < 5000
  ) |>
  mutate(Reef = factor(Reef)) |>
  droplevels()

old_results <- bind_rows(lapply(unique(old_data$horizon), function(horizon_used) {
  data_used <- droplevels(filter(old_data, horizon == horizon_used))
  fit <- glmmTMB(
    total_cots ~ perc_pos_reef + offset(log(total_bottom)) + (1 | Reef),
    family = nbinom2(), data = data_used
  )
  grid <- tibble(
    perc_pos_reef = 0:100,
    total_bottom = median(data_used$total_bottom),
    Reef = factor(levels(data_used$Reef)[1], levels = levels(data_used$Reef))
  )
  conditional <- exp(predict(fit, newdata = grid, type = "link", re.form = NA)) /
    grid$total_bottom
  reef_sd <- sqrt(as.numeric(VarCorr(fit)$cond$Reef[1, 1]))
  marginal <- conditional * exp(0.5 * reef_sd^2)
  bind_rows(lapply(c("conditional_typical_reef", "marginal_new_reef"), function(scale_name) {
    curve <- if (scale_name == "conditional_typical_reef") conditional else marginal
    tibble(
      horizon = horizon_used,
      prediction_scale = scale_name,
      edna_pct = 0:100,
      estimate = curve,
      n = nrow(data_used),
      reef_sd = reef_sd,
      crossing_002 = find_crossing(curve, 0.02),
      crossing_004 = find_crossing(curve, 0.04),
      crossing_008 = find_crossing(curve, 0.08)
    )
  }))
}))
write.csv(
  old_results,
  file.path(output_dir, "old_glmm_curve_scale_audit.csv"),
  row.names = FALSE
)

print(new_glmm_results |>
  distinct(model, prediction_scale, AIC, reef_sd, crossing_002, crossing_004))
print(old_results |>
  distinct(horizon, prediction_scale, n, reef_sd, crossing_002, crossing_004))
