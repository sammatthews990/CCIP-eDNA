suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(posterior)
  library(tidyr)
  library(ggplot2)
})

project_root <- "."
output_dir <- file.path(project_root, "analysis", "site", "output")
fit <- readRDS(file.path(output_dir, "site_visit_brms_final.rds"))
data <- readRDS(file.path(output_dir, "site_visit_model_data.rds"))
scale_spec <- readRDS(file.path(output_dir, "site_visit_scale_spec.rds"))

z_value <- function(raw, variable, transform = identity) {
  transformed <- transform(raw)
  (transformed - scale_spec[[variable]]["center"]) /
    scale_spec[[variable]]["scale"]
}
effort <- median(data$bottom_time)

draws <- as_draws_df(fit)
total_variance <- draws$sd_Reef__Intercept^2 +
  draws$sd_site_id__Intercept^2 + draws$sd_edna_campaign_id__Intercept^2

get_curve_data <- function(lag_days) {
  lag_z_val <- z_value(lag_days, "lag_model", log1p)
  
  df_grid <- tibble(
    Reef = rep(as.character(data$Reef[[1]]), 101),
    site_id = rep(as.character(data$site_id[[1]]), 101),
    edna_campaign_id = rep(as.character(data$edna_campaign_id[[1]]), 101),
    bottom_time = effort,
    log_effort = log(effort),
    distance_z = 0,
    lag_z = lag_z_val,
    edna_pct = 0:100,
    edna_prop_z = z_value((0:100) / 100, "edna_prop_model")
  )
  
  cond_counts <- posterior_epred(fit, newdata = df_grid, re_formula = NA)
  marg_cpue <- sweep(cond_counts, 1, exp(0.5 * total_variance), "*") / effort
  
  # 80% interval (quantiles 0.10 and 0.90)
  med <- apply(marg_cpue, 2, median)
  lo80 <- apply(marg_cpue, 2, quantile, probs = 0.10)
  hi80 <- apply(marg_cpue, 2, quantile, probs = 0.90)
  
  tibble(
    edna_pct = 0:100,
    lag_label = paste0(lag_days, " days"),
    expected = med,
    lower80 = lo80,
    upper80 = hi80
  )
}

res_60 <- get_curve_data(60)
res_270 <- get_curve_data(270)
combined <- bind_rows(res_60, res_270)

print(head(combined))
print(tail(combined))
