library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(purrr)
library(tidyr)
library(sf)
devtools::load_all("reefDNA")

cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull.dat <- read_excel(cull_file, sheet = "Cull")
edna.dat <- read_excel(edna_file, sheet = "eDNA_data_ALL")

cull <- cull.dat %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

edna_agg <- edna.dat %>%
    filter(!is.na(Year)) %>%
    rename(Collection.organisation = `Collection organisation`) %>%
    mutate(
        Reef = ReefName,
        Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
        date_edna = as.Date(Date),
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    arrange(Reef, Year, date_edna) %>%
    group_by(Reef, Collection.org, Year) %>%
    mutate(
        grp = cumsum(
            if_else(
                is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > 7,
                1L, 0L
            )
        )
    ) %>%
    ungroup() %>%
    group_by(Reef, Collection.org, Year, grp) %>%
    summarise(
        date_edna = min(date_edna),
        conc_mean = mean(Conc_mean, na.rm = TRUE),
        perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        .groups   = "drop"
    )

reef_prior_3m  <- get_prior_cohort(cull, edna_agg, 91,  "0-3 Months", min_days = 0)
reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "3-6 Months", min_days = 92)
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "6-12 Months", min_days = 184)

dat_glmm <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(
        conc_t    = log1p(conc_mean_reef),
        combo_index = perc_pos_reef * conc_t,
        obs_cpue  = total_cots / total_bottom,
        horizon   = factor(horizon, levels = c("0-3 Months", "3-6 Months", "6-12 Months"))
    ) %>%
    filter(
        !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
        !is.na(perc_pos_reef), !is.na(conc_t), !is.na(Reef),
        conc_mean_reef < 5000
    ) %>%
    mutate(Reef = factor(Reef)) %>%
    droplevels()

cat("=== 1. GLMM COMBO INDEX FITTING ===\n")
for (hor in levels(dat_glmm$horizon)) {
    dat_sub <- dat_glmm %>% filter(horizon == hor)
    m_combo <- glmmTMB(total_cots ~ combo_index + offset(log(total_bottom)) + (1 | Reef), family = nbinom2(), data = dat_sub)
    c_cpue <- fitted(m_combo) / dat_sub$total_bottom
    rmse_combo <- sqrt(mean((dat_sub$obs_cpue - c_cpue)^2))
    r2_combo   <- cor(dat_sub$obs_cpue, c_cpue)^2
    s_combo    <- summary(m_combo)$coefficients$cond["combo_index", ]
    cat(sprintf("Horizon: %s | N=%d | Beta=%.4f ± %.4f | RMSE=%.4f | R2=%.3f\n",
                hor, nrow(dat_sub), s_combo["Estimate"], s_combo["Std. Error"], rmse_combo, r2_combo))
}

cat("\n=== 2. CV THRESHOLD SEARCH FOR COMBO INDEX ===\n")
dat_evt <- dat_glmm %>%
    filter(horizon == "0-3 Months") %>%
    group_by(Reef, Collection.org, Year) %>%
    summarise(
        counts = sum(total_cots, na.rm = TRUE),
        total_bottom = sum(total_bottom, na.rm = TRUE),
        perc_pos = mean(perc_pos_reef, na.rm = TRUE),
        conc_mean_reef = mean(conc_mean_reef, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        cpue = counts / total_bottom,
        conc_t = log1p(conc_mean_reef),
        combo_index = perc_pos * conc_t
    )

combo_grid <- seq(min(dat_evt$combo_index), quantile(dat_evt$combo_index, 0.99), length.out = 50)
cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)

metrics_for_var <- function(p_star, cpue_thr, dat, var_col = "combo_index") {
    classified <- dat %>%
        mutate(
            pred = if_else(.data[[var_col]] >= p_star, "Pred +", "Pred -"),
            actual = if_else(cpue >= cpue_thr, "Actual +", "Actual -")
        )
    tp <- sum(classified$pred == "Pred +" & classified$actual == "Actual +")
    fp <- sum(classified$pred == "Pred +" & classified$actual == "Actual -")
    fn <- sum(classified$pred == "Pred -" & classified$actual == "Actual +")
    tn <- sum(classified$pred == "Pred -" & classified$actual == "Actual -")
    acc <- (tp + tn) / nrow(classified)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    rec <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((prec + rec) > 0) 2 * (prec * rec) / (prec + rec) else 0
    c(F1 = f1, Accuracy = acc, Precision = prec, Recall = rec)
}

run_cv_var <- function(dat_evt, cpue_thr, var_col, thr_grid, k = 5, R = 100) {
    dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
    strat_ok <- length(unique(dat2$actual)) > 1
    folds <- if (strat_ok) rsample::vfold_cv(dat2, v = k, repeats = R, strata = actual) else rsample::vfold_cv(dat2, v = k, repeats = R)
    
    cv_long <- folds %>%
        mutate(assess = map(splits, rsample::assessment)) %>%
        dplyr::select(id, id2, assess) %>%
        tidyr::expand_grid(thresh = thr_grid) %>%
        mutate(
            out = pmap(list(thresh, assess), \(p, d) metrics_for_var(p, cpue_thr, d, var_col)),
            F1 = vapply(out, `[[`, numeric(1), "F1"),
            Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
            Precision = vapply(out, `[[`, numeric(1), "Precision"),
            Recall = vapply(out, `[[`, numeric(1), "Recall")
        ) %>%
        dplyr::select(-out, -assess)
    
    cv_sum <- cv_long %>%
        group_by(thresh) %>%
        summarise(F1_mean = mean(F1, na.rm = TRUE), F1_sd = sd(F1, na.rm = TRUE), .groups = "drop")
    
    best <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
    p_star <- best$thresh
    
    perf_split <- folds %>%
        mutate(
            assess = map(splits, rsample::assessment),
            m = map(assess, ~ {
                out <- metrics_for_var(p_star, cpue_thr, .x, var_col)
                tibble(F1 = out["F1"], Accuracy = out["Accuracy"], Precision = out["Precision"], Recall = out["Recall"])
            })
        ) %>%
        dplyr::select(m) %>%
        unnest(m)
    
    perf_sum <- perf_split %>%
        summarise(
            cpue_thr = cpue_thr, thr_star = p_star,
            F1_mean = mean(F1), Prec_mean = mean(Precision), Rec_mean = mean(Recall), Acc_mean = mean(Accuracy),
            .groups = "drop"
        )
    list(perf = perf_sum)
}

res_combo <- map(cpue_levels, function(cp) {
    run_cv_var(dat_evt, cp, "combo_index", combo_grid)
})

summary_combo <- map_dfr(compact(res_combo), "perf")
print(summary_combo %>% dplyr::select(cpue_thr, thr_star, F1_mean, Prec_mean, Rec_mean, Acc_mean))

cat("\n=== 3. FALSE POSITIVES & FALSE NEGATIVES DISTANCE ANALYSIS ===\n")

# Extract eDNA sample-level coordinates and nearest cull distance
edna_coords <- edna.dat %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Year)) %>%
    mutate(date_edna = as.Date(Date)) %>%
    group_by(ReefName, Site_name, Year) %>%
    summarise(
        Lat = mean(Lat),
        Long = mean(Long),
        date_edna = min(date_edna),
        conc_mean = mean(as.numeric(Conc_mean), na.rm = TRUE),
        perc_pos = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        .groups = "drop"
    ) %>%
    rename(Reef = ReefName)

cull_coords <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

edna_sf <- st_as_sf(edna_coords, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)

# Calculate distance to nearest cull event for each eDNA site
dist_matrix <- st_distance(edna_sf, cull_coords)
edna_coords$min_cull_dist_m <- apply(dist_matrix, 1, min)

# Aggregate to event level including mean distance
event_dist <- edna_coords %>%
    group_by(Reef, Year) %>%
    summarise(
        mean_cull_dist_m = mean(min_cull_dist_m),
        min_cull_dist_m  = min(min_cull_dist_m),
        .groups = "drop"
    )

# Join with dat_evt at CPUE = 0.04
cpue_star <- 0.04
opt_perc <- 50
opt_conc_t <- log1p(10.5)

dat_eval <- dat_evt %>%
    left_join(event_dist, by = c("Reef", "Year")) %>%
    mutate(
        actual_pos = cpue >= cpue_star,
        pred_pos   = (perc_pos >= opt_perc) & (conc_t >= opt_conc_t),
        class = case_when(
            actual_pos & pred_pos   ~ "True Positive (TP)",
            !actual_pos & !pred_pos ~ "True Negative (TN)",
            !actual_pos & pred_pos  ~ "False Positive (FP)",
            actual_pos & !pred_pos  ~ "False Negative (FN)"
        )
    )

cat("\nSummary of Distances (meters) by Classification Class:\n")
dat_eval %>%
    group_by(class) %>%
    summarise(
        N = n(),
        Mean_Dist_m = mean(mean_cull_dist_m, na.rm = TRUE),
        Median_Dist_m = median(mean_cull_dist_m, na.rm = TRUE),
        Min_Dist_m = min(mean_cull_dist_m, na.rm = TRUE),
        Max_Dist_m = max(mean_cull_dist_m, na.rm = TRUE)
    ) %>%
    print()

cat("\nDetail Table of False Positives & False Negatives:\n")
fp_fn_table <- dat_eval %>%
    filter(class %in% c("False Positive (FP)", "False Negative (FN)")) %>%
    dplyr::select(Reef, Year, class, perc_pos, conc_mean_reef, cpue, mean_cull_dist_m, min_cull_dist_m) %>%
    arrange(class, desc(cpue))

print(fp_fn_table)
