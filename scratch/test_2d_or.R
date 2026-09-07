library(dplyr)
library(readxl)
library(ggplot2)
library(purrr)
library(tidyr)
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

reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "6 Months", min_days = 0)

dat_glmm <- reef_prior_6m %>%
    mutate(
        conc_t    = log1p(conc_mean_reef),
        obs_cpue  = total_cots / total_bottom
    ) %>%
    filter(
        !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
        !is.na(perc_pos_reef), !is.na(conc_t), !is.na(Reef),
        conc_mean_reef < 5000
    ) %>%
    mutate(Reef = factor(Reef)) %>%
    droplevels()

dat_evt <- dat_glmm %>%
    group_by(Reef, Collection.org, Year) %>%
    summarise(
        counts = sum(total_cots, na.rm = TRUE),
        total_bottom = sum(total_bottom, na.rm = TRUE),
        perc_pos = mean(perc_pos_reef, na.rm = TRUE),
        conc_mean_reef = mean(conc_mean_reef, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        reef = Reef,
        cpue = counts / total_bottom,
        conc_t = log1p(conc_mean_reef)
    )

cpue_levels <- c(0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15)
conc_t_max <- quantile(dat_evt$conc_t, 0.99, na.rm = TRUE)
conc_grid_t <- seq(min(dat_evt$conc_t, na.rm = TRUE), conc_t_max, length.out = 30)

cat("=== COMPUTING 2D AND RULE ===\n")
results_2d_and <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "and"))
summary_and <- map_dfr(compact(results_2d_and), "perf")

cat("=== COMPUTING 2D OR RULE ===\n")
results_2d_or <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "or"))
summary_or <- map_dfr(compact(results_2d_or), "perf")

cat("\nAND Rule Performance:\n")
print(summary_and %>% dplyr::select(cpue_thr, perc_star, conc_mean_star, F1_mean, Prec_mean, Rec_mean))

cat("\nOR Rule Performance:\n")
print(summary_or %>% dplyr::select(cpue_thr, perc_star, conc_mean_star, F1_mean, Prec_mean, Rec_mean))
