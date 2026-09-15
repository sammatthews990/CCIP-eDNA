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

reef_prior_3m <- get_prior_cohort(cull, edna_agg, 91, "3 Months")
reef_prior_6m <- get_prior_cohort(cull, edna_agg, 183, "6 Months")
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "12 Months")

reef_before <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(horizon = factor(horizon, levels = c("3 Months", "6 Months", "12 Months")))

dat_glmm <- reef_before %>%
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
    filter(horizon == "6 Months") %>%
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

perc_grid <- seq(0, 100, by = 1)

for (cp in c(0.02, 0.04)) {
    cat("\n==========================================\n")
    cat("EVALUATING CPUE TARGET =", cp, "\n")
    cat("==========================================\n")
    grid_eval <- tibble(perc_thresh = perc_grid) %>%
        mutate(
            out = map(perc_thresh, ~ metrics_for(.x, cp, dat_evt)),
            F1 = vapply(out, `[[`, numeric(1), "F1"),
            Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
            Precision = vapply(out, `[[`, numeric(1), "Precision"),
            Recall = vapply(out, `[[`, numeric(1), "Recall")
        )

    best_f1 <- grid_eval %>% slice_max(F1, n = 1, with_ties = FALSE)
    cat("F1 Optimal Threshold:\n")
    print(best_f1)

    max_rec <- max(grid_eval$Recall)
    best_rec_min_fn <- grid_eval %>% filter(Recall >= 0.98) %>% slice_max(perc_thresh, n = 1, with_ties = FALSE)
    cat("Highest threshold with Recall >= 98% (Minimizing FN):\n")
    print(best_rec_min_fn)

    cm_f1 <- make_confusion(dat_evt, perc_thresh = best_f1$perc_thresh, cpue_thresh = cp)
    cat("\nConfusion Matrix (F1 Optimised, %pos >=", best_f1$perc_thresh, "):\n")
    print(cm_f1$cm)

    cm_rec <- make_confusion(dat_evt, perc_thresh = best_rec_min_fn$perc_thresh, cpue_thresh = cp)
    cat("\nConfusion Matrix (Recall Optimised, %pos >=", best_rec_min_fn$perc_thresh, "):\n")
    print(cm_rec$cm)
}
