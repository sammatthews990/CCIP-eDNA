library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(glmmTMB)
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

reef_prior_3m  <- get_prior_cohort(cull, edna_agg, 91,  "0-3 Months", min_days = 0)
reef_prior_6m  <- get_prior_cohort(cull, edna_agg, 183, "3-6 Months", min_days = 92)
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "6-12 Months", min_days = 184)

dat_glmm <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(
        conc_t    = log1p(conc_mean_reef),
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

horizons <- levels(dat_glmm$horizon)
nd_perc_list <- list()
nd_conc_list <- list()
stats_perc_list <- list()
stats_conc_list <- list()

for (hor in horizons) {
    dat_sub <- dat_glmm %>% filter(horizon == hor)
    
    m_perc <- glmmTMB(total_cots ~ perc_pos_reef + offset(log(total_bottom)) + (1 | Reef), family = nbinom2(), data = dat_sub)
    fit_p <- predict(m_perc, newdata = tibble(perc_pos_reef = seq(min(dat_sub$perc_pos_reef), max(dat_sub$perc_pos_reef), length.out = 100), total_bottom = median(dat_sub$total_bottom, na.rm = TRUE), Reef = NA), type = "link", se.fit = TRUE, re.form = NA)
    
    nd_perc_list[[hor]] <- tibble(
        horizon = hor,
        perc_pos_reef = seq(min(dat_sub$perc_pos_reef), max(dat_sub$perc_pos_reef), length.out = 100),
        fit_cpue = exp(fit_p$fit) / median(dat_sub$total_bottom, na.rm = TRUE),
        lcl_cpue = exp(fit_p$fit - 1.96 * fit_p$se.fit) / median(dat_sub$total_bottom, na.rm = TRUE),
        ucl_cpue = exp(fit_p$fit + 1.96 * fit_p$se.fit) / median(dat_sub$total_bottom, na.rm = TRUE)
    )
    
    p_cpue <- fitted(m_perc) / dat_sub$total_bottom
    rmse_p <- sqrt(mean((dat_sub$obs_cpue - p_cpue)^2))
    r2_p   <- cor(dat_sub$obs_cpue, p_cpue)^2
    s_p    <- summary(m_perc)$coefficients$cond["perc_pos_reef", ]
    beta_p <- s_p["Estimate"]
    se_p   <- s_p["Std. Error"]
    
    stats_perc_list[[hor]] <- tibble(
        horizon = hor,
        label = sprintf("N = %d\nEffect (β) = %.4f ± %.4f\nRMSE = %.4f | R² = %.3f", nrow(dat_sub), beta_p, se_p, rmse_p, r2_p),
        x = min(dat_sub$perc_pos_reef) + 0.02 * (max(dat_sub$perc_pos_reef) - min(dat_sub$perc_pos_reef)),
        y = 0.25
    )
    
    m_conc <- glmmTMB(total_cots ~ conc_t + offset(log(total_bottom)) + (1 | Reef), family = nbinom2(), data = dat_sub)
    fit_c <- predict(m_conc, newdata = tibble(conc_t = seq(min(dat_sub$conc_t), max(dat_sub$conc_t), length.out = 100), total_bottom = median(dat_sub$total_bottom, na.rm = TRUE), Reef = NA), type = "link", se.fit = TRUE, re.form = NA)
    
    nd_conc_list[[hor]] <- tibble(
        horizon = hor,
        conc_t = seq(min(dat_sub$conc_t), max(dat_sub$conc_t), length.out = 100),
        fit_cpue = exp(fit_c$fit) / median(dat_sub$total_bottom, na.rm = TRUE),
        lcl_cpue = exp(fit_c$fit - 1.96 * fit_c$se.fit) / median(dat_sub$total_bottom, na.rm = TRUE),
        ucl_cpue = exp(fit_c$fit + 1.96 * fit_c$se.fit) / median(dat_sub$total_bottom, na.rm = TRUE)
    )
    
    c_cpue <- fitted(m_conc) / dat_sub$total_bottom
    rmse_c <- sqrt(mean((dat_sub$obs_cpue - c_cpue)^2))
    r2_c   <- cor(dat_sub$obs_cpue, c_cpue)^2
    s_c    <- summary(m_conc)$coefficients$cond["conc_t", ]
    beta_c <- s_c["Estimate"]
    se_c   <- s_c["Std. Error"]
    
    stats_conc_list[[hor]] <- tibble(
        horizon = hor,
        label = sprintf("N = %d\nEffect (β) = %.3f ± %.3f\nRMSE = %.4f | R² = %.3f", nrow(dat_sub), beta_c, se_c, rmse_c, r2_c),
        x = min(dat_sub$conc_t) + 0.02 * (max(dat_sub$conc_t) - min(dat_sub$conc_t)),
        y = 0.25
    )
}

nd_perc <- bind_rows(nd_perc_list) %>% mutate(horizon = factor(horizon, levels = horizons))
nd_conc <- bind_rows(nd_conc_list) %>% mutate(horizon = factor(horizon, levels = horizons))
df_stats_perc <- bind_rows(stats_perc_list) %>% mutate(horizon = factor(horizon, levels = horizons))
df_stats_conc <- bind_rows(stats_conc_list) %>% mutate(horizon = factor(horizon, levels = horizons))

window_colors <- c(
    "0-3 Months"  = "#009E73",
    "3-6 Months"  = "#E69F00",
    "6-12 Months" = "#D55E00"
)

p_perc <- ggplot(dat_glmm, aes(x = perc_pos_reef, y = obs_cpue)) +
    geom_point(aes(color = horizon), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_perc, aes(x = perc_pos_reef, ymin = lcl_cpue, ymax = ucl_cpue, fill = horizon), alpha = 0.25, inherit.aes = FALSE) +
    geom_line(data = nd_perc, aes(x = perc_pos_reef, y = fit_cpue, color = horizon), linewidth = 1, inherit.aes = FALSE) +
    geom_text(data = df_stats_perc, aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 3.6, fontface = "bold", color = "black", inherit.aes = FALSE) +
    facet_wrap(~horizon, ncol = 1, dir = "v") +
    scale_color_manual(values = window_colors) +
    scale_fill_manual(values = window_colors) +
    coord_cartesian(ylim = c(0, 0.28)) +
    labs(x = "% eDNA Positive", y = "CPUE", title = "A. % eDNA Positive vs CPUE") +
    theme_bw(base_family = "Helvetica") +
    theme(
        legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = 11)
    )

p_conc <- ggplot(dat_glmm, aes(x = conc_t, y = obs_cpue)) +
    geom_point(aes(color = horizon), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_conc, aes(x = conc_t, ymin = lcl_cpue, ymax = ucl_cpue, fill = horizon), alpha = 0.25, inherit.aes = FALSE) +
    geom_line(data = nd_conc, aes(x = conc_t, y = fit_cpue, color = horizon), linewidth = 1, inherit.aes = FALSE) +
    geom_text(data = df_stats_conc, aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 3.6, fontface = "bold", color = "black", inherit.aes = FALSE) +
    facet_wrap(~horizon, ncol = 1, dir = "v") +
    scale_color_manual(values = window_colors) +
    scale_fill_manual(values = window_colors) +
    coord_cartesian(ylim = c(0, 0.28)) +
    labs(x = "log1p(Mean Concentration)", y = "CPUE", title = "B. log1p(Mean Concentration) vs CPUE") +
    theme_bw(base_family = "Helvetica") +
    theme(
        legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = 11)
    )

p_6panel <- p_perc | p_conc
ggsave("scratch/test_ribbon_coord_cartesian.png", p_6panel, width = 9.5, height = 11, dpi = 300)
cat("Ribbon fix plot saved successfully!\n")
