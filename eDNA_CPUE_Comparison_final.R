library(dplyr)
library(ggplot2)
library(readxl)
library(fuzzyjoin)
library(broom)
library(MASS)
library(glmmTMB)
library(mgcv)
library(tidyr)
library(patchwork)
library(rsample)
library(purrr)
library(tibble)
library(viridis)
library(scales)

# Setup a clean plotting theme
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
# Formatting labels for linear models
lm_label <- function(model, digits_r2 = 2, digits_p = 3, digits_rmse = 2) {
    g <- glance(model)
    r2 <- g$r.squared
    p <- g$p.value
    if (is.na(p)) {
        fs <- summary(model)$fstatistic
        p <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
    }
    rmse <- sqrt(mean(residuals(model)^2))
    p_txt <- if (p < 0.001) "< 0.001" else formatC(p, format = "f", digits = digits_p)
    sprintf(paste0("R^2 = %.", digits_r2, "f\n", "p = %s\n", "RMSE = %.", digits_rmse, "f"), r2, p_txt, rmse)
}

# Confusion matrices and Evaluation Metrics
make_confusion <- function(dat, perc_thresh, cpue_thresh) {
    classified <- dat %>%
        mutate(
            pred   = if_else(perc_pos >= perc_thresh, "Pred +", "Pred -"),
            actual = if_else(cpue >= cpue_thresh, "Actual +", "Actual -")
        )

    cm <- classified %>%
        count(actual, pred, .drop = FALSE) %>%
        complete(
            actual = c("Actual +", "Actual -"),
            pred = c("Pred +", "Pred -"),
            fill = list(n = 0)
        ) %>%
        group_by(actual) %>%
        mutate(row_prop = n / sum(n)) %>%
        ungroup() %>%
        mutate(label = paste0(n, "\n", scales::percent(row_prop, accuracy = 0.1)))

    TP <- cm$n[cm$actual == "Actual +" & cm$pred == "Pred +"]
    FP <- cm$n[cm$actual == "Actual -" & cm$pred == "Pred +"]
    TN <- cm$n[cm$actual == "Actual -" & cm$pred == "Pred -"]
    FN <- cm$n[cm$actual == "Actual +" & cm$pred == "Pred -"]

    acc <- (TP + TN) / (TP + TN + FP + FN + 1e-9)
    sens <- TP / (TP + FN + 1e-9) # Recall / TPR
    spec <- TN / (TN + FP + 1e-9)
    prec <- TP / (TP + FP + 1e-9) # PPV
    f1 <- 2 * prec * sens / (prec + sens + 1e-9)

    list(
        cm = cm,
        metrics = tibble(
            perc_thresh = perc_thresh, cpue_thresh = cpue_thresh, Accuracy = acc,
            Sensitivity_TPR = sens, Specificity_TNR = spec, Precision_PPV = prec, F1 = f1
        )
    )
}

plot_confusion <- function(cm_obj, title = "Confusion matrix") {
    ggplot(cm_obj$cm, aes(x = pred, y = actual, fill = row_prop)) +
        geom_tile(color = "grey85", linewidth = 0.4) +
        geom_text(aes(label = label), size = 5) +
        scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
        labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = title) +
        coord_fixed()
}

metrics_for <- function(p_thr, c_thr, x) {
    pred <- x$perc_pos >= p_thr
    actual <- x$cpue >= c_thr
    TP <- sum(pred & actual)
    FP <- sum(pred & !actual)
    TN <- sum(!pred & !actual)
    FN <- sum(!pred & actual)
    prec <- TP / (TP + FP + 1e-9)
    rec <- TP / (TP + FN + 1e-9)
    f1 <- 2 * prec * rec / (prec + rec + 1e-9)
    acc <- (TP + TN) / (TP + TN + FP + FN + 1e-9)
    c(F1 = f1, Accuracy = acc, Precision = prec, Recall = rec)
}
# Load Cull data
cull.dat <- read_excel("data/260201_COTS-Cull-Data-Ewels.xlsx", sheet = "Cull")

# Load eDNA data
edna.dat <- read_excel("data/eDNA data_ALL_20260225.xlsx", sheet = "eDNA_data_ALL")

# Clean Cull Data
cull <- cull.dat %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

# Aggregate eDNA data
edna_agg <- edna.dat %>%
    filter(!is.na(Year)) %>%
    rename(Collection.organisation = `Collection organisation`) %>%
    mutate(
        Reef = ReefName,
        Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
        date_edna = as.Date(Date), # Excel supplies native dates
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    arrange(Reef, Year, date_edna) %>%
    group_by(Reef, Collection.org, Year) %>%
    # cluster rows into groups where consecutive samples are ≤ 7 days apart
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
win_days <- 183

same6_any <- fuzzy_inner_join(
    cull, edna_agg,
    by = c("Reef" = "Reef", "date_cull" = "date_edna"),
    match_fun = list(
        `==`,
        function(cull_date, edna_date) {
            abs(as.numeric(cull_date - edna_date)) <= win_days
        }
    )
) %>%
    mutate(diff_days = as.numeric(date_cull - date_edna), Reef = Reef.x)

reef_any <- same6_any %>%
    group_by(Reef, Collection.org) %>%
    summarise(
        total_cots      = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
        total_bottom    = sum(Bottomtime, na.rm = TRUE),
        cpue_reef       = total_cots / total_bottom,
        perc_pos_reef   = mean(perc_pos, na.rm = TRUE),
        conc_mean_reef  = mean(conc_mean, na.rm = TRUE),
        n_matches       = n(),
        .groups         = "drop"
    )
get_prior_cohort <- function(culls, ednas, win_days, label) {
    res <- fuzzy_inner_join(
        culls, ednas,
        by = c("Reef" = "Reef", "date_cull" = "date_edna"),
        match_fun = list(
            `==`,
            function(cull_date, edna_date) {
                diff <- as.numeric(cull_date - edna_date)
                diff >= 0 & diff <= win_days
            }
        )
    ) %>% mutate(diff_days = as.numeric(date_cull - date_edna), Reef = Reef.x)

    res %>%
        group_by(Reef, Collection.org, Year, grp) %>%
        summarise(
            horizon         = label,
            total_cots      = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
            total_bottom    = sum(Bottomtime, na.rm = TRUE),
            cpue_reef       = total_cots / total_bottom,
            perc_pos_reef   = mean(perc_pos, na.rm = TRUE),
            conc_mean_reef  = mean(conc_mean, na.rm = TRUE),
            n_matches       = n(),
            .groups         = "drop"
        )
}

reef_prior_3m <- get_prior_cohort(cull, edna_agg, 91, "3 Months")
reef_prior_6m <- get_prior_cohort(cull, edna_agg, 183, "6 Months")
reef_prior_12m <- get_prior_cohort(cull, edna_agg, 365, "12 Months")

reef_before <- bind_rows(reef_prior_3m, reef_prior_6m, reef_prior_12m) %>%
    mutate(horizon = factor(horizon, levels = c("3 Months", "6 Months", "12 Months")))
p1 <- ggplot(reef_before, aes(x = perc_pos_reef, y = cpue_reef, color = horizon)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    labs(x = "% eDNA positive", y = "CPUE", title = "Prior eDNA (% Pos) vs CPUE") +
    facet_wrap(~horizon)

p2 <- ggplot(reef_before, aes(x = log1p(conc_mean_reef), y = cpue_reef, color = horizon)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    labs(x = "log1p(Mean Concentration)", y = "CPUE", title = "Prior eDNA (Conc) vs CPUE") +
    facet_wrap(~horizon)

p1 / p2
# Filter outliers systematically at modeling onset
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

effort_ref <- median(dat_glmm$total_bottom, na.rm = TRUE)
fit_horizon_glmm <- function(hor_val) {
    dat_sub <- dat_glmm %>% filter(horizon == hor_val)
    if (nrow(dat_sub) == 0) {
        return(NULL)
    }

    m_nb <- glmmTMB(
        total_cots ~ perc_pos_reef + offset(log(total_bottom)) + (1 | Reef),
        family = nbinom2(),
        data = dat_sub,
        control = glmmTMBControl(optCtrl = list(iter.max = 1e3, eval.max = 1e3))
    )

    nd <- tibble(
        perc_pos_reef = seq(min(dat_sub$perc_pos_reef), max(dat_sub$perc_pos_reef), length.out = 100),
        total_bottom = median(dat_sub$total_bottom, na.rm = TRUE),
        Reef = NA,
        horizon = hor_val
    )

    pred <- predict(m_nb, newdata = nd, type = "link", se.fit = TRUE, re.form = NA)
    nd %>%
        mutate(
            fit_cpue  = exp(pred$fit) / total_bottom,
            lcl_cpue  = exp(pred$fit - 1.96 * pred$se.fit) / total_bottom,
            ucl_cpue  = exp(pred$fit + 1.96 * pred$se.fit) / total_bottom
        )
}

horizons <- levels(dat_glmm$horizon)
nd_all <- map_dfr(horizons, fit_horizon_glmm) %>%
    mutate(horizon = factor(horizon, levels = horizons))

ggplot(dat_glmm, aes(x = perc_pos_reef, y = obs_cpue)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_all, aes(x = perc_pos_reef, ymin = lcl_cpue, ymax = ucl_cpue, fill = horizon), alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = nd_all, aes(x = perc_pos_reef, y = fit_cpue, color = horizon), linewidth = 1, inherit.aes = FALSE) +
    facet_wrap(~horizon) +
    labs(x = "% eDNA positive", y = "CPUE", title = "NegBin GLM Offset Modeling: % Pos Degradation") +
    theme(legend.position = "none")
m_nb_base <- glmmTMB(
    total_cots ~ perc_pos_reef + offset(log(total_bottom)) + (1 | Reef),
    family = nbinom2(),
    data = dat_glmm
)

m_nb_both <- glmmTMB(
    total_cots ~ perc_pos_reef * conc_t + offset(log(total_bottom)) + (1 | Reef),
    family = nbinom2(),
    data = dat_glmm
)

# Interaction Anova Validation
anova(m_nb_base, m_nb_both)
g_perc <- gam(
    total_cots ~ s(perc_pos_reef, k = 3) + s(Reef, bs = "re") + offset(log(total_bottom)),
    family = nb(), method = "REML",
    data = dat_glmm
)

nd_gam <- tibble(
    perc_pos_reef = seq(min(dat_glmm$perc_pos_reef), max(dat_glmm$perc_pos_reef), length.out = 200),
    total_bottom = effort_ref,
    Reef = dat_glmm$Reef[1]
)

pp_perc <- predict(g_perc, newdata = nd_gam, type = "link", se.fit = TRUE, exclude = "s(Reef)")
nd_gam <- nd_gam %>%
    mutate(
        fit_cpue  = exp(pp_perc$fit) / total_bottom,
        lcl_cpue  = exp(pp_perc$fit - 1.96 * pp_perc$se.fit) / total_bottom,
        ucl_cpue  = exp(pp_perc$fit + 1.96 * pp_perc$se.fit) / total_bottom
    )

ggplot(dat_glmm, aes(x = perc_pos_reef, y = obs_cpue)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0.04, linetype = "dashed", color = "purple") +
    geom_ribbon(data = nd_gam, aes(x = perc_pos_reef, ymin = lcl_cpue, ymax = ucl_cpue), alpha = 0.3, fill = "skyblue", inherit.aes = FALSE) +
    geom_line(data = nd_gam, aes(x = perc_pos_reef, y = fit_cpue), linewidth = 1, color = "black", inherit.aes = FALSE) +
    labs(title = "GAM Spline NegBin Estimate: % Pos")
# Set CPUE Management Threshold
cpue_thresh <- 0.04

dat_evt <- dat_glmm %>%
    rename(reef = Reef, counts = total_cots, cpue = obs_cpue, perc_pos = perc_pos_reef) %>%
    dplyr::select(reef, Collection.org, perc_pos, cpue)

res_all <- make_confusion(dat_evt, perc_thresh = 48, cpue_thresh = cpue_thresh)
plot_confusion(res_all, title = sprintf("Full Data Matrix (eDNA ≥ 48%%, CPUE ≥ %.3f)", cpue_thresh))
# Ensure the evt data contains the required columns including conc_t
dat_evt <- dat_glmm %>%
    rename(reef = Reef, counts = total_cots, cpue = obs_cpue, perc_pos = perc_pos_reef) %>%
    mutate(conc_t = log1p(conc_mean_reef)) %>%
    dplyr::select(reef, Collection.org, perc_pos, conc_t, conc_mean_reef, cpue)

# 1) Evaluate on a grid of %pos AND cpue thresholds (Heatmap)
perc_grid <- seq(0, 100, by = 2)
cpue_grid <- seq(0, quantile(dat_evt$cpue, 0.99, na.rm = TRUE), length.out = 40)

grid <- tidyr::expand_grid(perc_thresh = perc_grid, cpue_thresh = cpue_grid) %>%
    mutate(
        out = purrr::pmap(., ~ metrics_for(..1, ..2, dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1"),
        Accuracy = vapply(out, `[[`, numeric(1), "Accuracy")
    ) %>%
    dplyr::select(-out)

# F1 Heatmap
p_hm <- ggplot(grid, aes(x = perc_thresh, y = cpue_thresh, fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 1, by = 0.2), colour = "black", linewidth = 0.3, alpha = 0.7) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = "eDNA % positive threshold", y = "CPUE threshold", fill = "F1", title = "F1 across %pos thresholds with iso-F1 contours")

# 2) Cross validate across multiple CPUE levels
cpue_levels <- c(0.01, 0.02, 0.04, 0.06, 0.08, 0.10)
k <- 5
R <- 100 # Reduced from 200 for rendering speed
set.seed(1)

run_cv_for_cpue <- function(dat_evt, cpue_thr, perc_grid, k = 5, R = 100) {
    dat2 <- dat_evt %>% mutate(actual = if (cpue_thr == 0) (cpue > 0) else (cpue >= cpue_thr))
    strat_ok <- length(unique(dat2$actual)) > 1
    folds <- if (strat_ok) vfold_cv(dat2, v = k, repeats = R, strata = actual) else vfold_cv(dat2, v = k, repeats = R)

    cv_long <- folds %>%
        mutate(assess = map(splits, assessment)) %>%
        dplyr::select(id, id2, assess) %>%
        tidyr::expand_grid(perc_thresh = perc_grid) %>%
        mutate(
            out = pmap(list(perc_thresh, assess), \(p, d) metrics_for(p, cpue_thr, d)),
            F1 = vapply(out, `[[`, numeric(1), "F1"),
            Accuracy = vapply(out, `[[`, numeric(1), "Accuracy"),
            Precision = vapply(out, `[[`, numeric(1), "Precision"),
            Recall = vapply(out, `[[`, numeric(1), "Recall")
        ) %>%
        dplyr::select(-out, -assess)

    cv_sum <- cv_long %>%
        group_by(perc_thresh) %>%
        summarise(
            F1_mean = mean(F1, na.rm = TRUE), F1_sd = sd(F1, na.rm = TRUE), F1_se = F1_sd / sqrt(n()), .groups = "drop"
        )

    best <- cv_sum %>% slice_max(F1_mean, n = 1, with_ties = FALSE)
    p_star <- best$perc_thresh

    perf_split <- folds %>%
        mutate(
            assess = map(splits, assessment),
            m = map(assess, ~ {
                out <- metrics_for(p_star, cpue_thr, .x)
                tibble(F1 = out["F1"], Accuracy = out["Accuracy"], Precision = out["Precision"], Recall = out["Recall"])
            })
        ) %>%
        dplyr::select(m) %>%
        unnest(m)

    perf_sum <- perf_split %>%
        summarise(
            cpue_thr = cpue_thr, perc_star = p_star,
            F1_mean = mean(F1), F1_sd = sd(F1), F1_se = F1_sd / sqrt(n()),
            Acc_mean = mean(Accuracy), Acc_sd = sd(Accuracy), Acc_se = Acc_sd / sqrt(n()),
            Prec_mean = mean(Precision), Prec_sd = sd(Precision), Prec_se = Prec_sd / sqrt(n()),
            Rec_mean = mean(Recall), Rec_sd = sd(Recall), Rec_se = Rec_sd / sqrt(n()),
            .groups = "drop"
        )

    cm_full <- make_confusion(dat2 %>% dplyr::select(perc_pos, cpue), p_star, cpue_thr)

    list(cpue_thr = cpue_thr, perc_star = p_star, cv_sum = cv_sum, perf_sum = perf_sum, cm_full = cm_full)
}

results <- map(cpue_levels, ~ run_cv_for_cpue(dat_evt, .x, perc_grid, k = k, R = R))

summary_table <- map_dfr(results, "perf_sum") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(cpue_thr, perc_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

# Show the optimal threshold summary
summary_table

# Plot multi-panel confusion matrices
cm_all <- map_dfr(results, function(r) {
    r$cm_full$cm %>% mutate(cpue_thr = r$cpue_thr, perc_star = r$perc_star)
}) %>% mutate(panel = paste0("CPUE ≥ ", cpue_thr, "\n%pos* = ", perc_star))

p_cm <- ggplot(cm_all, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
    labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = "Confusion matrices with CV-tuned thresholds") +
    coord_fixed()

p_hm
p_cm
cpue_target <- 0.02

# Ensure we retain horizon data during formatting
dat_evt_hor <- dat_glmm %>%
    rename(reef = Reef, counts = total_cots, cpue = obs_cpue, perc_pos = perc_pos_reef) %>%
    dplyr::select(reef, Collection.org, perc_pos, cpue, horizon)

run_cv_for_horizon <- function(dat_full, hor_val, cpue_thr, perc_grid, k = 5, R = 100) {
    dat_sub <- dat_full %>% filter(horizon == hor_val)
    res <- run_cv_for_cpue(dat_sub, cpue_thr, perc_grid, k = k, R = R)
    res$perf_sum <- res$perf_sum %>% mutate(horizon = hor_val)
    res$cm_full$cm <- res$cm_full$cm %>% mutate(horizon = hor_val)
    res
}

horizons <- levels(dat_evt_hor$horizon)
results_hor <- map(horizons, ~ run_cv_for_horizon(dat_evt_hor, .x, cpue_target, perc_grid, k = k, R = R))

summary_hor <- map_dfr(results_hor, "perf_sum") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(horizon, cpue_thr, perc_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

summary_hor

cm_all_hor <- map_dfr(results_hor, function(r) {
    r$cm_full$cm %>% mutate(cpue_thr = r$cpue_thr, perc_star = r$perc_star)
}) %>% mutate(panel = factor(paste0(horizon, "\n%pos* = ", perc_star), levels = paste0(levels(dat_evt_hor$horizon), "\n%pos* = ", map_dbl(results_hor, "perc_star"))))

p_cm_hor <- ggplot(cm_all_hor, aes(x = pred, y = actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
    labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = sprintf("Confusion Matrices across Horizons (CPUE ≥ %.2f)", cpue_target)) +
    coord_fixed()

p_cm_hor
metrics_for_both <- function(p_thr, c_thr, cpue_thr, rule = c("or", "and"), x) {
    rule <- match.arg(rule)
    pred <- if (rule == "or") {
        (x$perc_pos >= p_thr) | (x$conc_t >= c_thr)
    } else {
        (x$perc_pos >= p_thr) & (x$conc_t >= c_thr)
    }
    actual <- x$cpue >= cpue_thr
    TP <- sum(pred & actual)
    FP <- sum(pred & !actual)
    TN <- sum(!pred & !actual)
    FN <- sum(!pred & actual)
    prec <- TP / (TP + FP + 1e-9)
    rec <- TP / (TP + FN + 1e-9)
    f1 <- 2 * prec * rec / (prec + rec + 1e-9)
    acc <- (TP + TN) / (TP + TN + FP + FN + 1e-9)
    c(F1 = f1, Accuracy = acc, Precision = prec, Recall = rec)
}

# Tune on transformed concentration to avoid skewness clustering 0s all the way across
conc_t_max <- quantile(dat_evt$conc_t, 0.99, na.rm = TRUE)
conc_grid_t <- seq(min(dat_evt$conc_t, na.rm = TRUE), conc_t_max, length.out = 30)

# Evaluate OR and AND grids over 2D thresholds for a baseline cpue
cpue_star <- 0.04

grid_both_and <- tidyr::expand_grid(perc_thresh = seq(0, 100, by = 2), conc_thresh_t = conc_grid_t) %>%
    mutate(
        out = purrr::pmap(., ~ metrics_for_both(..1, ..2, cpue_star, "and", dat_evt)),
        F1 = vapply(out, `[[`, numeric(1), "F1")
    )

p_2dhm <- ggplot(grid_both_and, aes(x = perc_thresh, y = expm1(conc_thresh_t), fill = F1)) +
    geom_raster() +
    geom_contour(aes(z = F1), breaks = seq(0.2, 0.9, by = 0.2), colour = "black", linewidth = 0.3) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(x = "% eDNA Positive Threshold", y = "Concentration Threshold (pseudo-log)", title = sprintf("2D Boolean Search: AND Rule F1 (CPUE ≥ %.3f)", cpue_star))


# Full 2D Loop over CPUE levels
run_cv_2d <- function(dat, cpue_thr, perc_grid, conc_grid_t, rule = "and") {
    best_pair <- function(df, cpue_thr) {
        actual <- df$cpue >= cpue_thr
        P <- length(perc_grid)
        C <- length(conc_grid_t)
        perc_ge <- outer(df$perc_pos, perc_grid, `>=`)
        conc_ge <- outer(df$conc_t, conc_grid_t, `>=`)
        f1_mat <- matrix(NA_real_, nrow = P, ncol = C)
        for (pi in seq_len(P)) {
            pred_mat <- if (rule == "and") conc_ge & perc_ge[, pi] else conc_ge | perc_ge[, pi]
            TP <- colSums(pred_mat & actual)
            FP <- colSums(pred_mat & !actual)
            FN <- colSums(!pred_mat & actual)
            prec <- TP / (TP + FP + 1e-9)
            rec <- TP / (TP + FN + 1e-9)
            f1_mat[pi, ] <- 2 * prec * rec / (prec + rec + 1e-9)
        }
        max_f1 <- max(f1_mat, na.rm = TRUE)
        idx <- which(f1_mat == max_f1, arr.ind = TRUE)
        idx <- idx[order(idx[, 1], idx[, 2], decreasing = TRUE), , drop = FALSE][1, ]
        list(perc_star = perc_grid[idx[1]], conc_t_star = conc_grid_t[idx[2]])
    }

    dat2 <- dat %>% mutate(actual = cpue >= cpue_thr)
    if (length(unique(dat2$actual)) < 2) {
        return(NULL)
    }

    folds <- vfold_cv(dat2, v = k, repeats = R, strata = actual)

    split_res <- map_dfr(folds$splits, function(spl) {
        tr <- analysis(spl)
        te <- assessment(spl)
        best <- best_pair(tr, cpue_thr)
        met <- metrics_for_both(best$perc_star, best$conc_t_star, cpue_thr, rule, te)
        tibble(perc_star = best$perc_star, conc_t_star = best$conc_t_star, conc_mean_star = pmax(expm1(best$conc_t_star), 0)) %>%
            bind_cols(tibble(F1 = met["F1"], Accuracy = met["Accuracy"], Precision = met["Precision"], Recall = met["Recall"]))
    })

    thr <- split_res %>% summarise(perc_star = median(perc_star), conc_t_star = median(conc_t_star), conc_mean_star = median(conc_mean_star), .groups = "drop")

    perf <- split_res %>%
        summarise(
            cpue_thr = cpue_thr, perc_star = thr$perc_star, conc_mean_star = thr$conc_mean_star,
            F1_mean = mean(F1), F1_se = sd(F1) / sqrt(n()),
            Acc_mean = mean(Accuracy), Acc_se = sd(Accuracy) / sqrt(n()),
            Prec_mean = mean(Precision), Prec_se = sd(Precision) / sqrt(n()),
            Rec_mean = mean(Recall), Rec_se = sd(Recall) / sqrt(n()), .groups = "drop"
        )

    list(thr = thr, perf = perf)
}

results_2d <- map(cpue_levels, ~ run_cv_2d(dat_evt, .x, seq(0, 100, by = 5), conc_grid_t, rule = "and")) # %>% discard(is.null)

summary_table_2d <- map_dfr(results_2d, "perf") %>%
    mutate(
        F1_mean_se   = sprintf("%.3f ± %.3f", F1_mean, F1_se),
        Acc_mean_se  = sprintf("%.3f ± %.3f", Acc_mean, Acc_se),
        Prec_mean_se = sprintf("%.3f ± %.3f", Prec_mean, Prec_se),
        Rec_mean_se  = sprintf("%.3f ± %.3f", Rec_mean, Rec_se)
    ) %>%
    dplyr::select(cpue_thr, perc_star, conc_mean_star, F1_mean_se, Acc_mean_se, Prec_mean_se, Rec_mean_se)

# Show operational 2D thresholds
summary_table_2d

# Plot 2D Multi-CM
make_confusion_2d <- function(dat, perc_star, conc_t_star, cpue_thr, rule = "and") {
    classified <- dat %>%
        mutate(
            pred = if_else(perc_pos >= perc_star & conc_t >= conc_t_star, "Pred +", "Pred -"),
            actual = if_else(cpue >= cpue_thr, "Actual +", "Actual -")
        )
    cm <- classified %>%
        count(actual, pred, .drop = FALSE) %>%
        complete(actual = c("Actual +", "Actual -"), pred = c("Pred +", "Pred -"), fill = list(n = 0)) %>%
        group_by(actual) %>%
        mutate(row_prop = n / sum(n)) %>%
        ungroup() %>%
        mutate(label = paste0(n, "\n", scales::percent(row_prop, accuracy = 0.1)))
    list(cm = cm)
}

thr_df <- map_dfr(results_2d, "thr") %>% mutate(cpue_thr = sapply(results_2d, function(x) x$perf$cpue_thr))
cm_all_2d <- pmap_dfr(list(thr_df$cpue_thr, thr_df$perc_star, thr_df$conc_t_star), function(cp, p, c) {
    obj <- make_confusion_2d(dat_evt, p, c, cp)
    obj$cm %>% mutate(panel = paste0("CPUE ≥ ", cp, "\n%pos*=", round(p, 1), "  conc*≈", round(pmax(expm1(c), 0), 1)))
})

p_2dcm <- ggplot(cm_all_2d, aes(pred, actual, fill = row_prop)) +
    geom_tile(color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4, color = "black") +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, limits = c(0, 1), name = "Row Rate") +
    labs(x = "Prediction (2D rule)", y = "Ground truth (CPUE)", title = "2D AND Rule Confusion Matrices") +
    coord_fixed()

p_2dhm
p_2dcm


# Safely load ggplot and export defined variables
library(ggplot2)
dir.create("plots", showWarnings=FALSE)
if(exists("p_hm")) ggsave("plots/temporal_degradation_GLM.png", plot = p_hm, width=10, height=8, dpi=300)
if(exists("p_cm")) ggsave("plots/multi_horizon_CPUE.png", plot = p_cm, width=10, height=8, dpi=300)
if(exists("p_cm_hor")) ggsave("plots/cv_horizons_02_CPUE.png", plot = p_cm_hor, width=10, height=4, dpi=300)
if(exists("p_2dhm")) ggsave("plots/2D_AND_F1_Heatmap.png", plot = p_2dhm, width=10, height=8, dpi=300)
if(exists("p_2dcm")) ggsave("plots/2D_AND_Multi_CM.png", plot = p_2dcm, width=10, height=8, dpi=300)
