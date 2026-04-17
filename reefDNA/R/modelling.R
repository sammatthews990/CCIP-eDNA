#' Fit Horizon GLMM
#'
#' @param hor_val Horizon Value ("3 Months", etc)
#' @param dat_glmm Dataset containing modelling inputs.
#' @return Predictions tibble
#' @export
#' @import dplyr tidyr glmmTMB
#' @importFrom stats median pf predict residuals sd
fit_horizon_glmm <- function(hor_val, dat_glmm) {
    dat_sub <- dat_glmm %>% filter(horizon == hor_val)
    if (nrow(dat_sub) == 0) return(NULL)

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

#' Run CV for CPUE threshold
#'
#' @param dat_evt Event dataset
#' @param cpue_thr CPUE Threshold
#' @param perc_grid Grid
#' @param k Folds
#' @param R Repeats
#' @return List of performance
#' @export
#' @import dplyr rsample purrr
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

#' Run CV for Horizon
#'
#' @param dat_full Full dataset
#' @param hor_val Horizon
#' @param cpue_thr Target cpue
#' @param perc_grid Grid
#' @param k folds
#' @param R repeats
#' @return List of performance
#' @export
#' @import dplyr tidyr
run_cv_for_horizon <- function(dat_full, hor_val, cpue_thr, perc_grid, k = 5, R = 100) {
    dat_sub <- dat_full %>% filter(horizon == hor_val)
    res <- run_cv_for_cpue(dat_sub, cpue_thr, perc_grid, k = k, R = R)
    res$perf_sum <- res$perf_sum %>% mutate(horizon = hor_val)
    res$cm_full$cm <- res$cm_full$cm %>% mutate(horizon = hor_val)
    res
}

#' Run 2D Cross Validation
#'
#' @param dat Dataset
#' @param cpue_thr CPUE target
#' @param perc_grid Grid 1
#' @param conc_grid_t Grid 2
#' @param rule Boolean rule
#' @param k Folds
#' @param R repeats
#' @return list
#' @export
#' @import dplyr tidyr rsample purrr
run_cv_2d <- function(dat, cpue_thr, perc_grid, conc_grid_t, rule = "and", k = 5, R = 100) {
    dat2 <- dat %>% mutate(actual = cpue >= cpue_thr)
    if (length(unique(dat2$actual)) < 2) return(NULL)

    folds <- vfold_cv(dat2, v = k, repeats = R, strata = actual)

    split_res <- map_dfr(folds$splits, function(spl) {
        tr <- analysis(spl)
        te <- assessment(spl)
        best <- best_pair(tr, cpue_thr, perc_grid, conc_grid_t, rule)
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

#' Get Prior Cohort
#'
#' @param culls Cull dataset
#' @param ednas eDNA dataset
#' @param win_days Window
#' @param label Label string
#' @return Data frame
#' @export
#' @import dplyr tidyr fuzzyjoin
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

    # Define group variables based on availability
    group_vars <- intersect(c("Reef", "Collection.org", "Year", "grp"), names(res))

    res %>%
        group_by(across(all_of(group_vars))) %>%
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

#' Evaluate Unified Spatial
#'
#' @param df Spatial data frame
#' @param cp CPUE target
#' @param perc_grid Grid array
#' @return list
#' @export
#' @import dplyr tidyr
evaluate_unified_spatial <- function(df, cp, perc_grid) {
    cms_all <- list()
    for (r in levels(df$radius)) {
        sub_df <- df %>% filter(radius == r)
        if (nrow(sub_df) < 5) next

        # Calculate Single Global Threshold for the combined Radius
        out <- run_cv_for_cpue(sub_df, cp, perc_grid, k = 5, R = 50)
        if (is.null(out)) next

        opt_thresh <- out$perc_star

        # Apply that same unified threshold across both Regimes individually
        for (reg in unique(sub_df$regime)) {
            reg_df <- sub_df %>% filter(regime == reg)
            if (nrow(reg_df) > 0) {
                cm_reg <- make_confusion(reg_df, perc_thresh = opt_thresh, cpue_thresh = cp)
                cm_reg$cm <- cm_reg$cm %>% mutate(radius = r, regime = reg, cpue_thr = cp, perc_star = opt_thresh)
                cms_all[[paste0(r, "_", reg)]] <- cm_reg$cm
            }
        }
    }

    if (length(cms_all) > 0) {
        df_cm <- bind_rows(cms_all)
        lbl_order <- df_cm %>%
            dplyr::select(radius, perc_star) %>%
            distinct() %>%
            arrange(radius) %>%
            mutate(lbl = paste0(radius, "\nGlobal %pos* = ", perc_star)) %>%
            pull(lbl)

        df_cm <- df_cm %>%
            mutate(radius_lbl = factor(paste0(radius, "\nGlobal %pos* = ", perc_star), levels = lbl_order))

        return(list(cm = df_cm))
    }
    return(NULL)
}
