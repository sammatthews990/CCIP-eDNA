#' Evaluate Performance Metrics
#'
#' @param p_thr Numeric percentage positive threshold
#' @param c_thr Numeric CPUE threshold
#' @param x Dataframe with `perc_pos` and `cpue` columns
#' @return Named numeric vector containing F1, Accuracy, Precision, and Recall.
#' @export
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

#' Metrics for Both
#'
#' @param p_thr Numeric percentage positive threshold
#' @param c_thr Numeric Concentration threshold
#' @param cpue_thr CPUE target
#' @param rule Boolean logic rule: "or" or "and"
#' @param x Dataframe with `perc_pos`, `conc_t` and `cpue`
#' @return Numeric performance metrics
#' @export
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

#' Create Confusion Matrix Object
#'
#' @param dat Data containing `perc_pos` and `cpue`.
#' @param perc_thresh Percentage positive threshold.
#' @param cpue_thresh CPUE target threshold.
#' @return A list containing the confusion matrix data frame (`cm`) and computed metrics.
#' @export
#' @import dplyr tidyr
#' @importFrom scales percent
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

#' Make 2D Confusion Matrix
#'
#' @param dat Dataframe
#' @param perc_star Chosen perc pos threshold
#' @param conc_t_star Chosen conc threshold
#' @param cpue_thr CPUE threshold
#' @param rule "and" or "or"
#' @return List with cm
#' @export
#' @import dplyr tidyr
#' @importFrom scales percent
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

#' Best Pair Threshold
#'
#' @param df Data
#' @param cpue_thr Target CPUE
#' @param perc_grid Grid
#' @param conc_grid_t Grid
#' @param rule logic
#' @return list
#' @export
best_pair <- function(df, cpue_thr, perc_grid, conc_grid_t, rule = "and") {
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

#' Get Threshold wrapper
#'
#' @param cpue Cpue value
#' @return String label for threshold
#' @export
get_thr <- function(cpue) {
    if (cpue == 0.04) "0.04 (Borderline/Needed)" else "0.02 (Wasted/Borderline)"
}
