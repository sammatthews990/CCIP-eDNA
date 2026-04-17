#' Format Linear Model Label for Plotting
#'
#' @param model A fitted linear model object.
#' @param digits_r2 Digits for R squared
#' @param digits_p Digits for p-value
#' @param digits_rmse Digits for RMSE
#' @return Formatted string label
#' @export
#' @import broom
lm_label <- function(model, digits_r2 = 2, digits_p = 3, digits_rmse = 2) {
    g <- broom::glance(model)
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

#' Plot Confusion Matrix
#'
#' @param cm_obj Confusion matrix list object
#' @param title Plot title
#' @return ggplot object
#' @export
#' @import ggplot2
plot_confusion <- function(cm_obj, title = "Confusion matrix") {
    ggplot(cm_obj$cm, aes(x = pred, y = actual, fill = row_prop)) +
        geom_tile(color = "grey85", linewidth = 0.4) +
        geom_text(aes(label = label), size = 5) +
        scale_fill_viridis_c(option = "mako", direction = -1, begin = 0.25, name = "Count Rate") +
        labs(x = "Prediction from eDNA (% pos)", y = "Ground truth from CPUE", title = title) +
        coord_fixed()
}

#' Build Leaflet Popup Html
#'
#' @param reef Reef name
#' @param year Year
#' @param edna_val eDNA metric
#' @param scar_val Scar metric
#' @param cpue_val CPUE metric
#' @param truth Logical true outcome
#' @param conf_class Confusion classification string
#' @param method Description of prediction method
#' @return String HTML representation
#' @export
build_popup <- function(reef, year, edna_val, scar_val, cpue_val, truth, conf_class, method) {
    paste0(
        "<b>Reef: ", reef, " (", year, ")</b><br>",
        "<i>Method: ", method, "</i><hr style='margin:4px 0'>",
        "<b>eDNA % Positive: </b>", round(edna_val, 1), "%<br>",
        "<b>Manta Scars: </b>", round(scar_val, 1), "%<br>",
        "<b>Initial CPUE: </b>", round(cpue_val, 3), " COTS/min<br>",
        "<hr style='margin:4px 0'>",
        "Truth Required Culling: <b>", truth, "</b><br>",
        "Outcome: <b>", conf_class, "</b>"
    )
}

#' Make Leaflet Popup Html for Droppoints
#'
#' @param reef Reef Name
#' @param date_edna Date
#' @param perc_pos Percent positive
#' @param conc_mean Mean concentration
#' @param n_samples Sample count
#' @return String HTML
#' @export
make_popup <- function(reef, date_edna, perc_pos, conc_mean, n_samples) {
    paste0(
        "<b>Reef: ", reef, "</b><br>",
        "Survey Date: ", format(as.Date(date_edna), "%d %b %Y"), "<br>",
        "<b>% Positive: ", round(perc_pos, 1), "%</b><br>",
        "<b>Mean Conc: ", round(conc_mean, 2), " copies/rxn</b><br>",
        "Samples: ", n_samples
    )
}
