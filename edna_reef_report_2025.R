# =============================================================================
# edna_reef_report_2025.R
#
# Purpose: Reef eDNA survey status report for surveys since January 2025.
#   1. Re-creates edna_agg (with sample date) from raw eDNA data
#   2. Derives 2D AND-rule CPUE thresholds (perc_pos AND conc_mean) for four
#      CPUE categories, mirroring the logic in eDNA_CPUE_Comparison_final.R
#   3. Classifies each survey event into four expected CPUE categories:
#        < 0.01  |  0.01–0.04  |  0.04–0.08  |  > 0.08  COTS/min
#   4. Filters to surveys since 2025-01-01
#   5. Joins to the DRAFT Target Reef List 2026-27
#   6. Reports reefs NOT on the target list
#   7. Produces an interactive leaflet map (saved as edna_reef_report_2025.html)
#
# Dependencies: dplyr, readxl, tidyr, purrr, stringr, fuzzyjoin, leaflet,
#               htmlwidgets, scales
# =============================================================================

library(dplyr)
library(readxl)
library(tidyr)
library(purrr)
library(stringr)
library(fuzzyjoin)
library(leaflet)
library(htmlwidgets)
library(scales)

# =============================================================================
# 1.  Load raw data
# =============================================================================

edna.dat <- read_excel("data/eDNA data_ALL_20260225.xlsx", sheet = "eDNA_data_ALL")
cull.dat <- read_excel("data/260201_COTS-Cull-Data-Ewels.xlsx", sheet = "Cull")

# =============================================================================
# 2.  Re-create edna_agg (mirrors eDNA_CPUE_Comparison_final.R lines 107-135)
#     — includes date_edna for every aggregated event
# =============================================================================

edna_agg <- edna.dat |>
    filter(!is.na(Year)) |>
    rename(Collection.organisation = `Collection organisation`) |>
    mutate(
        Reef           = ReefName,
        Collection.org = if_else(Collection.organisation == "AIMS", "AIMS", "Other"),
        date_edna      = as.Date(Date),
        Conc_mean      = as.numeric(Conc_mean),
        Lat            = as.numeric(Lat),
        Long           = as.numeric(Long)
    ) |>
    arrange(Reef, Year, date_edna) |>
    group_by(Reef, Collection.org, Year) |>
    mutate(
        grp = cumsum(
            if_else(is.na(lag(date_edna)) | as.numeric(date_edna - lag(date_edna)) > 7,
                1L, 0L
            )
        )
    ) |>
    ungroup() |>
    group_by(Reef, Collection.org, Year, grp) |>
    summarise(
        date_edna = min(date_edna),
        conc_mean = mean(Conc_mean, na.rm = TRUE),
        perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        # Mean coordinates from the raw sample locations for this survey event
        lat_edna  = mean(Lat, na.rm = TRUE),
        lon_edna  = mean(Long, na.rm = TRUE),
        .groups   = "drop"
    )

# =============================================================================
# 3.  Derive 2D AND-rule thresholds from the full matched cull × eDNA dataset
#
#     This re-creates the dat_glmm pipeline from the final script (condensed),
#     then finds optimal (perc_star, conc_t_star) for each CPUE boundary using
#     the same best_pair() matrix search — no slow CV, uses the full dataset.
#     The four category boundaries are: 0.01, 0.04, 0.08 COTS/min
# =============================================================================

# --- 3a. Build cull & 6-month matched dataset (same as final script) ---------
win_days <- 183 # 6-month window

cull <- cull.dat |>
    rename(Reef = ReefName) |>
    mutate(date_cull = as.Date(SurveyDate))

same6_any <- fuzzy_inner_join(
    cull, edna_agg,
    by = c("Reef" = "Reef", "date_cull" = "date_edna"),
    match_fun = list(`==`, function(cd, ed) abs(as.numeric(cd - ed)) <= win_days)
) |>
    mutate(diff_days = as.numeric(date_cull - date_edna), Reef = Reef.x)

dat_glmm <- same6_any |>
    group_by(Reef, Collection.org, Year, grp) |>
    summarise(
        total_cots = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
        total_bottom = sum(Bottomtime, na.rm = TRUE),
        perc_pos_reef = mean(perc_pos, na.rm = TRUE),
        conc_mean_reef = mean(conc_mean, na.rm = TRUE),
        .groups = "drop"
    ) |>
    mutate(
        obs_cpue = total_cots / total_bottom,
        conc_t   = log1p(conc_mean_reef)
    ) |>
    filter(
        !is.na(total_cots), !is.na(total_bottom), total_bottom > 0,
        !is.na(perc_pos_reef), !is.na(conc_t),
        conc_mean_reef < 5000 # remove extreme outliers
    )

# --- 3b. best_pair() — identical to the function in the final script ---------
#         Returns (perc_star, conc_t_star) for a given CPUE threshold

perc_grid_thr <- seq(0, 100, by = 5)
conc_grid_t <- seq(min(dat_glmm$conc_t, na.rm = TRUE),
    quantile(dat_glmm$conc_t, 0.99, na.rm = TRUE),
    length.out = 30
)

best_pair <- function(df, cpue_thr) {
    actual <- df$obs_cpue >= cpue_thr
    if (length(unique(actual)) < 2) {
        return(list(perc_star = NA, conc_t_star = NA))
    }
    P <- length(perc_grid_thr)
    C <- length(conc_grid_t)
    perc_ge <- outer(df$perc_pos_reef, perc_grid_thr, `>=`)
    conc_ge <- outer(df$conc_t, conc_grid_t, `>=`)
    f1_mat <- matrix(NA_real_, nrow = P, ncol = C)
    for (pi in seq_len(P)) {
        pred_mat <- perc_ge[, pi] & conc_ge
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
    list(
        perc_star = perc_grid_thr[idx[1]],
        conc_t_star = conc_grid_t[idx[2]]
    )
}

# --- 3c. Compute thresholds for each CPUE category boundary -----------------
cpue_boundaries <- c(0.01, 0.04, 0.08)

thresh_2d <- map_dfr(cpue_boundaries, function(cthr) {
    bp <- best_pair(dat_glmm, cthr)
    tibble(
        cpue_thr       = cthr,
        perc_star      = bp$perc_star,
        conc_t_star    = bp$conc_t_star,
        conc_mean_star = pmax(expm1(bp$conc_t_star), 0)
    )
})

message("\n── 2D AND-rule thresholds (full-data optimum) ──")
print(thresh_2d)

# --- 3d. Extract the three boundary pairs ------------------------------------
get_thr <- function(cpue) {
    thresh_2d |>
        filter(cpue_thr == cpue) |>
        slice(1)
}

thr_01 <- get_thr(0.01) # boundary between <0.01 and 0.01-0.04
thr_04 <- get_thr(0.04) # boundary between 0.01-0.04 and 0.04-0.08
thr_08 <- get_thr(0.08) # boundary between 0.04-0.08 and >0.08

# =============================================================================
# 4.  Apply 4-category CPUE classification to ALL edna_agg rows
#     (AND rule: perc_pos AND conc_mean must both exceed the threshold)
# =============================================================================

edna_agg_classified <- edna_agg |>
    mutate(
        conc_t = log1p(conc_mean),
        cpue_category = case_when(
            perc_pos >= thr_08$perc_star & conc_t >= thr_08$conc_t_star ~ "> 0.08 CPUE",
            perc_pos >= thr_04$perc_star & conc_t >= thr_04$conc_t_star ~ "0.04–0.08 CPUE",
            perc_pos >= thr_01$perc_star & conc_t >= thr_01$conc_t_star ~ "0.01–0.04 CPUE",
            TRUE ~ "< 0.01 CPUE"
        ),
        cpue_category = factor(cpue_category,
            levels = c(
                "< 0.01 CPUE", "0.01–0.04 CPUE",
                "0.04–0.08 CPUE", "> 0.08 CPUE"
            )
        )
    )

# =============================================================================
# 5.  Filter to surveys since 1 January 2025
# =============================================================================

edna_2025 <- edna_agg_classified |>
    filter(date_edna >= as.Date("2025-01-01"))

message(sprintf(
    "\n── %d survey events since 2025-01-01 across %d reefs ──",
    nrow(edna_2025), length(unique(edna_2025$Reef))
))

# =============================================================================
# 6.  Load Target Reef List and extract coordinates
# =============================================================================

target_list <- read.csv("data/DRAFT Target Reef List 2026-27.csv",
    stringsAsFactors = FALSE,
    check.names      = TRUE # default: spaces → dots; 'Area (ha)' → 'Area..ha.'
) |>
    rename(
        Management_Region = Management.Region,
        Area_ha           = Area..ha.,
        TUMRA             = TUMRA
    ) |>
    mutate(ReefID = as.character(ReefID))

# --- 6a. Extract ReefID from the LAST parenthetical in eDNA reef names -------
#         e.g.  "Lizard Island Reef (North West) (14-116a)"  →  "14-116a"
#         Uses str_extract_all to grab all parenthetical groups, then takes last
edna_2025 <- edna_2025 |>
    mutate(
        ReefID = sapply(
            str_extract_all(Reef, "(?<=\\()[^)]+(?=\\))"),
            function(x) if (length(x) == 0) NA_character_ else x[length(x)]
        )
    )

# --- 6b. Join eDNA 2025 to Target Reef List by ReefID -----------------------
edna_target <- inner_join(edna_2025, target_list,
    by = "ReefID", suffix = c("", ".trl")
) |>
    mutate(in_target_list = TRUE)

edna_not_target <- anti_join(edna_2025, target_list,
    by = "ReefID"
) |>
    mutate(
        in_target_list    = FALSE,
        # Use coordinates from the raw eDNA sample locations
        Longitude         = lon_edna,
        Latitude          = lat_edna,
        Management_Region = "Not on target list",
        TUMRA             = NA_character_
    )

# Target reefs with NO eDNA data since Jan 2025
target_no_edna <- anti_join(target_list,
    edna_2025 |> filter(!is.na(ReefID)),
    by = "ReefID"
) |>
    mutate(
        cpue_category = NA_character_,
        date_edna = NA,
        perc_pos = NA_real_,
        conc_mean = NA_real_,
        n_samples = NA_integer_,
        in_target_list = TRUE
    )

# =============================================================================
# 7.  Table of reefs surveyed but NOT in the target list
# =============================================================================

message("\n── Reefs with eDNA data (Jan 2025+) NOT in 2026-27 target list ──")
non_target_table <- edna_not_target |>
    select(Reef, date_edna, perc_pos, conc_mean, n_samples, cpue_category) |>
    arrange(date_edna)
print(non_target_table, n = Inf)

# =============================================================================
# 8.  Leaflet interactive map
# =============================================================================

# Colour palette for the 4 CPUE categories
cpue_levels <- c("< 0.01 CPUE", "0.01–0.04 CPUE", "0.04–0.08 CPUE", "> 0.08 CPUE")
cpue_colours <- c("#3288BD", "#FEE08B", "#FC8D59", "#D53E4F") # blue→yellow→orange→red
pal_cpue <- colorFactor(palette = cpue_colours, levels = cpue_levels, ordered = TRUE)

# --- Helper: build popup HTML -------------------------------------------------
make_popup <- function(reef, date_edna, perc_pos, conc_mean, n_samples,
                       cpue_category, region = "", tumra = "") {
    paste0(
        "<b>", reef, "</b><br>",
        ifelse(!is.na(region) & region != "", paste0("Region: ", region, "<br>"), ""),
        ifelse(!is.na(tumra) & tumra != "", paste0("TUMRA: ", tumra, "<br>"), ""),
        "Survey date: ", format(date_edna, "%d %b %Y"), "<br>",
        "% eDNA positive: ", round(perc_pos, 1), "%<br>",
        "Mean concentration: ", round(conc_mean, 2), " copies/rxn<br>",
        "Samples: ", n_samples, "<br>",
        "<b>Expected CPUE: ", cpue_category, "</b>"
    )
}

# --- 8a. eDNA target reefs (have lat/lon from target list) -------------------
target_popups <- unname(with(
    edna_target,
    mapply(make_popup, ReefName, date_edna, perc_pos, conc_mean,
        n_samples, as.character(cpue_category), Management_Region, TUMRA,
        SIMPLIFY = TRUE
    )
))

# --- 8b. Target reefs with NO eDNA data (grey) --------------------------------
# Filter to only rows with coordinates (must match what addCircleMarkers receives)
target_no_edna_map <- target_no_edna |> filter(!is.na(Longitude))

no_edna_popups <- with(
    target_no_edna_map,
    paste0(
        "<b>", ReefName, "</b><br>",
        "Region: ", Management_Region, "<br>",
        ifelse(!is.na(TUMRA) & TUMRA != "", paste0("TUMRA: ", TUMRA, "<br>"), ""),
        "<i>No eDNA data since Jan 2025</i>"
    )
)

# --- 8c. Build map -----------------------------------------------------------
# Build popup and data for non-target reefs (use eDNA-derived coords)
edna_not_target_map <- edna_not_target |> filter(!is.na(Longitude))

non_target_popups <- unname(with(
    edna_not_target_map,
    mapply(make_popup, Reef, date_edna, perc_pos, conc_mean,
        n_samples, as.character(cpue_category),
        rep("", nrow(edna_not_target_map)),   # no region
        rep("", nrow(edna_not_target_map)),   # no TUMRA
        SIMPLIFY = TRUE
    )
))

# Map marker colours to CPUE categories for non-target star markers
# Use the integer factor level to avoid any en-dash encoding mismatch
edna_not_target_map <- edna_not_target_map |>
    mutate(
        marker_col = case_when(
            as.integer(cpue_category) == 1L ~ "blue",      # < 0.01
            as.integer(cpue_category) == 2L ~ "orange",    # 0.01-0.04
            as.integer(cpue_category) == 3L ~ "red",       # 0.04-0.08
            as.integer(cpue_category) == 4L ~ "darkred",   # > 0.08
            TRUE                            ~ "gray"
        )
    )

non_target_icons <- awesomeIcons(
    icon        = "star",
    iconColor   = "white",
    library     = "fa",
    markerColor = edna_not_target_map$marker_col
)

map <- leaflet() |>
    addProviderTiles(
        providers$Esri.WorldImagery,
        group   = "ESRI Imagery",
        options = providerTileOptions(maxZoom = 18)
    ) |>
    setView(lng = 147, lat = -18, zoom = 6) |>

    # Layer 1: Target reefs – NO eDNA data (light grey, small circles)
    addCircleMarkers(
        data        = target_no_edna_map,
        lng         = ~Longitude,
        lat         = ~Latitude,
        radius      = 5,
        color       = "#ffffff",
        fillColor   = "#cccccc",
        fillOpacity = 0.6,
        weight      = 1,
        popup       = no_edna_popups,
        group       = "Target reefs (no eDNA since Jan 2025)"
    ) |>

    # Layer 2: Target reefs WITH eDNA data (coloured circles by CPUE category)
    addCircleMarkers(
        data        = edna_target,
        lng         = ~Longitude,
        lat         = ~Latitude,
        radius      = 9,
        color       = "#ffffff",
        fillColor   = ~pal_cpue(cpue_category),
        fillOpacity = 0.9,
        weight      = 1.5,
        popup       = target_popups,
        group       = "Target reefs (eDNA surveyed Jan 2025+)"
    ) |>

    # Layer 3: Non-target reefs WITH eDNA data (star pin markers, CPUE-coloured)
    addAwesomeMarkers(
        data   = edna_not_target_map,
        lng    = ~Longitude,
        lat    = ~Latitude,
        icon   = non_target_icons,
        popup  = non_target_popups,
        group  = "Non-target reefs (eDNA surveyed Jan 2025+)"
    ) |>

    # Legend
    addLegend(
        position = "bottomright",
        colors   = c(cpue_colours, "#cccccc"),
        labels   = c(cpue_levels, "No eDNA (target reef)"),
        title    = "Expected CPUE<br><small>Circles = target | Stars = non-target</small>",
        opacity  = 1
    ) |>

    # Layer control
    addLayersControl(
        overlayGroups = c(
            "Target reefs (eDNA surveyed Jan 2025+)",
            "Target reefs (no eDNA since Jan 2025)",
            "Non-target reefs (eDNA surveyed Jan 2025+)"
        ),
        options = layersControlOptions(collapsed = FALSE)
    )

# =============================================================================
# 9.  Save map
# =============================================================================

saveWidget(map, file = "edna_reef_report_2025.html", selfcontained = FALSE)
message("\n── Map saved to: edna_reef_report_2025.html ──")

# =============================================================================
# 10. Summary counts
# =============================================================================

message("\n── Summary ──")
message(sprintf("  Target reefs with eDNA data (Jan 2025+): %d", nrow(edna_target)))
message(sprintf("  Target reefs WITHOUT eDNA data:           %d", nrow(target_no_edna)))
message(sprintf("  Non-target reefs with eDNA data:          %d", nrow(edna_not_target)))

message("\n── CPUE category breakdown (target reefs with eDNA) ──")
print(edna_target |> count(cpue_category))

# =============================================================================
# 11. Export combined reef CSV
# =============================================================================

# Standardise columns across all three groups and bind
common_cols <- c("Reef", "ReefID", "Longitude", "Latitude",
                 "date_edna", "perc_pos", "conc_mean", "n_samples",
                 "cpue_category", "category")

reef_csv <- bind_rows(
    edna_target |>
        mutate(Reef = ReefName, category = "Target reef with eDNA") |>
        select(any_of(common_cols)),
    target_no_edna |>
        mutate(Reef = ReefName, category = "Target reef (no eDNA)") |>
        select(any_of(common_cols)),
    edna_not_target |>
        mutate(category = "Non-target reef with eDNA") |>
        select(any_of(common_cols))
) |>
    arrange(category, Reef)

write.csv(reef_csv, "edna_reef_report_2025.csv", row.names = FALSE)
message(sprintf("\n── CSV saved: edna_reef_report_2025.csv (%d rows) ──", nrow(reef_csv)))

# Clean up helper objects
rm(same6_any, cull.dat, edna.dat)
