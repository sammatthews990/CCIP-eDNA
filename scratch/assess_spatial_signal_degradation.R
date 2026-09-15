# Script to quantitatively assess spatial signal degradation across buffer radii
# Computes F1, Precision, Recall, AUC, GLMM log-likelihood, and Marginal R2 across 200m - 5000m radii

library(dplyr)
library(ggplot2)
library(readxl)
library(purrr)
library(tidyr)
library(sf)
library(glmmTMB)
library(rsample)
library(patchwork)
library(pROC)
library(tibble)

devtools::load_all("reefDNA")

# Load data
cull_file <- "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file <- "data/eDNA data_ALL_20260528.xlsx"

cull <- read_excel(cull_file, sheet = "Cull") %>%
    rename(Reef = ReefName) %>%
    mutate(date_cull = as.Date(SurveyDate))

edna_site <- read_excel(edna_file, sheet = "eDNA_data_ALL") %>%
    filter(!is.na(Lat), !is.na(Long), !is.na(Year)) %>%
    rename(Collection.organisation = `Collection organisation`) %>%
    mutate(
        Reef = ReefName,
        Collection.org = ifelse(Collection.organisation == "AIMS", "AIMS", "Other"),
        date_edna = as.Date(Date),
        Conc_mean = as.numeric(Conc_mean)
    ) %>%
    group_by(Reef, Site_name, Collection.org, Year) %>%
    summarise(
        Lat = mean(Lat, na.rm = TRUE),
        Long = mean(Long, na.rm = TRUE),
        date_edna = min(date_edna),
        conc_mean = mean(Conc_mean, na.rm = TRUE),
        perc_pos  = mean(LOD_sample_positive, na.rm = TRUE) * 100,
        n_samples = n(),
        .groups   = "drop"
    )

edna_sf <- st_as_sf(edna_site, coords = c("Long", "Lat"), crs = 4326) %>% st_transform(3112)
cull_sf <- cull %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(3112)

buffer_radii <- c(200, 500, 1000, 1500, 2000, 3000, 5000)
cpue_target <- 0.04
set.seed(42)

degradation_results <- map_dfr(buffer_radii, function(r) {
    edna_buf <- st_buffer(edna_sf, dist = r)
    intersect_cull <- st_join(edna_buf, cull_sf, join = st_intersects) %>% filter(!is.na(date_cull))
    
    valid_encounters <- intersect_cull %>%
        mutate(diff_days = as.numeric(date_cull - date_edna)) %>%
        filter(diff_days >= 0 & diff_days <= 183)
    
    site_cpue <- valid_encounters %>%
        st_drop_geometry() %>%
        rename(any_of(c(Reef = "Reef.x", Year = "Year.x"))) %>%
        group_by(Reef, Site_name, Collection.org, Year, perc_pos, conc_mean) %>%
        summarise(
            counts = sum(Cohort1 + Cohort2 + Cohort3 + Cohort4, na.rm = TRUE),
            total_bottom = sum(Bottomtime, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        filter(total_bottom > 0) %>%
        mutate(
            cpue = counts / total_bottom,
            actual = cpue >= cpue_target
        )
    
    N <- nrow(site_cpue)
    N_pos <- sum(site_cpue$actual)
    
    # Tune optimal F1 for this radius
    p_grid <- seq(0, 100, by = 1)
    grid_eval <- map_dfr(p_grid, function(p) {
        m <- metrics_for(p, cpue_target, site_cpue)
        tibble(
            p_star = p, 
            F1 = unname(m["F1"]), 
            Precision = unname(m["Precision"]), 
            Recall = unname(m["Recall"]), 
            Accuracy = unname(m["Accuracy"]), 
            FPR = unname(m["FPR"])
        )
    })
    
    best_row <- grid_eval %>% slice_max(F1, n = 1, with_ties = FALSE)
    
    # Calculate AUC-ROC
    roc_obj <- tryCatch(pROC::roc(site_cpue$actual, site_cpue$perc_pos, quiet = TRUE), error = function(e) NULL)
    auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA
    
    # GLMM fit effect size and AIC
    glmm_fit <- tryCatch(
        glmmTMB(counts ~ perc_pos + offset(log(total_bottom)) + (1|Reef), data = site_cpue, family = nbinom2),
        error = function(e) NULL
    )
    
    beta_perc <- if (!is.null(glmm_fit)) fixef(glmm_fit)$cond["perc_pos"] else NA
    p_val     <- if (!is.null(glmm_fit)) summary(glmm_fit)$coefficients$cond["perc_pos", "Pr(>|z|)"] else NA
    aic_val   <- if (!is.null(glmm_fit)) AIC(glmm_fit) else NA
    
    tibble(
        radius_m = r,
        N = N,
        N_pos = N_pos,
        Best_P_Star = best_row$p_star,
        Max_F1 = best_row$F1,
        Precision = best_row$Precision,
        Recall = best_row$Recall,
        Accuracy = best_row$Accuracy,
        FPR = best_row$FPR,
        AUC = auc_val,
        Beta_perc_pos = beta_perc,
        P_value = p_val,
        AIC = aic_val
    )
})

cat("=== SPATIAL SIGNAL DECAY ANALYSIS SUMMARY ===\n")
print(degradation_results)

# Calculate relative decay from peak performance
peak_f1 <- max(degradation_results$Max_F1, na.rm = TRUE)
peak_auc <- max(degradation_results$AUC, na.rm = TRUE)

degradation_results <- degradation_results %>%
    mutate(
        F1_pct_peak = Max_F1 / peak_f1 * 100,
        AUC_pct_peak = AUC / peak_auc * 100
    )

# Save metrics table
dir.create("plots", showWarnings = FALSE)
write.csv(degradation_results, "plots/spatial_signal_degradation_metrics.csv", row.names = FALSE)

# Generate diagnostic plot
p1 <- ggplot(degradation_results, aes(x = radius_m, y = Max_F1)) +
    geom_line(color = "#2b5c8f", linewidth = 1.2) +
    geom_point(color = "#2b5c8f", size = 3) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.9) +
    annotate("text", x = 1100, y = min(degradation_results$Max_F1, na.rm=TRUE) + 0.02, label = "Optimal Threshold (1000m)", color = "red", hjust = 0, fontface = "bold") +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "A. Max F1 Score vs Buffer Radius", x = "Spatial Buffer Radius", y = "Max F1 Score") +
    theme_bw(base_family = "Helvetica")

p2 <- ggplot(degradation_results, aes(x = radius_m, y = AUC)) +
    geom_line(color = "#2a9d8f", linewidth = 1.2) +
    geom_point(color = "#2a9d8f", size = 3) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.9) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "B. Area Under ROC Curve (AUC) vs Buffer Radius", x = "Spatial Buffer Radius", y = "AUC-ROC") +
    theme_bw(base_family = "Helvetica")

p3 <- ggplot(degradation_results, aes(x = radius_m, y = Beta_perc_pos)) +
    geom_line(color = "#e76f51", linewidth = 1.2) +
    geom_point(color = "#e76f51", size = 3) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.9) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "C. GLMM Effect Size (Beta) vs Buffer Radius", x = "Spatial Buffer Radius", y = "Beta (%pos coefficient)") +
    theme_bw(base_family = "Helvetica")

p4 <- ggplot(degradation_results, aes(x = radius_m, y = Precision)) +
    geom_line(color = "#e63946", linewidth = 1.2) +
    geom_point(color = "#e63946", size = 3) +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 0.9) +
    scale_x_continuous(breaks = buffer_radii, labels = paste0(buffer_radii, "m")) +
    labs(title = "D. Precision vs Buffer Radius", x = "Spatial Buffer Radius", y = "Precision") +
    theme_bw(base_family = "Helvetica")

combined_decay <- (p1 + p2) / (p3 + p4)
ggsave("plots/spatial_signal_degradation_quantified.png", combined_decay, width = 11, height = 8.5, dpi = 300)
cat("Saved diagnostic decay plot to plots/spatial_signal_degradation_quantified.png!\n")
