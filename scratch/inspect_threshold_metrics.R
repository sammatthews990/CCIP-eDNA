suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rsample)
  library(glmmTMB)
  library(pROC)
})

# Extract dat_evt without running the full script each time
# We can read dat_evt logic directly:
dat_glmm <- readRDS("dat_glmm.rds") # if exists, or recreate dat_evt

if (!file.exists("dat_evt.rds")) {
  source("eDNA_CPUE_Comparison_final.R", local = FALSE)
  saveRDS(dat_evt, "dat_evt.rds")
} else {
  dat_evt <- readRDS("dat_evt.rds")
}

cpue_levels <- c(0.02, 0.04, 0.08)
perc_grid <- seq(0, 100, by = 2)

inspect_cpue <- function(cpue_thr) {
  cat("==========================================\n")
  cat(sprintf("CPUE THRESHOLD: %.2f\n", cpue_thr))
  cat("==========================================\n")
  
  grid_res <- map_dfr(perc_grid, function(p) {
    m <- metrics_for(p, cpue_thr, dat_evt)
    prec <- m["Precision"]
    rec <- m["Recall"]
    f1 <- m["F1"]
    f2 <- if (is.na(prec) || is.na(rec) || (4 * prec + rec == 0)) 0 else (5 * prec * rec / (4 * prec + rec))
    f1_5 <- if (is.na(prec) || is.na(rec) || (2.25 * prec + rec == 0)) 0 else (3.25 * prec * rec / (2.25 * prec + rec))
    
    tibble(
      perc_thresh = p,
      Accuracy = m["Accuracy"],
      Precision = prec,
      Recall = rec,
      F1 = f1,
      F1.5 = f1_5,
      F2 = f2
    )
  })
  
  print(grid_res %>% filter(perc_thresh %in% seq(10, 90, by = 5)), n = 30)
}

walk(cpue_levels, inspect_cpue)
