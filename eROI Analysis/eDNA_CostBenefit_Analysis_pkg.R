# eDNA Cost-Benefit Analysis for COTS Culling (eROI Analysis)
# This script uses the reefDNA package to process historical COTS Cull and Manta Tow data.

library(dplyr)
pkg_path <- if(dir.exists("reefDNA")) "reefDNA" else "../reefDNA"
devtools::load_all(pkg_path)
library(readr)

# Set base path (assuming this script is run from the project root or eROI folder)
root_dir <- if (dir.exists("data")) "." else ".."

cull_file <- file.path(root_dir, "data/260201_COTS-Cull-Data-Ewels.xlsx")
manta_file <- file.path(root_dir, "data/COTS Program  Manta Tow Data-2026-02-04.csv")
edna_file <- file.path(root_dir, "data/eDNA data_ALL_20260225.xlsx")

# 1. Load and prepare Cull Data using reefDNA
cull_data <- prepare_cull_data(cull_file)
cull_summary <- aggregate_cull_data(cull_data)
first_day_summary <- get_scout_metrics(cull_data)

# 2. Load and prepare Manta Tow Data using reefDNA
manta_data <- prepare_manta_data(manta_file)
manta_summary <- summarize_manta_data(manta_data)

# 3. Load and prepare eDNA Data using reefDNA
edna_summary <- summarize_edna_data(edna_file)

# 4. Join ALL data together using reefDNA
combined_data <- combine_reef_data(cull_summary, first_day_summary, manta_summary, edna_summary)

# 5. Extract Full Matched Dataset (eDNA + Cull both exist)
full_matched_data <- combined_data %>%
  filter(!is.na(eDNA_Pct_Positive))

# Keep the top 20 logic for the original deep dive section
top_20_annual <- combined_data %>%
  group_by(Year) %>%
  slice_max(order_by = Wasted_Diver_Hours, n = 20, with_ties = FALSE) %>%
  ungroup()

# Save datasets
write_csv(combined_data, file.path(root_dir, "eROI Analysis", "combined_reef_analysis_v2.csv"))
write_csv(top_20_annual, file.path(root_dir, "eROI Analysis", "top_20_annual_wasted_v2.csv"))
write_csv(full_matched_data, file.path(root_dir, "eROI Analysis", "full_matched_data_v2.csv"))
