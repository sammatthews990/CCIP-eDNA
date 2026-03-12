# Save plots
library(ggplot2)
dir.create("plots", showWarnings = FALSE)
ggsave("plots/temporal_degradation_GLM.png", plot = p_hm, width = 10, height = 8, dpi = 300)
ggsave("plots/multi_horizon_CPUE.png", plot = p_cm, width = 10, height = 8, dpi = 300)
ggsave("plots/2D_AND_F1_Heatmap.png", plot = p_2dhm, width = 10, height = 8, dpi = 300)
ggsave("plots/2D_AND_Multi_CM.png", plot = p_2dcm, width = 10, height = 8, dpi = 300)
