import codecs

def extract():
    with codecs.open('eDNA_CPUE_Comparison.qmd', 'r', 'utf-8') as f:
        lines = f.readlines()

    in_chunk = False
    r_code = []

    for line in lines:
        if line.startswith('```'):
            if '{r' in line or line.strip() == '```r':
                in_chunk = True
            else:
                in_chunk = False
        elif in_chunk:
            r_code.append(line)

    code = ''.join(r_code)

    code += '\n\n# Safely load ggplot and export defined variables\n'
    code += 'library(ggplot2)\n'
    code += 'dir.create("plots", showWarnings=FALSE)\n'
    code += 'if(exists("p_hm")) ggsave("plots/temporal_degradation_GLM.png", plot = p_hm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_cm")) ggsave("plots/multi_horizon_CPUE.png", plot = p_cm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_cm_hor")) ggsave("plots/cv_horizons_02_CPUE.png", plot = p_cm_hor, width=10, height=4, dpi=300)\n'
    code += 'if(exists("p_2dhm")) ggsave("plots/2D_AND_F1_Heatmap.png", plot = p_2dhm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_2dcm")) ggsave("plots/2D_AND_Multi_CM.png", plot = p_2dcm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_spatial_04")) ggsave("plots/spatial_buffer_CM_04.png", plot = p_spatial_04, width=12, height=8, dpi=300)\n'
    code += 'if(exists("p_spatial_02")) ggsave("plots/spatial_buffer_CM_02.png", plot = p_spatial_02, width=12, height=8, dpi=300)\n'
    code += 'if(exists("df_metrics")) write.csv(df_metrics, "plots/spatial_metrics.csv", row.names=FALSE)\n'

    with codecs.open("eDNA_CPUE_Comparison_final.R", 'w', 'utf-8') as f:
        f.write(code)

extract()
