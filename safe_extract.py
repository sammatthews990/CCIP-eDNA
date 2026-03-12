import re
import codecs

def extract_r_code(qmd_path, out_path):
    with codecs.open(qmd_path, 'r', 'utf-8') as f:
        content = f.read()

    # Safely match ```{r ...} or ```r up until the closing ```
    matches = re.findall(r'```(?:\{r[^\}]*\}|r)\n([\s\S]*?)\n```', content)

    code = '\n\n'.join(matches)

    code += '\n\n# Safely load ggplot and export defined variables\n'
    code += 'library(ggplot2)\n'
    code += 'dir.create("plots", showWarnings=FALSE)\n'
    code += 'if(exists("p_hm")) ggsave("plots/temporal_degradation_GLM.png", plot = p_hm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_cm")) ggsave("plots/multi_horizon_CPUE.png", plot = p_cm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_2dhm")) ggsave("plots/2D_AND_F1_Heatmap.png", plot = p_2dhm, width=10, height=8, dpi=300)\n'
    code += 'if(exists("p_2dcm")) ggsave("plots/2D_AND_Multi_CM.png", plot = p_2dcm, width=10, height=8, dpi=300)\n'

    with codecs.open(out_path, 'w', 'utf-8') as f:
        f.write(code)

extract_r_code("eDNA_CPUE_Comparison.qmd", "eDNA_CPUE_Comparison_final.R")
