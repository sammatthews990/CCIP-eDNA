import re
import codecs

with codecs.open('eDNA_CPUE_Comparison.qmd', 'r', 'utf-8') as f:
    rmd = f.read()

chunks = re.findall(r'```(?:\{r.*?\}|r)\n(.*?)\n```', rmd, re.DOTALL)
code = '\n\n'.join(chunks)

code += '\n\n# Save plots\n'
code += 'dir.create("plots", showWarnings=FALSE)\n'
code += 'ggsave("plots/p_hm.png", p_hm, width=10, height=8, dpi=300)\n'
code += 'ggsave("plots/p_cm.png", p_cm, width=10, height=8, dpi=300)\n'
code += 'ggsave("plots/p_2dhm.png", p_2dhm, width=10, height=8, dpi=300)\n'
code += 'ggsave("plots/p_2dcm.png", p_2dcm, width=10, height=8, dpi=300)\n'

with codecs.open('eDNA_CPUE_Comparison.R', 'w', 'utf-8') as f:
    f.write(code)
