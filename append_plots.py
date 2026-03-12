import codecs

append_code = """

# Save plots
dir.create("plots", showWarnings=FALSE)
ggsave("plots/p_hm.png", p_hm, width=10, height=8, dpi=300)
ggsave("plots/p_cm.png", p_cm, width=10, height=8, dpi=300)
ggsave("plots/p_2dhm.png", p_2dhm, width=10, height=8, dpi=300)
ggsave("plots/p_2dcm.png", p_2dcm, width=10, height=8, dpi=300)
"""

with codecs.open('eDNA_CPUE_Comparison.R', 'a', 'utf-8') as f:
    f.write(append_code)
