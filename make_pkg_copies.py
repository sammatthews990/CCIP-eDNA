import os
import re

files_to_refactor = [
    ("eDNA_CPUE_Comparison.qmd", "eDNA_CPUE_Comparison_pkg.qmd"),
    ("eDNA_CPUE_Comparison_final.R", "eDNA_CPUE_Comparison_final_pkg.R"),
    ("edna_reef_report_2025.R", "edna_reef_report_2025_pkg.R"),
    ("eROI Analysis/eDNA_CostBenefit_Report.qmd", "eROI Analysis/eDNA_CostBenefit_Report_pkg.qmd")
]

functions_to_remove = [
    "lm_label <- function", "make_confusion <- function", "plot_confusion <- function",
    "metrics_for <- function", "get_prior_cohort <- function", "fit_horizon_glmm <- function",
    "run_cv_for_cpue <- function", "run_cv_for_horizon <- function", "metrics_for_both <- function",
    "run_cv_2d <- function", "make_confusion_2d <- function", "evaluate_unified_spatial <- function",
    "best_pair <- function", "get_thr <- function", "build_popup <- function", "make_popup <- function"
]

for orig_rel, new_rel in files_to_refactor:
    orig = os.path.join(r"c:\Users\smatthew\Documents\GitKraken\CCIP-eDNA", orig_rel)
    new_f = os.path.join(r"c:\Users\smatthew\Documents\GitKraken\CCIP-eDNA", new_rel)
    
    if not os.path.exists(orig): 
        print(f"Skipping {orig}, does not exist")
        continue

    with open(orig, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    out = []
    skip = False
    brace_count = 0
    added_lib = False
    nested_braces = 0

    for line in lines:
        stripped = line.lstrip()
        
        # Check if we should start skipping
        if not skip:
            is_func = False
            for fn in functions_to_remove:
                if stripped.startswith(fn):
                    is_func = True
                    skip = True
                    # Initialize brace count for this line
                    # Functions defined like 'name <- function(...) {'
                    brace_count = line.count('{') - line.count('}')
                    break
            
            if is_func:
                continue
                
            # Add package load on the first standard package found
            if "library(dplyr)" in line and not added_lib:
                out.append(line)
                out.append("devtools::load_all('reefDNA')\n")
                added_lib = True
                continue
                
            # Provide an alternative place to add library if dplyr isnt loaded independently
            if "library(" in line and not added_lib:
                out.append(line)
                out.append("devtools::load_all('reefDNA')\n")
                added_lib = True
                continue
                
            out.append(line)
        else:
            # We are skipping inside a function block
            brace_count += line.count('{') - line.count('}')
            
            # Special case for run_cv_2d which contains best_pair
            # best_pair might finish and bring brace_count to 0 momentarily if they are sequentially declared? No, usually nested.
            # But the logic is: if brace_count <= 0, we've exited the outer function
            if brace_count <= 0:
                skip = False
                
    with open(new_f, 'w', encoding='utf-8') as f:
        f.writelines(out)
        
    print(f"Processed {orig} -> {new_f}")
