"""Run the full pipeline in order (from the repository root):  python scripts/run_all.py

Steps 01-02 take ~45 min on a 12-core / 32 GB workstation; 05 (sub-clustering + GSEA) ~25 min.
Re-plot Figure 2 from cached results with:  python scripts/05_fig2_deep_dives.py --replot
"""
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).parent
STEPS = ["01_qc_doublets.py", "02_integrate_annotate.py", "03_fig1_global_landscape.py",
         "04_pseudobulk_de.py", "05_fig2_deep_dives.py", "06_fig3_convergence.py",
         "07_fig4_cell_communication.py", "08_fig5_target_prioritization.py", "09_qc_figures.py", "10_donor_covariates.py"]

if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    for s in STEPS:
        print(f"===== {s}", flush=True)
        sys.argv = [s]
        runpy.run_path(str(HERE / s), run_name="__main__")
