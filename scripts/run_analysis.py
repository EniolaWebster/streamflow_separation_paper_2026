#!/usr/bin/env python3
"""
Single entry point to reproduce the paper's streamflow-partitioning analysis:

  1. Run PyBFS baseflow separation for all 50 sites over the 2018 evaluation
     year (main_baseflow_2018_all_sites.py) -> baseflow_only_2018_all_sites_copy.csv
  2. Score PyBFS against the modified strict baseflow reference for all sites
     (compute_baseflow_skill_all_sites.py) -> nrmse_summary_all_sites.csv

Both write into outputs/. Run from anywhere:
    python scripts/run_analysis.py

See README.md for details on each file and output.
"""
import runpy
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


def run(script_name, argv):
    sys.argv = [str(SCRIPT_DIR / script_name)] + argv
    runpy.run_path(str(SCRIPT_DIR / script_name), run_name="__main__")


def main():
    extra_args = sys.argv[1:]  # forwarded to both steps, e.g. --skip_missing

    print("=" * 70)
    print("STEP 1/2: PyBFS baseflow separation for all 50 sites (2018)")
    print("=" * 70)
    run("main_baseflow_2018_all_sites.py", extra_args)

    print("\n" + "=" * 70)
    print("STEP 2/2: Scoring PyBFS against modified strict baseflow reference")
    print("=" * 70)
    run("compute_baseflow_skill_all_sites.py", extra_args)

    print("\nDone. Outputs are in outputs/; see README.md for details.")


if __name__ == "__main__":
    main()
