#!/usr/bin/env python3
"""
Score PyBFS baseflow separation against the modified strict baseflow reference
(Xie et al., 2020, as modified in Section 2.4 of the paper) for all sites, and
report NRMSE (Eq. 50) for the 2018 validation year.

This reproduces the labeled-baseflow comparison behind Figures 12-14 and the
NRMSE values reported in Section 3.1, using pybfs.modified_strict_baseflow()
and pybfs.separation_skill() (pybfs/skill.py).
"""
import os
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import pybfs

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = REPO_ROOT / "outputs"


from main_baseflow_2018_all_sites import pick_site_col


def main():
    parser = argparse.ArgumentParser(
        description="Compute NRMSE between PyBFS baseflow and the modified strict "
                     "baseflow reference for all sites, for the 2018 validation year."
    )
    parser.add_argument("--params", default=str(DATA_DIR / "bfs_params_python_all_sites.csv"))
    parser.add_argument("--streamflow", default=str(DATA_DIR / "streamflow_with_date.csv"))
    parser.add_argument("--outdir", default=str(OUT_DIR))
    parser.add_argument("--outfile", default="nrmse_summary_all_sites.csv")
    parser.add_argument("--start", default="2018-01-01", help="Start of validation window (inclusive)")
    parser.add_argument("--end", default="2018-12-31", help="End of validation window (inclusive)")
    parser.add_argument("--date_col", default="Date")
    parser.add_argument("--quantile", type=float, default=0.8,
                         help="Quantile used by modified_strict_baseflow (default 0.8, per paper)")
    parser.add_argument("--skip_missing", action="store_true",
                         help="Skip sites missing in streamflow columns or failing BFS")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    bfs_params = pd.read_csv(args.params, encoding="utf-8-sig")
    site_col = pick_site_col(bfs_params)
    site_ids = bfs_params[site_col].astype(str).str.strip().unique().tolist()
    print(f"Loaded {len(site_ids)} sites from: {args.params} (site col='{site_col}')")

    qdv = pd.read_csv(args.streamflow, encoding="utf-8-sig")
    if args.date_col not in qdv.columns:
        raise KeyError(f"'{args.date_col}' not found in {args.streamflow}. Columns: {list(qdv.columns)}")
    qdv[args.date_col] = pd.to_datetime(qdv[args.date_col], errors="coerce")
    qdv.columns = qdv.columns.astype(str)

    start = pd.to_datetime(args.start)
    end = pd.to_datetime(args.end)

    # Restrict to the validation window *before* running PyBFS, matching
    # main_baseflow_2018_all_sites.py. Running on the full multi-year record
    # instead changes the model's internal storage state (no cold start at
    # the window's beginning) and gives different, less accurate results.
    qdv = qdv[(qdv[args.date_col] >= start) & (qdv[args.date_col] <= end)].copy()
    streamflow_cols = set(qdv.columns)

    rows = []
    for idx, site_id in enumerate(site_ids, start=1):
        site_id_str = str(site_id).strip()

        if site_id_str not in streamflow_cols:
            msg = f"Site {site_id_str} not found as a column in streamflow.csv"
            if args.skip_missing:
                print("SKIP -", msg)
                rows.append({"SiteID": site_id_str, "NRMSE": np.nan, "n_labeled_days": 0, "Reason": "Missing streamflow column"})
                continue
            raise KeyError(msg + " (use --skip_missing to continue)")

        print(f"[{idx}/{len(site_ids)}] Scoring site {site_id_str}")
        try:
            sf = pd.DataFrame({
                "Date": qdv[args.date_col].values,
                "Streamflow": pd.to_numeric(qdv[site_id_str], errors="coerce").values,
            })

            try:
                site_int = int(float(site_id_str))
                basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params, site_int)
            except Exception:
                basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params, site_id_str)

            lb, x1, wb, por = basin_char[1], basin_char[2], basin_char[3], basin_char[4]
            alpha, beta, ks, kb, kz = gw_hyd

            SBT = pybfs.base_table(lb, x1, wb, beta, kb, sf, por)
            skill_df, _ = pybfs.separation_skill(sf, SBT, basin_char, gw_hyd, flow, quantile=args.quantile)

            window = (skill_df["Date"] >= start) & (skill_df["Date"] <= end)
            mask = window & skill_df["BF_strict"] & skill_df["RES"].notna()
            n_days = int(mask.sum())

            if n_days == 0:
                rows.append({"SiteID": site_id_str, "NRMSE": np.nan, "n_labeled_days": 0, "Reason": "No labeled baseflow days in window"})
                print(f"WARN - {site_id_str}: no labeled baseflow days in {args.start}..{args.end}")
                continue

            mean_q = skill_df.loc[mask, "Q"].mean()
            rmse = np.sqrt((skill_df.loc[mask, "RES"] ** 2).mean())
            nrmse = rmse / mean_q if mean_q > 0 else np.nan

            rows.append({"SiteID": site_id_str, "NRMSE": nrmse, "n_labeled_days": n_days, "Reason": ""})
            print(f"OK   - {site_id_str}: NRMSE={nrmse:.3f} (n={n_days})")

        except Exception as e:
            print(f"FAIL - {site_id_str}: {e}")
            rows.append({"SiteID": site_id_str, "NRMSE": np.nan, "n_labeled_days": 0, "Reason": str(e)})
            if not args.skip_missing:
                raise

    out_path = os.path.join(args.outdir, args.outfile)
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print("\nSaved NRMSE summary:", out_path)


if __name__ == "__main__":
    main()
