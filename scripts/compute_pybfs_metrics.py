#!/usr/bin/env python3
"""Compute the baseflow-evaluation metrics reported in the PyBFS paper.

Per-gauge NRMSE, KGE with its components, percent bias and the CV of labeled
baseflow; regional and per-state summaries; Pearson, Spearman and partial
correlations of the metrics with catchment size and flow magnitude; recession
behaviour on strict baseflow days.

Inputs are read from the repository's data/ folder:
    streamflow_with_date.csv          observed daily streamflow, 2013-2018 (m3/day)
    baseflow_only_2018_all_sites.csv  PyBFS simulated baseflow, 2018 (m3/day)
    siteinfo_paper.csv                drainage area used in the paper (m2)
    site_metadata_nwis.csv            station name, state code, NWIS drainage area (sq mi)

No PyBFS runs are needed. Labeled baseflow points are derived with the modified
strict baseflow algorithm applied to the 2018 record.

Usage:
    python scripts/compute_pybfs_metrics.py               # data/ -> results/
    python scripts/compute_pybfs_metrics.py --fetch-nwis  # re-download site_metadata_nwis.csv first

Requires numpy, pandas and scipy.
"""
import argparse
import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = REPO_ROOT / "results"

INPUT_FILES = ["streamflow_with_date.csv", "baseflow_only_2018_all_sites.csv", "siteinfo_paper.csv"]
NWIS_SITE_URL = "https://waterservices.usgs.gov/nwis/site/?format=rdb&siteOutput=expanded&sites="
NWIS_FILE = "site_metadata_nwis.csv"

STUDY_YEAR = 2018
LABEL_QUANTILE = 0.8
M3D_PER_CMS = 86400.0
M2_PER_SQMI = 2_589_988.110336
MIN_LABELED = 30       # flag gauges with fewer labeled points
LOW_CV = 0.05          # flag gauges whose labeled baseflow is nearly flat (KGE is unreliable there)

STATE_FIPS = {"01": "AL", "12": "FL", "13": "GA", "28": "MS"}
# NWIS lists 02228500 (North Prong St. Marys River at Moniac, GA) under Florida; the paper counts it
# as Georgia.
STATE_OVERRIDES = {"2228500": "GA"}
TEST_SITES = {"2422500": "AL", "2312200": "FL", "2388975": "GA", "2472000": "MS"}

SUMMARY_METRICS = ["nrmse", "kge", "kge_r", "kge_alpha", "kge_beta", "pbias_pct", "bias_share_of_mse",
                   "cv_labeled_bf"]
CORR_METRICS = ["nrmse", "kge", "kge_r", "kge_alpha", "kge_beta"]
SIZE_COVARIATES = ["drainage_area_sqmi", "drainage_area_nwis_sqmi", "mean_streamflow_2018_cms"]


# ------------------------------------------------------------------------------------ inputs

def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def parse_rdb(text):
    """Parse a USGS RDB table: '#' comment lines, a header line, then a column-format line."""
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[2:])), sep="\t", dtype=str)


def fetch_nwis_metadata():
    """Download station name, state code and drainage area from NWIS for every site in siteinfo_paper.csv."""
    print("Downloading site metadata from USGS NWIS")
    sites = pd.read_csv(DATA_DIR / "siteinfo_paper.csv", dtype=str)["site_no"]
    rdb = fetch(NWIS_SITE_URL + ",".join(s.zfill(8) for s in sites)).decode("utf-8")
    parse_rdb(rdb)[["site_no", "station_nm", "state_cd", "drain_area_va"]].to_csv(DATA_DIR / NWIS_FILE,
                                                                                  index=False)


def load_inputs():
    for name in INPUT_FILES + [NWIS_FILE]:
        if not (DATA_DIR / name).exists():
            hint = " (run with --fetch-nwis to download it)" if name == NWIS_FILE else ""
            raise SystemExit(f"Missing input file: {DATA_DIR / name}{hint}")
    flow = pd.read_csv(DATA_DIR / "streamflow_with_date.csv")
    flow.index = pd.to_datetime(flow.pop("Date"), format="%m/%d/%Y")
    baseflow = pd.read_csv(DATA_DIR / "baseflow_only_2018_all_sites.csv")
    baseflow.index = pd.to_datetime(baseflow.pop("Date"), format="%Y-%m-%d")
    siteinfo = pd.read_csv(DATA_DIR / "siteinfo_paper.csv", dtype={"site_no": str}).set_index("site_no")
    nwis = pd.read_csv(DATA_DIR / NWIS_FILE, dtype=str)
    nwis.index = nwis["site_no"].str.lstrip("0")

    sites = list(siteinfo.index)
    missing = [s for s in sites if s not in flow or s not in baseflow or s not in nwis.index]
    if missing:
        raise ValueError(f"Sites missing from an input file: {missing}")
    year = flow.loc[flow.index.year == STUDY_YEAR, sites]
    if len(year) != len(baseflow) or not (year.index == baseflow.index).all():
        raise ValueError(f"Streamflow and baseflow dates differ for {STUDY_YEAR}")
    if year.isna().any().any() or baseflow[sites].isna().any().any():
        raise ValueError(f"Missing values in the {STUDY_YEAR} streamflow or baseflow")
    return flow[sites], baseflow[sites], siteinfo, nwis


# ---------------------------------------------------------------------------------- labeling

def strict_baseflow(Q, quantile=LABEL_QUANTILE):
    """Strict baseflow days: the modified strict baseflow rules of the paper's Methods.

    Rule 3 removes 4 days before and 5 days after each high-flow day (flow at or
    above the given quantile) on which -dQ/dt is increasing.
    """
    dQ = (Q[2:] - Q[:-2]) / 2

    # 1. Days with dQ/dt >= 0.
    wet1 = np.concatenate([[True], dQ >= 0, [True]])

    # 2. Two days before and three days after each run of such days.
    idx_first = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == 1)[0] + 1
    idx_last = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == -1)[0]
    idx_before = np.repeat([idx_first], 2) - np.tile(range(1, 3), idx_first.shape)
    idx_next = np.repeat([idx_last], 3) + np.tile(range(1, 4), idx_last.shape)
    idx_remove = np.concatenate([idx_before, idx_next])
    wet2 = np.full(Q.shape, False)
    wet2[idx_remove.clip(min=0, max=Q.shape[0] - 1)] = True

    # 3. Window around high-flow recession days (flow at or above the quantile).
    wet3_core = np.concatenate([[True], dQ[1:] - dQ[:-1] < 0, [True, True]])
    idx3_all = np.where(wet3_core)[0]
    idx3 = idx3_all[Q[idx3_all] >= np.quantile(Q, quantile)]
    idx3_buffer = np.repeat([idx3], 10) + np.tile(range(-4, 6), idx3.shape)
    wet3 = np.full(Q.shape, False)
    wet3[idx3_buffer.clip(min=0, max=Q.shape[0] - 1)] = True

    # 4. Days followed by a day with a larger -dQ/dt.
    wet4 = np.concatenate([[True], dQ[1:] - dQ[:-1] < 0, [True, True]])

    return ~(wet1 | wet2 | wet3 | wet4)


def label_baseflow(Q):
    """Labeled days: strict days plus days at or below the mean strict-day flow (Eq. 23)."""
    strict = strict_baseflow(Q)
    return strict | (Q <= Q[strict].mean()), strict


# ----------------------------------------------------------------------------------- metrics

def error_metrics(obs, sim):
    """NRMSE (Eq. 24), original KGE (Eq. 25) with its components, and percent bias (Eq. 26).

    Percent bias is positive when PyBFS overestimates the labeled baseflow.
    bias_share_of_mse is mean error squared over mean squared error: the part of
    NRMSE that is bias rather than scatter.
    """
    error = sim - obs
    r = np.corrcoef(sim, obs)[0, 1]
    alpha = sim.std(ddof=1) / obs.std(ddof=1)
    beta = sim.mean() / obs.mean()
    rmse = np.sqrt(np.mean(error ** 2))
    return {
        "nrmse": rmse / obs.mean(),
        "kge": 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2),
        "kge_r": r,
        "kge_alpha": alpha,
        "kge_beta": beta,
        "pbias_pct": 100 * error.sum() / obs.sum(),
        "bias_share_of_mse": error.mean() ** 2 / np.mean(error ** 2),
        "rmse_cms": rmse / M3D_PER_CMS,
    }


def nrmse_pbias(obs, sim):
    if len(obs) < 2:
        return np.nan, np.nan
    return (np.sqrt(np.mean((sim - obs) ** 2)) / obs.mean(),
            100 * (sim.sum() - obs.sum()) / obs.sum())


def recession_exponent(flow, days):
    """Exponent b of -dQ/dt = a Q^b, by OLS in log space over the given days.

    dQ/dt is the central difference used by the labeling rules. Days on which the
    series is not receding are dropped; the dropped fraction is returned.
    """
    idx = np.where(days)[0]
    idx = idx[(idx > 0) & (idx < len(flow) - 1)]
    rate = -(flow[idx + 1] - flow[idx - 1]) / 2
    receding = rate > 0
    if receding.sum() < 3:
        return np.nan, int(receding.sum()), 1 - receding.mean()
    b, _ = np.polyfit(np.log10(flow[idx][receding]), np.log10(rate[receding]), 1)
    return b, int(receding.sum()), 1 - receding.mean()


def evaluate_gauges(flow, baseflow, siteinfo, nwis):
    year = flow[flow.index.year == STUDY_YEAR]
    gauges, recession, pairs = [], [], []
    for site in siteinfo.index:
        Q = year[site].to_numpy(float)
        B = baseflow[site].to_numpy(float)
        labeled, strict = label_baseflow(Q)
        threshold_only = labeled & ~strict
        obs, sim = Q[labeled], B[labeled]
        record = flow[site].dropna()
        state_cd = nwis.at[site, "state_cd"]
        state = STATE_OVERRIDES.get(site, STATE_FIPS.get(state_cd, state_cd))
        area = siteinfo.at[site, "AREA"] / M2_PER_SQMI
        area_nwis = float(nwis.at[site, "drain_area_va"])

        gauges.append({
            "gauge_id": site,
            "state": state,
            "station_nm": nwis.at[site, "station_nm"],
            "drainage_area_sqmi": area,
            "mean_streamflow_2018_cms": Q.mean() / M3D_PER_CMS,
            "n_labeled": int(labeled.sum()),
            "mean_labeled_bf_cms": obs.mean() / M3D_PER_CMS,
            "sd_labeled_bf_cms": obs.std(ddof=1) / M3D_PER_CMS,
            **error_metrics(obs, sim),
            "cv_labeled_bf": obs.std(ddof=1) / obs.mean(),
            "cv_streamflow_2018": Q.std(ddof=1) / Q.mean(),
            "cv_streamflow_2013_2018": record.std(ddof=1) / record.mean(),
            "mean_streamflow_2013_2018_cms": record.mean() / M3D_PER_CMS,
            "n_days_2013_2018": len(record),
            "n_strict": int(strict.sum()),
            "drainage_area_nwis_sqmi": area_nwis,
        })

        nrmse_strict, pbias_strict = nrmse_pbias(Q[strict], B[strict])
        nrmse_flat, pbias_flat = nrmse_pbias(Q[threshold_only], B[threshold_only])
        b_obs, n_obs, _ = recession_exponent(Q, strict)
        b_sim, n_sim, not_receding = recession_exponent(B, strict)
        recession.append({
            "gauge_id": site,
            "state": state,
            "n_strict": int(strict.sum()),
            "n_threshold_only": int(threshold_only.sum()),
            "nrmse_strict": nrmse_strict,
            "nrmse_threshold_only": nrmse_flat,
            "pbias_strict_pct": pbias_strict,
            "pbias_threshold_only_pct": pbias_flat,
            "recession_b_observed": b_obs,
            "recession_b_pybfs": b_sim,
            "n_fit_observed": n_obs,
            "n_fit_pybfs": n_sim,
            "frac_pybfs_not_receding": not_receding,
        })

        pairs.append(pd.DataFrame({
            "gauge_id": site,
            "date": year.index[labeled].strftime("%Y-%m-%d"),
            "labeled_bf_cms": obs / M3D_PER_CMS,
            "pybfs_bf_cms": sim / M3D_PER_CMS,
            "strict": strict[labeled],
        }))

    return pd.DataFrame(gauges), pd.DataFrame(recession), pd.concat(pairs, ignore_index=True)


# ------------------------------------------------------------------------------- summaries

def quartiles(values):
    v = pd.Series(values).dropna()
    return {"n": len(v), "median": v.median(), "q25": v.quantile(0.25), "q75": v.quantile(0.75),
            "min": v.min(), "max": v.max()}


def region_summary(gauges):
    return pd.DataFrame([{"metric": m, **quartiles(gauges[m])} for m in SUMMARY_METRICS])


def state_summary(gauges):
    rows = []
    for state, g in [*gauges.groupby("state"), ("All", gauges)]:
        rows.append({
            "state": state,
            "n_gauges": len(g),
            "nrmse_median": g["nrmse"].median(),
            "nrmse_q25": g["nrmse"].quantile(0.25),
            "nrmse_q75": g["nrmse"].quantile(0.75),
            "kge_median": g["kge"].median(),
            "kge_q25": g["kge"].quantile(0.25),
            "kge_q75": g["kge"].quantile(0.75),
            "pbias_median_pct": g["pbias_pct"].median(),
            "pct_nrmse_lt_0.5": 100 * (g["nrmse"] < 0.5).mean(),
            "pct_kge_gt_0.5": 100 * (g["kge"] > 0.5).mean(),
        })
    return pd.DataFrame(rows)


def fisher_ci(coef, n, variance_factor=1.0, n_controls=0):
    """95% CI by Fisher z; variance_factor 1.06 for Spearman (Fieller et al., 1957)."""
    half = stats.norm.ppf(0.975) * np.sqrt(variance_factor / (n - 3 - n_controls))
    z = np.arctanh(coef)
    return np.tanh(z - half), np.tanh(z + half)


def partial_spearman(x, y, z):
    """Spearman partial correlation of x and y controlling for z, t-test on n - 3 df."""
    rx, ry, rz = (pd.Series(v).rank().to_numpy() for v in (x, y, z))
    rxy, rxz, ryz = (np.corrcoef(a, b)[0, 1] for a, b in ((rx, ry), (rx, rz), (ry, rz)))
    coef = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    df = len(rx) - 3
    t = coef * np.sqrt(df / (1 - coef ** 2))
    return coef, 2 * stats.t.sf(abs(t), df)


def detectable_r(n, alpha=0.05, power=0.8):
    """Smallest |r| detected with the given power by a two-sided test (Fisher z approximation)."""
    return np.tanh((stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)) / np.sqrt(n - 3))


def correlation_table(gauges):
    n = len(gauges)
    rows = []
    for covariate in SIZE_COVARIATES + ["cv_labeled_bf"]:
        x = gauges[covariate].to_numpy(float)
        for metric in CORR_METRICS:
            y = gauges[metric].to_numpy(float)
            for method, result, factor in [
                ("pearson", stats.pearsonr(x, y), 1.0),
                ("pearson_log10x", stats.pearsonr(np.log10(x), y), 1.0),
                ("spearman", stats.spearmanr(x, y), 1.06),
            ]:
                coef = float(result.statistic)
                low, high = fisher_ci(coef, n, factor)
                rows.append({"covariate": covariate, "metric": metric, "method": method, "control": "",
                             "n": n, "coef": coef, "p": float(result.pvalue),
                             "ci95_low": low, "ci95_high": high})
    for covariate in SIZE_COVARIATES:
        for metric in ["nrmse", "kge"]:
            coef, p = partial_spearman(gauges[covariate], gauges[metric], gauges["cv_labeled_bf"])
            low, high = fisher_ci(coef, n, 1.06, n_controls=1)
            rows.append({"covariate": covariate, "metric": metric, "method": "spearman_partial",
                         "control": "cv_labeled_bf", "n": n, "coef": coef, "p": p,
                         "ci95_low": low, "ci95_high": high})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------- report

def fmt(value, digits=3):
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    if isinstance(value, (float, np.floating)):
        return "" if np.isnan(value) else f"{value:.{digits}f}"
    return str(value)


def md_table(df, digits=3):
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "|".join("---" for _ in df.columns) + "|"]
    lines += ["| " + " | ".join(fmt(v, digits) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def write_report(gauges, recession, regional, states, corr, path):
    n = len(gauges)

    def corr_rows(covariates, metrics, methods, control=""):
        sel = corr[corr.covariate.isin(covariates) & corr.metric.isin(metrics)
                   & corr.method.isin(methods) & (corr.control == control)]
        return sel[["covariate", "metric", "method", "coef", "p", "ci95_low", "ci95_high"]]

    out = [
        "# PyBFS evaluation metrics",
        "",
        f"Generated by compute_pybfs_metrics.py from the files in data/. {n} gauges. "
        f"Labeled baseflow derived on the {STUDY_YEAR} record; all metrics at labeled points "
        "unless stated; no new PyBFS runs. Percent bias is positive when PyBFS overestimates.",
        "",
        "## 1. Regional distribution",
        "",
        md_table(regional),
        "",
        "## 2. By state",
        "",
        "02228500 (North Prong St. Marys River at Moniac, GA) is counted as Georgia.",
        "",
        md_table(states),
        "",
        "## 3. Catchment scale and flow magnitude (Section 3.1.3)",
        "",
        "`drainage_area_sqmi` is the area used in the paper; `drainage_area_nwis_sqmi` is the NWIS value. "
        "`pearson_log10x` is Pearson on log10 of the covariate.",
        "",
        md_table(corr_rows(SIZE_COVARIATES, ["nrmse", "kge"], ["pearson", "pearson_log10x", "spearman"])),
        "",
        f"Smallest |r| detectable with 80% power at n = {n}, two-sided alpha = 0.05 "
        f"(Fisher z approximation): {detectable_r(n):.3f}.",
        "",
        "## 4. CV of labeled baseflow",
        "",
        "The CV of streamflow at the labeled timestamps equals the CV of labeled baseflow, because the "
        "labeled values are the observed flows; the contrast is with the CV of the whole record.",
        "",
        md_table(pd.DataFrame([
            {"metric": "CV labeled baseflow", **quartiles(gauges.cv_labeled_bf)},
            {"metric": f"CV streamflow {STUDY_YEAR}", **quartiles(gauges.cv_streamflow_2018)},
            {"metric": "CV streamflow 2013-2018", **quartiles(gauges.cv_streamflow_2013_2018)},
        ])),
        "",
        "Spearman correlation of each metric with CV of labeled baseflow:",
        "",
        md_table(corr_rows(["cv_labeled_bf"], CORR_METRICS, ["spearman"])),
        "",
        "Spearman partial correlations controlling for CV of labeled baseflow:",
        "",
        md_table(corr_rows(SIZE_COVARIATES, ["nrmse", "kge"], ["spearman_partial"], "cv_labeled_bf")),
        "",
        "## 5. KGE decomposition at the test sites",
        "",
        md_table(gauges[gauges.gauge_id.isin(TEST_SITES)][
            ["gauge_id", "state", "n_labeled", "nrmse", "kge", "kge_r", "kge_alpha", "kge_beta",
             "pbias_pct", "cv_labeled_bf"]]),
        "",
        "## 6. Recession behaviour",
        "",
        "Strict days are the recession points that pass rules 1-4; threshold-only days are the flat "
        "low-flow points added by Eq. 23. `recession_b` is the exponent b in -dQ/dt = a Q^b, fitted by "
        "OLS in log space over strict days (central-difference derivative), for observed streamflow and "
        "for PyBFS baseflow on the same days. b = 1 is a linear reservoir.",
        "",
        md_table(recession[recession.gauge_id.isin(TEST_SITES)].drop(columns=["n_fit_observed"])),
        "",
        "Regional medians:",
        "",
        md_table(pd.DataFrame([{"metric": c, **quartiles(recession[c])} for c in [
            "nrmse_strict", "nrmse_threshold_only", "pbias_strict_pct", "pbias_threshold_only_pct",
            "recession_b_observed", "recession_b_pybfs", "frac_pybfs_not_receding"]])),
        "",
        "## 7. Data flags",
        "",
        f"Gauges with fewer than {MIN_LABELED} labeled points: "
        f"{', '.join(gauges.gauge_id[gauges.n_labeled < MIN_LABELED]) or 'none'} "
        f"(minimum {gauges.n_labeled.min()}).",
        "",
        f"Gauges with CV of labeled baseflow below {LOW_CV} (KGE is unreliable on a near-flat series):",
        "",
        md_table(gauges[gauges.cv_labeled_bf < LOW_CV][
            ["gauge_id", "station_nm", "n_labeled", "cv_labeled_bf", "nrmse", "kge", "kge_r"]]),
        "",
        "Gauges with missing days in 2013-2018 (affects only the full-record columns): "
        + (", ".join(f"{r.gauge_id} ({r.n_days_2013_2018} of 2191 days)"
                     for r in gauges[gauges.n_days_2013_2018 < 2191].itertuples()) or "none") + ".",
        "",
    ]
    path.write_text("\n".join(out), encoding="utf-8")


# ------------------------------------------------------------------------------------ main

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fetch-nwis", action="store_true",
                        help=f"download {NWIS_FILE} from USGS NWIS before computing (needs network access)")
    parser.add_argument("--datadir",
                        help="read inputs from somewhere other than the repository's data/ folder")
    parser.add_argument("--outdir",
                        help="write tables somewhere other than the repository's results/ folder")
    args = parser.parse_args()

    global OUT_DIR, DATA_DIR
    if args.datadir:
        DATA_DIR = Path(args.datadir).resolve()
        if not DATA_DIR.exists():
            raise SystemExit(f"No such data directory: {DATA_DIR}")
    if args.outdir:
        OUT_DIR = Path(args.outdir).resolve()

    if args.fetch_nwis:
        fetch_nwis_metadata()
    flow, baseflow, siteinfo, nwis = load_inputs()
    gauges, recession, pairs = evaluate_gauges(flow, baseflow, siteinfo, nwis)
    regional = region_summary(gauges)
    states = state_summary(gauges)
    corr = correlation_table(gauges)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tables = {
        "gauge_metrics.csv": gauges,
        "labeled_pairs.csv": pairs,
        "summary_region.csv": regional,
        "summary_by_state.csv": states,
        "correlations.csv": corr,
        "recession_metrics.csv": recession,
    }
    for name, table in tables.items():
        table.to_csv(OUT_DIR / name, index=False, float_format="%.6g")
    write_report(gauges, recession, regional, states, corr, OUT_DIR / "metrics_report.md")
    print(f"Wrote {len(tables)} tables and metrics_report.md to {OUT_DIR}")


if __name__ == "__main__":
    main()
