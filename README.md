# Streamflow Separation Paper — Data and Code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22780128.svg)](https://doi.org/10.5281/zenodo.22780128)

Input data, model parameters, analysis scripts, result tables and figures to
reproduce the streamflow separation results reported in:

> Webster-Esho, E., Konrad, C. P., Talukdar, J., Aghababaei, A., Van der
> Heijden, R., Li, X., Williams, G. P., Jones, N. L., Rizzo, D. M., & Clement,
> T. P. *Development and Testing of a Process-Based Streamflow Separation
> Model for Baseflow Simulation and Short-Term Low-Flow Forecasts.*
> (manuscript under review)

The study applies [PyBFS](https://github.com/BYU-Hydroinformatics/pybfs) to 50
unregulated USGS gauges across Alabama, Florida, Georgia, and Mississippi for
the 2018 evaluation year, and compares the simulated baseflow to a "labeled"
reference baseflow derived from the modified strict baseflow algorithm
(Section 2.4 of the paper).

The PyBFS model code lives in its own repository; this repository holds only
the data and scripts specific to the paper.

## Citation

This repository is archived on Zenodo. To cite the data and code:

> Webster-Esho, E., Konrad, C. P., Talukdar, J., Aghababaei, A., Van der
> Heijden, R., Li, X., Williams, G. P., Jones, N. L., Rizzo, D. M., & Clement,
> T. P. (2026). *Data and analysis code for "Development and Testing of a
> Process-Based Streamflow Separation Model for Baseflow Simulation and
> Short-Term Low-Flow Forecasts"* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.22780128

The DOI above always resolves to the most recent version. To cite a specific
one, use its own DOI: 10.5281/zenodo.22797807 for v1.1.0 (current) or
10.5281/zenodo.22780129 for v1.0.0.

## Quick start

Install PyBFS and run the analysis:

```bash
pip install pybfs
python scripts/run_analysis.py --skip_missing
```

Both steps write into `outputs/`, which is not tracked. The full run over all
50 sites takes a few minutes.

To regenerate the evaluation metrics and the data figures:

```bash
python scripts/compute_pybfs_metrics.py   # data/    -> results/
python scripts/make_figures.py            # results/ -> figures/
```

`results/` and `figures/` are tracked, so these two commands overwrite the
committed copies; the metrics script takes a few seconds and needs no PyBFS
run. Environment: Python 3.x with numpy, pandas and scipy; matplotlib for the
figures; the pybfs package for `fig14` and for `run_analysis.py`.

## Files

### Inputs — `data/`

| File | Description |
|---|---|
| `streamflow_with_date.csv` | Daily streamflow (m³/day) for 50 USGS gauges, 2013–2018, wide format (`Date` as `M/D/YYYY` + one column per site number). Site numbers are written without the leading zero of the USGS form (`2422500` is USGS 02422500). Missing days are blank. |
| `siteinfo_paper.csv` | Site metadata: USGS site number (`site_no`) and drainage area (`AREA`, m²). |
| `bfs_params_python_all_sites.csv` | Calibrated PyBFS parameters per site (`Lb`, `X1`, `Wb`, `POR`, `ALPHA`, `BETA`, `Ks`, `Kb`, `Kz`, `Qthresh`, `Rs`, `Rb1`, `Rb2`, `Prec`, `Frac4Rise`), estimated with the parameter-estimation tool described in Section 2.2.4. Also carries `site_no`, `AREA` (m²), and the `Error` and `BFF` values reported by the calibration. |
| `baseflow_only_2018_all_sites.csv` | Reference output: PyBFS-simulated baseflow for 2018 (m³/day, `Date` as `YYYY-MM-DD`), one column per site — the data used in the paper. |
| `site_metadata_nwis.csv` | USGS NWIS site metadata for the 50 gauges: `site_no` (8-digit USGS form), `station_nm`, `state_cd` (FIPS state code) and `drain_area_va` (NWIS drainage area, square miles). Used by `compute_pybfs_metrics.py` for station names, state grouping and the NWIS drainage-area column. |

### Scripts — `scripts/`

| File | Description |
|---|---|
| `run_analysis.py` | Entry point. Runs both steps below in order. |
| `main_baseflow_2018_all_sites.py` | Step 1: runs PyBFS (`pybfs.base_table` + `pybfs.bfs`) for all 50 sites over the 2018 evaluation window and writes the simulated baseflow, reproducing the model output shown as the green line/points in Figure 11. |
| `compute_baseflow_skill_all_sites.py` | Step 2: derives the labeled reference baseflow per site using `pybfs.modified_strict_baseflow()` (Section 2.4), scores the PyBFS baseflow against it with `pybfs.separation_skill()`, and reports NRMSE (Eq. 24) for 2018. This reproduces the comparison behind Figures 11–13. |
| `compute_pybfs_metrics.py` | Computes the evaluation metrics reported in the paper from `data/` alone (no PyBFS run): per-gauge NRMSE, KGE and its components, percent bias, CV of labeled baseflow, regional and per-state summaries, correlations with catchment size and flow magnitude, and recession behaviour. Writes the tables in `results/`. Requires numpy, pandas and scipy. `--fetch-nwis` re-downloads `data/site_metadata_nwis.csv` from NWIS. |
| `make_figures.py` | Draws the data figures in `figures/` from `data/` and `results/`: the four test-site hydrographs with labeled and simulated baseflow (`fig11`), model performance against drainage area and mean streamflow (`fig13`), and the short-term low-flow forecasts (`fig14`). `--figure` selects one figure or `all`. Requires matplotlib; `fig14` also requires pybfs. |

The first three scripts take `--skip_missing` to continue past a site that
fails (for example a missing streamflow column) instead of stopping, along
with `--start`, `--end`, and `--quantile` flags. Run any script with `--help`
for details.

### Outputs — `outputs/`

| File | Description |
|---|---|
| `baseflow_only_2018_all_sites_copy.csv` | From step 1. Compare it against `data/baseflow_only_2018_all_sites.csv` to confirm the pipeline reproduces the paper's output; the two should match exactly. |
| `nrmse_summary_all_sites.csv` | From step 2. One row per site with NRMSE and the number of labeled baseflow days used, reproducing the data behind Figures 11–13. |

### Results — `results/`

Written by `compute_pybfs_metrics.py`. All metrics are computed at the labeled
baseflow days of 2018 unless stated; flows are in m³/s (`cms`).

| File | Description |
|---|---|
| `gauge_metrics.csv` | One row per gauge; columns defined below. |
| `labeled_pairs.csv` | One row per gauge and labeled day: `gauge_id`, `date`, `labeled_bf_cms` (observed streamflow on that day, the reference baseflow), `pybfs_bf_cms` (PyBFS baseflow on that day), `strict` (True if the day passes the strict-baseflow rules, False if it was added by the threshold rule). |
| `summary_region.csv` | Median, quartiles, minimum and maximum of each metric over all 50 gauges. |
| `summary_by_state.csv` | Per-state and overall medians and quartiles of NRMSE, KGE and percent bias, and the share of gauges with NRMSE < 0.5 and KGE > 0.5. |
| `correlations.csv` | Pearson, Pearson on log10 of the covariate, Spearman, and Spearman partial correlations (controlling for CV of labeled baseflow) between each metric and drainage area or mean streamflow, with p-values and 95% confidence intervals. |
| `recession_metrics.csv` | Per gauge: NRMSE and percent bias split into strict days and threshold-only days, and the recession exponent b of -dQ/dt = a·Q^b fitted over strict days for observed streamflow (`recession_b_observed`) and PyBFS baseflow (`recession_b_pybfs`). |

Columns of `gauge_metrics.csv`:

| Column | Definition |
|---|---|
| `gauge_id` | USGS site number without the leading zero. |
| `state` | State used for the per-state summaries (02228500 is counted as GA). |
| `station_nm` | NWIS station name. |
| `drainage_area_sqmi` | Drainage area used in the paper (`AREA` from `siteinfo_paper.csv`, converted to square miles). |
| `drainage_area_nwis_sqmi` | Drainage area reported by NWIS (`drain_area_va`), square miles. |
| `mean_streamflow_2018_cms` | Mean observed streamflow in 2018. |
| `mean_streamflow_2013_2018_cms`, `cv_streamflow_2013_2018`, `n_days_2013_2018` | Mean, coefficient of variation and number of days with data over the full 2013–2018 record. |
| `cv_streamflow_2018` | Coefficient of variation of the 2018 streamflow. |
| `n_labeled`, `n_strict` | Number of labeled baseflow days, and how many of them are strict-baseflow days. |
| `mean_labeled_bf_cms`, `sd_labeled_bf_cms`, `cv_labeled_bf` | Mean, standard deviation and coefficient of variation of the labeled baseflow. |
| `nrmse` | Root-mean-square error of PyBFS baseflow against labeled baseflow, divided by the mean labeled baseflow. |
| `rmse_cms` | The same root-mean-square error, in m³/s. |
| `kge`, `kge_r`, `kge_alpha`, `kge_beta` | Kling–Gupta efficiency and its components: correlation, ratio of standard deviations (simulated/labeled), ratio of means (simulated/labeled). |
| `pbias_pct` | Percent bias; positive when PyBFS overestimates the labeled baseflow. |
| `bias_share_of_mse` | Squared mean error divided by mean squared error: the fraction of the error that is bias rather than scatter. |

### Figures — `figures/`

Written by `make_figures.py`; PNG at 400 dpi and vector PDF of each figure.

| File | Description |
|---|---|
| `fig11.png`, `fig11.pdf` | Four test sites (02422500 AL, 02312200 FL, 02388975 GA, 02472000 MS): observed streamflow, PyBFS baseflow, labeled baseflow points and the simulated value on those days; linear and logarithmic scales. |
| `fig13.png`, `fig13.pdf` | NRMSE and KGE against drainage area and mean 2018 streamflow for all 50 gauges. |
| `fig14.png`, `fig14.pdf` | Short-term low-flow forecasts at the four test sites: the 2018 hydrograph with the forecast window shaded, and the forecast baseflow against observed streamflow over that window. |
| `fig14_windows.csv` | The forecast window chosen for each test site (`start`, `end`, number of rise days, mean flow). The selection rule is stated in the docstring of `make_figures.py`. |
| `fig14_forecast_series.csv` | Daily observed streamflow and forecast baseflow (m³/s) over each forecast window. |

## Notes

- The `modified_strict_baseflow` algorithm and the RMSE/NRMSE scoring logic
  live in the PyBFS package (`pybfs/skill.py`), so they stay in sync with the
  rest of the model code. An earlier single-site draft of the algorithm,
  adapted from a Colab notebook, is superseded by
  `pybfs.modified_strict_baseflow()`.
- `compute_baseflow_skill_all_sites.py` restricts the streamflow record to the
  2018 evaluation window before running PyBFS, matching
  `main_baseflow_2018_all_sites.py`, so both scripts evaluate the same period.
- `compute_pybfs_metrics.py` carries its own copy of the labeling rules so
  that the metrics can be computed without PyBFS installed; it applies them to
  the 2018 record.
- PyBFS version: to reproduce the simulated baseflow and the forecasts with
  the same model code, use PyBFS at commit `2efaca8` of
  https://github.com/BYU-Hydroinformatics/pybfs. The parameter-estimation
  routine described in Section 2.2.4 (`pybfs.calibrate_beta1.bfs_calibrate_beta1`)
  is on the `calibration-beta1` branch of that repository (commit `3cc4044`).
- The manuscript and response-to-reviewer documents are kept in the working
  folder but kept out of version control by `.gitignore`.
