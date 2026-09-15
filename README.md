# Streamflow Separation Paper — Data and Code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22780129.svg)](https://doi.org/10.5281/zenodo.22780129)

Input data, model parameters, and analysis scripts to reproduce the
streamflow separation results reported in:

> Webster-Esho, E., Konrad, C. P., Talukdar, J., Aghababaei, A., Van der
> Heijden, R., Li, X., Williams, G. P., Jones, N. L., Rizzo, D. M., & Clement,
> T. P. *Development and Testing of a Process-Based Streamflow Separation
> Model for Baseflow Simulation and Short-Term Low-Flow Forecasts.*
> (manuscript under review)

The study applies [PyBFS](https://github.com/BYU-Hydroinformatics/pybfs) to 50
unregulated USGS gauges across Alabama, Florida, Georgia, and Mississippi for
the 2018 validation year, and compares the simulated baseflow to a "labeled"
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
> Short-Term Low-Flow Forecasts"* (Version 1.0.0) [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.22780129

## Quick start

Install PyBFS and run the analysis:

```bash
pip install pybfs
python scripts/run_analysis.py --skip_missing
```

Both steps write into `outputs/`, which is not tracked. The full run over all
50 sites takes a few minutes.

## Files

### Inputs — `data/`

| File | Description |
|---|---|
| `streamflow_with_date.csv` | Daily streamflow (m³/day) for 50 USGS gauges, 2013–2018, wide format (`Date` + one column per site number). |
| `siteinfo_paper.csv` | Site metadata: USGS site number and drainage area (`AREA`). |
| `bfs_params_python_all_sites.csv` | Calibrated PyBFS parameters per site (`Lb`, `X1`, `Wb`, `POR`, `ALPHA`, `BETA`, `Ks`, `Kb`, `Kz`, `Qthresh`, `Rs`, `Rb1`, `Rb2`, `Prec`, `Frac4Rise`), estimated with the parameter-estimation tool described in Section 2.2.4. |
| `baseflow_only_2018_all_sites.csv` | Reference output: PyBFS-simulated baseflow for 2018, one column per site — the data used in the paper. |

### Scripts — `scripts/`

| File | Description |
|---|---|
| `run_analysis.py` | Entry point. Runs both steps below in order. |
| `main_baseflow_2018_all_sites.py` | Step 1: runs PyBFS (`pybfs.base_table` + `pybfs.bfs`) for all 50 sites over the 2018 validation window and writes the simulated baseflow, reproducing the model output shown as the green line/points in Figure 12. |
| `compute_baseflow_skill_all_sites.py` | Step 2: derives the labeled reference baseflow per site using `pybfs.modified_strict_baseflow()` (Section 2.4), scores the PyBFS baseflow against it with `pybfs.separation_skill()`, and reports NRMSE (Eq. 50) for 2018. This reproduces the comparison behind Figures 12–14. |

Both scripts take `--skip_missing` to continue past a site that fails (for
example a missing streamflow column) instead of stopping, along with
`--start`, `--end`, and `--quantile` flags. Run either with `--help` for
details.

### Outputs — `outputs/`

| File | Description |
|---|---|
| `baseflow_only_2018_all_sites_copy.csv` | From step 1. Compare it against `data/baseflow_only_2018_all_sites.csv` to confirm the pipeline reproduces the paper's output; the two should match exactly. |
| `nrmse_summary_all_sites.csv` | From step 2. One row per site with NRMSE and the number of labeled baseflow days used, reproducing the data behind Figures 12–14. |

## Notes

- The `modified_strict_baseflow` algorithm and the RMSE/NRMSE scoring logic
  live in the PyBFS package (`pybfs/skill.py`), so they stay in sync with the
  rest of the model code. An earlier single-site draft of the algorithm,
  adapted from a Colab notebook, is superseded by
  `pybfs.modified_strict_baseflow()`.
- `compute_baseflow_skill_all_sites.py` restricts the streamflow record to the
  2018 validation window before running PyBFS, matching
  `main_baseflow_2018_all_sites.py`, so both scripts evaluate the same period.
- The manuscript and response-to-reviewer documents are kept in the working
  folder but excluded from version control by `.gitignore`.
