#!/usr/bin/env python3
"""Draw the paper's data figures from the files in data/ and results/.

  fig11  Four test sites: observed streamflow (black), PyBFS baseflow (green), labeled baseflow points
         (blue), and the simulated value on those same days (green points). Left column linear, right
         column logarithmic.
  fig13  Model performance (NRMSE, KGE) against drainage area and mean 2018 streamflow, all 50 gauges.
  fig14  Short-term low-flow forecasts at the four test sites. Left column marks the forecast window on
         the 2018 hydrograph; right column compares the forecast baseflow against observed streamflow
         over that window.

FORECAST WINDOW SELECTION (fig14). Each window is chosen by rule: among all 60-day windows starting
between 1 June and 1 November 2018, take the one with the fewest "rise days", a rise day being one where
streamflow increases by more than 5% over the previous day; ties break toward the lower mean flow. That
selects the longest effectively rain-free recession in the dry half of the year, which is the condition
the forecast routine assumes. The chosen windows are written to figures/fig14_windows.csv and the
forecast series to figures/fig14_forecast_series.csv.

Run scripts/compute_pybfs_metrics.py first: fig11 and fig13 read results/gauge_metrics.csv and
results/labeled_pairs.csv. fig11 and fig13 need numpy, pandas and matplotlib; fig14 also needs the
pybfs package.

Usage:
    python scripts/make_figures.py --figure fig11
    python scripts/make_figures.py --figure all
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
CMS = 86400.0
YEAR = 2018

TEST_SITES = [("2422500", "Alabama"), ("2312200", "Florida"),
              ("2388975", "Georgia"), ("2472000", "Mississippi")]

# Forecast window rule, see module docstring.
WINDOW_DAYS = 60
WINDOW_EARLIEST = f"{YEAR}-06-01"
WINDOW_LATEST = f"{YEAR}-11-01"
RISE_THRESHOLD = 0.05

STREAM_COLOR = "black"
BASEFLOW_COLOR = "#1a7f37"
LABELED_COLOR = "#1f6feb"
WINDOW_COLOR = "#d1242f"


def use_paper_style():
    """Matplotlib's standard look, as used for the manuscript's figures: DejaVu Sans, a framed legend,
    a light dashed grid, and all four spines."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.linewidth": 0.9,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.grid": True,
        "grid.linestyle": "--",
        "grid.linewidth": 0.5,
        "grid.alpha": 0.45,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.3",
        "lines.linewidth": 1.0,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
    })


def panel_label(ax, letter):
    ax.text(0.012, 0.94, f"({letter})", transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="top", ha="left")


def month_axis(ax):
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())


def load(data_dir, results_dir):
    flow = pd.read_csv(data_dir / "streamflow_with_date.csv")
    flow.index = pd.to_datetime(flow.pop("Date"), format="%m/%d/%Y")
    baseflow = pd.read_csv(data_dir / "baseflow_only_2018_all_sites.csv")
    baseflow.index = pd.to_datetime(baseflow.pop("Date"), format="%Y-%m-%d")
    pairs = pd.read_csv(results_dir / "labeled_pairs.csv", dtype={"gauge_id": str}, parse_dates=["date"])
    gauges = pd.read_csv(results_dir / "gauge_metrics.csv", dtype={"gauge_id": str})
    return flow[flow.index.year == YEAR] / CMS, baseflow / CMS, pairs, gauges


def save(fig, path):
    pdf = path.parent / f"{path.stem}.pdf"
    fig.savefig(path)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"wrote {path} and {pdf.name}")


# ------------------------------------------------------------------------------------- figure 11

def figure_11(flow, baseflow, pairs, gauges, path):
    fig, axes = plt.subplots(4, 2, figsize=(7.5, 8.6), sharex="col")
    letters = iter("abcdefgh")
    ordered = [(row, col) for col in (0, 1) for row in range(4)]
    labels = {pos: next(letters) for pos in ordered}

    for row, (site, state) in enumerate(TEST_SITES):
        labeled = pairs[pairs.gauge_id == site]
        for col, log in enumerate([False, True]):
            ax = axes[row, col]
            ax.plot(flow.index, flow[site], color=STREAM_COLOR, linewidth=0.7, label="Observed streamflow")
            ax.plot(baseflow.index, baseflow[site], color=BASEFLOW_COLOR, linewidth=1.2,
                    label="PyBFS baseflow")
            ax.plot(labeled.date, labeled.labeled_bf_cms, "o", color=LABELED_COLOR, markersize=2.6,
                    linestyle="none", label="Labeled baseflow")
            ax.plot(labeled.date, labeled.pybfs_bf_cms, "o", color=BASEFLOW_COLOR, markersize=2.6,
                    linestyle="none", label="Simulated, labeled days")
            if log:
                ax.set_yscale("log")
            panel_label(ax, labels[(row, col)])
            ax.set_ylabel("Flow (cms)")
            month_axis(ax)
            ax.margins(x=0.01)
            if row == 0:
                ax.set_title(f"{'Linear' if not log else 'Logarithmic'} scale", fontsize=9, pad=6)
            ax.annotate(f"{site} - {state}", xy=(0.115, 0.94), xycoords="axes fraction",
                        ha="left", va="top", fontsize=8)
            row_metrics = gauges.loc[gauges.gauge_id == site].iloc[0]
            ax.annotate(f"NRMSE: {row_metrics.nrmse:.2f}\nKGE: {row_metrics.kge:.2f}",
                        xy=(0.985, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=7.5,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.3",
                                  linewidth=0.7))

    handles, labs = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labs, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.012))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    save(fig, path)


# ------------------------------------------------------------------------------------- figure 14

def pick_window(series):
    """Fewest rise days in a 60-day window; ties break toward lower mean flow. See module docstring."""
    rises = (series.pct_change() > RISE_THRESHOLD).astype(int)
    best = None
    starts = series.loc[WINDOW_EARLIEST:WINDOW_LATEST].index
    for start in starts:
        end = start + pd.Timedelta(days=WINDOW_DAYS - 1)
        if end > series.index[-1]:
            break
        window = series.loc[start:end]
        score = (int(rises.loc[start:end].sum()), float(window.mean()))
        if best is None or score < best[0]:
            best = (score, start, end)
    return best[1], best[2], best[0][0]


def figure_14(flow, data_dir, path):
    import pybfs  # only needed here

    params = pd.read_csv(data_dir / "bfs_params_python_all_sites.csv")
    fig, axes = plt.subplots(4, 2, figsize=(7.5, 8.6))
    letters = iter("abcdefgh")
    ordered = [(row, col) for col in (0, 1) for row in range(4)]
    labels = {pos: next(letters) for pos in ordered}
    chosen = []
    series_out = []

    for row, (site, state) in enumerate(TEST_SITES):
        series = flow[site].dropna()
        start, end, rise_days = pick_window(series)
        chosen.append({"site": site, "state": state, "start": start.date(), "end": end.date(),
                       "rise_days": rise_days, "mean_flow_cms": round(float(series.loc[start:end].mean()), 3)})

        basin_char, gw_hyd, flow_params = pybfs.get_values_for_site(params, int(site))
        train = pd.DataFrame({"Date": series.loc[:start - pd.Timedelta(days=1)].index,
                              "Streamflow": series.loc[:start - pd.Timedelta(days=1)].values * CMS})
        table = pybfs.base_table(basin_char[1], basin_char[2], basin_char[3],
                                 gw_hyd[1], gw_hyd[3], train, basin_char[4])
        trained = pybfs.bfs(train, table, basin_char, gw_hyd, flow_params)
        state_cols = ["X", "Zb.L", "Zs.L", "StBase", "StSur", "SurfaceFlow", "Baseflow", "Rech"]
        initial = tuple(trained.iloc[-1][state_cols])

        dates = pd.date_range(start, end, freq="D")
        forecast = pybfs.forecast(pd.DataFrame({"date": dates, "streamflow": np.nan}),
                                  table, basin_char, gw_hyd, flow_params, initial)

        series_out.append(pd.DataFrame({
            "site": site, "state": state, "date": forecast["Date"].values,
            "observed_streamflow_cms": series.loc[start:end].values,
            "forecast_baseflow_cms": forecast["Baseflow"].values / CMS}))

        left = axes[row, 0]
        left.plot(series.index, series.values, color=STREAM_COLOR, linewidth=0.7,
                  label="Observed streamflow")
        left.axvspan(start, end, color=WINDOW_COLOR, alpha=0.16, linewidth=0,
                     label="Forecast window")
        left.set_yscale("log")
        left.set_ylabel("Flow (cms)")
        month_axis(left)
        left.margins(x=0.01)
        panel_label(left, labels[(row, 0)])
        left.annotate(f"{site} - {state}", xy=(0.115, 0.94), xycoords="axes fraction",
                      ha="left", va="top", fontsize=8)

        right = axes[row, 1]
        right.plot(series.loc[start:end].index, series.loc[start:end].values, color=LABELED_COLOR,
                   linewidth=1.0, label="Observed streamflow")
        right.plot(forecast["Date"], forecast["Baseflow"] / CMS, color=BASEFLOW_COLOR,
                   linewidth=1.3, linestyle="--", label="Forecast baseflow")
        right.set_ylabel("Flow (cms)")
        right.xaxis.set_major_locator(mdates.DayLocator(interval=15))
        right.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        right.margins(x=0.02)
        panel_label(right, labels[(row, 1)])
        right.annotate(f"{start:%d %b} to {end:%d %b %Y}", xy=(0.115, 0.94), xycoords="axes fraction",
                       ha="left", va="top", fontsize=8)
        if row == 0:
            left.set_title("2018 hydrograph, forecast window shaded", fontsize=9, pad=6)
            right.set_title("Forecast period", fontsize=9, pad=6)

    handles = (axes[0, 0].get_legend_handles_labels()[0] + axes[0, 1].get_legend_handles_labels()[0][1:])
    labs = (axes[0, 0].get_legend_handles_labels()[1] + axes[0, 1].get_legend_handles_labels()[1][1:])
    fig.legend(handles, labs, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.012))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    save(fig, path)

    table_out = pd.DataFrame(chosen)
    table_out.to_csv(path.parent / "fig14_windows.csv", index=False)
    print(table_out.to_string(index=False))
    pd.concat(series_out, ignore_index=True).to_csv(path.parent / "fig14_forecast_series.csv", index=False)


# ------------------------------------------------------------------------------------- figure 13

def figure_13(gauges, path):
    """Viridis points coloured by the plotted metric, a red least-squares line, and the Pearson
    coefficient in a boxed annotation, on linear axes."""
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.6))
    specs = [("drainage_area_sqmi", "nrmse", "Drainage Area (square miles)", "NRMSE"),
             ("mean_streamflow_2018_cms", "nrmse", "Mean Streamflow (cms)", "NRMSE"),
             ("drainage_area_sqmi", "kge", "Drainage Area (square miles)", "KGE"),
             ("mean_streamflow_2018_cms", "kge", "Mean Streamflow (cms)", "KGE")]

    for ax, (x, y, xlabel, ylabel), letter in zip(axes.ravel(), specs, "abcd"):
        xv = gauges[x].to_numpy(float)
        yv = gauges[y].to_numpy(float)
        points = ax.scatter(xv, yv, c=yv, cmap="viridis", s=48, edgecolor="white", linewidth=0.4,
                            zorder=3)
        fig.colorbar(points, ax=ax, fraction=0.046, pad=0.02)

        slope, intercept = np.polyfit(xv, yv, 1)
        line = np.array([xv.min(), xv.max()])
        ax.plot(line, slope * line + intercept, color="red", linewidth=1.8, zorder=2)

        r = float(np.corrcoef(xv, yv)[0, 1])
        ax.annotate(f"Pearson r: {r:.3f}", xy=(0.97, 0.06), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.3",
                              linewidth=0.8))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        panel_label(ax, letter)

    fig.tight_layout()
    save(fig, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figure", choices=["fig11", "fig13", "fig14", "all"],
                        help="which figure to draw (all if omitted)")
    parser.add_argument("--data", help="input folder, relative to the repository root (data)")
    parser.add_argument("--results",
                        help="folder holding the metrics tables, relative to the repository root (results)")
    parser.add_argument("--outdir", help="output folder, relative to the repository root (figures)")
    args = parser.parse_args()

    use_paper_style()
    data_dir = REPO_ROOT / (args.data or "data")
    results_dir = REPO_ROOT / (args.results or "results")
    figures = REPO_ROOT / (args.outdir or "figures")
    figures.mkdir(parents=True, exist_ok=True)

    flow, baseflow, pairs, gauges = load(data_dir, results_dir)
    figure = args.figure or "all"
    wanted = ["fig11", "fig13", "fig14"] if figure == "all" else [figure]
    for name in wanted:
        path = figures / f"{name}.png"
        if name == "fig11":
            figure_11(flow, baseflow, pairs, gauges, path)
        elif name == "fig13":
            figure_13(gauges, path)
        else:
            figure_14(flow, data_dir, path)


if __name__ == "__main__":
    main()
