"""Plot year-five NPZD output against published monthly observational climatologies."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

STATIONS = {
    "biotrans": (2, "BIOTRANS · 47°N, 20°W", {}),
    "india": (1, "India · 60°N, 20°W", {"mz": 0.0}),
    "papa": (4, "Papa · 50°N, 145°W", {"Vp0": 1.25, "alpha": 0.075, "Imax": 1.25}),
    "kerfix": (
        3,
        "KERFIX · 50.67°S, 68.42°E",
        {"Vp0": 1.25, "alpha": 0.075, "Imax": 2.0},
    ),
}
MONTH_EDGES = np.r_[0, np.cumsum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])]
MONTH_CENTRES = (MONTH_EDGES[:-1] + MONTH_EDGES[1:]) / 2
OBSERVATIONS = ROOT / "examples/data/phydra_stations_verification.csv"


def to_chlorophyll(phytoplankton: np.ndarray, c_to_chl: float) -> np.ndarray:
    """Convert mmol N m^-3 to mg chlorophyll m^-3 using EMPOWER's fixed ratios."""
    return phytoplankton * 6.625 * 12.0 / c_to_chl


def load_observations(path: Path = OBSERVATIONS) -> pd.DataFrame:
    """Read the authors' monthly table, preserving missing satellite observations."""
    raw = pd.read_csv(path, header=[0, 1], index_col=0)
    if not np.array_equal(raw.index, np.arange(12)):
        raise ValueError("Expected January–December observations indexed 0–11")
    return pd.concat(
        [
            pd.DataFrame(
                {
                    "station": station,
                    "month": np.arange(1, 13),
                    "observed_chlorophyll": raw[station]["Chla"].to_numpy(),
                    "observed_nitrate": raw[station]["N"].to_numpy(),
                }
            )
            for station in STATIONS
        ],
        ignore_index=True,
    )


def simulate_stations() -> pd.DataFrame:
    """Run five years at each station with the parameters published in Table 3."""
    from pynpzd import NPZDModel

    results = []
    for station, (number, _, overrides) in STATIONS.items():
        model = NPZDModel(ROOT / "examples/empower_observations.toml")
        model.settings["flag_stn"] = number
        model.params.update(overrides)
        model.initialize_model()
        time, state = model.run()
        assert np.all(np.isfinite(state)) and np.all(state >= 0), station
        final_year = time >= (model.settings["nyears"] - 1) * 365
        day = time[final_year] - (model.settings["nyears"] - 1) * 365
        frame = pd.DataFrame(state[final_year], columns=["P", "N", "Z", "D"])
        frame["station"] = station
        frame["day"] = day
        frame["month"] = np.searchsorted(MONTH_EDGES, day, side="right")
        frame["chlorophyll"] = to_chlorophyll(
            frame.P.to_numpy(), model.params["CtoChl"]
        )
        results.append(frame)
    return pd.concat(results, ignore_index=True)


def monthly_comparison(model: pd.DataFrame, observed: pd.DataFrame) -> pd.DataFrame:
    """Pair monthly model means with observations; missing values remain missing."""
    means = model.groupby(["station", "month"])[["chlorophyll", "N"]].mean()
    means = means.rename(
        columns={"chlorophyll": "model_chlorophyll", "N": "model_nitrate"}
    )
    return observed.merge(
        means, on=["station", "month"], how="left", validate="one_to_one"
    )


def comparison_metrics(monthly: pd.DataFrame) -> pd.DataFrame:
    """Calculate RMSE and model-minus-observation bias across available months."""
    rows = []
    for station, group in monthly.groupby("station", sort=False):
        for variable, units in [
            ("chlorophyll", "mg chlorophyll m^-3"),
            ("nitrate", "mmol N m^-3"),
        ]:
            pairs = group[[f"model_{variable}", f"observed_{variable}"]].dropna()
            error = pairs.iloc[:, 0] - pairs.iloc[:, 1]
            rows.append(
                {
                    "station": station,
                    "variable": variable,
                    "months": len(pairs),
                    "rmse": np.sqrt(np.mean(error**2)),
                    "bias": np.mean(error),
                    "units": units,
                }
            )
    return pd.DataFrame(rows)


def plot_comparison(
    model: pd.DataFrame, monthly: pd.DataFrame, metrics: pd.DataFrame, output: Path
) -> None:
    """Save a four-station, two-variable comparison as PNG and vector PDF."""
    ink, muted, model_colour, observed_colour = "#222222", "#6b6b6b", "#0072B2", "#D55E00"
    style = {
        "font.size": 9,
        "axes.edgecolor": muted,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": muted,
        "ytick.color": muted,
        "xtick.labelcolor": ink,
        "ytick.labelcolor": ink,
        "pdf.fonttype": 42,
    }
    columns = [
        ("chlorophyll", "chlorophyll", "Chlorophyll (mg m$^{-3}$)"),
        ("nitrate", "N", "Nitrate (mmol N m$^{-3}$)"),
    ]
    with plt.rc_context(style):
        fig, axes = plt.subplots(4, 2, figsize=(9, 9), sharex=True, sharey="col")
        fig.subplots_adjust(
            left=0.07, right=0.98, top=0.80, bottom=0.10, hspace=0.45, wspace=0.14
        )
        fig.text(0.07, 0.965, "Simulated and observed seasonal cycles",
                 fontsize=15, weight="bold", color=ink)
        fig.text(0.07, 0.937, "Year 5 of 5 · published EMPOWER station parameters · "
                 "not fitted to these observations", fontsize=9.5, color=muted)
        for row, (station, (_, label, _)) in enumerate(STATIONS.items()):
            simulated = model[model.station == station]
            paired = monthly[monthly.station == station]
            for col, (variable, field, _) in enumerate(columns):
                ax = axes[row, col]
                ax.plot(simulated.day, simulated[field], color=model_colour, lw=1.8,
                        label="Model, daily", zorder=2)
                ax.scatter(MONTH_CENTRES, paired[f"model_{variable}"], marker="s",
                           s=16, facecolors="white", edgecolors=model_colour,
                           linewidths=1, label="Model, monthly mean", zorder=3)
                ax.scatter(MONTH_CENTRES, paired[f"observed_{variable}"], s=24,
                           color=observed_colour, edgecolors="white", linewidths=0.6,
                           label="Observed, monthly climatology", zorder=4)
                metric = metrics[
                    (metrics.station == station) & (metrics.variable == variable)
                ].iloc[0]
                ax.set_title(label, loc="left", fontsize=9.5, weight="bold", color=ink)
                ax.set_title(f"RMSE {metric.rmse:.2f} · {metric.months}/12 months",
                             loc="right", fontsize=8, color=muted)
                ax.grid(axis="y", color="#e8e8e8", lw=0.7)
                ax.set_axisbelow(True)
                ax.tick_params(length=3)
        for col, (_, _, header) in enumerate(columns):
            ticks = axes[0, col].get_yticks()
            data_max = max(ax.dataLim.ymax for ax in axes[:, col])
            # End each shared axis on a labelled tick so no column looks cut off.
            axes[0, col].set_ylim(0, np.ceil(data_max / (ticks[1] - ticks[0])) * (ticks[1] - ticks[0]))
            fig.text(axes[0, col].get_position().x0, 0.845, header,
                     fontsize=11, weight="bold", color=ink)
        axes[-1, 0].set_xlim(0, 365)
        axes[-1, 0].set_xticks(MONTH_CENTRES, list("JFMAMJJASOND"))
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.06, 0.92),
                   ncol=3, frameon=False, fontsize=9, handletextpad=0.4,
                   columnspacing=1.6)
        fig.text(
            0.07,
            0.03,
            "Observations: MODIS-Aqua chlorophyll and WOA2018 mixed-layer nitrate via Phydra "
            "(Post et al., 2024); blank months have no satellite data.\nRMSE uses paired "
            "monthly means; no observational uncertainty was supplied. Updated climatologies, "
            "so not an exact reproduction of the 2015 figures.",
            fontsize=7.5,
            color=muted,
            linespacing=1.5,
        )
        for extension in ("png", "pdf"):
            fig.savefig(
                output / f"empower_observations.{extension}", dpi=300, facecolor="white"
            )
        plt.close(fig)


def main() -> None:
    """Run the example from the repository root and write its figure and source tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "data/examples/empower_observations"
    )
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    observations = load_observations()
    model = simulate_stations()
    monthly = monthly_comparison(model, observations)
    metrics = comparison_metrics(monthly)
    model.to_csv(output / "model_last_year.csv", index=False)
    monthly.to_csv(output / "monthly_comparison.csv", index=False)
    metrics.to_csv(output / "metrics.csv", index=False)
    plot_comparison(model, monthly, metrics, output)
    print(metrics.to_string(index=False))
    print(f"Figure and data saved to {output}")


if __name__ == "__main__":
    main()
