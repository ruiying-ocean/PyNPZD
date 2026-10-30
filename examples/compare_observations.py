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
    from npzd import NPZDModel

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
    """Save an eight-panel comparison as PNG and vector PDF."""
    style = {
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(4, 2, figsize=(10, 11), sharex=True, sharey="col")
        fig.subplots_adjust(
            left=0.09, right=0.97, bottom=0.11, top=0.86, hspace=0.39, wspace=0.24
        )
        fig.suptitle("NPZD and observed seasonal cycles", fontsize=19, y=0.975)
        fig.text(
            0.5,
            0.941,
            "Year 5 · Published EMPOWER station parameters · No fitting to these observations",
            ha="center",
            fontsize=10,
        )
        for row, (station, (_, label, _)) in enumerate(STATIONS.items()):
            simulated = model[model.station == station]
            paired = monthly[monthly.station == station]
            for col, (variable, field, ylabel) in enumerate(
                [
                    ("chlorophyll", "chlorophyll", "Chlorophyll (mg m$^{-3}$)"),
                    ("nitrate", "N", "Nitrate (mmol N m$^{-3}$)"),
                ]
            ):
                ax = axes[row, col]
                ax.plot(
                    simulated.day,
                    simulated[field],
                    color="#0072B2",
                    lw=1.8,
                    label="Model: seasonal curve",
                    zorder=2,
                )
                ax.scatter(
                    MONTH_CENTRES,
                    paired[f"model_{variable}"],
                    marker="s",
                    s=24,
                    facecolors="white",
                    edgecolors="#0072B2",
                    linewidths=1.1,
                    label="Model: monthly mean",
                    zorder=4,
                )
                ax.scatter(
                    MONTH_CENTRES,
                    paired[f"observed_{variable}"],
                    s=25,
                    color="#222222",
                    marker="o",
                    label="Observed monthly climatology",
                    zorder=5,
                )
                metric = metrics[
                    (metrics.station == station) & (metrics.variable == variable)
                ].iloc[0]
                ax.text(
                    0.97,
                    0.93,
                    f"RMSE {metric.rmse:.2f} · {metric.months}/12 months",
                    ha="right",
                    va="top",
                    transform=ax.transAxes,
                    fontsize=9,
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.9,
                        "pad": 2,
                    },
                )
                ax.set_title(f"{chr(97 + 2 * row + col)}) {label}", loc="left", pad=8)
                ax.set_ylabel(ylabel)
                ax.set_xlim(0, 365)
                ax.set_xticks(MONTH_CENTRES, list("JFMAMJJASOND"))
                ax.grid(axis="y", color="#dddddd", lw=0.6)
                ax.set_axisbelow(True)
        for col in range(2):
            axes[0, col].set_ylim(0, axes[0, col].get_ylim()[1] * 1.12)
            axes[-1, col].set_xlabel("Month")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.921),
            ncol=3,
            frameon=False,
            fontsize=9,
        )
        fig.text(
            0.09,
            0.052,
            "Observations: MODIS-Aqua chlorophyll and WOA2018 mixed-layer nitrate, distributed with Phydra (Post et al., 2024).\n"
            "Missing satellite months remain blank. RMSE uses paired monthly means; no observational uncertainty was supplied.\n"
            "Updated climatologies and original EMPOWER forcing: this is not an exact reproduction of the 2015 figures.",
            fontsize=8,
            va="center",
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
