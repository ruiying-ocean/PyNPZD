"""Plot the year-five seasonal cycle at BIOTRANS with its mixed-layer forcing."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MONTH_EDGES = np.r_[0, np.cumsum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])]
MONTH_CENTRES = (MONTH_EDGES[:-1] + MONTH_EDGES[1:]) / 2
# Okabe-Ito colours, each paired with a distinct line style.
SERIES = [
    ("P", "Phytoplankton (P)", "#009E73", "-"),
    ("Z", "Zooplankton (Z)", "#D55E00", "--"),
    ("D", "Detritus (D)", "#CC79A7", ":"),
]


def simulate():
    """Run the README example: five years at BIOTRANS with EMPOWER-1.0 inputs."""
    from pynpzd import NPZDModel

    model = NPZDModel(ROOT / "examples/empower_observations.toml")
    time, state = model.run()
    assert np.all(np.isfinite(state)) and np.all(state >= 0)
    start = (model.settings["nyears"] - 1) * 365
    final_year = time >= start
    return model.MLD[:365], time[final_year] - start, state[final_year]


def plot(mld, day, state, output: Path) -> None:
    style = {
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(
            3, 1, figsize=(7.5, 7.5), sharex=True, layout="constrained"
        )
        fig.suptitle("Seasonal cycle at BIOTRANS (47°N, 20°W)", fontsize=14)

        axes[0].plot(np.arange(365), mld, color="#555555", lw=1.8)
        axes[0].set_ylim(np.max(mld) * 1.05, 0)  # depth increases downward
        axes[0].set_ylabel("Mixed-layer\ndepth (m)")
        axes[0].set_title("a) Forcing: monthly climatology, interpolated daily", loc="left")

        axes[1].plot(day, state[:, 1], color="#0072B2", lw=1.8)
        axes[1].set_ylim(0, np.max(state[:, 1]) * 1.1)
        axes[1].set_ylabel("Nitrate\n(mmol N m$^{-3}$)")
        axes[1].set_title("b) Nutrient (N)", loc="left")

        for column, label, colour, line in SERIES:
            index = "PNZD".index(column)
            axes[2].plot(day, state[:, index], color=colour, ls=line, lw=1.8, label=label)
        axes[2].set_ylim(0, np.max(state[:, [0, 2, 3]]) * 1.15)
        axes[2].set_ylabel("Concentration\n(mmol N m$^{-3}$)")
        axes[2].set_title("c) Plankton and detritus", loc="left")
        axes[2].legend(frameon=False, ncol=3, loc="upper right", fontsize=9)

        for ax in axes:
            ax.grid(axis="y", color="#dddddd", lw=0.6)
            ax.set_axisbelow(True)
        axes[-1].set_xlim(0, 365)
        axes[-1].set_xticks(MONTH_CENTRES, list("JFMAMJJASOND"))
        axes[-1].set_xlabel("Month (year 5 of 5; EMPOWER-1.0 baseline inputs)")
        fig.savefig(output / "station_cycle.png", dpi=200, facecolor="white")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/figures")
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    plot(*simulate(), output)
    print(f"Figure saved to {output / 'station_cycle.png'}")


if __name__ == "__main__":
    main()
