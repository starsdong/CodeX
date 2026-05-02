#!/usr/bin/env python3
"""Plot total bbbar production cross sections versus collision energy."""

from __future__ import annotations

import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MPL_CACHE = Path("/private/tmp/codex-matplotlib-cache")
MPL_CACHE.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DATA_PATH = ROOT / "data" / "bbbar_total_cross_sections.csv"
FIG_DIR = ROOT / "fig"


def read_points(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def main() -> None:
    points = read_points(DATA_PATH)
    FIG_DIR.mkdir(exist_ok=True)

    styles = {
        "pp": {
            "marker": "o",
            "color": "#1f77b4",
            "label": "pp",
            "zorder": 3,
        },
        "pA": {
            "marker": "s",
            "color": "#d95f02",
            "label": "pA, divided by A",
            "zorder": 4,
        },
    }

    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)

    for system in ("pp", "pA"):
        subset = [row for row in points if row["system"] == system]
        x = [as_float(row, "sqrts_GeV") for row in subset]
        y = [as_float(row, "sigma_ub_per_nn") for row in subset]
        yerr_low = [as_float(row, "err_low_ub") for row in subset]
        yerr_high = [as_float(row, "err_high_ub") for row in subset]
        style = styles[system]
        ax.errorbar(
            x,
            y,
            yerr=[yerr_low, yerr_high],
            fmt=style["marker"],
            ms=6,
            lw=1.2,
            capsize=2.8,
            color=style["color"],
            markeredgecolor="black",
            markeredgewidth=0.35,
            label=style["label"],
            zorder=style["zorder"],
        )

    # A subtle guide-to-the-eye helps read the steep energy dependence without
    # implying a fit model.
    guide = sorted(
        (
            (as_float(row, "sqrts_GeV"), as_float(row, "sigma_ub_per_nn"))
            for row in points
            if row["system"] == "pp"
            and not row["label"].startswith("PHENIX pp 200 GeV e-h")
        ),
        key=lambda item: item[0],
    )
    ax.plot(
        [item[0] for item in guide],
        [item[1] for item in guide],
        color="#1f77b4",
        alpha=0.23,
        lw=1.6,
        zorder=1,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(30, 20000)
    ax.set_ylim(0.002, 1000)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ or $\sqrt{s}$ (GeV)")
    ax.set_ylabel(r"$\sigma_{b\bar{b}}$ ($\mu$b per NN)")
    ax.set_title(r"Total $b\bar{b}$ production cross section")
    ax.grid(True, which="major", color="#dddddd", lw=0.8)
    ax.grid(True, which="minor", color="#eeeeee", lw=0.45, alpha=0.7)
    ax.legend(frameon=False, loc="upper left")

    tick_values = [40, 200, 510, 2760, 5020, 7000, 13000]
    tick_labels = ["40", "200", "510", "2.76k", "5.02k", "7k", "13k"]
    ax.set_xticks(tick_values, tick_labels)

    ax.text(
        0.98,
        0.03,
        "pA points are normalized per target nucleon; errors are quadrature sums.",
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="#444444",
        transform=ax.transAxes,
    )

    out_png = FIG_DIR / "bbbar_total_cross_section_vs_energy.png"
    out_pdf = FIG_DIR / "bbbar_total_cross_section_vs_energy.pdf"
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()
