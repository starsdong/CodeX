#!/usr/bin/env python3
"""Plot published Lambda global-polarization points versus collision energy."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "lambda_global_polarization.csv"
FIG_DIR = ROOT / "figures"
MPL_DIR = ROOT / ".matplotlib-cache"

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
MPL_DIR.mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt


STYLE = {
    "HADES": {"marker": "s", "color": "#c44536", "label": "HADES"},
    "STAR": {"marker": "o", "color": "#1976a2", "label": "STAR"},
    "ALICE": {"marker": "^", "color": "#4c8f45", "label": "ALICE"},
}


def as_float(value: str) -> float:
    if value == "":
        return 0.0
    return float(value)


def load_main_points() -> list[dict[str, object]]:
    points: list[dict[str, object]] = []
    with DATA.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["particle"] != "Lambda" or row["main_plot"].upper() != "TRUE":
                continue
            stat = as_float(row["stat_err_percent"])
            sys_minus = as_float(row["sys_err_minus_percent"])
            sys_plus = as_float(row["sys_err_plus_percent"])
            points.append(
                {
                    "experiment": row["experiment"],
                    "system": row["system"],
                    "sqrt_snn": as_float(row["sqrt_snn_GeV"]),
                    "p": as_float(row["p_h_percent"]),
                    "err_minus": math.hypot(stat, sys_minus),
                    "err_plus": math.hypot(stat, sys_plus),
                }
            )
    return sorted(points, key=lambda item: item["sqrt_snn"])


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    points = load_main_points()

    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Sans",
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)

    for experiment, style in STYLE.items():
        selected = [p for p in points if p["experiment"] == experiment]
        if not selected:
            continue
        x = [p["sqrt_snn"] for p in selected]
        y = [p["p"] for p in selected]
        yerr = [
            [p["err_minus"] for p in selected],
            [p["err_plus"] for p in selected],
        ]
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=style["marker"],
            ms=6.5,
            color=style["color"],
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markeredgewidth=0.7,
            ecolor=style["color"],
            elinewidth=1.0,
            capsize=2.5,
            label=style["label"],
            zorder=3,
        )

    # A light guide-to-eye through the selected Lambda points.
    ax.plot(
        [p["sqrt_snn"] for p in points],
        [p["p"] for p in points],
        color="0.55",
        lw=0.9,
        alpha=0.55,
        zorder=1,
    )

    ax.axhline(0, color="0.2", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(1.9, 7000)
    ax.set_ylim(-0.8, 9.7)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$P_\Lambda$ (%)")
    ax.set_title(r"Global $\Lambda$ Polarization in Heavy-Ion Collisions")
    ax.grid(True, which="major", axis="both", color="0.86", lw=0.7)
    ax.grid(True, which="minor", axis="x", color="0.92", lw=0.5)
    ax.legend(frameon=False, loc="upper right")
    png = FIG_DIR / "lambda_global_polarization_vs_energy.png"
    pdf = FIG_DIR / "lambda_global_polarization_vs_energy.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
