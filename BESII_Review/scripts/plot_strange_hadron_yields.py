#!/usr/bin/env python3
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = Path("data/strange_hadron_yields_vs_energy.csv")
OUT_PNG = Path("data/strange_hadron_yields_vs_energy.png")

PARTICLE_ORDER = ["K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"]
STYLE = {
    "K+": {"marker": "o", "color": "tab:blue"},
    "K-": {"marker": "s", "color": "tab:cyan"},
    "Ks0": {"marker": "D", "color": "tab:green"},
    "Lambda": {"marker": "^", "color": "tab:red"},
    "Lambda_bar": {"marker": "v", "color": "tab:orange"},
    "Xi": {"marker": "P", "color": "tab:purple"},
    "Xi_bar": {"marker": "X", "color": "tab:pink"},
    "phi": {"marker": "*", "color": "tab:brown"},
}


def to_float(text):
    try:
        return float(text)
    except Exception:
        return None


def combined_error(row):
    stat = to_float(row.get("stat_err")) or 0.0
    sys = to_float(row.get("sys_err")) or 0.0
    return math.hypot(stat, sys) if (stat or sys) else 0.0


def main():
    rows = list(csv.DictReader(SRC.open(encoding="utf-8")))
    fig, ax = plt.subplots(figsize=(9.2, 6.0))

    for particle in PARTICLE_ORDER:
        pts = [r for r in rows if r["particle"] == particle]
        if not pts:
            continue
        pts.sort(key=lambda r: to_float(r["energy_GeV"]) or 0.0)
        xs = [to_float(r["energy_GeV"]) for r in pts]
        ys = [to_float(r["value"]) for r in pts]
        yerr = [combined_error(r) for r in pts]
        style = STYLE[particle]
        ax.errorbar(
            xs,
            ys,
            yerr=yerr,
            marker=style["marker"],
            ms=5.0,
            lw=1.4,
            capsize=2,
            color=style["color"],
            label=particle,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title("Extracted Strange-Hadron Yields vs Energy")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=9, ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=240)
    plt.close(fig)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
