#!/usr/bin/env python3
"""Plot particle-to-pion yield ratios with retained GC and SCE THERMUS curves."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INPUT_CSV = Path("data/thermus_sce_particle_set_comparison.csv")
OUT_CSV = Path("data/particle_to_pi_ratios_with_thermus.csv")
OUT_PNG = Path("data/particle_to_pi_ratios_with_thermus.png")
OUT_PDF = Path("data/particle_to_pi_ratios_with_thermus.pdf")

NUMERATOR_ORDER = [
    "K+",
    "K-",
    "p",
    "pbar",
    "Ks0",
    "Lambda",
    "Lambda_bar",
    "Xi",
    "Xi_bar",
    "phi",
]

STYLE = {
    "K+": {"marker": "^", "color": "tab:green"},
    "K-": {"marker": "v", "color": "tab:red"},
    "p": {"marker": "D", "color": "tab:orange"},
    "pbar": {"marker": "P", "color": "tab:brown"},
    "Ks0": {"marker": "X", "color": "tab:olive"},
    "Lambda": {"marker": "<", "color": "tab:pink"},
    "Lambda_bar": {"marker": ">", "color": "tab:purple"},
    "Xi": {"marker": "h", "color": "tab:gray"},
    "Xi_bar": {"marker": "8", "color": "tab:gray"},
    "phi": {"marker": "*", "color": "black"},
}

MODEL_CURVE_MIN_ENERGY = {
    "gc_ratio": {"pbar": 7.7, "Lambda_bar": 7.7, "Xi_bar": 7.7},
    "sce_ratio": {"pbar": 7.7, "Lambda_bar": 7.7, "Xi_bar": 7.7},
}


def to_float(text: str | None) -> float | None:
    if text is None or text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def average_pair(a: float, b: float) -> float:
    return 0.5 * (a + b)


def average_pair_err(a_err: float | None, b_err: float | None) -> float | None:
    a_term = a_err or 0.0
    b_term = b_err or 0.0
    if a_term == 0.0 and b_term == 0.0:
        return None
    # Use fully correlated propagation for the pion average because
    # the uncertainties are systematics-dominated.
    return 0.5 * (a_term + b_term)


def ratio_and_err(num: float | None, num_err: float | None, den: float | None, den_err: float | None) -> tuple[float | None, float | None]:
    if num is None or den is None or den <= 0.0:
        return None, None
    ratio = num / den
    if ratio <= 0.0:
        return None, None
    rel2 = 0.0
    has_err = False
    if num_err is not None and num > 0.0:
        rel2 += (num_err / num) ** 2
        has_err = True
    if den_err is not None and den > 0.0:
        rel2 += (den_err / den) ** 2
        has_err = True
    return ratio, ratio * math.sqrt(rel2) if has_err else None


def load_rows():
    by_energy: dict[float, dict[str, dict[str, float | str | None]]] = {}
    with INPUT_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            energy = float(row["energy_GeV"])
            particle = row["particle"]
            by_energy.setdefault(energy, {})[particle] = {
                "data_yield": to_float(row["data_yield"]),
                "data_err": to_float(row["data_err"]),
                "gc_yield": to_float(row["gc_effective_yield"]),
                "sce_yield": to_float(row["sce_blended_yield"]),
            }
    return by_energy


def build_ratio_rows():
    ratio_rows: list[dict[str, str]] = []
    by_energy = load_rows()
    for energy in sorted(by_energy):
        entries = by_energy[energy]
        pi_plus = entries.get("pi+")
        pi_minus = entries.get("pi-")
        if pi_plus is None or pi_minus is None:
            continue

        data_pi = average_pair(to_float(pi_plus["data_yield"]), to_float(pi_minus["data_yield"])) if (
            pi_plus["data_yield"] is not None and pi_minus["data_yield"] is not None
        ) else None
        data_pi_err = average_pair_err(to_float(pi_plus["data_err"]), to_float(pi_minus["data_err"]))

        gc_pi = average_pair(pi_plus["gc_yield"], pi_minus["gc_yield"]) if (
            pi_plus["gc_yield"] is not None and pi_minus["gc_yield"] is not None
        ) else None
        sce_pi = average_pair(pi_plus["sce_yield"], pi_minus["sce_yield"]) if (
            pi_plus["sce_yield"] is not None and pi_minus["sce_yield"] is not None
        ) else None

        for particle in NUMERATOR_ORDER:
            entry = entries.get(particle)
            if entry is None:
                continue
            data_ratio, data_ratio_err = ratio_and_err(
                entry["data_yield"], entry["data_err"], data_pi, data_pi_err
            )
            gc_ratio, _ = ratio_and_err(entry["gc_yield"], None, gc_pi, None)
            sce_ratio, _ = ratio_and_err(entry["sce_yield"], None, sce_pi, None)
            ratio_rows.append(
                {
                    "energy_GeV": f"{energy:g}",
                    "particle": particle,
                    "pi_avg_data_yield": f"{data_pi:.10g}" if data_pi is not None else "",
                    "pi_avg_data_err": f"{data_pi_err:.10g}" if data_pi_err is not None else "",
                    "data_ratio": f"{data_ratio:.10g}" if data_ratio is not None else "",
                    "data_ratio_err": f"{data_ratio_err:.10g}" if data_ratio_err is not None else "",
                    "gc_ratio": f"{gc_ratio:.10g}" if gc_ratio is not None else "",
                    "sce_ratio": f"{sce_ratio:.10g}" if sce_ratio is not None else "",
                }
            )
    return ratio_rows


def write_csv(rows: list[dict[str, str]]) -> None:
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "pi_avg_data_yield",
                "pi_avg_data_err",
                "data_ratio",
                "data_ratio_err",
                "gc_ratio",
                "sce_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 6.4))

    for particle in NUMERATOR_ORDER:
        pts = [r for r in rows if r["particle"] == particle]
        if not pts:
            continue
        pts.sort(key=lambda r: float(r["energy_GeV"]))
        style = STYLE[particle]

        data_pts = [r for r in pts if r["data_ratio"]]
        if data_pts:
            ax.errorbar(
                [float(r["energy_GeV"]) for r in data_pts],
                [float(r["data_ratio"]) for r in data_pts],
                yerr=[float(r["data_ratio_err"]) if r["data_ratio_err"] else 0.0 for r in data_pts],
                linestyle="none",
                marker=style["marker"],
                color=style["color"],
                ms=4.8,
                lw=0.0,
                elinewidth=1.2,
                capsize=2,
            )

        for field, line_style in (("gc_ratio", "-"), ("sce_ratio", "-.")):
            min_energy = MODEL_CURVE_MIN_ENERGY.get(field, {}).get(particle, 0.0)
            model_pts = [r for r in pts if r[field] and float(r["energy_GeV"]) >= min_energy]
            if len(model_pts) >= 2:
                ax.plot(
                    [float(r["energy_GeV"]) for r in model_pts],
                    [float(r[field]) for r in model_pts],
                    line_style,
                    lw=1.7,
                    color=style["color"],
                    alpha=0.95,
                )
            elif len(model_pts) == 1:
                ax.plot(
                    [float(model_pts[0]["energy_GeV"])],
                    [float(model_pts[0][field])],
                    marker="_",
                    ms=10,
                    mew=1.7,
                    linestyle="none",
                    color=style["color"],
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2.5, 300.0)
    ax.set_ylim(1.0e-4, 2.0)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy \; / \; [0.5(\pi^{+}+\pi^{-})]$")
    ax.set_title("Particle-to-Pion Yield Ratios vs Energy with THERMUS", fontsize=12)
    ax.grid(True, which="both", ls="--", alpha=0.28)

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in NUMERATOR_ORDER
    ]
    model_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured ratio"),
        Line2D([0], [0], color="black", lw=1.7, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color="black", lw=1.7, ls="-.", label="THERMUS SCE (endpoint-interpolated)"),
    ]

    leg_particles = ax.legend(
        handles=particle_handles,
        loc="lower right",
        fontsize=8,
        ncol=2,
        frameon=False,
        title="Particle",
        title_fontsize=9,
    )
    ax.add_artist(leg_particles)
    ax.legend(
        handles=model_handles,
        loc="upper left",
        fontsize=8,
        frameon=False,
        title="Model",
        title_fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=240)
    fig.savefig(OUT_PDF)
    plt.close(fig)


def main() -> None:
    rows = build_ratio_rows()
    write_csv(rows)
    make_plot(rows)
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
