#!/usr/bin/env python3
"""Build a comparison of measured yields with THERMUS GC and retained SCE curves."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RAW_CSV = Path("data/first_group_dn_dy_vs_energy.csv")
STRANGE_MERGED_CSV = Path("data/strange_hadron_yields_with_thermus.csv")
GC_LOW_DIR = Path("data/thermus_gc_from_sce_effective/prediction_points")
GC_FIT_DIR = Path("data/thermus_fit_points")
SCE_BLEND_DIR = Path("data/thermus_sce_blended/prediction_points")

OUT_PNG = Path("data/strange_hadron_model_comparison_triptych.png")
OUT_PDF = Path("data/strange_hadron_model_comparison_triptych.pdf")

PARTICLE_ORDER = [
    "pi+",
    "pi-",
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
    "pi+": {"marker": "o", "color": "tab:blue"},
    "pi-": {"marker": "s", "color": "tab:cyan"},
    "K+": {"marker": "^", "color": "tab:green"},
    "K-": {"marker": "v", "color": "tab:red"},
    "p": {"marker": "D", "color": "tab:orange"},
    "pbar": {"marker": "P", "color": "tab:brown"},
    "Ks0": {"marker": "X", "color": "tab:olive"},
    "Lambda": {"marker": "<", "color": "tab:pink"},
    "Lambda_bar": {"marker": ">", "color": "tab:purple"},
    "Xi": {"marker": "h", "color": "tab:gray"},
    "Xi_bar": {"marker": "8", "color": "dimgray"},
    "phi": {"marker": "*", "color": "black"},
}


def to_float(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def combined_error(stat: str | None, sys: str | None) -> float | None:
    stat_val = to_float(stat) or 0.0
    sys_val = to_float(sys) or 0.0
    if stat_val == 0.0 and sys_val == 0.0:
        return None
    return math.hypot(stat_val, sys_val)


def normalize_particle(name: str) -> str | None:
    if name in PARTICLE_ORDER:
        return name
    return None


def load_measured_rows() -> dict[tuple[float, str], dict[str, str]]:
    merged: dict[tuple[float, str], dict[str, str]] = {}

    with RAW_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle"])
            if particle is None:
                continue
            err = combined_error(row.get("stat_err"), row.get("sys_err"))
            key = (float(row["energy_GeV"]), particle)
            merged[key] = {
                "energy_GeV": f"{key[0]:g}",
                "particle": particle,
                "data_yield": row["value"],
                "data_err": f"{err:.10g}" if err is not None else "",
                "data_source": row.get("source", ""),
                "gc_model_yield": "",
                "sce_blended_yield": "",
                "note": row.get("note", ""),
            }

    # Supplement strange-hadron measurements that were carried from THERMUS fit-point
    # tables into the merged strange-hadron CSV.
    with STRANGE_MERGED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle"])
            if particle is None or not row.get("data_yield"):
                continue
            key = (float(row["energy_GeV"]), particle)
            if key in merged:
                continue
            merged[key] = {
                "energy_GeV": f"{key[0]:g}",
                "particle": particle,
                "data_yield": row["data_yield"],
                "data_err": row.get("data_err", ""),
                "data_source": row.get("data_source", ""),
                "gc_model_yield": "",
                "sce_blended_yield": "",
                "note": row.get("note", ""),
            }

    return merged


def attach_fit_point_rows(merged: dict[tuple[float, str], dict[str, str]], model_dir: Path, column: str) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = normalize_particle(row["particle_name"])
                if particle is None:
                    continue
                key = (float(row["energy_GeV"]), particle)
                if key not in merged:
                    merged[key] = {
                        "energy_GeV": f"{key[0]:g}",
                        "particle": particle,
                        "data_yield": row.get("data_yield", ""),
                        "data_err": row.get("data_err", ""),
                        "data_source": str(path),
                        "gc_model_yield": "",
                        "sce_blended_yield": "",
                        "note": "model-only point",
                    }
                merged[key][column] = row.get("model_yield", "")


def attach_prediction_rows(merged: dict[tuple[float, str], dict[str, str]], model_dir: Path, column: str) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = normalize_particle(row["particle_name"])
                if particle is None:
                    continue
                key = (float(row["energy_GeV"]), particle)
                if key not in merged:
                    merged[key] = {
                        "energy_GeV": f"{key[0]:g}",
                        "particle": particle,
                        "data_yield": "",
                        "data_err": "",
                        "data_source": "",
                        "gc_model_yield": "",
                        "sce_blended_yield": "",
                        "note": "model-only point",
                    }
                merged[key][column] = row.get("model_yield", "")


def draw_panel(ax: plt.Axes, rows: list[dict[str, str]], title: str, model_field: str, line_style: str) -> None:
    for particle in PARTICLE_ORDER:
        pts = [r for r in rows if r["particle"] == particle]
        if not pts:
            continue
        pts.sort(key=lambda r: float(r["energy_GeV"]))
        style = STYLE[particle]

        data_pts = [r for r in pts if r["data_yield"]]
        if data_pts:
            ax.errorbar(
                [float(r["energy_GeV"]) for r in data_pts],
                [float(r["data_yield"]) for r in data_pts],
                yerr=[float(r["data_err"]) if r["data_err"] else 0.0 for r in data_pts],
                marker=style["marker"],
                ms=4.2,
                lw=0.0,
                capsize=2,
                linestyle="none",
                color=style["color"],
                alpha=0.95,
            )

        model_pts = [r for r in pts if r[model_field]]
        if len(model_pts) >= 2:
            ax.plot(
                [float(r["energy_GeV"]) for r in model_pts],
                [float(r[model_field]) for r in model_pts],
                line_style,
                lw=1.7,
                color=style["color"],
                alpha=0.95,
            )
        elif len(model_pts) == 1:
            ax.plot(
                [float(model_pts[0]["energy_GeV"])],
                [float(model_pts[0][model_field])],
                marker="_",
                ms=10,
                mew=1.7,
                linestyle="none",
                color=style["color"],
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title, fontsize=11)
    ax.grid(True, which="both", ls="--", alpha=0.28)
    ax.set_xlim(2.5, 300.0)
    ax.set_ylim(5.0e-7, 1.0e3)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")


def main() -> None:
    merged = load_measured_rows()
    attach_prediction_rows(merged, GC_LOW_DIR, "gc_model_yield")
    attach_fit_point_rows(merged, GC_FIT_DIR, "gc_model_yield")
    attach_prediction_rows(merged, SCE_BLEND_DIR, "sce_blended_yield")

    rows = sorted(
        merged.values(),
        key=lambda r: (float(r["energy_GeV"]), PARTICLE_ORDER.index(r["particle"])),
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 6.1), sharey=True)
    draw_panel(
        axes[0],
        rows,
        "Measured + THERMUS GC\n(using effective SCE trace)",
        "gc_model_yield",
        "-",
    )
    draw_panel(
        axes[1],
        rows,
        "Measured + THERMUS SCE\n(endpoint-interpolated 3.0 to 7.7 GeV)",
        "sce_blended_yield",
        "-.",
    )
    axes[0].set_ylabel(r"$dN/dy$")

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in PARTICLE_ORDER
    ]
    model_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured yield"),
        Line2D([0], [0], color="black", lw=1.7, ls="-", label="THERMUS GC (effective trace)"),
        Line2D(
            [0], [0], color="black", lw=1.7, ls="-.",
            label="THERMUS SCE endpoint-interpolated"
        ),
    ]

    leg_particles = axes[0].legend(
        handles=particle_handles,
        loc="upper left",
        fontsize=8,
        ncol=2,
        frameon=False,
        title="Particle",
        title_fontsize=9,
    )
    axes[0].add_artist(leg_particles)
    axes[1].legend(
        handles=model_handles,
        loc="lower right",
        fontsize=8,
        frameon=False,
        title="Model",
        title_fontsize=9,
    )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=240)
    fig.savefig(OUT_PDF)
    plt.close(fig)

    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
