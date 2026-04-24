#!/usr/bin/env python3
"""Overlay measured yields with retained GC and low-energy SCE THERMUS calculations."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MEASURED_CSV = Path("data/strange_hadron_yields_vs_energy.csv")
RAW_MEASURED_CSV = Path("data/first_group_dn_dy_vs_energy.csv")
GC_LOW_DIR = Path("data/thermus_gc_from_sce_effective/prediction_points")
GC_FIT_FULL_DIR = Path("data/thermus_fit_predictions_full/prediction_points")
SCE_BLEND_CSV = Path("data/thermus_sce_blended/strange_hadron_yields_with_sce.csv")
OUT_CSV = Path("data/strange_hadron_yields_gc_sce_overlay.csv")
OUT_PNG = Path("data/strange_hadron_yields_gc_sce_overlay.png")
OUT_PDF = Path("data/strange_hadron_yields_gc_sce_overlay.pdf")
ZOOM_OUT_PNG = Path("data/strange_hadron_yields_gc_sce_overlay_zoom.png")
ZOOM_OUT_PDF = Path("data/strange_hadron_yields_gc_sce_overlay_zoom.pdf")

FULL_PARTICLE_ORDER = ["pi+", "pi-", "p", "K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"]
ZOOM_PARTICLE_ORDER = ["p", "K+", "K-", "Ks0", "Lambda", "Xi", "phi"]
STYLE = {
    "pi+": {"marker": "o", "color": "tab:blue"},
    "pi-": {"marker": "s", "color": "tab:cyan"},
    "p": {"marker": "D", "color": "tab:orange"},
    "K+": {"marker": "^", "color": "tab:green"},
    "K-": {"marker": "v", "color": "tab:red"},
    "Ks0": {"marker": "X", "color": "tab:olive"},
    "Lambda": {"marker": "<", "color": "tab:pink"},
    "Lambda_bar": {"marker": "v", "color": "tab:orange"},
    "Xi": {"marker": "P", "color": "tab:purple"},
    "Xi_bar": {"marker": "X", "color": "tab:pink"},
    "phi": {"marker": "*", "color": "tab:brown"},
}


def empty_row(energy: float, particle: str) -> dict[str, str]:
    return {
        "energy_GeV": f"{energy:g}",
        "particle": particle,
        "data_yield": "",
        "data_err": "",
        "gc_model_yield": "",
        "sce_blended_model_yield": "",
    }


def to_float(text: str | None) -> float | None:
    if text is None or text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def combined_error(stat: str | None, sys: str | None) -> str:
    stat_val = to_float(stat) or 0.0
    sys_val = to_float(sys) or 0.0
    if stat_val == 0.0 and sys_val == 0.0:
        return ""
    return f"{math.hypot(stat_val, sys_val):.10g}"


def load_measured() -> dict[tuple[float, str], dict[str, str]]:
    rows: dict[tuple[float, str], dict[str, str]] = {}
    with RAW_MEASURED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in FULL_PARTICLE_ORDER:
                continue
            key = (float(row["energy_GeV"]), particle)
            rows[key] = empty_row(key[0], particle)
            rows[key]["data_yield"] = row["value"]
            rows[key]["data_err"] = combined_error(row.get("stat_err"), row.get("sys_err"))

    with MEASURED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in FULL_PARTICLE_ORDER:
                continue
            key = (float(row["energy_GeV"]), particle)
            if key in rows:
                continue
            rows[key] = empty_row(key[0], particle)
            rows[key]["data_yield"] = row["value"]
            rows[key]["data_err"] = combined_error(row.get("stat_err"), row.get("sys_err"))
    return rows


def load_zoom_measured() -> dict[tuple[float, str], dict[str, str]]:
    rows: dict[tuple[float, str], dict[str, str]] = {}
    with RAW_MEASURED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in ZOOM_PARTICLE_ORDER:
                continue
            key = (float(row["energy_GeV"]), particle)
            rows[key] = empty_row(key[0], particle)
            rows[key]["data_yield"] = row["value"]
            rows[key]["data_err"] = combined_error(row.get("stat_err"), row.get("sys_err"))
    return rows


def attach_model(
    rows: dict[tuple[float, str], dict[str, str]],
    path: Path,
    column: str,
    model_field: str,
    allowed_particles: set[str],
) -> None:
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in allowed_particles:
                continue
            key = (float(row["energy_GeV"]), particle)
            if key not in rows:
                rows[key] = empty_row(key[0], particle)
            rows[key][column] = row[model_field]


def attach_prediction_rows(
    rows: dict[tuple[float, str], dict[str, str]],
    model_dir: Path,
    column: str,
    allowed_particles: set[str],
) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = row["particle_name"]
                if particle not in allowed_particles:
                    continue
                key = (float(row["energy_GeV"]), particle)
                if key not in rows:
                    rows[key] = empty_row(key[0], particle)
                rows[key][column] = row["model_yield"]


def write_csv(rows: list[dict[str, str]]) -> None:
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "data_yield",
                "data_err",
                "gc_model_yield",
                "sce_blended_model_yield",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def make_plot(
    rows: list[dict[str, str]],
    particle_order: list[str],
    out_png: Path,
    out_pdf: Path,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    title: str = "",
    curve_min_energy_by_particle: dict[str, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 6.4))
    curve_min_energy_by_particle = curve_min_energy_by_particle or {}

    for particle in particle_order:
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
                linestyle="none",
                marker=style["marker"],
                color=style["color"],
                ms=5.0,
                lw=0.0,
                elinewidth=1.2,
                capsize=2,
            )

        min_energy = curve_min_energy_by_particle.get(particle, 0.0)
        gc_pts = [r for r in pts if r["gc_model_yield"] and float(r["energy_GeV"]) >= min_energy]
        if len(gc_pts) >= 2:
            ax.plot(
                [float(r["energy_GeV"]) for r in gc_pts],
                [float(r["gc_model_yield"]) for r in gc_pts],
                "-",
                lw=1.8,
                color=style["color"],
            )

        sce_blended_pts = [
            r for r in pts if r["sce_blended_model_yield"] and float(r["energy_GeV"]) >= min_energy
        ]
        if len(sce_blended_pts) >= 2:
            ax.plot(
                [float(r["energy_GeV"]) for r in sce_blended_pts],
                [float(r["sce_blended_model_yield"]) for r in sce_blended_pts],
                "-.",
                lw=1.8,
                color=style["color"],
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title(title)
    ax.grid(True, which="both", ls="--", alpha=0.35)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in particle_order
        if any(r["particle"] == p for r in rows)
    ]
    model_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured yield"),
        Line2D([0], [0], color="black", lw=1.8, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color="black", lw=1.8, ls="-.", label="THERMUS SCE (endpoint-interpolated)"),
    ]
    leg1 = ax.legend(handles=particle_handles, fontsize=9, ncol=2, loc="lower right", title="Particle")
    ax.add_artist(leg1)
    ax.legend(handles=model_handles, fontsize=9, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    rows = load_measured()
    attach_prediction_rows(rows, GC_LOW_DIR, "gc_model_yield", set(FULL_PARTICLE_ORDER))
    attach_prediction_rows(rows, GC_FIT_FULL_DIR, "gc_model_yield", set(FULL_PARTICLE_ORDER))
    attach_model(rows, SCE_BLEND_CSV, "sce_blended_model_yield", "sce_model_yield", set(FULL_PARTICLE_ORDER))
    ordered = sorted(
        rows.values(),
        key=lambda r: (float(r["energy_GeV"]), FULL_PARTICLE_ORDER.index(r["particle"])),
    )
    write_csv(ordered)

    zoom_rows = load_zoom_measured()
    attach_prediction_rows(zoom_rows, GC_LOW_DIR, "gc_model_yield", set(ZOOM_PARTICLE_ORDER))
    attach_prediction_rows(zoom_rows, GC_FIT_FULL_DIR, "gc_model_yield", set(ZOOM_PARTICLE_ORDER))
    attach_model(
        zoom_rows,
        SCE_BLEND_CSV,
        "sce_blended_model_yield",
        "sce_model_yield",
        set(ZOOM_PARTICLE_ORDER),
    )
    zoom_ordered = sorted(
        zoom_rows.values(),
        key=lambda r: (float(r["energy_GeV"]), ZOOM_PARTICLE_ORDER.index(r["particle"])),
    )

    make_plot(
        ordered,
        FULL_PARTICLE_ORDER,
        OUT_PNG,
        OUT_PDF,
        ylim=(3.0e-3, 1.0e3),
        title="Strange-Hadron Yields: Measured, GC Effective Trace, and SCE",
        curve_min_energy_by_particle={"pbar": 7.7, "Lambda_bar": 7.7, "Xi_bar": 7.7},
    )
    make_plot(
        zoom_ordered,
        ZOOM_PARTICLE_ORDER,
        ZOOM_OUT_PNG,
        ZOOM_OUT_PDF,
        xlim=(2.0, 10.0),
        ylim=(1.0e-3, 1.0e2),
        title="Low-Energy Identified-Hadron Yields: Measured, GC Effective Trace, and SCE",
    )
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")
    print(f"wrote {ZOOM_OUT_PNG}")
    print(f"wrote {ZOOM_OUT_PDF}")


if __name__ == "__main__":
    main()
