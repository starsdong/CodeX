#!/usr/bin/env python3
"""Overlay measured strange-hadron yields with GC and low-energy SCE THERMUS calculations."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MEASURED_CSV = Path("data/strange_hadron_yields_vs_energy.csv")
GC_LOW_DIR = Path("data/thermus_gc_from_sce_effective/prediction_points")
GC_FIT_DIR = Path("data/thermus_fit_points")
SCE_BLEND_CSV = Path("data/thermus_sce_blended/strange_hadron_yields_with_sce.csv")
OUT_CSV = Path("data/strange_hadron_yields_gc_sce_overlay.csv")
OUT_PNG = Path("data/strange_hadron_yields_gc_sce_overlay.png")
OUT_PDF = Path("data/strange_hadron_yields_gc_sce_overlay.pdf")
ZOOM_OUT_PNG = Path("data/strange_hadron_yields_gc_sce_overlay_zoom.png")
ZOOM_OUT_PDF = Path("data/strange_hadron_yields_gc_sce_overlay_zoom.pdf")

PARTICLE_ORDER = ["K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"]
ANTI_PARTICLES = {"K-", "Lambda_bar", "Xi_bar"}
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


def load_measured() -> dict[tuple[float, str], dict[str, str]]:
    rows: dict[tuple[float, str], dict[str, str]] = {}
    with MEASURED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (float(row["energy_GeV"]), row["particle"])
            rows[key] = {
                "energy_GeV": f"{key[0]:g}",
                "particle": row["particle"],
                "data_yield": row["value"],
                "data_err": row["stat_err"],
                "gc_model_yield": "",
                "sce_blended_model_yield": "",
            }
    return rows


def attach_model(rows: dict[tuple[float, str], dict[str, str]], path: Path, column: str, model_field: str) -> None:
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            key = (float(row["energy_GeV"]), particle)
            if key not in rows:
                rows[key] = {
                    "energy_GeV": f"{key[0]:g}",
                    "particle": particle,
                    "data_yield": "",
                    "data_err": "",
                    "gc_model_yield": "",
                    "sce_blended_model_yield": "",
                }
            rows[key][column] = row[model_field]


def attach_prediction_rows(rows: dict[tuple[float, str], dict[str, str]], model_dir: Path, column: str) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = row["particle_name"]
                if particle not in PARTICLE_ORDER:
                    continue
                key = (float(row["energy_GeV"]), particle)
                if key not in rows:
                    rows[key] = {
                        "energy_GeV": f"{key[0]:g}",
                        "particle": particle,
                        "data_yield": "",
                        "data_err": "",
                        "gc_model_yield": "",
                        "sce_blended_model_yield": "",
                    }
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
    out_png: Path,
    out_pdf: Path,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    title_suffix: str = "",
    suppress_curve_particles: set[str] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 6.4))
    suppress_curve_particles = suppress_curve_particles or set()

    for particle in PARTICLE_ORDER:
        pts = [r for r in rows if r["particle"] == particle]
        if not pts:
            continue
        pts.sort(key=lambda r: float(r["energy_GeV"]))
        style = STYLE[particle]

        data_pts = [r for r in pts if r["data_yield"]]
        if data_pts:
            ax.plot(
                [float(r["energy_GeV"]) for r in data_pts],
                [float(r["data_yield"]) for r in data_pts],
                linestyle="none",
                marker=style["marker"],
                color=style["color"],
                ms=5.0,
            )

        gc_pts = [r for r in pts if r["gc_model_yield"]]
        if len(gc_pts) >= 2 and particle not in suppress_curve_particles:
            ax.plot(
                [float(r["energy_GeV"]) for r in gc_pts],
                [float(r["gc_model_yield"]) for r in gc_pts],
                "-",
                lw=1.8,
                color=style["color"],
            )

        sce_blended_pts = [r for r in pts if r["sce_blended_model_yield"]]
        if len(sce_blended_pts) >= 2 and particle not in suppress_curve_particles:
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
    ax.set_title("Strange-Hadron Yields: Measured, GC Effective Trace, and SCE" + title_suffix)
    ax.grid(True, which="both", ls="--", alpha=0.35)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in PARTICLE_ORDER
        if any(r["particle"] == p for r in rows)
    ]
    model_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured yield"),
        Line2D([0], [0], color="black", lw=1.8, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color="black", lw=1.8, ls="-.", label="THERMUS SCE (endpoint-interpolated)"),
    ]
    leg1 = ax.legend(handles=particle_handles, fontsize=9, ncol=2, loc="upper left", title="Particle")
    ax.add_artist(leg1)
    ax.legend(handles=model_handles, fontsize=9, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    rows = load_measured()
    attach_prediction_rows(rows, GC_LOW_DIR, "gc_model_yield")
    attach_prediction_rows(rows, GC_FIT_DIR, "gc_model_yield")
    attach_model(rows, SCE_BLEND_CSV, "sce_blended_model_yield", "sce_model_yield")
    ordered = sorted(
        rows.values(),
        key=lambda r: (float(r["energy_GeV"]), PARTICLE_ORDER.index(r["particle"])),
    )
    write_csv(ordered)
    make_plot(ordered, OUT_PNG, OUT_PDF)
    make_plot(
        ordered,
        ZOOM_OUT_PNG,
        ZOOM_OUT_PDF,
        xlim=(2.0, 10.0),
        ylim=(1.0e-3, 1.0e2),
        title_suffix=" (3.0 to 7.7 GeV Zoom)",
        suppress_curve_particles=ANTI_PARTICLES,
    )
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")
    print(f"wrote {ZOOM_OUT_PNG}")
    print(f"wrote {ZOOM_OUT_PDF}")


if __name__ == "__main__":
    main()
