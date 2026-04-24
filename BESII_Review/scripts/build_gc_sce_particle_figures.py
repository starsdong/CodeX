#!/usr/bin/env python3
"""Build low-energy per-particle figures with measured data and retained THERMUS GC/SCE curves."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RAW_CSV = Path("data/first_group_dn_dy_vs_energy.csv")
GC_LOW_DIR = Path("data/thermus_gc_from_sce_effective/prediction_points")
GC_FIT_DIR = Path("data/thermus_fit_predictions_full/prediction_points")
SCE_BLEND_DIR = Path("data/thermus_sce_blended/prediction_points")

OUT_DIR = Path("data/gc_sce_particle_figures")
MANIFEST = OUT_DIR / "manifest.csv"
POINTS_CSV = OUT_DIR / "particle_model_points.csv"

PARTICLE_ORDER = ["p", "K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"]
STYLE = {
    "p": {"marker": "D", "color": "tab:orange"},
    "K+": {"marker": "o", "color": "tab:blue"},
    "K-": {"marker": "s", "color": "tab:cyan"},
    "Ks0": {"marker": "D", "color": "tab:green"},
    "Lambda": {"marker": "^", "color": "tab:red"},
    "Lambda_bar": {"marker": "v", "color": "tab:orange"},
    "Xi": {"marker": "P", "color": "tab:purple"},
    "Xi_bar": {"marker": "X", "color": "tab:pink"},
    "phi": {"marker": "*", "color": "tab:brown"},
}
CURVE_MIN_ENERGY = {"Lambda_bar": 7.7, "Xi_bar": 7.7}

XMIN = 2.0
XMAX = 10.0


def to_float(text: str | None) -> float | None:
    if text is None or text == "":
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


def empty_row(energy: float, particle: str) -> dict[str, str]:
    return {
        "energy_GeV": f"{energy:g}",
        "particle": particle,
        "data_yield": "",
        "data_err": "",
        "gc_model_yield": "",
        "sce_blended_yield": "",
        "data_source": "",
        "note": "",
    }


def load_measured_rows() -> dict[tuple[float, str], dict[str, str]]:
    merged: dict[tuple[float, str], dict[str, str]] = {}
    with RAW_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in PARTICLE_ORDER:
                continue
            energy = float(row["energy_GeV"])
            key = (energy, particle)
            err = combined_error(row.get("stat_err"), row.get("sys_err"))
            out = empty_row(energy, particle)
            out["data_yield"] = row["value"]
            out["data_err"] = f"{err:.10g}" if err is not None else ""
            out["data_source"] = row.get("source", "")
            out["note"] = row.get("note", "")
            merged[key] = out
    return merged


def attach_prediction_rows(
    merged: dict[tuple[float, str], dict[str, str]],
    model_dir: Path,
    column: str,
) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = row["particle_name"]
                if particle not in PARTICLE_ORDER:
                    continue
                energy = float(row["energy_GeV"])
                key = (energy, particle)
                if key not in merged:
                    merged[key] = empty_row(energy, particle)
                merged[key][column] = row["model_yield"]


def attach_fit_point_rows(
    merged: dict[tuple[float, str], dict[str, str]],
    model_dir: Path,
    column: str,
) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = row["particle_name"]
                if particle not in PARTICLE_ORDER:
                    continue
                energy = float(row["energy_GeV"])
                key = (energy, particle)
                if key not in merged:
                    merged[key] = empty_row(energy, particle)
                merged[key][column] = row["model_yield"]


def write_points_csv(rows: list[dict[str, str]]) -> None:
    with POINTS_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "data_yield",
                "data_err",
                "gc_model_yield",
                "sce_blended_yield",
                "data_source",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def nice_limits(values: list[float]) -> tuple[float, float]:
    positive = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positive:
        return 1.0e-3, 1.0
    ymin = min(positive)
    ymax = max(positive)
    lower = 10 ** math.floor(math.log10(ymin) - 0.2)
    upper = 10 ** math.ceil(math.log10(ymax) + 0.2)
    if lower == upper:
        upper *= 10.0
    return lower, upper


def build_particle_figure(rows: list[dict[str, str]], particle: str) -> tuple[Path, Path]:
    pts = [r for r in rows if r["particle"] == particle and XMIN <= float(r["energy_GeV"]) <= XMAX]
    pts.sort(key=lambda r: float(r["energy_GeV"]))
    style = STYLE[particle]

    fig, ax = plt.subplots(figsize=(7.6, 5.4))

    data_pts = [r for r in pts if r["data_yield"]]
    if data_pts:
        ax.errorbar(
            [float(r["energy_GeV"]) for r in data_pts],
            [float(r["data_yield"]) for r in data_pts],
            yerr=[float(r["data_err"]) if r["data_err"] else 0.0 for r in data_pts],
            marker=style["marker"],
            ms=5.5,
            lw=0.0,
            capsize=2,
            linestyle="none",
            color=style["color"],
            label="Measured yield",
        )

    model_specs = [
        ("gc_model_yield", "-", "THERMUS GC (effective trace)"),
        ("sce_blended_yield", "-.", "THERMUS SCE (endpoint-interpolated)"),
    ]
    y_values: list[float] = [float(r["data_yield"]) for r in data_pts if r["data_yield"]]
    for field, linestyle, label in model_specs:
        min_energy = CURVE_MIN_ENERGY.get(particle, 0.0)
        model_pts = [r for r in pts if r[field] and float(r["energy_GeV"]) >= min_energy]
        if len(model_pts) >= 2:
            xs = [float(r["energy_GeV"]) for r in model_pts]
            ys = [float(r[field]) for r in model_pts]
            y_values.extend(ys)
            ax.plot(xs, ys, linestyle, lw=1.9, color=style["color"], label=label)
        elif len(model_pts) == 1:
            y_values.append(float(model_pts[0][field]))
            ax.plot(
                [float(model_pts[0]["energy_GeV"])],
                [float(model_pts[0][field])],
                marker="_",
                ms=10,
                mew=1.8,
                linestyle="none",
                color=style["color"],
                label=label,
            )

    ymin, ymax = nice_limits(y_values)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(max(3.0e-3, ymin), ymax)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title(f"{particle} yield vs energy")
    ax.grid(True, which="both", ls="--", alpha=0.35)

    handles = [
        Line2D([0], [0], marker=style["marker"], color=style["color"], lw=0, markersize=6, label="Measured yield"),
        Line2D([0], [0], color=style["color"], lw=1.9, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color=style["color"], lw=1.9, ls="-.", label="THERMUS SCE (endpoint-interpolated)"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=8, frameon=False)

    fig.tight_layout()
    png_path = OUT_DIR / f"{particle}_yield_models_2to10GeV.png"
    pdf_path = OUT_DIR / f"{particle}_yield_models_2to10GeV.pdf"
    fig.savefig(png_path, dpi=240)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def write_manifest(entries: list[dict[str, str]]) -> None:
    with MANIFEST.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["particle", "png_file", "pdf_file", "x_range", "y_scale", "notes"],
        )
        writer.writeheader()
        writer.writerows(entries)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged = load_measured_rows()
    attach_prediction_rows(merged, GC_LOW_DIR, "gc_model_yield")
    attach_fit_point_rows(merged, GC_FIT_DIR, "gc_model_yield")
    attach_prediction_rows(merged, SCE_BLEND_DIR, "sce_blended_yield")

    rows = sorted(
        merged.values(),
        key=lambda r: (PARTICLE_ORDER.index(r["particle"]), float(r["energy_GeV"])),
    )
    write_points_csv(rows)

    manifest_rows = []
    for particle in PARTICLE_ORDER:
        png_path, pdf_path = build_particle_figure(rows, particle)
        manifest_rows.append(
            {
                "particle": particle,
                "png_file": str(png_path),
                "pdf_file": str(pdf_path),
                "x_range": f"{XMIN:g}-{XMAX:g} GeV",
                "y_scale": "log auto",
                "notes": "Measured data with GC effective-trace and SCE endpoint-interpolated curves",
            }
        )

    write_manifest(manifest_rows)
    print(f"wrote {POINTS_CSV}")
    print(f"wrote {MANIFEST}")
    print(f"wrote {len(manifest_rows)} particle figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
