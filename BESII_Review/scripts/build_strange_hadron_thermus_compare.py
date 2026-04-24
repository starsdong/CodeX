#!/usr/bin/env python3
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MEASURED_CSV = Path("data/strange_hadron_yields_vs_energy.csv")
THERMUS_DIR = Path("data/thermus_fit_points")
THERMUS_GC_EFFECTIVE_DIR = Path("data/thermus_gc_from_sce_effective/prediction_points")
THERMUS_SCE_BLEND_CSV = Path("data/thermus_sce_blended/strange_hadron_yields_with_sce.csv")
OUT_CSV = Path("data/strange_hadron_yields_with_thermus.csv")
OUT_PNG = Path("data/strange_hadron_yields_with_thermus.png")
OUT_PDF = Path("data/strange_hadron_yields_with_thermus.pdf")

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
THERMUS_PARTICLES = {"K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"}


def to_float(text):
    try:
        return float(text)
    except Exception:
        return None


def normalize_particle(name):
    mapping = {
        "Lambda_bar": "Lambda_bar",
        "Xi_bar": "Xi_bar",
        "Ks0": "Ks0",
        "phi": "phi",
        "Lambda": "Lambda",
        "Xi": "Xi",
        "K+": "K+",
        "K-": "K-",
    }
    return mapping.get(name)


def combined_error(stat, sys):
    stat_val = to_float(stat) or 0.0
    sys_val = to_float(sys) or 0.0
    return math.hypot(stat_val, sys_val) if (stat_val or sys_val) else 0.0


def load_measured_rows():
    merged = {}
    with MEASURED_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (float(row["energy_GeV"]), row["particle"])
            merged[key] = {
                "energy_GeV": f"{float(row['energy_GeV']):g}",
                "particle": row["particle"],
                "data_yield": row["value"],
                "data_err": f"{combined_error(row.get('stat_err'), row.get('sys_err')):.10g}",
                "gc_effective_model_yield": "",
                "sce_blended_model_yield": "",
                "data_source": row["source"],
                "model_source": "",
                "note": row.get("note", ""),
            }
    return merged


def attach_thermus_fit_points(merged, path, column, model_source):
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle_name"])
            if particle not in THERMUS_PARTICLES:
                continue
            energy = float(row["energy_GeV"])
            key = (energy, particle)
            if key not in merged:
                merged[key] = {
                    "energy_GeV": f"{energy:g}",
                    "particle": particle,
                    "data_yield": row["data_yield"],
                    "data_err": row["data_err"],
                    "gc_effective_model_yield": "",
                    "sce_blended_model_yield": "",
                    "data_source": f"{path}",
                    "model_source": "",
                    "note": "measured yield carried from THERMUS fit-point table",
                }
            merged[key][column] = row["model_yield"]
            merged[key]["model_source"] = model_source


def attach_thermus_prediction_points(merged, path, column, model_source):
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle_name"])
            if particle not in THERMUS_PARTICLES:
                continue
            energy = float(row["energy_GeV"])
            key = (energy, particle)
            if key not in merged:
                merged[key] = {
                    "energy_GeV": f"{energy:g}",
                    "particle": particle,
                    "data_yield": "",
                    "data_err": "",
                    "gc_effective_model_yield": "",
                    "sce_blended_model_yield": "",
                    "data_source": "",
                    "model_source": "",
                    "note": "model-only point",
                }
            merged[key][column] = row["model_yield"]
            merged[key]["model_source"] = model_source


def attach_sce_csv(merged, path, column, model_source):
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle"])
            if particle not in THERMUS_PARTICLES:
                continue
            energy = float(row["energy_GeV"])
            key = (energy, particle)
            if key not in merged:
                merged[key] = {
                    "energy_GeV": f"{energy:g}",
                    "particle": particle,
                    "data_yield": row.get("data_yield", ""),
                    "data_err": row.get("data_err", ""),
                    "gc_effective_model_yield": "",
                    "sce_blended_model_yield": "",
                    "data_source": row.get("data_source", ""),
                    "model_source": "",
                    "note": row.get("note", "model-only point"),
                }
            merged[key][column] = row["sce_model_yield"]
            merged[key]["model_source"] = model_source


def write_csv(rows):
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "data_yield",
                "data_err",
                "gc_effective_model_yield",
                "sce_blended_model_yield",
                "data_source",
                "model_source",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows):
    fig, ax = plt.subplots(figsize=(9.6, 6.2))

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
                yerr=[float(r["data_err"]) for r in data_pts],
                marker=style["marker"],
                ms=5.0,
                lw=0.0,
                capsize=2,
                linestyle="none",
                color=style["color"],
                label=particle,
            )

        model_specs = [
            ("gc_effective_model_yield", "-"),
            ("sce_blended_model_yield", "-."),
        ]
        for field, linestyle in model_specs:
            model_pts = [r for r in pts if r[field]]
            if len(model_pts) >= 2:
                ax.plot(
                    [float(r["energy_GeV"]) for r in model_pts],
                    [float(r[field]) for r in model_pts],
                    linestyle,
                    lw=1.6,
                    color=style["color"],
                    alpha=0.95,
                )
            elif len(model_pts) == 1:
                ax.plot(
                    [float(model_pts[0]["energy_GeV"])],
                    [float(model_pts[0][field])],
                    marker="_",
                    ms=10,
                    mew=1.8,
                    linestyle="none",
                    color=style["color"],
                    alpha=0.95,
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title("Strange-Hadron Yields vs Energy with THERMUS")
    ax.grid(True, which="both", ls="--", alpha=0.35)

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in PARTICLE_ORDER
        if any(r["particle"] == p for r in rows)
    ]
    style_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured yield"),
        Line2D([0], [0], color="black", lw=1.6, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color="black", lw=1.6, ls="-.", label="THERMUS SCE (endpoint-interpolated)"),
    ]
    leg1 = ax.legend(handles=particle_handles, fontsize=9, ncol=2, loc="upper left", title="Particle")
    ax.add_artist(leg1)
    ax.legend(handles=style_handles, fontsize=9, loc="lower right")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=240)
    fig.savefig(OUT_PDF)
    plt.close(fig)


def main():
    merged = load_measured_rows()
    for path in sorted(THERMUS_GC_EFFECTIVE_DIR.glob("sqrts_*GeV_points.csv")):
        attach_thermus_prediction_points(merged, path, "gc_effective_model_yield", "THERMUS GC effective trace")
    for path in sorted(THERMUS_DIR.glob("sqrts_*GeV_points.csv")):
        attach_thermus_fit_points(merged, path, "gc_effective_model_yield", "THERMUS GC fit")
    attach_sce_csv(merged, THERMUS_SCE_BLEND_CSV, "sce_blended_model_yield", "THERMUS SCE endpoint-interpolated")

    rows = sorted(
        merged.values(),
        key=lambda r: (float(r["energy_GeV"]), PARTICLE_ORDER.index(r["particle"])),
    )
    write_csv(rows)
    make_plot(rows)
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
