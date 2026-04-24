#!/usr/bin/env python3
"""Plot measured strange-hadron yields with low-energy THERMUS SCE predictions."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MEASURED_CSV = Path("data/strange_hadron_yields_with_thermus.csv")

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


def normalize_particle(name: str) -> str | None:
    mapping = {
        "K+": "K+",
        "K-": "K-",
        "Ks0": "Ks0",
        "Lambda": "Lambda",
        "Lambda_bar": "Lambda_bar",
        "Xi": "Xi",
        "Xi_bar": "Xi_bar",
        "phi": "phi",
    }
    return mapping.get(name)


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot low-energy SCE predictions against measured yields.")
    parser.add_argument("--model-dir", default="data/thermus_sce_low_energy/prediction_points")
    parser.add_argument("--outdir", default="data/thermus_sce_low_energy")
    parser.add_argument("--title", default="Strange-Hadron Yields vs Energy with THERMUS SCE")
    parser.add_argument("--model-label", default="THERMUS SCE")
    parser.add_argument("--measured-csv", default=str(MEASURED_CSV))
    return parser.parse_args()


def load_rows(measured_csv: Path) -> dict[tuple[float, str], dict[str, str]]:
    merged: dict[tuple[float, str], dict[str, str]] = {}
    with measured_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = set(reader.fieldnames or [])
        for row in reader:
            particle = normalize_particle(row["particle"])
            if particle is None:
                continue
            key = (float(row["energy_GeV"]), particle)
            if "data_yield" in fieldnames:
                data_yield = row.get("data_yield", "")
                data_err = row.get("data_err", "")
                data_source = row.get("data_source", "")
            else:
                data_yield = row.get("value", "")
                err = combined_error(row.get("stat_err"), row.get("sys_err"))
                data_err = f"{err:.10g}" if err is not None else ""
                data_source = row.get("source", "")
            merged[key] = {
                "energy_GeV": f"{key[0]:g}",
                "particle": particle,
                "data_yield": data_yield,
                "data_err": data_err,
                "gc_effective_model_yield": row.get("gc_effective_model_yield", row.get("gc_model_yield", "")),
                "sce_model_yield": "",
                "data_source": data_source,
                "model_source": "",
                "note": row.get("note", ""),
            }
    return merged


def attach_model_rows(merged: dict[tuple[float, str], dict[str, str]], model_dir: Path, model_label: str) -> None:
    for path in sorted(model_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = normalize_particle(row["particle_name"])
                if particle is None:
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
                        "sce_model_yield": "",
                        "data_source": "",
                        "model_source": "",
                        "note": "model-only point",
                    }
                merged[key]["sce_model_yield"] = row["model_yield"]
                merged[key]["model_source"] = model_label


def write_csv(rows: list[dict[str, str]], out_csv: Path) -> None:
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "data_yield",
                "data_err",
                "gc_effective_model_yield",
                "sce_model_yield",
                "data_source",
                "model_source",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows: list[dict[str, str]], out_png: Path, out_pdf: Path, title: str, model_label: str) -> None:
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
            )

        model_specs = [
            ("gc_effective_model_yield", "-"),
            ("sce_model_yield", "-."),
        ]
        for field, linestyle in model_specs:
            model_pts = [r for r in pts if r[field]]
            if len(model_pts) >= 2:
                ax.plot(
                    [float(r["energy_GeV"]) for r in model_pts],
                    [float(r[field]) for r in model_pts],
                    linestyle,
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
                    mew=1.8,
                    linestyle="none",
                    color=style["color"],
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title(title)
    ax.grid(True, which="both", ls="--", alpha=0.35)

    particle_handles = [
        Line2D([0], [0], marker=STYLE[p]["marker"], color=STYLE[p]["color"], lw=0, markersize=6, label=p)
        for p in PARTICLE_ORDER
        if any(r["particle"] == p for r in rows)
    ]
    style_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label="Measured yield"),
        Line2D([0], [0], color="black", lw=1.7, ls="-", label="THERMUS GC (effective trace)"),
        Line2D([0], [0], color="black", lw=1.7, ls="-.", label=model_label),
    ]
    leg1 = ax.legend(handles=particle_handles, fontsize=9, ncol=2, loc="upper left", title="Particle")
    ax.add_artist(leg1)
    ax.legend(handles=style_handles, fontsize=9, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / "strange_hadron_yields_with_sce.csv"
    out_png = outdir / "strange_hadron_yields_with_sce.png"
    out_pdf = outdir / "strange_hadron_yields_with_sce.pdf"
    merged = load_rows(Path(args.measured_csv))
    attach_model_rows(merged, Path(args.model_dir), args.model_label)
    rows = sorted(
        merged.values(),
        key=lambda r: (float(r["energy_GeV"]), PARTICLE_ORDER.index(r["particle"])),
    )
    write_csv(rows, out_csv)
    make_plot(rows, out_png, out_pdf, args.title, args.model_label)
    print(f"wrote {out_csv}")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
