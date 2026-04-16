#!/usr/bin/env python3
"""Make per-energy plots comparing experimental yields to THERMUS fitted yields."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_points(path: Path):
    rows = []
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    rows.sort(key=lambda r: int(r["pdg_id"]))
    return rows


def load_fit_summary(path: Path):
    by_energy = {}
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            e = round(float(r["energy_GeV"]), 6)
            by_energy[e] = {
                "T_MeV": float(r["T_MeV"]),
                "T_err_MeV": float(r["T_err_MeV"]),
                "muB_MeV": float(r["muB_MeV"]),
                "muB_err_MeV": float(r["muB_err_MeV"]),
                "muS_MeV": float(r["muS_MeV"]),
                "muS_err_MeV": float(r["muS_err_MeV"]),
                "gammaS": float(r["gammaS"]),
                "gammaS_err": float(r["gammaS_err"]),
                "R_fm": float(r["R_fm"]),
                "R_err_fm": float(r["R_err_fm"]),
                "dVdy_fm3": float(r["dVdy_fm3"]),
                "dVdy_err_fm3": float(r["dVdy_err_fm3"]),
                "chi2": float(r["chi2"]),
                "ndf": int(float(r["ndf"])),
                "chi2_ndf": float(r["chi2_ndf"]),
            }
    return by_energy


def make_plot(rows, out_png: Path, fit_info: dict | None = None):
    if not rows:
        return

    energy = float(rows[0]["energy_GeV"])
    labels = [r["particle_name"] for r in rows]
    data = [float(r["data_yield"]) for r in rows]
    err = [float(r["data_err"]) for r in rows]
    model = [float(r["model_yield"]) for r in rows]
    residual = [((d - m) / e) if e > 0.0 else float("nan") for d, m, e in zip(data, model, err)]

    x = list(range(len(labels)))
    fig, (ax, axr) = plt.subplots(
        2,
        1,
        figsize=(max(10, 0.62 * len(labels)), 7.2),
        constrained_layout=True,
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.3], "hspace": 0.05},
    )

    ax.errorbar(x, data, yerr=err, fmt="o", ms=5, capsize=2, lw=1.0, color="C0", label="Data")
    halfw = 0.35
    for i, yi in enumerate(model):
        label = "THERMUS fit" if i == 0 else None
        ax.hlines(yi, i - halfw, i + halfw, colors="C1", lw=2.4, label=label)

    ax.set_yscale("log")
    ax.set_ylabel(r"$dN/dy$")
    ax.set_title(rf"$\sqrt{{s_{{NN}}}}={energy:g}$ GeV: data vs THERMUS")
    ax.grid(True, axis="y", which="both", ls="--", alpha=0.35)
    ax.legend()

    if fit_info is not None:
        txt = (
            rf"$T={fit_info['T_MeV']:.2f}\pm{fit_info['T_err_MeV']:.2f}\ \mathrm{{MeV}}$" + "\n"
            rf"$\mu_B={fit_info['muB_MeV']:.2f}\pm{fit_info['muB_err_MeV']:.2f}\ \mathrm{{MeV}}$" + "\n"
            rf"$\mu_S={fit_info['muS_MeV']:.2f}\pm{fit_info['muS_err_MeV']:.2f}\ \mathrm{{MeV}}$" + "\n"
            rf"$\gamma_S={fit_info['gammaS']:.3f}\pm{fit_info['gammaS_err']:.3f}$" + "\n"
            rf"$R={fit_info['R_fm']:.3f}\pm{fit_info['R_err_fm']:.3f}\ \mathrm{{fm}}$" + "\n"
            rf"$dV/dy={fit_info['dVdy_fm3']:.1f}\pm{fit_info['dVdy_err_fm3']:.1f}\ \mathrm{{fm^3}}$" + "\n"
            rf"$\chi^2={fit_info['chi2']:.2f}$" + "\n"
            rf"$\mathrm{{ndf}}={fit_info['ndf']}$" + "\n"
            rf"$\chi^2/\mathrm{{ndf}}={fit_info['chi2_ndf']:.2f}$"
        )
        ax.text(
            0.015,
            0.98,
            txt,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
        )

    axr.axhline(0.0, color="black", lw=1.0)
    axr.axhline(1.0, color="gray", lw=0.8, ls="--")
    axr.axhline(-1.0, color="gray", lw=0.8, ls="--")
    axr.plot(x, residual, "o", ms=4.5, color="C2")
    axr.set_ylabel(r"$(D-M)/\sigma$")
    axr.set_xlabel("Particle")
    axr.grid(True, axis="y", ls="--", alpha=0.35)
    axr.set_xticks(x)
    axr.set_xticklabels(labels, rotation=50, ha="right")

    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot per-energy yield comparison figures for THERMUS.")
    parser.add_argument("--points-dir", default="data/thermus_fit_points")
    parser.add_argument("--outdir", default="data/thermus_figures/yield_comparisons")
    parser.add_argument("--fit-summary", default="data/thermus_fit_results.csv")
    args = parser.parse_args()

    points_dir = Path(args.points_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fit_map = load_fit_summary(Path(args.fit_summary))

    files = sorted(points_dir.glob("sqrts_*GeV_points.csv"))
    if not files:
        raise RuntimeError(f"No points files found in {points_dir}")

    for f in files:
        rows = load_points(f)
        tag = f.stem.replace("_points", "")
        out_png = outdir / f"{tag}_yield_compare.png"
        energy = round(float(rows[0]["energy_GeV"]), 6) if rows else None
        make_plot(rows, out_png, fit_map.get(energy) if energy is not None else None)

    print(f"wrote {len(files)} yield-comparison figures to {outdir}")


if __name__ == "__main__":
    main()
