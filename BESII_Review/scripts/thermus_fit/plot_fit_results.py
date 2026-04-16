#!/usr/bin/env python3
"""Plot THERMUS fit parameters vs energy and T vs muB."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_rows(path: Path):
    rows = []
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    rows.sort(key=lambda r: float(r["energy_GeV"]))
    return rows


def fval(row: dict, key: str) -> float:
    return float(row[key])


def plot_vs_energy(rows, y_key, y_err_key, ylabel, out_png: Path):
    x = [fval(r, "energy_GeV") for r in rows]
    y = [fval(r, y_key) for r in rows]
    yerr = [fval(r, y_err_key) for r in rows]

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.errorbar(x, y, yerr=yerr, fmt="o-", ms=4, lw=1.3, capsize=2)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def plot_t_vs_mub(rows, out_png: Path):
    x = [fval(r, "muB_MeV") for r in rows]
    y = [fval(r, "T_MeV") for r in rows]
    xerr = [fval(r, "muB_err_MeV") for r in rows]
    yerr = [fval(r, "T_err_MeV") for r in rows]
    en = [fval(r, "energy_GeV") for r in rows]

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o", ms=4, capsize=2, lw=1.0)

    for xi, yi, ei in zip(x, y, en):
        ax.annotate(f"{ei:g}", (xi, yi), textcoords="offset points", xytext=(5, 4), fontsize=8)

    ax.set_xlabel(r"$\mu_B$ (MeV)")
    ax.set_ylabel(r"$T$ (MeV)")
    ax.grid(True, ls="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot THERMUS fit results.")
    parser.add_argument("--input", default="data/thermus_fit_results.csv")
    parser.add_argument("--outdir", default="data/thermus_figures")
    args = parser.parse_args()

    in_csv = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(in_csv)
    if not rows:
        raise RuntimeError(f"No data rows found in {in_csv}")

    plot_vs_energy(rows, "T_MeV", "T_err_MeV", r"$T$ (MeV)", outdir / "T_vs_energy.png")
    plot_vs_energy(rows, "muB_MeV", "muB_err_MeV", r"$\mu_B$ (MeV)", outdir / "muB_vs_energy.png")
    plot_vs_energy(rows, "muS_MeV", "muS_err_MeV", r"$\mu_S$ (MeV)", outdir / "muS_vs_energy.png")
    plot_vs_energy(rows, "gammaS", "gammaS_err", r"$\gamma_S$", outdir / "gammaS_vs_energy.png")
    plot_vs_energy(rows, "dVdy_fm3", "dVdy_err_fm3", r"$dV/dy$ (fm$^3$)", outdir / "dVdy_vs_energy.png")
    plot_t_vs_mub(rows, outdir / "T_vs_muB.png")

    print(f"wrote figures to {outdir}")


if __name__ == "__main__":
    main()
