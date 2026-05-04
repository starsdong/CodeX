#!/usr/bin/env python3
"""Regenerate the 3 GeV phi/K- comparison plot from packaged TSV tables."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"


def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    urqmd = read_tsv(DATA / "phi_k_centrality3_besii_matched_vs_npart_3gev_10k.tsv")
    besii = read_tsv(DATA / "besii_measured_phi_k_ratio_vs_npart_3gev.tsv")
    thermus = read_tsv(DATA / "thermus_phi_k_ratio_rc_scan_3gev.tsv")

    ux = np.array([float(r["npart_mean"]) for r in urqmd])
    uxerr = np.array([float(r["npart_rms"]) for r in urqmd])
    uy = np.array([float(r["ratio"]) for r in urqmd])
    uyerr = np.array([float(r["stat_err"]) for r in urqmd])

    mx = np.array([float(r["npart_mean"]) for r in besii])
    mxerr = np.array([float(r["npart_rms"]) for r in besii])
    my = np.array([float(r["ratio"]) for r in besii])
    myerr = np.array([float(r["total_err"]) for r in besii])

    plt.rcParams.update({
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "legend.fontsize": 9.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 160,
    })

    fig, ax = plt.subplots(figsize=(8.1, 5.35))
    ax.errorbar(
        ux, uy, yerr=uyerr, xerr=uxerr,
        fmt="o", color="black", ecolor="0.40", elinewidth=1.15,
        capsize=4, markersize=6.5, label="UrQMD 4.0, 10k events", zorder=5,
    )
    ax.errorbar(
        mx, my, yerr=myerr, xerr=mxerr,
        fmt="s", color="#b2182b", markerfacecolor="white", markeredgewidth=1.8,
        ecolor="#b2182b", elinewidth=1.15, capsize=4, markersize=7.0,
        label="STAR BESII data, mid-y dN/dy ratio", zorder=6,
    )

    label_offsets = {"0-10%": (-50, 0.016), "10-40%": (-53, 0.016), "40-60%": (9, -0.020)}
    for row in besii:
        dx, dy = label_offsets[row["centrality"]]
        ax.text(
            float(row["npart_mean"]) + dx,
            float(row["ratio"]) + dy,
            row["centrality"],
            fontsize=10.2,
            color="#7f1d1d",
        )

    colors = {2: "#9b2226", 3: "#2a9d8f", 4: "#005f73", 6: "#7a5195"}
    linestyles = {2: "-", 3: (0, (5, 2)), 4: "--", 6: "-."}
    for row in thermus:
        rc = int(float(row["Rc_fm"]))
        ratio = float(row["ratio"])
        ax.axhline(
            ratio,
            color=colors[rc],
            linestyle=linestyles[rc],
            linewidth=1.8,
            label=fr"THERMUS SCE $R_c={rc}$ fm: {ratio:.3f}",
            zorder=2,
        )

    ax.set_title(r"Au+Au $\sqrt{s_{NN}}=3$ GeV")
    ax.set_xlabel(r"$\langle N_{\mathrm{part}}\rangle$")
    ax.set_ylabel(r"$\phi/K^-$")
    ax.set_xlim(0, 380)
    ymax = max(np.max(uy + uyerr), np.max(my + myerr), max(float(r["ratio"]) for r in thermus)) * 1.12
    ax.set_ylim(0, max(0.44, ymax))
    ax.grid(True, which="major", axis="both", alpha=0.22, linewidth=0.8)
    ax.legend(loc="upper right", frameon=True, framealpha=0.96, borderpad=0.6)
    ax.text(
        0.02,
        0.025,
        r"THERMUS: BESII Review 3 GeV SCE fit, $T=81.66$ MeV, $\mu_B=751.42$ MeV, $R=8.08$ fm",
        transform=ax.transAxes,
        fontsize=9.4,
        color="0.25",
        va="bottom",
    )
    fig.tight_layout()

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "phi_k_centrality3_vs_npart_3gev_10k.png", dpi=240)
    fig.savefig(FIGURES / "phi_k_centrality3_vs_npart_3gev_10k.pdf")


if __name__ == "__main__":
    main()

