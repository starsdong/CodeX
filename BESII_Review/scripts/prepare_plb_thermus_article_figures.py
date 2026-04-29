#!/usr/bin/env python3
import csv
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUT_DIR = Path("manuscripts/plb_thermus_star/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

YIELD_PANEL_INPUTS = [
    (
        3.0,
        "SCE",
        Path("data/thermus_sce_3gev/fit_points/sqrts_3GeV_points.csv"),
    ),
    (
        7.7,
        "GCE",
        Path("data/thermus_fit_points/sqrts_7p7GeV_points.csv"),
    ),
    (
        200.0,
        "GCE",
        Path("data/thermus_fit_points/sqrts_200p0GeV_points.csv"),
    ),
]

PARTICLE_ORDER = [
    "pi+",
    "pi-",
    "K+",
    "K-",
    "Ks0",
    "p",
    "pbar",
    "Lambda",
    "Lambda_bar",
    "Xi",
    "Xi_bar",
    "phi",
]

PARTICLE_LABELS = {
    "pi+": r"$\pi^+$",
    "pi-": r"$\pi^-$",
    "K+": r"$K^+$",
    "K-": r"$K^-$",
    "Ks0": r"$K^0_S$",
    "p": r"$p$",
    "pbar": r"$\bar{p}$",
    "Lambda": r"$\Lambda$",
    "Lambda_bar": r"$\bar{\Lambda}$",
    "Xi": r"$\Xi^-$",
    "Xi_bar": r"$\bar{\Xi}^+$",
    "phi": r"$\phi$",
}


def read_rows(path):
    with Path(path).open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row, key):
    text = row.get(key, "")
    return float(text) if text else 0.0


def maybe_f(row, key):
    text = row.get(key, "")
    return float(text) if text else None


def sort_yield_rows(rows):
    order = {name: i for i, name in enumerate(PARTICLE_ORDER)}
    return sorted(rows, key=lambda r: order.get(r["particle_name"], len(order)))


def fit_info_by_energy(rows):
    return {round(f(row, "energy_GeV"), 6): row for row in rows}


def draw_parameter_panel(
    ax,
    gce,
    sce_row,
    y_key,
    y_err_key,
    ylabel,
    *,
    gce_default=None,
    gce_err_default=None,
    sce_default=None,
    sce_err_default=None,
    ylim=None,
    yscale=None,
    annotate=None,
    legend=False,
):
    def collect(rows, default, err_default):
        x, y, yerr = [], [], []
        for row in rows:
            value = maybe_f(row, y_key)
            if value is None:
                value = default
            if value is None:
                continue
            err = maybe_f(row, y_err_key) if y_err_key else None
            if err is None:
                err = err_default
            x.append(f(row, "energy_GeV"))
            y.append(value)
            yerr.append(err if err is not None else 0.0)
        return x, y, yerr

    e_gce, y_gce, yerr_gce = collect(gce, gce_default, gce_err_default)
    e_sce, y_sce, yerr_sce = collect([sce_row], sce_default, sce_err_default)

    if e_gce:
        ax.errorbar(
            e_gce,
            y_gce,
            yerr=yerr_gce,
            fmt="o",
            color="#1f77b4",
            capsize=2,
            ms=4.6,
            lw=1.0,
            label="GCE, 7.7-200 GeV" if legend else None,
        )
    if e_sce:
        ax.errorbar(
            e_sce,
            y_sce,
            yerr=yerr_sce,
            fmt="s",
            color="#d62728",
            capsize=2,
            ms=4.8,
            lw=1.0,
            label="SCE, 3 GeV" if legend else None,
        )
    ax.set_xscale("log")
    if yscale:
        ax.set_yscale(yscale)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls="--", alpha=0.35)
    if annotate:
        ax.text(
            0.04,
            0.92,
            annotate,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7.4,
            color="0.35",
        )
    if legend:
        ax.legend(loc="best", fontsize=7.6, frameon=False)


def plot_freezeout_energy(gce, sce_row):
    fig, axes = plt.subplots(2, 4, figsize=(12.8, 6.1), sharex=True)
    panels = [
        {
            "y_key": "T_MeV",
            "y_err_key": "T_err_MeV",
            "ylabel": r"$T_{\rm ch}$ (MeV)",
            "ylim": (60, 178),
            "legend": True,
        },
        {
            "y_key": "muB_MeV",
            "y_err_key": "muB_err_MeV",
            "ylabel": r"$\mu_B$ (MeV)",
            "ylim": (0, 810),
        },
        {
            "y_key": "muS_MeV",
            "y_err_key": "muS_err_MeV",
            "ylabel": r"$\mu_S$ (MeV)",
            "ylim": (0, 112),
            "annotate": "GCE only",
        },
        {
            "y_key": "muQ_MeV",
            "y_err_key": "muQ_err_MeV",
            "ylabel": r"$\mu_Q$ (MeV)",
            "gce_default": 0.0,
            "gce_err_default": 0.0,
            "ylim": (-30, 8),
            "annotate": r"GCE fixed at $0$",
        },
        {
            "y_key": "gammaS",
            "y_err_key": "gammaS_err",
            "ylabel": r"$\gamma_S$",
            "ylim": (0.94, 1.06),
            "annotate": "fixed",
        },
        {
            "y_key": "R_fm",
            "y_err_key": "R_err_fm",
            "ylabel": r"$R$ (fm)",
            "ylim": (4.3, 9.6),
        },
        {
            "y_key": "Rc_fm",
            "y_err_key": "Rc_err_fm",
            "ylabel": r"$R_C$ (fm)",
            "ylim": (0, 5.3),
            "annotate": "SCE only",
        },
        {
            "y_key": "dVdy_fm3",
            "y_err_key": "dVdy_err_fm3",
            "ylabel": r"$dV/dy$ (fm$^3$)",
            "yscale": "log",
            "ylim": (4.0e2, 4.0e3),
        },
    ]
    for ax, spec in zip(axes.flat, panels):
        draw_parameter_panel(ax, gce, sce_row, **spec)

    fig.tight_layout(w_pad=1.0, h_pad=1.0)
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"freezeout_summary.{ext}"
        fig.savefig(out, dpi=260)
        print(f"wrote {out}")
    plt.close(fig)


def plot_t_vs_mub(gce, sce_row):
    mub_gce = [f(r, "muB_MeV") for r in gce]
    muberr_gce = [f(r, "muB_err_MeV") for r in gce]
    t_gce = [f(r, "T_MeV") for r in gce]
    terr_gce = [f(r, "T_err_MeV") for r in gce]
    e_gce = [f(r, "energy_GeV") for r in gce]
    chi_gce = [f(r, "chi2_ndf") for r in gce]

    mub_sce = f(sce_row, "muB_MeV")
    muberr_sce = f(sce_row, "muB_err_MeV")
    t_sce = f(sce_row, "T_MeV")
    terr_sce = f(sce_row, "T_err_MeV")
    chi_sce = f(sce_row, "chi2_ndf")

    fig, ax = plt.subplots(figsize=(5.4, 4.1))
    ax.errorbar(mub_gce, t_gce, xerr=muberr_gce, yerr=terr_gce, fmt="o", color="#1f77b4", capsize=2, label="GCE")
    ax.errorbar([mub_sce], [t_sce], xerr=[muberr_sce], yerr=[terr_sce], fmt="s", color="#d62728", capsize=2, label="SCE")

    offsets = {
        200.0: (8, 8),
        62.4: (8, -12),
        39.0: (8, 5),
        27.0: (8, 5),
        19.6: (8, 5),
        11.5: (8, 5),
        7.7: (8, 5),
    }
    for x, y, energy, chi in zip(mub_gce, t_gce, e_gce, chi_gce):
        ax.annotate(
            f"{energy:g}",
            (x, y),
            xytext=offsets.get(energy, (5, 4)),
            textcoords="offset points",
            fontsize=8,
        )
    ax.annotate(
        f"3\n$\\chi^2/\\mathrm{{ndf}}={chi_sce:.1f}$",
        (mub_sce, t_sce),
        xytext=(-8, 6),
        textcoords="offset points",
        ha="right",
        fontsize=8,
    )
    ax.set_xlabel(r"$\mu_B$ (MeV)")
    ax.set_ylabel(r"$T_{\rm ch}$ (MeV)")
    ax.set_xlim(0, 820)
    ax.set_ylim(60, 178)
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    fig.tight_layout()

    for ext in ("pdf", "png"):
        out = OUT_DIR / f"freezeout_t_vs_mub.{ext}"
        fig.savefig(out, dpi=260)
        print(f"wrote {out}")
    plt.close(fig)


def plot_selected_yield_comparisons(gce_by_energy, sce_row):
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(12.6, 5.8),
        sharey="row",
        gridspec_kw={"height_ratios": [3.0, 1.15], "hspace": 0.05},
    )
    fit_lookup = dict(gce_by_energy)
    fit_lookup[3.0] = sce_row

    for col, (energy, ensemble, path) in enumerate(YIELD_PANEL_INPUTS):
        ax = axes[0, col]
        axr = axes[1, col]
        rows = sort_yield_rows(read_rows(path))
        labels = [PARTICLE_LABELS.get(r["particle_name"], r["particle_name"]) for r in rows]
        data = [f(r, "data_yield") for r in rows]
        err = [f(r, "data_err") for r in rows]
        model = [f(r, "model_yield") for r in rows]
        residual = [((d - m) / e) if e > 0.0 else 0.0 for d, m, e in zip(data, model, err)]
        x = list(range(len(rows)))

        ax.errorbar(x, data, yerr=err, fmt="o", ms=4.2, capsize=2, lw=1.0, color="#1f77b4", label="Data")
        for i, yval in enumerate(model):
            ax.hlines(yval, i - 0.33, i + 0.33, colors="#ff7f0e", lw=2.2, label="THERMUS" if i == 0 else None)

        info = fit_lookup[round(energy, 6)]
        ax.text(
            0.03,
            0.95,
            (
                rf"{ensemble}, $\sqrt{{s_{{NN}}}}={energy:g}$ GeV"
                "\n"
                rf"$T={f(info, 'T_MeV'):.1f}$ MeV, $\mu_B={f(info, 'muB_MeV'):.0f}$ MeV"
                "\n"
                rf"$\chi^2/\mathrm{{ndf}}={f(info, 'chi2_ndf'):.2f}$"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.75", "pad": 2},
        )
        ax.set_yscale("log")
        ax.set_ylim(3.0e-3, 5.0e2)
        ax.set_xticks(x)
        ax.set_xticklabels([])
        ax.grid(True, axis="y", which="both", ls="--", alpha=0.32)

        axr.axhline(0.0, color="black", lw=1.0)
        axr.axhline(1.0, color="0.55", lw=0.8, ls="--")
        axr.axhline(-1.0, color="0.55", lw=0.8, ls="--")
        axr.plot(x, residual, "o", ms=4.0, color="#2ca02c")
        axr.set_ylim(-3.2, 3.2)
        axr.set_xticks(x)
        axr.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
        axr.grid(True, axis="y", ls="--", alpha=0.32)
        axr.set_xlabel("Particle")

    axes[0, 0].set_ylabel(r"$dN/dy$")
    axes[1, 0].set_ylabel(r"$(D-M)/\sigma$")
    axes[0, -1].legend(loc="upper right", fontsize=8, frameon=False)
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.22, top=0.98, wspace=0.08, hspace=0.06)
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"yield_comparison_selected_energies.{ext}"
        fig.savefig(out, dpi=260)
        print(f"wrote {out}")
    plt.close(fig)


def main():
    gce = read_rows("data/thermus_fit_results.csv")
    sce = read_rows("data/thermus_sce_3gev/fit_summary.csv")
    sce_row = sce[0]
    gce_by_energy = fit_info_by_energy(gce)

    plot_freezeout_energy(gce, sce_row)
    plot_t_vs_mub(gce, sce_row)
    plot_selected_yield_comparisons(gce_by_energy, sce_row)


if __name__ == "__main__":
    main()
