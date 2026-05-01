#!/usr/bin/env python3
"""Plot compiled heavy-ion flow excitation functions."""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / ".matplotlib-cache"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def fnum(value: str) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def total_error(row: dict[str, str]) -> tuple[float | None, float | None]:
    stat_plus = fnum(row.get("stat_err_plus", ""))
    stat_minus = fnum(row.get("stat_err_minus", ""))
    syst_plus = fnum(row.get("syst_err_plus", ""))
    syst_minus = fnum(row.get("syst_err_minus", ""))

    plus_terms = [x for x in (stat_plus, syst_plus) if x is not None]
    minus_terms = [x for x in (stat_minus, syst_minus) if x is not None]
    plus = math.sqrt(sum(x * x for x in plus_terms)) if plus_terms else None
    minus = math.sqrt(sum(x * x for x in minus_terms)) if minus_terms else None
    return minus, plus


def yerr_for(rows: list[dict[str, str]]) -> list[list[float]] | None:
    lower: list[float] = []
    upper: list[float] = []
    has_any = False
    for row in rows:
        err_minus, err_plus = total_error(row)
        if err_minus is None and err_plus is None:
            lower.append(0.0)
            upper.append(0.0)
            continue
        has_any = True
        lower.append(err_minus or err_plus or 0.0)
        upper.append(err_plus or err_minus or 0.0)
    return [lower, upper] if has_any else None


def setup_energy_axis(ax: plt.Axes, include_sis: bool = False) -> None:
    ax.set_xscale("log")
    if include_sis:
        ticks = [2, 2.4, 3, 4.5, 7.7, 10, 20, 39, 62.4, 100, 200]
        ax.set_xlim(1.75, 260.0)
    else:
        ticks = [5, 7.7, 10, 20, 39, 62.4, 100, 200]
        ax.set_xlim(4.0, 260.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick:g}" for tick in ticks])
    ax.set_xlabel(r"$\sqrt{s_{NN}}$ (GeV)")
    ax.grid(True, which="major", color="#d9d9d9", linewidth=0.8)
    ax.grid(True, which="minor", color="#eeeeee", linewidth=0.5)
    ax.tick_params(which="both", direction="in", top=True, right=True)


def save(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_dv1dy() -> None:
    rows = read_rows(DATA / "flow_dv1dy_vs_energy.csv")

    marker_by_particle = {
        "proton": "o",
        "net-proton": "D",
        "antiproton": "s",
        "pi+": "^",
        "pi-": "v",
    }
    color_by_particle = {
        "proton": "#1f77b4",
        "net-proton": "#111111",
        "antiproton": "#d62728",
        "pi+": "#2ca02c",
        "pi-": "#9467bd",
    }

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    by_particle: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_particle[row["particle"]].append(row)

    order = ["proton", "net-proton", "antiproton", "pi+", "pi-"]
    for particle in order:
        particle_rows = sorted(
            by_particle[particle],
            key=lambda item: (item["experiment"] != "STAR", fnum(item["sqrt_s_NN_GeV"]) or 0),
        )
        if not particle_rows:
            continue

        star_rows = [row for row in particle_rows if row["experiment"] == "STAR"]
        comparison_rows = [row for row in particle_rows if row["experiment"] != "STAR"]

        for subset, label, alpha, fillstyle in (
            (star_rows, particle, 1.0, "full"),
            (comparison_rows, f"{particle} older data", 0.78, "none"),
        ):
            if not subset:
                continue
            x = [fnum(row["sqrt_s_NN_GeV"]) for row in subset]
            y = [fnum(row["value"]) for row in subset]
            ax.errorbar(
                x,
                y,
                yerr=yerr_for(subset),
                marker=marker_by_particle.get(particle, "o"),
                color=color_by_particle.get(particle, "#333333"),
                markerfacecolor=(
                    color_by_particle.get(particle, "#333333")
                    if fillstyle == "full"
                    else "white"
                ),
                markeredgecolor=color_by_particle.get(particle, "#333333"),
                linestyle="-",
                linewidth=1.2 if star_rows and subset is star_rows else 0.0,
                elinewidth=1.0,
                capsize=2.5,
                markersize=6.0,
                alpha=alpha,
                label=label,
            )

    ax.axhline(0, color="#4d4d4d", linewidth=0.9)
    setup_energy_axis(ax)
    ax.set_ylabel(r"$dv_1/dy|_{y=0}$")
    ax.set_ylim(-0.09, 0.10)
    ax.set_title(r"Directed-flow slope excitation function", pad=10)

    handles, labels = ax.get_legend_handles_labels()
    kept = []
    seen = set()
    for handle, label in zip(handles, labels):
        if label.endswith("older data"):
            continue
        if label not in seen:
            kept.append((handle, label))
            seen.add(label)
    older_handle = Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor="white",
        markeredgecolor="#666666",
        color="#666666",
        label="E895/NA49 comparison",
    )
    ax.legend(
        [h for h, _ in kept] + [older_handle],
        [label for _, label in kept] + ["E895/NA49 comparison"],
        frameon=False,
        fontsize=9,
        ncol=2,
        loc="upper right",
    )
    ax.text(
        0.03,
        0.04,
        "10-40% centrality; total uncertainties shown where available",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#555555",
    )
    fig.tight_layout()
    save(fig, "dv1dy_vs_energy")


def plot_v2_cross_experiment() -> None:
    rows = read_rows(DATA / "flow_v2_cross_experiment_vs_energy.csv")
    star_rows = read_rows(DATA / "flow_v2_star_bes_eta0_vs_energy.csv")

    color_by_experiment = {
        "E877": "#8c564b",
        "NA49": "#1f77b4",
        "CERES": "#ff7f0e",
        "STAR": "#111111",
        "PHENIX": "#d62728",
        "PHOBOS": "#2ca02c",
    }
    marker_by_experiment = {
        "E877": "P",
        "NA49": "o",
        "CERES": "s",
        "STAR": "D",
        "PHENIX": "^",
        "PHOBOS": "v",
    }

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    by_experiment: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_experiment[row["experiment"]].append(row)

    for experiment in ["E877", "NA49", "CERES", "STAR", "PHENIX", "PHOBOS"]:
        subset = sorted(
            by_experiment.get(experiment, []),
            key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0,
        )
        if not subset:
            continue
        ax.errorbar(
            [fnum(row["sqrt_s_NN_GeV"]) for row in subset],
            [fnum(row["value"]) for row in subset],
            yerr=yerr_for(subset),
            marker=marker_by_experiment[experiment],
            color=color_by_experiment[experiment],
            markerfacecolor=color_by_experiment[experiment],
            linestyle="none",
            elinewidth=1.0,
            capsize=2.5,
            markersize=6.0,
            label=experiment,
        )

    star_rows_sorted = sorted(star_rows, key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0)
    ax.errorbar(
        [fnum(row["sqrt_s_NN_GeV"]) for row in star_rows_sorted],
        [fnum(row["value"]) for row in star_rows_sorted],
        yerr=yerr_for(star_rows_sorted),
        marker="o",
        color="#111111",
        markerfacecolor="white",
        markeredgecolor="#111111",
        linestyle="-",
        linewidth=1.2,
        elinewidth=1.0,
        capsize=2.5,
        markersize=5.2,
        label=r"STAR BES $\eta\approx0$",
    )

    setup_energy_axis(ax)
    ax.set_ylabel(r"$v_2$")
    ax.set_ylim(0.0, 0.07)
    ax.set_title(r"Elliptic-flow excitation function", pad=10)
    ax.legend(frameon=False, fontsize=9, ncol=2, loc="upper left")
    ax.text(
        0.03,
        0.04,
        "Mid-central, near midrapidity; cuts differ among experiments",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#555555",
    )
    fig.tight_layout()
    save(fig, "v2_vs_energy_cross_experiment")


def plot_v2_star_only() -> None:
    rows = sorted(
        read_rows(DATA / "flow_v2_star_bes_eta0_vs_energy.csv"),
        key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0,
    )

    fig, ax = plt.subplots(figsize=(6.8, 4.7))
    ax.errorbar(
        [fnum(row["sqrt_s_NN_GeV"]) for row in rows],
        [fnum(row["value"]) for row in rows],
        yerr=yerr_for(rows),
        marker="o",
        color="#111111",
        markerfacecolor="#ffffff",
        markeredgecolor="#111111",
        linestyle="-",
        linewidth=1.3,
        elinewidth=1.0,
        capsize=2.5,
        markersize=6.0,
    )
    setup_energy_axis(ax)
    ax.set_ylabel(r"$v_2\{EP\}$ at $\eta\approx0$")
    ax.set_ylim(0.035, 0.064)
    ax.set_title(r"STAR BES charged-hadron elliptic flow", pad=10)
    ax.text(
        0.03,
        0.05,
        r"Au+Au 10-40%, $0.2<p_T<2.0$ GeV/$c$",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#555555",
    )
    fig.tight_layout()
    save(fig, "v2_vs_energy_star_bes_eta0")


def plot_dv1dy_low_energy_expanded() -> None:
    rows = read_rows(DATA / "flow_dv1dy_vs_energy.csv")
    rows.extend(read_rows(DATA / "flow_dv1dy_low_energy_expanded.csv"))

    styles = {
        "proton": ("#1f77b4", "o"),
        "net-proton": ("#111111", "D"),
        "antiproton": ("#d62728", "s"),
        "pi+": ("#2ca02c", "^"),
        "pi-": ("#9467bd", "v"),
        "Lambda": ("#ff7f0e", "P"),
        "K0S": ("#17becf", "X"),
    }

    fig, ax = plt.subplots(figsize=(7.6, 5.1))
    by_particle: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["particle"] in styles:
            by_particle[row["particle"]].append(row)

    for particle in ["proton", "net-proton", "antiproton", "pi+", "pi-", "Lambda", "K0S"]:
        subset = sorted(by_particle.get(particle, []), key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0)
        if not subset:
            continue
        color, marker = styles[particle]
        star_like = [row for row in subset if row["experiment"] in {"STAR", "STAR FXT"}]
        other = [row for row in subset if row["experiment"] not in {"STAR", "STAR FXT"}]
        for rows_part, face, label in (
            (star_like, color, particle),
            (other, "white", f"{particle} comparison"),
        ):
            if not rows_part:
                continue
            ax.errorbar(
                [fnum(row["sqrt_s_NN_GeV"]) for row in rows_part],
                [fnum(row["value"]) for row in rows_part],
                yerr=yerr_for(rows_part),
                marker=marker,
                color=color,
                markerfacecolor=face,
                markeredgecolor=color,
                linestyle="-" if rows_part is star_like and len(rows_part) > 1 else "none",
                linewidth=1.1,
                elinewidth=1.0,
                capsize=2.4,
                markersize=5.8,
                label=label,
            )

    ax.axhline(0, color="#4d4d4d", linewidth=0.9)
    setup_energy_axis(ax, include_sis=True)
    ax.set_ylabel(r"$dv_1/dy|_{y=0}$")
    ax.set_ylim(-0.12, 0.43)
    ax.set_title(r"Directed-flow slope with low-energy additions", pad=10)

    handles, labels = ax.get_legend_handles_labels()
    unique = []
    seen = set()
    for handle, label in zip(handles, labels):
        if label in seen:
            continue
        seen.add(label)
        unique.append((handle, label))
    ax.legend([h for h, _ in unique], [label for _, label in unique], frameon=False, fontsize=8.3, ncol=2, loc="upper right")
    ax.text(
        0.03,
        0.04,
        "Low-energy points use HADES/STAR FXT compilations; cuts are not identical across experiments",
        transform=ax.transAxes,
        fontsize=8.2,
        color="#555555",
    )
    fig.tight_layout()
    save(fig, "dv1dy_vs_energy_low_energy_expanded")


def plot_v2_low_energy_expanded() -> None:
    rows = read_rows(DATA / "flow_v2_low_energy_expanded.csv")
    styles = {
        "FOPI": ("#1f77b4", "o"),
        "EOS/E895/E877": ("#8c564b", "P"),
        "STAR FXT": ("#111111", "D"),
        "CERES": ("#ff7f0e", "s"),
        "NA49": ("#17becf", "X"),
        "STAR": ("#2ca02c", "^"),
        "PHENIX": ("#d62728", "v"),
        "PHOBOS": ("#9467bd", "<"),
    }

    fig, ax = plt.subplots(figsize=(7.6, 5.1))
    by_experiment: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_experiment[row["experiment"]].append(row)

    for experiment in ["FOPI", "EOS/E895/E877", "STAR FXT", "CERES", "NA49", "STAR", "PHENIX", "PHOBOS"]:
        subset = sorted(by_experiment.get(experiment, []), key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0)
        if not subset:
            continue
        color, marker = styles[experiment]
        ax.errorbar(
            [fnum(row["sqrt_s_NN_GeV"]) for row in subset],
            [fnum(row["value"]) for row in subset],
            yerr=yerr_for(subset),
            marker=marker,
            color=color,
            markerfacecolor=color if experiment not in {"STAR FXT"} else "white",
            markeredgecolor=color,
            linestyle="-" if experiment in {"FOPI", "EOS/E895/E877", "STAR"} else "none",
            linewidth=1.1,
            elinewidth=1.0,
            capsize=2.4,
            markersize=5.6,
            label=experiment,
        )

    ax.axhline(0, color="#4d4d4d", linewidth=0.9)
    setup_energy_axis(ax, include_sis=True)
    ax.set_ylabel(r"$v_2$")
    ax.set_ylim(-0.095, 0.085)
    ax.set_title(r"Expanded elliptic-flow excitation function", pad=10)
    ax.legend(frameon=False, fontsize=8.3, ncol=2, loc="lower right")
    ax.text(
        0.03,
        0.95,
        "Compiled from STAR FXT Fig. 14 HEPData; centrality/species differ across data sets",
        transform=ax.transAxes,
        fontsize=8.2,
        color="#555555",
        va="top",
    )
    fig.tight_layout()
    save(fig, "v2_vs_energy_low_energy_expanded")


def plot_combined_summary() -> None:
    dv1_rows = read_rows(DATA / "flow_dv1dy_vs_energy.csv")
    v2_rows = read_rows(DATA / "flow_v2_star_bes_eta0_vs_energy.csv")
    selected_particles = ["proton", "net-proton", "pi+", "pi-"]
    color_by_particle = {
        "proton": "#1f77b4",
        "net-proton": "#111111",
        "pi+": "#2ca02c",
        "pi-": "#9467bd",
    }
    marker_by_particle = {
        "proton": "o",
        "net-proton": "D",
        "pi+": "^",
        "pi-": "v",
    }

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), constrained_layout=True)
    ax = axes[0]
    for particle in selected_particles:
        subset = sorted(
            [
                row
                for row in dv1_rows
                if row["particle"] == particle and row["experiment"] == "STAR"
            ],
            key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0,
        )
        ax.errorbar(
            [fnum(row["sqrt_s_NN_GeV"]) for row in subset],
            [fnum(row["value"]) for row in subset],
            yerr=yerr_for(subset),
            marker=marker_by_particle[particle],
            color=color_by_particle[particle],
            linestyle="-",
            linewidth=1.2,
            elinewidth=1.0,
            capsize=2.2,
            markersize=5.4,
            label=particle,
        )
    ax.axhline(0, color="#4d4d4d", linewidth=0.9)
    setup_energy_axis(ax)
    ax.set_ylabel(r"$dv_1/dy|_{y=0}$")
    ax.set_ylim(-0.09, 0.08)
    ax.set_title(r"STAR $dv_1/dy$")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")

    ax = axes[1]
    v2_rows = sorted(v2_rows, key=lambda item: fnum(item["sqrt_s_NN_GeV"]) or 0)
    ax.errorbar(
        [fnum(row["sqrt_s_NN_GeV"]) for row in v2_rows],
        [fnum(row["value"]) for row in v2_rows],
        yerr=yerr_for(v2_rows),
        marker="o",
        color="#111111",
        markerfacecolor="white",
        markeredgecolor="#111111",
        linestyle="-",
        linewidth=1.3,
        elinewidth=1.0,
        capsize=2.2,
        markersize=5.4,
    )
    setup_energy_axis(ax)
    ax.set_ylabel(r"$v_2\{EP\}$ at $\eta\approx0$")
    ax.set_ylim(0.035, 0.064)
    ax.set_title(r"STAR charged-hadron $v_2$")
    fig.suptitle("Flow excitation functions", y=1.03, fontsize=14)
    save(fig, "flow_excitation_summary")


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    plot_dv1dy()
    plot_v2_cross_experiment()
    plot_v2_star_only()
    plot_dv1dy_low_energy_expanded()
    plot_v2_low_energy_expanded()
    plot_combined_summary()
    print(f"Wrote figures to {FIGURES}")


if __name__ == "__main__":
    main()
