#!/usr/bin/env python3
"""Compare two-track DCA observable definitions on a reconstructed track tree."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import fast_vertex_mc as fvmc

SCENARIO_NAMES = {
    0: "silicon_only",
    1: "ideal",
    2: "distorted",
}


@dataclass
class RunningMoments:
    nbins: int

    def __post_init__(self) -> None:
        self.count = np.zeros(self.nbins, dtype=np.int64)
        self.sum = np.zeros(self.nbins, dtype=np.float64)
        self.sum_sq = np.zeros(self.nbins, dtype=np.float64)

    def add(self, idx: np.ndarray, values: np.ndarray) -> None:
        valid = (idx >= 0) & (idx < self.nbins) & np.isfinite(values)
        if not np.any(valid):
            return
        np.add.at(self.count, idx[valid], 1)
        np.add.at(self.sum, idx[valid], values[valid])
        np.add.at(self.sum_sq, idx[valid], values[valid] * values[valid])

    def mean(self) -> np.ndarray:
        out = np.full(self.nbins, np.nan, dtype=np.float64)
        mask = self.count > 0
        out[mask] = self.sum[mask] / self.count[mask]
        return out

    def std(self) -> np.ndarray:
        out = np.full(self.nbins, np.nan, dtype=np.float64)
        mask = self.count > 1
        mean = self.mean()
        variance = np.zeros(self.nbins, dtype=np.float64)
        variance[mask] = self.sum_sq[mask] / self.count[mask] - mean[mask] * mean[mask]
        out[mask] = np.sqrt(np.maximum(variance[mask], 0.0))
        return out

    def mean_error(self) -> np.ndarray:
        out = np.full(self.nbins, np.nan, dtype=np.float64)
        mask = self.count > 1
        out[mask] = self.std()[mask] / np.sqrt(self.count[mask])
        return out

    def std_error(self) -> np.ndarray:
        out = np.full(self.nbins, np.nan, dtype=np.float64)
        mask = self.count > 1
        out[mask] = self.std()[mask] / np.sqrt(2.0 * (self.count[mask] - 1.0))
        return out


@dataclass
class MethodProfiles:
    edges: np.ndarray

    def __post_init__(self) -> None:
        nbins = len(self.edges) - 1
        self.current_xy = RunningMoments(nbins)
        self.current_z = RunningMoments(nbins)
        self.alt_xy = RunningMoments(nbins)
        self.alt_z = RunningMoments(nbins)
        self.diff_xy = RunningMoments(nbins)
        self.diff_z = RunningMoments(nbins)

    def add(
        self,
        x: np.ndarray,
        current_xy: np.ndarray,
        current_z: np.ndarray,
        alt_xy: np.ndarray,
        alt_z: np.ndarray,
    ) -> None:
        idx = np.searchsorted(self.edges, x, side="right") - 1
        if len(self.edges) > 1:
            idx = np.where(np.isclose(x, self.edges[-1]), len(self.edges) - 2, idx)
        self.current_xy.add(idx, current_xy)
        self.current_z.add(idx, current_z)
        self.alt_xy.add(idx, alt_xy)
        self.alt_z.add(idx, alt_z)
        self.diff_xy.add(idx, alt_xy - current_xy)
        self.diff_z.add(idx, alt_z - current_z)

    def to_jsonable(self) -> dict[str, Any]:
        centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        return {
            "bin_edges": self.edges.tolist(),
            "bin_centers": centers.tolist(),
            "count": self.current_xy.count.tolist(),
            "current_mean_xy_cm": self.current_xy.mean().tolist(),
            "current_mean_xy_err_cm": self.current_xy.mean_error().tolist(),
            "current_sigma_xy_cm": self.current_xy.std().tolist(),
            "current_sigma_xy_err_cm": self.current_xy.std_error().tolist(),
            "current_mean_z_cm": self.current_z.mean().tolist(),
            "current_mean_z_err_cm": self.current_z.mean_error().tolist(),
            "current_sigma_z_cm": self.current_z.std().tolist(),
            "current_sigma_z_err_cm": self.current_z.std_error().tolist(),
            "alt_mean_xy_cm": self.alt_xy.mean().tolist(),
            "alt_mean_xy_err_cm": self.alt_xy.mean_error().tolist(),
            "alt_sigma_xy_cm": self.alt_xy.std().tolist(),
            "alt_sigma_xy_err_cm": self.alt_xy.std_error().tolist(),
            "alt_mean_z_cm": self.alt_z.mean().tolist(),
            "alt_mean_z_err_cm": self.alt_z.mean_error().tolist(),
            "alt_sigma_z_cm": self.alt_z.std().tolist(),
            "alt_sigma_z_err_cm": self.alt_z.std_error().tolist(),
            "delta_mean_xy_cm": self.diff_xy.mean().tolist(),
            "delta_mean_xy_err_cm": self.diff_xy.mean_error().tolist(),
            "delta_sigma_xy_cm": self.diff_xy.std().tolist(),
            "delta_sigma_xy_err_cm": self.diff_xy.std_error().tolist(),
            "delta_mean_z_cm": self.diff_z.mean().tolist(),
            "delta_mean_z_err_cm": self.diff_z.mean_error().tolist(),
            "delta_sigma_z_cm": self.diff_z.std().tolist(),
            "delta_sigma_z_err_cm": self.diff_z.std_error().tolist(),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True, help="Input reconstructed tree .npz.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for comparison outputs.")
    parser.add_argument("--pair-dp-cut", type=float, default=0.1, help="Maximum |p1-p2| cut in GeV/c.")
    parser.add_argument("--scenario-codes", type=int, nargs="+", default=[1, 2], help="Scenario codes to analyze.")
    parser.add_argument("--max-events", type=int, default=None, help="Optional cap on event entries.")
    return parser.parse_args()


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return float("nan")
    x_std = float(np.std(x))
    y_std = float(np.std(y))
    if x_std <= 0.0 or y_std <= 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def render_profiles(summary: dict[str, Any], outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "mathtext.fontset": "dejavusans",
        }
    )

    payload = summary["pavg_profiles"]
    centers = np.asarray(payload["bin_centers"], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), dpi=180, sharex=True)

    axes[0, 0].errorbar(
        centers,
        payload["current_sigma_xy_cm"],
        yerr=payload["current_sigma_xy_err_cm"],
        color="#0072B2",
        marker="o",
        ms=3.0,
        capsize=2.0,
        label="Current",
    )
    axes[0, 0].errorbar(
        centers,
        payload["alt_sigma_xy_cm"],
        yerr=payload["alt_sigma_xy_err_cm"],
        color="#D55E00",
        marker="s",
        ms=3.0,
        capsize=2.0,
        label="Alternate",
    )
    axes[0, 1].errorbar(
        centers,
        payload["current_sigma_z_cm"],
        yerr=payload["current_sigma_z_err_cm"],
        color="#0072B2",
        marker="o",
        ms=3.0,
        capsize=2.0,
        label="Current",
    )
    axes[0, 1].errorbar(
        centers,
        payload["alt_sigma_z_cm"],
        yerr=payload["alt_sigma_z_err_cm"],
        color="#D55E00",
        marker="s",
        ms=3.0,
        capsize=2.0,
        label="Alternate",
    )
    axes[1, 0].errorbar(
        centers,
        payload["delta_mean_xy_cm"],
        yerr=payload["delta_mean_xy_err_cm"],
        color="#009E73",
        marker="o",
        ms=3.0,
        capsize=2.0,
    )
    axes[1, 1].errorbar(
        centers,
        payload["delta_mean_z_cm"],
        yerr=payload["delta_mean_z_err_cm"],
        color="#009E73",
        marker="o",
        ms=3.0,
        capsize=2.0,
    )

    axes[0, 0].set_ylabel(r"$\sigma$(2-track DCA$_{xy}$) [cm]")
    axes[0, 1].set_ylabel(r"$\sigma$(2-track DCA$_{z}$) [cm]")
    axes[1, 0].set_ylabel(r"Alt - current mean DCA$_{xy}$ [cm]")
    axes[1, 1].set_ylabel(r"Alt - current mean DCA$_{z}$ [cm]")
    axes[1, 0].set_xlabel(r"$p_{\mathrm{avg}}$ [GeV/$c$]")
    axes[1, 1].set_xlabel(r"$p_{\mathrm{avg}}$ [GeV/$c$]")
    axes[0, 0].legend(loc="best")
    axes[0, 1].legend(loc="best")
    axes[0, 0].set_title("Resolution comparison")
    axes[0, 1].set_title("Resolution comparison")
    axes[1, 0].set_title("Mean shift between methods")
    axes[1, 1].set_title("Mean shift between methods")
    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
        ax.axhline(0.0, color="0.7", lw=0.8)
    fig.tight_layout()
    fig.savefig(outdir / "compare_two_track_methods_vs_pavg.png", bbox_inches="tight")
    fig.savefig(outdir / "compare_two_track_methods_vs_pavg.pdf", bbox_inches="tight")
    plt.close(fig)


def analyze_scenario(
    data: Any,
    scenario_code: int,
    pair_dp_cut: float,
    max_events: int | None,
) -> dict[str, Any]:
    event_id = data["event_id"]
    event_scenario = data["event_scenario"]
    event_track_start = data["event_track_start"]
    event_track_stop = data["event_track_stop"]
    track_p = data["track_p"]
    track_eta = data["track_eta"]
    track_phi = data["track_phi"]
    track_pca_x = data["track_pca_x"]
    track_pca_y = data["track_pca_y"]
    track_pca_z = data["track_pca_z"]
    track_dir_x = data["track_dir_x"]
    track_dir_y = data["track_dir_y"]
    track_dir_z = data["track_dir_z"]

    pavg_profiles = MethodProfiles(np.linspace(0.0, 3.0, 61))
    all_current_xy: list[np.ndarray] = []
    all_current_z: list[np.ndarray] = []
    all_alt_xy: list[np.ndarray] = []
    all_alt_z: list[np.ndarray] = []
    all_pavg: list[np.ndarray] = []
    total_events = 0
    used_events = 0
    total_pairs = 0

    event_limit = event_id.shape[0] if max_events is None else min(int(max_events), int(event_id.shape[0]))
    for event_index in range(event_limit):
        total_events += 1
        if int(event_scenario[event_index]) != int(scenario_code):
            continue
        start = int(event_track_start[event_index])
        stop = int(event_track_stop[event_index])
        if stop - start < 2:
            continue

        p = np.asarray(track_p[start:stop], dtype=float)
        pairs = fvmc.selected_pair_indices(p, float(pair_dp_cut))
        if len(pairs) == 0:
            continue

        pair_idx = np.asarray(pairs, dtype=np.int64)
        first = pair_idx[:, 0]
        second = pair_idx[:, 1]
        pca = np.column_stack((track_pca_x[start:stop], track_pca_y[start:stop], track_pca_z[start:stop])).astype(float)
        direction = np.column_stack((track_dir_x[start:stop], track_dir_y[start:stop], track_dir_z[start:stop])).astype(float)
        momentum = direction * p[:, None]

        pca_1 = pca[first]
        pca_2 = pca[second]
        mom_1 = momentum[first]
        mom_2 = momentum[second]
        current_xy, current_z, parallel_current = fvmc.signed_pair_dca_batch(pca_1, mom_1, pca_2, mom_2)
        alt_xy, alt_z, parallel_alt = fvmc.pair_cross_projection_dca_batch(pca_1, mom_1, pca_2, mom_2)
        valid = ~(parallel_current | parallel_alt)
        if not np.any(valid):
            continue

        current_xy = current_xy[valid]
        current_z = current_z[valid]
        alt_xy = alt_xy[valid]
        alt_z = alt_z[valid]
        pavg = 0.5 * (p[first][valid] + p[second][valid])

        pavg_profiles.add(pavg, current_xy, current_z, alt_xy, alt_z)
        all_current_xy.append(current_xy)
        all_current_z.append(current_z)
        all_alt_xy.append(alt_xy)
        all_alt_z.append(alt_z)
        all_pavg.append(pavg)
        total_pairs += int(current_xy.shape[0])
        used_events += 1

    current_xy = np.concatenate(all_current_xy) if all_current_xy else np.empty(0, dtype=float)
    current_z = np.concatenate(all_current_z) if all_current_z else np.empty(0, dtype=float)
    alt_xy = np.concatenate(all_alt_xy) if all_alt_xy else np.empty(0, dtype=float)
    alt_z = np.concatenate(all_alt_z) if all_alt_z else np.empty(0, dtype=float)
    pavg = np.concatenate(all_pavg) if all_pavg else np.empty(0, dtype=float)

    return {
        "scenario_code": int(scenario_code),
        "scenario_name": SCENARIO_NAMES.get(int(scenario_code), f"scenario_{scenario_code}"),
        "pair_dp_cut": float(pair_dp_cut),
        "max_events": max_events,
        "events_seen": total_events,
        "events_with_pairs": used_events,
        "pairs_used": int(total_pairs),
        "overall": {
            "current_xy_mean_cm": float(np.mean(current_xy)) if current_xy.size else float("nan"),
            "current_xy_sigma_cm": float(np.std(current_xy)) if current_xy.size else float("nan"),
            "current_z_mean_cm": float(np.mean(current_z)) if current_z.size else float("nan"),
            "current_z_sigma_cm": float(np.std(current_z)) if current_z.size else float("nan"),
            "alt_xy_mean_cm": float(np.mean(alt_xy)) if alt_xy.size else float("nan"),
            "alt_xy_sigma_cm": float(np.std(alt_xy)) if alt_xy.size else float("nan"),
            "alt_z_mean_cm": float(np.mean(alt_z)) if alt_z.size else float("nan"),
            "alt_z_sigma_cm": float(np.std(alt_z)) if alt_z.size else float("nan"),
            "xy_corr": pearson_corr(current_xy, alt_xy),
            "z_corr": pearson_corr(current_z, alt_z),
            "abs_xy_corr": pearson_corr(np.abs(current_xy), np.abs(alt_xy)),
            "abs_z_corr": pearson_corr(np.abs(current_z), np.abs(alt_z)),
            "sign_same_xy_fraction": float(np.mean(np.sign(current_xy) == np.sign(alt_xy))) if current_xy.size else float("nan"),
            "sign_same_z_fraction": float(np.mean(np.sign(current_z) == np.sign(alt_z))) if current_z.size else float("nan"),
            "mean_abs_diff_xy_cm": float(np.mean(np.abs(alt_xy - current_xy))) if current_xy.size else float("nan"),
            "mean_abs_diff_z_cm": float(np.mean(np.abs(alt_z - current_z))) if current_z.size else float("nan"),
            "pavg_mean_gev": float(np.mean(pavg)) if pavg.size else float("nan"),
        },
        "pavg_profiles": pavg_profiles.to_jsonable(),
    }


def main() -> None:
    args = parse_args()
    tree_path = args.tree.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with np.load(tree_path) as data:
        available = {int(code) for code in np.unique(data["event_scenario"])}
        scenario_codes = [int(code) for code in args.scenario_codes if int(code) in available]
        manifest: dict[str, Any] = {
            "tree": str(tree_path),
            "pair_dp_cut": float(args.pair_dp_cut),
            "max_events": args.max_events,
            "scenarios": {},
        }
        for scenario_code in scenario_codes:
            scenario_summary = analyze_scenario(data, scenario_code, args.pair_dp_cut, args.max_events)
            scenario_outdir = outdir / scenario_summary["scenario_name"]
            scenario_outdir.mkdir(parents=True, exist_ok=True)
            (scenario_outdir / "two_track_method_comparison.json").write_text(
                json.dumps(scenario_summary, indent=2) + "\n",
                encoding="utf-8",
            )
            render_profiles(scenario_summary, scenario_outdir)
            manifest["scenarios"][scenario_summary["scenario_name"]] = str(
                (scenario_outdir / "two_track_method_comparison.json").resolve()
            )

    (outdir / "compare_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
