#!/usr/bin/env python3
"""Reconstruct tracks from the baseline hit sample and export tree-like arrays."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import analyze_pointing_recursive_kf as aprk
import fast_hit_mc as fhmc


SCENARIO_ALIASES = {
    "ideal": "ideal",
    "with_tpc": "ideal",
    "distorted": "distorted",
    "with_tpc_distorted": "distorted",
    "silicon_only": "silicon_only",
}

SCENARIO_CODES = {
    "silicon_only": 0,
    "ideal": 1,
    "distorted": 2,
}

TRACKING_METHODS = ("inward_kf", "straight_line")
DEFAULT_TRACKING_METHOD = "inward_kf"
DISTORTION_MODELS = ("recursive_legacy", "recursive_legacy_zsmear003", "fast_hit_mc")
DEFAULT_DISTORTION_MODEL = "recursive_legacy"


@dataclass(frozen=True)
class SampleMetadata:
    sample_dir: Path
    layers: list[dict[str, Any]]

    @property
    def tpc_layer_ids(self) -> set[int]:
        return {idx for idx, layer in enumerate(self.layers) if layer["detector"] == "TPC"}

    @property
    def mvtx_layer_ids(self) -> set[int]:
        return {idx for idx, layer in enumerate(self.layers) if layer["detector"] == "MVTX"}

    @property
    def intt_layer_ids(self) -> set[int]:
        return {idx for idx, layer in enumerate(self.layers) if layer["detector"] == "INTT"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-dir",
        type=Path,
        required=True,
        help="Directory containing baseline_*.npz shards and metadata.json.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output .npz file for the reconstructed event/track tree.",
    )
    parser.add_argument(
        "--root-output",
        type=Path,
        default=None,
        help="Optional ROOT file output with a flat per-track TTree.",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=["silicon_only", "ideal", "distorted"],
        help="Scenario list: silicon_only, ideal, distorted.",
    )
    parser.add_argument(
        "--tracking-method",
        choices=TRACKING_METHODS,
        default=DEFAULT_TRACKING_METHOD,
        help="Track fit / PCA reconstruction method.",
    )
    parser.add_argument(
        "--distortion-model",
        choices=DISTORTION_MODELS,
        default=DEFAULT_DISTORTION_MODEL,
        help="TPC distortion model used for the distorted scenario.",
    )
    parser.add_argument(
        "--distortion-seed",
        type=int,
        default=12345,
        help="Random seed for stochastic distortion models.",
    )
    parser.add_argument(
        "--max-shards",
        type=int,
        default=None,
        help="Optional cap on processed baseline shards for test runs.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional cap on processed source events across all shards.",
    )
    return parser.parse_args()


def canonicalize_scenarios(names: list[str]) -> list[str]:
    canonical: list[str] = []
    seen: set[str] = set()
    for name in names:
        if name not in SCENARIO_ALIASES:
            allowed = ", ".join(sorted(SCENARIO_ALIASES))
            raise ValueError(f"unknown scenario '{name}', allowed: {allowed}")
        resolved = SCENARIO_ALIASES[name]
        if resolved not in seen:
            seen.add(resolved)
            canonical.append(resolved)
    return canonical


def load_sample_metadata(sample_dir: Path) -> SampleMetadata:
    metadata_path = sample_dir / "metadata.json"
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    return SampleMetadata(sample_dir=sample_dir.resolve(), layers=payload["layers"])


def baseline_shards(sample_dir: Path) -> list[Path]:
    shards = sorted(sample_dir.glob("baseline_*.npz"))
    if not shards:
        raise FileNotFoundError(f"no baseline shards found in {sample_dir}")
    return shards


def apply_tpc_distortion_model(
    hits: np.ndarray,
    layer_ids: np.ndarray,
    tpc_layer_ids: set[int],
    distortion_model: str = DEFAULT_DISTORTION_MODEL,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    if distortion_model == "recursive_legacy":
        return aprk.apply_offline_tpc_distortion(hits, layer_ids)
    if distortion_model == "recursive_legacy_zsmear003":
        out = np.array(hits, copy=True)
        tpc_mask = np.array([int(layer_id) in tpc_layer_ids for layer_id in layer_ids], dtype=bool)
        if not np.any(tpc_mask):
            return out
        if rng is None:
            rng = np.random.default_rng(12345)
        out[tpc_mask, 1] += 0.2
        pos_z = tpc_mask & (out[:, 2] >= 0.0)
        neg_z = tpc_mask & (out[:, 2] < 0.0)
        out[pos_z, 2] += -0.07 + rng.normal(0.0, 0.03, size=int(np.count_nonzero(pos_z)))
        out[neg_z, 2] += 0.07 + rng.normal(0.0, 0.03, size=int(np.count_nonzero(neg_z)))
        return out
    if distortion_model == "fast_hit_mc":
        out = np.array(hits, copy=True)
        for idx, layer_id in enumerate(layer_ids):
            if int(layer_id) in tpc_layer_ids:
                out[idx] = fhmc.apply_tpc_distortion(out[idx])
        return out
    raise ValueError(f"unsupported distortion model: {distortion_model}")


def scenario_hits(
    layer_ids: np.ndarray,
    hits: np.ndarray,
    scenario: str,
    tpc_layer_ids: set[int],
    distortion_model: str = DEFAULT_DISTORTION_MODEL,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if scenario == "ideal":
        return np.array(hits, copy=True), np.array(layer_ids, copy=True)
    if scenario == "distorted":
        out = apply_tpc_distortion_model(
            hits,
            layer_ids,
            tpc_layer_ids,
            distortion_model=distortion_model,
            rng=rng,
        )
        return out, np.array(layer_ids, copy=True)
    if scenario == "silicon_only":
        keep = np.array([int(layer_id) not in tpc_layer_ids for layer_id in layer_ids], dtype=bool)
        return np.array(hits[keep], copy=True), np.array(layer_ids[keep], copy=True)
    raise ValueError(f"unsupported scenario: {scenario}")


def outward_direction(direction: np.ndarray, hits: np.ndarray) -> np.ndarray:
    oriented = np.array(direction, copy=True)
    if hits.shape[0] >= 2:
        radial_step = hits[-1] - hits[0]
        if float(np.dot(oriented, radial_step)) < 0.0:
            oriented = -oriented
    return oriented / np.linalg.norm(oriented)


def track_slices(track_n_hits: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    offsets = np.zeros(track_n_hits.shape[0] + 1, dtype=np.int64)
    if track_n_hits.size > 0:
        offsets[1:] = np.cumsum(track_n_hits.astype(np.int64))
    return offsets[:-1], offsets[1:]


def metadata_layer_array(metadata: SampleMetadata, key: str, default: float = 0.0) -> np.ndarray:
    return np.asarray([float(layer.get(key, default)) for layer in metadata.layers], dtype=float)


def reconstruct_track_parameters(
    hits: np.ndarray,
    layer_ids: np.ndarray,
    momentum: float,
    eta: float,
    sigma_rphi_by_layer: np.ndarray,
    sigma_z_by_layer: np.ndarray,
    thickness_x0_by_layer: np.ndarray,
    tracking_method: str = DEFAULT_TRACKING_METHOD,
) -> tuple[np.ndarray, np.ndarray]:
    if tracking_method == "straight_line":
        point, direction = fhmc.fit_line(hits)
        direction = outward_direction(direction, hits)
        pca = fhmc.closest_point_on_line_to_vertex(point, direction, np.zeros(3, dtype=float))
        return pca, direction
    if tracking_method == "inward_kf":
        return aprk.inward_recursive_pca_direction(
            hits,
            layer_ids,
            momentum,
            eta,
            sigma_rphi_by_layer,
            sigma_z_by_layer,
            thickness_x0_by_layer,
        )
    raise ValueError(f"unsupported tracking method: {tracking_method}")


def export_tracking_tree(
    sample_dir: Path,
    output_path: Path,
    scenarios: list[str],
    tracking_method: str = DEFAULT_TRACKING_METHOD,
    distortion_model: str = DEFAULT_DISTORTION_MODEL,
    distortion_seed: int = 12345,
    max_shards: int | None = None,
    max_events: int | None = None,
) -> dict[str, Any]:
    metadata = load_sample_metadata(sample_dir)
    tpc_layer_ids = metadata.tpc_layer_ids
    sigma_rphi_by_layer = metadata_layer_array(metadata, "sigma_rphi_um") * 1.0e-4
    sigma_z_by_layer = metadata_layer_array(metadata, "sigma_z_um") * 1.0e-4
    thickness_x0_by_layer = metadata_layer_array(metadata, "thickness_x0")

    event_id: list[int] = []
    event_truth_vz: list[float] = []
    event_scenario: list[int] = []
    event_n_tracks: list[int] = []
    event_track_start: list[int] = []
    event_track_stop: list[int] = []

    track_event_id: list[int] = []
    track_scenario: list[int] = []
    track_index: list[int] = []
    track_n_hits: list[int] = []
    track_n_mvtx: list[int] = []
    track_n_intt: list[int] = []
    track_n_tpc: list[int] = []
    track_p: list[float] = []
    track_eta: list[float] = []
    track_phi: list[float] = []
    track_pca_x: list[float] = []
    track_pca_y: list[float] = []
    track_pca_z: list[float] = []
    track_dir_x: list[float] = []
    track_dir_y: list[float] = []
    track_dir_z: list[float] = []

    total_input_tracks = 0
    mvtx_layer_ids = metadata.mvtx_layer_ids
    intt_layer_ids = metadata.intt_layer_ids
    distortion_rng = np.random.default_rng(distortion_seed)

    shards = baseline_shards(sample_dir)
    if max_shards is not None:
        shards = shards[:max_shards]
    processed_source_events = 0

    for shard in shards:
        if max_events is not None and processed_source_events >= max_events:
            break
        with np.load(shard) as data:
            shard_event_id = data["track_event_id"].astype(np.int64)
            shard_track_id = data["track_id"].astype(np.int64)
            shard_track_n_hits = data["track_n_hits"].astype(np.int64)
            shard_track_truth_vz = data["track_vertex_z"].astype(np.float64) if "track_vertex_z" in data.files else np.zeros(shard_track_id.shape[0], dtype=np.float64)
            shard_hit_track_id = data["hit_track_id"].astype(np.int64)
            shard_hit_layer_id = data["hit_layer_id"].astype(np.uint8)
            shard_hits = np.column_stack((data["hit_x"], data["hit_y"], data["hit_z"])).astype(float)

            total_input_tracks += int(shard_track_id.shape[0])
            starts, stops = track_slices(shard_track_n_hits)
            track_lookup = {int(track_id): idx for idx, track_id in enumerate(shard_track_id)}

            if shard_hit_track_id.size:
                unique_hit_ids = set(np.unique(shard_hit_track_id).tolist())
                missing = [int(track_id) for track_id in shard_track_id if int(track_id) not in unique_hit_ids and int(shard_track_n_hits[track_lookup[int(track_id)]]) > 0]
                if missing:
                    raise ValueError(f"missing hit payload for track ids: {missing[:5]}")

            unique_events = np.unique(shard_event_id)
            if max_events is not None:
                remaining_events = max_events - processed_source_events
                if remaining_events <= 0:
                    break
                unique_events = unique_events[:remaining_events]
            for scenario in scenarios:
                scenario_code = SCENARIO_CODES[scenario]
                for current_event_id in unique_events:
                    current_event_id = int(current_event_id)
                    selected = np.nonzero(shard_event_id == current_event_id)[0]
                    event_start = len(track_event_id)
                    kept_in_event = 0
                    truth_vz = float(shard_track_truth_vz[selected[0]]) if selected.size else 0.0

                    for local_track_pos in selected:
                        start = int(starts[local_track_pos])
                        stop = int(stops[local_track_pos])
                        layer_ids = shard_hit_layer_id[start:stop]
                        hits = shard_hits[start:stop]
                        momentum = float(data["track_p"][local_track_pos])
                        eta = float(data["track_eta"][local_track_pos])
                        phi = float(data["track_phi"][local_track_pos])
                        scenario_specific_hits, scenario_layer_ids = scenario_hits(
                            layer_ids,
                            hits,
                            scenario,
                            tpc_layer_ids,
                            distortion_model=distortion_model,
                            rng=distortion_rng,
                        )
                        if scenario_specific_hits.shape[0] < 2:
                            continue
                        pca, direction = reconstruct_track_parameters(
                            scenario_specific_hits,
                            scenario_layer_ids,
                            momentum=momentum,
                            eta=eta,
                            sigma_rphi_by_layer=sigma_rphi_by_layer,
                            sigma_z_by_layer=sigma_z_by_layer,
                            thickness_x0_by_layer=thickness_x0_by_layer,
                            tracking_method=tracking_method,
                        )
                        n_mvtx = sum(int(layer_id) in mvtx_layer_ids for layer_id in scenario_layer_ids)
                        n_intt = sum(int(layer_id) in intt_layer_ids for layer_id in scenario_layer_ids)
                        n_tpc = sum(int(layer_id) in tpc_layer_ids for layer_id in scenario_layer_ids)

                        track_event_id.append(current_event_id)
                        track_scenario.append(scenario_code)
                        track_index.append(kept_in_event)
                        track_n_hits.append(int(scenario_specific_hits.shape[0]))
                        track_n_mvtx.append(int(n_mvtx))
                        track_n_intt.append(int(n_intt))
                        track_n_tpc.append(int(n_tpc))
                        track_p.append(momentum)
                        track_eta.append(eta)
                        track_phi.append(phi)
                        track_pca_x.append(float(pca[0]))
                        track_pca_y.append(float(pca[1]))
                        track_pca_z.append(float(pca[2]))
                        track_dir_x.append(float(direction[0]))
                        track_dir_y.append(float(direction[1]))
                        track_dir_z.append(float(direction[2]))
                        kept_in_event += 1

                    event_id.append(current_event_id)
                    event_truth_vz.append(truth_vz)
                    event_scenario.append(scenario_code)
                    event_n_tracks.append(kept_in_event)
                    event_track_start.append(event_start)
                    event_track_stop.append(len(track_event_id))
            processed_source_events += int(unique_events.shape[0])

    arrays = {
        "event_id": np.asarray(event_id, dtype=np.int32),
        "event_truth_vz": np.asarray(event_truth_vz, dtype=np.float32),
        "event_scenario": np.asarray(event_scenario, dtype=np.int8),
        "event_n_tracks": np.asarray(event_n_tracks, dtype=np.int32),
        "event_track_start": np.asarray(event_track_start, dtype=np.int64),
        "event_track_stop": np.asarray(event_track_stop, dtype=np.int64),
        "track_event_id": np.asarray(track_event_id, dtype=np.int32),
        "track_scenario": np.asarray(track_scenario, dtype=np.int8),
        "track_index": np.asarray(track_index, dtype=np.int32),
        "track_n_hits": np.asarray(track_n_hits, dtype=np.int16),
        "track_n_mvtx": np.asarray(track_n_mvtx, dtype=np.int16),
        "track_n_intt": np.asarray(track_n_intt, dtype=np.int16),
        "track_n_tpc": np.asarray(track_n_tpc, dtype=np.int16),
        "track_p": np.asarray(track_p, dtype=np.float32),
        "track_eta": np.asarray(track_eta, dtype=np.float32),
        "track_phi": np.asarray(track_phi, dtype=np.float32),
        "track_pca_x": np.asarray(track_pca_x, dtype=np.float32),
        "track_pca_y": np.asarray(track_pca_y, dtype=np.float32),
        "track_pca_z": np.asarray(track_pca_z, dtype=np.float32),
        "track_dir_x": np.asarray(track_dir_x, dtype=np.float32),
        "track_dir_y": np.asarray(track_dir_y, dtype=np.float32),
        "track_dir_z": np.asarray(track_dir_z, dtype=np.float32),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **arrays)

    summary = {
        "sample_dir": str(metadata.sample_dir),
        "output": str(output_path.resolve()),
        "scenarios": scenarios,
        "tracking_method": tracking_method,
        "distortion_model": distortion_model,
        "distortion_seed": distortion_seed,
        "n_shards_processed": len(shards),
        "n_source_events_processed": processed_source_events,
        "scenario_codes": SCENARIO_CODES,
        "n_events": int(arrays["event_id"].shape[0]),
        "n_tracks": int(arrays["track_event_id"].shape[0]),
        "n_unique_input_tracks": total_input_tracks,
        "fields": list(arrays),
    }
    summary_path = output_path.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def write_event_track_text(npz_path: Path, text_path: Path) -> None:
    with np.load(npz_path) as arrays:
        event_n_tracks = arrays["event_n_tracks"]
        event_ids = arrays["event_id"]
        event_truth_vz = arrays["event_truth_vz"] if "event_truth_vz" in arrays.files else np.zeros(arrays["event_id"].shape[0], dtype=np.float32)
        event_scenarios = arrays["event_scenario"]
        event_track_start = arrays["event_track_start"]
        event_track_stop = arrays["event_track_stop"]
        track_n_hits = arrays["track_n_hits"]
        track_n_mvtx = arrays["track_n_mvtx"]
        track_n_intt = arrays["track_n_intt"]
        track_n_tpc = arrays["track_n_tpc"]
        track_p = arrays["track_p"]
        track_eta = arrays["track_eta"]
        track_phi = arrays["track_phi"]
        track_pca_x = arrays["track_pca_x"]
        track_pca_y = arrays["track_pca_y"]
        track_pca_z = arrays["track_pca_z"]
        track_dir_x = arrays["track_dir_x"]
        track_dir_y = arrays["track_dir_y"]
        track_dir_z = arrays["track_dir_z"]

        lines = [
            "event_id\tevent_truth_vz\tscenario\tevent_n_tracks\ttrack_n_hits\ttrack_n_mvtx\ttrack_n_intt\ttrack_n_tpc\ttrack_p\ttrack_eta\ttrack_phi\tpca_x\tpca_y\tpca_z\tdir_x\tdir_y\tdir_z"
        ]
        for event_idx in range(event_ids.shape[0]):
            start = int(event_track_start[event_idx])
            stop = int(event_track_stop[event_idx])
            n_tracks = int(event_n_tracks[event_idx])
            lines.append(
                "\t".join(
                    [
                        str(int(event_ids[event_idx])),
                        f"{float(event_truth_vz[event_idx]):.9g}",
                        str(int(event_scenarios[event_idx])),
                        str(n_tracks),
                        ";".join(str(int(value)) for value in track_n_hits[start:stop]),
                        ";".join(str(int(value)) for value in track_n_mvtx[start:stop]),
                        ";".join(str(int(value)) for value in track_n_intt[start:stop]),
                        ";".join(str(int(value)) for value in track_n_tpc[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_p[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_eta[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_phi[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_pca_x[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_pca_y[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_pca_z[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_dir_x[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_dir_y[start:stop]),
                        ";".join(f"{float(value):.9g}" for value in track_dir_z[start:stop]),
                    ]
                )
            )
        text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_root_tree(npz_path: Path, root_output: Path) -> None:
    macro_path = Path(__file__).with_name("write_tracking_tree_root.C").resolve()
    if not macro_path.exists():
        raise FileNotFoundError(f"missing ROOT writer macro: {macro_path}")
    root_output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tracking_tree_root_") as tmpdir:
        text_path = Path(tmpdir) / "events.tsv"
        write_event_track_text(npz_path, text_path)
        command = [
            "root",
            "-l",
            "-b",
            "-q",
            f'{macro_path.as_posix()}("{text_path.as_posix()}","{root_output.resolve().as_posix()}")',
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(
                "ROOT writer failed\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )


def main() -> None:
    args = parse_args()
    scenarios = canonicalize_scenarios(args.scenarios)
    summary = export_tracking_tree(
        sample_dir=args.sample_dir.resolve(),
        output_path=args.output.resolve(),
        scenarios=scenarios,
        tracking_method=args.tracking_method,
        distortion_model=args.distortion_model,
        distortion_seed=args.distortion_seed,
        max_shards=args.max_shards,
        max_events=args.max_events,
    )
    if args.root_output is not None:
        write_root_tree(args.output.resolve(), args.root_output.resolve())
        summary["root_output"] = str(args.root_output.resolve())
        args.output.resolve().with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
