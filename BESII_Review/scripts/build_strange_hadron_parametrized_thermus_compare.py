#!/usr/bin/env python3
"""Merge strange-hadron data with parametrized THERMUS GCE predictions."""

from __future__ import annotations

import csv
import math
from pathlib import Path

CANONICAL_MEASURED_CSV = Path("data/first_group_dn_dy_vs_energy.csv")
MEASURED_CSV = Path("data/strange_hadron_yields_vs_energy.csv")
PARAM_GCE_DIR = Path("data/thermus_fit_predictions_parametrized/prediction_points")
VARIATION_INDEX = Path("data/thermus_fit_predictions_parametrized/uncertainty_variations.csv")
SCE_3GEV_POINTS = Path("data/thermus_sce_3gev/fit_points/sqrts_3GeV_points.csv")
SCE_VARIATION_INDEX = Path("data/thermus_sce_3gev/uncertainty_variations.csv")
OUT_CSV = Path("data/strange_hadron_yields_parametrized_thermus.csv")

PARTICLE_ORDER = ["K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi"]
THERMUS_PARTICLES = set(PARTICLE_ORDER)


def to_float(text: str | None) -> float | None:
    if text is None or text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def normalize_particle(name: str) -> str | None:
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


def combined_error(stat: str | None, sys: str | None) -> float:
    stat_val = to_float(stat) or 0.0
    sys_val = to_float(sys) or 0.0
    return math.hypot(stat_val, sys_val) if (stat_val or sys_val) else 0.0


def row_key(energy: float, particle: str) -> tuple[float, str]:
    return (round(float(energy), 6), particle)


def empty_row(energy: float, particle: str) -> dict[str, str]:
    return {
        "energy_GeV": f"{energy:g}",
        "particle": particle,
        "data_yield": "",
        "data_err": "",
        "gc_model_yield": "",
        "gc_model_err_low": "",
        "gc_model_err_high": "",
        "sce_blended_model_yield": "",
        "sce_blended_model_err_low": "",
        "sce_blended_model_err_high": "",
        "data_source": "",
        "model_source": "",
        "note": "",
    }


def merge_measured_csv(merged: dict[tuple[float, str], dict[str, str]], path: Path, overwrite: bool = False) -> None:
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = row["particle"]
            if particle not in THERMUS_PARTICLES:
                continue
            energy = float(row["energy_GeV"])
            if energy < 3.0 or energy > 200.0:
                continue
            key = row_key(energy, particle)
            if key in merged and not overwrite:
                continue
            merged[key] = empty_row(energy, particle)
            merged[key]["data_yield"] = row["value"]
            merged[key]["data_err"] = f"{combined_error(row.get('stat_err'), row.get('sys_err')):.10g}"
            merged[key]["data_source"] = row["source"]
            merged[key]["note"] = row.get("note", "")


def load_measured_rows() -> dict[tuple[float, str], dict[str, str]]:
    merged: dict[tuple[float, str], dict[str, str]] = {}
    merge_measured_csv(merged, CANONICAL_MEASURED_CSV, overwrite=True)
    merge_measured_csv(merged, MEASURED_CSV, overwrite=False)
    return merged


def attach_parametrized_predictions(merged: dict[tuple[float, str], dict[str, str]]) -> dict[tuple[float, str], float]:
    central_predictions: dict[tuple[float, str], float] = {}
    for path in sorted(PARAM_GCE_DIR.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = normalize_particle(row["particle_name"])
                if particle not in THERMUS_PARTICLES:
                    continue
                energy = float(row["energy_GeV"])
                if energy < 7.7 or energy > 200.0:
                    continue
                key = row_key(energy, particle)
                model_yield = to_float(row["model_yield"])
                if model_yield is None:
                    continue
                central_predictions[key] = model_yield
                if key not in merged:
                    merged[key] = empty_row(energy, particle)
                merged[key]["gc_model_yield"] = row["model_yield"]
                merged[key]["model_source"] = "THERMUS GCE parametrized fit"
                if not merged[key]["data_yield"] and not merged[key]["note"]:
                    merged[key]["note"] = "model-only point from smooth parametrized GCE trace"
    return central_predictions


def load_prediction_points(points_dir: Path) -> dict[tuple[float, str], float]:
    predictions: dict[tuple[float, str], float] = {}
    for path in sorted(points_dir.glob("sqrts_*GeV_points.csv")):
        with path.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                particle = normalize_particle(row["particle_name"])
                if particle not in THERMUS_PARTICLES:
                    continue
                value = to_float(row["model_yield"])
                if value is None:
                    continue
                energy = float(row["energy_GeV"])
                predictions[row_key(energy, particle)] = value
    return predictions


def compute_prediction_uncertainties(
    index_path: Path,
    central_predictions: dict[tuple[float, str], float],
) -> dict[tuple[float, str], tuple[float, float]]:
    if not index_path.exists():
        return {}

    shifts: dict[tuple[float, str], dict[str, list[float]]] = {}
    with index_path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            parameter = row["parameter"]
            points_dir = Path(row["points"])
            if not points_dir.exists():
                continue
            varied_predictions = load_prediction_points(points_dir)
            for key, central_value in central_predictions.items():
                varied_value = varied_predictions.get(key)
                if varied_value is None:
                    continue
                shifts.setdefault(key, {}).setdefault(parameter, []).append(varied_value - central_value)

    uncertainties: dict[tuple[float, str], tuple[float, float]] = {}
    for key, parameter_shifts in shifts.items():
        low2 = 0.0
        high2 = 0.0
        for values in parameter_shifts.values():
            low = max([0.0] + [-delta for delta in values])
            high = max([0.0] + [delta for delta in values])
            low2 += low * low
            high2 += high * high
        uncertainties[key] = (math.sqrt(low2), math.sqrt(high2))
    return uncertainties


def attach_prediction_uncertainties(
    merged: dict[tuple[float, str], dict[str, str]],
    central_predictions: dict[tuple[float, str], float],
) -> None:
    uncertainties = compute_prediction_uncertainties(VARIATION_INDEX, central_predictions)
    for key, (err_low, err_high) in uncertainties.items():
        row = merged.get(key)
        if row is None:
            continue
        row["gc_model_err_low"] = f"{err_low:.10g}"
        row["gc_model_err_high"] = f"{err_high:.10g}"


def attach_sce_fit_points(merged: dict[tuple[float, str], dict[str, str]]) -> None:
    central_predictions: dict[tuple[float, str], float] = {}
    with SCE_3GEV_POINTS.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            particle = normalize_particle(row["particle_name"])
            if particle not in THERMUS_PARTICLES:
                continue
            energy = float(row["energy_GeV"])
            key = row_key(energy, particle)
            model_yield = to_float(row["model_yield"])
            if model_yield is not None:
                central_predictions[key] = model_yield
            if key not in merged:
                merged[key] = empty_row(energy, particle)
                merged[key]["data_yield"] = row.get("data_yield", "")
                merged[key]["data_err"] = row.get("data_err", "")
                merged[key]["data_source"] = str(SCE_3GEV_POINTS)
            merged[key]["sce_blended_model_yield"] = row["model_yield"]
            merged[key]["model_source"] = "THERMUS SCE fit, 3 GeV"

    uncertainties = compute_prediction_uncertainties(SCE_VARIATION_INDEX, central_predictions)
    for key, (err_low, err_high) in uncertainties.items():
        row = merged.get(key)
        if row is None:
            continue
        row["sce_blended_model_err_low"] = f"{err_low:.10g}"
        row["sce_blended_model_err_high"] = f"{err_high:.10g}"


def main() -> None:
    merged = load_measured_rows()
    central_predictions = attach_parametrized_predictions(merged)
    attach_prediction_uncertainties(merged, central_predictions)
    attach_sce_fit_points(merged)
    rows = sorted(
        merged.values(),
        key=lambda r: (float(r["energy_GeV"]), PARTICLE_ORDER.index(r["particle"])),
    )
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "data_yield",
                "data_err",
                "gc_model_yield",
                "gc_model_err_low",
                "gc_model_err_high",
                "sce_blended_model_yield",
                "sce_blended_model_err_low",
                "sce_blended_model_err_high",
                "data_source",
                "model_source",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
