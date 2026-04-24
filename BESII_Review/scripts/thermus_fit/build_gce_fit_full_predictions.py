#!/usr/bin/env python3
"""Prepare full-particle THERMUS GC predictions from fitted high-energy parameters."""

from __future__ import annotations

import csv
from pathlib import Path

FIT_RESULTS = Path("data/thermus_fit_results.csv")
OUT_DIR = Path("data/thermus_fit_predictions_full")
INPUT_DIR = OUT_DIR / "prediction_inputs"
MANIFEST = OUT_DIR / "prediction_manifest.csv"
PARAMETERS = OUT_DIR / "parameter_grid.csv"

PREDICTION_PARTICLES = [
    (211, "pi+"),
    (-211, "pi-"),
    (2212, "p"),
    (-2212, "pbar"),
    (321, "K+"),
    (-321, "K-"),
    (310, "Ks0"),
    (3122, "Lambda"),
    (-3122, "Lambda_bar"),
    (3312, "Xi"),
    (-3312, "Xi_bar"),
    (333, "phi"),
]


def energy_tag(energy_gev: float) -> str:
    return f"{energy_gev:.1f}".replace(".", "p")


def write_prediction_input(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("# id<TAB>descriptor<TAB>value<TAB>error\n")
        for pdg, label in PREDICTION_PARTICLES:
            out.write(f"{pdg}\t{label}\t0.0\t0.0\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    parameter_rows = []
    manifest_rows = []

    with FIT_RESULTS.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            energy = float(row["energy_GeV"])
            prediction_file = INPUT_DIR / f"sqrts_{energy_tag(energy)}GeV.txt"
            write_prediction_input(prediction_file)
            note = "GC prediction using exact fitted high-energy THERMUS parameters"
            parameter_rows.append(
                {
                    "energy_GeV": f"{energy:g}",
                    "T_MeV": row["T_MeV"],
                    "muB_MeV": row["muB_MeV"],
                    "muS_MeV": row["muS_MeV"],
                    "muQ_MeV": "0.000000",
                    "gammaS": row["gammaS"],
                    "R_fm": row["R_fm"],
                    "note": note,
                }
            )
            manifest_rows.append(
                {
                    "energy_GeV": f"{energy:g}",
                    "prediction_file": str(prediction_file),
                    "T_MeV": row["T_MeV"],
                    "muB_MeV": row["muB_MeV"],
                    "muS_MeV": row["muS_MeV"],
                    "muQ_MeV": "0.000000",
                    "gammaS": row["gammaS"],
                    "R_fm": row["R_fm"],
                    "interpolation_note": note,
                }
            )

    with PARAMETERS.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["energy_GeV", "T_MeV", "muB_MeV", "muS_MeV", "muQ_MeV", "gammaS", "R_fm", "note"],
        )
        writer.writeheader()
        writer.writerows(parameter_rows)

    with MANIFEST.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "prediction_file",
                "T_MeV",
                "muB_MeV",
                "muS_MeV",
                "muQ_MeV",
                "gammaS",
                "R_fm",
                "interpolation_note",
            ],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"wrote {PARAMETERS}")
    print(f"wrote {MANIFEST}")
    print(f"wrote {len(manifest_rows)} prediction input files to {INPUT_DIR}")


if __name__ == "__main__":
    main()
