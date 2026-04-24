#!/usr/bin/env python3
"""Build the 7.7 GeV endpoint manifest for blended SCE yield interpolation."""

from __future__ import annotations

import csv
from pathlib import Path

OUT_DIR = Path("data/thermus_sce_blended")
INPUT_DIR = OUT_DIR / "endpoint_7p7" / "prediction_inputs"
MANIFEST = OUT_DIR / "endpoint_7p7" / "prediction_manifest.csv"
SCE_3GEV_SUMMARY = Path("data/thermus_sce_3gev/fit_summary.csv")
GC_RESULTS = Path("data/thermus_fit_results.csv")

PREDICTION_PARTICLES = [
    (211, "pi+"),
    (-211, "pi-"),
    (321, "K+"),
    (-321, "K-"),
    (310, "Ks0"),
    (2212, "p"),
    (-2212, "pbar"),
    (3122, "Lambda"),
    (-3122, "Lambda_bar"),
    (3312, "Xi"),
    (-3312, "Xi_bar"),
    (333, "phi"),
]


def read_sce_reference() -> dict[str, float]:
    with SCE_3GEV_SUMMARY.open(newline="", encoding="utf-8") as f:
        row = next(csv.DictReader(f))
    return {
        "Rc_fm": float(row["Rc_fm"]),
        "R_fm": float(row["R_fm"]),
        "gammaS": float(row["gammaS"]),
        "muQ_constraint_B2Q": float(row["muQ_constraint_B2Q"]),
    }


def read_gc_77() -> dict[str, float]:
    with GC_RESULTS.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if abs(float(row["energy_GeV"]) - 7.7) < 1e-9:
                return {
                    "T_MeV": float(row["T_MeV"]),
                    "muB_MeV": float(row["muB_MeV"]),
                    "R_fm": float(row["R_fm"]),
                }
    raise ValueError("7.7 GeV GC fit point not found in data/thermus_fit_results.csv")


def write_prediction_input(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("# id<TAB>descriptor<TAB>value<TAB>error\n")
        for pdg, label in PREDICTION_PARTICLES:
            out.write(f"{pdg}\t{label}\t0.0\t0.0\n")


def main() -> None:
    ref = read_sce_reference()
    gc77 = read_gc_77()
    rc_over_r = ref["Rc_fm"] / ref["R_fm"] if ref["R_fm"] else 0.0
    rc_77 = rc_over_r * gc77["R_fm"]

    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    pred_file = INPUT_DIR / "sqrts_7p7GeV.txt"
    write_prediction_input(pred_file)

    with MANIFEST.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "prediction_file",
                "T_MeV",
                "muB_MeV",
                "Rc_fm",
                "R_fm",
                "gammaS",
                "muQ_constraint_B2Q",
                "interpolation_note",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "energy_GeV": "7.7",
                "prediction_file": str(pred_file),
                "T_MeV": f"{gc77['T_MeV']:.6f}",
                "muB_MeV": f"{gc77['muB_MeV']:.6f}",
                "Rc_fm": f"{rc_77:.6f}",
                "R_fm": f"{gc77['R_fm']:.6f}",
                "gammaS": f"{ref['gammaS']:.6f}",
                "muQ_constraint_B2Q": f"{ref['muQ_constraint_B2Q']:.6f}",
                "interpolation_note": (
                    "7.7 GeV endpoint for blended SCE interpolation: "
                    "T and muB from the 7.7 GeV GC THERMUS fit; "
                    "R from the exact 7.7 GeV GC fit; Rc/R fixed to the 3 GeV SCE fit ratio"
                ),
            }
        )

    print(f"wrote {MANIFEST}")
    print(f"wrote {pred_file}")


if __name__ == "__main__":
    main()
