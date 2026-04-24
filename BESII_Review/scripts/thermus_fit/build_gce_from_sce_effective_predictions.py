#!/usr/bin/env python3
"""Prepare low-energy THERMUS GC inputs using the SCE effective endpoint trace."""

from __future__ import annotations

import csv
from pathlib import Path

OUT_DIR = Path("data/thermus_gc_from_sce_effective")
INPUT_DIR = OUT_DIR / "prediction_inputs"
MANIFEST = OUT_DIR / "prediction_manifest.csv"
PARAMETERS = OUT_DIR / "parameter_grid.csv"
GC_RESULTS = Path("data/thermus_fit_results.csv")
SCE_ENDPOINTS = Path("data/thermus_sce_blended/endpoint_parameters.csv")
SCE_WEIGHTS = Path("data/thermus_sce_blended/interpolation_grid.csv")

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


def average_muS_over_muB() -> float:
    ratios = []
    with GC_RESULTS.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            mu_b = float(row["muB_MeV"])
            mu_s = float(row["muS_MeV"])
            if mu_b > 0.0:
                ratios.append(mu_s / mu_b)
    return sum(ratios) / len(ratios)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def energy_map(rows: list[dict[str, str]]) -> dict[float, dict[str, str]]:
    return {float(row["energy_GeV"]): row for row in rows}


def write_prediction_input(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("# id<TAB>descriptor<TAB>value<TAB>error\n")
        for pdg, label in PREDICTION_PARTICLES:
            out.write(f"{pdg}\t{label}\t0.0\t0.0\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    mu_s_over_mu_b = average_muS_over_muB()
    endpoints = energy_map(read_csv_rows(SCE_ENDPOINTS))
    weights = energy_map(read_csv_rows(SCE_WEIGHTS))
    low = endpoints[3.0]
    high = endpoints[7.7]
    parameter_rows = []
    manifest_rows = []

    for energy in sorted(weights):
        weight_high = float(weights[energy]["weight_7p7GeV"])
        weight_low = 1.0 - weight_high
        t_mev = weight_low * float(low["T_MeV"]) + weight_high * float(high["T_MeV"])
        mub_mev = weight_low * float(low["muB_MeV"]) + weight_high * float(high["muB_MeV"])
        muq_mev = weight_low * float(low["muQ_MeV"]) + weight_high * float(high["muQ_MeV"])
        gamma_s = weight_low * float(low["gammaS"]) + weight_high * float(high["gammaS"])
        r_fm = weight_low * float(low["R_fm"]) + weight_high * float(high["R_fm"])
        mus_mev = mu_s_over_mu_b * mub_mev
        prediction_file = INPUT_DIR / f"sqrts_{energy_tag(energy)}GeV.txt"
        write_prediction_input(prediction_file)
        note = (
            "GC calculation using T, muB, muQ, and R from the SCE endpoint-interpolated effective trace; "
            "muS derived from the mean muS/muB ratio of existing THERMUS GC fits"
        )
        parameter_rows.append(
            {
                "energy_GeV": f"{energy:g}",
                "T_MeV": f"{t_mev:.6f}",
                "muB_MeV": f"{mub_mev:.6f}",
                "muS_MeV": f"{mus_mev:.6f}",
                "muQ_MeV": f"{muq_mev:.6f}",
                "gammaS": f"{gamma_s:.6f}",
                "R_fm": f"{r_fm:.6f}",
                "note": note,
            }
        )
        manifest_rows.append(
            {
                "energy_GeV": f"{energy:g}",
                "prediction_file": str(prediction_file),
                "T_MeV": f"{t_mev:.6f}",
                "muB_MeV": f"{mub_mev:.6f}",
                "muS_MeV": f"{mus_mev:.6f}",
                "muQ_MeV": f"{muq_mev:.6f}",
                "gammaS": f"{gamma_s:.6f}",
                "R_fm": f"{r_fm:.6f}",
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
