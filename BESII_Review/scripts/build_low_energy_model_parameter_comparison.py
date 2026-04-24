#!/usr/bin/env python3
"""Build a side-by-side low-energy parameter comparison table for the retained models."""

from __future__ import annotations

import csv
from pathlib import Path

GC_EFFECTIVE = Path("data/thermus_gc_from_sce_effective/parameter_grid.csv")
SCE_ENDPOINTS = Path("data/thermus_sce_blended/endpoint_parameters.csv")
SCE_WEIGHTS = Path("data/thermus_sce_blended/interpolation_grid.csv")
OUT = Path("data/low_energy_model_parameter_comparison.csv")

TARGET_ENERGIES = [3.0, 3.2, 3.5, 3.9, 4.5, 5.2, 7.7]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def energy_map(rows: list[dict[str, str]]) -> dict[float, dict[str, str]]:
    return {float(row["energy_GeV"]): row for row in rows}


def sce_effective_map() -> dict[float, dict[str, str]]:
    endpoints = energy_map(read_csv_rows(SCE_ENDPOINTS))
    weights = energy_map(read_csv_rows(SCE_WEIGHTS))
    low = endpoints[3.0]
    high = endpoints[7.7]
    out: dict[float, dict[str, str]] = {}
    for energy in TARGET_ENERGIES:
        weight_high = float(weights[energy]["weight_7p7GeV"])
        weight_low = 1.0 - weight_high
        out[energy] = {
            "sce_effective_T_MeV": f"{weight_low * float(low['T_MeV']) + weight_high * float(high['T_MeV']):.6f}",
            "sce_effective_muB_MeV": f"{weight_low * float(low['muB_MeV']) + weight_high * float(high['muB_MeV']):.6f}",
            "sce_effective_muQ_MeV": f"{weight_low * float(low['muQ_MeV']) + weight_high * float(high['muQ_MeV']):.6f}",
            "sce_effective_gammaS": f"{weight_low * float(low['gammaS']) + weight_high * float(high['gammaS']):.6f}",
            "sce_effective_Rc_fm": f"{weight_low * float(low['Rc_fm']) + weight_high * float(high['Rc_fm']):.6f}",
            "sce_effective_R_fm": f"{weight_low * float(low['R_fm']) + weight_high * float(high['R_fm']):.6f}",
        }
    return out


def main() -> None:
    gc_effective = energy_map(read_csv_rows(GC_EFFECTIVE))
    sce_effective = sce_effective_map()

    with OUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "gc_effective_T_MeV",
                "gc_effective_muB_MeV",
                "gc_effective_muS_MeV",
                "gc_effective_muQ_MeV",
                "gc_effective_gammaS",
                "gc_effective_R_fm",
                "sce_effective_T_MeV",
                "sce_effective_muB_MeV",
                "sce_effective_muQ_MeV",
                "sce_effective_gammaS",
                "sce_effective_Rc_fm",
                "sce_effective_R_fm",
            ],
        )
        writer.writeheader()
        for energy in TARGET_ENERGIES:
            writer.writerow(
                {
                    "energy_GeV": f"{energy:g}",
                    "gc_effective_T_MeV": gc_effective[energy]["T_MeV"],
                    "gc_effective_muB_MeV": gc_effective[energy]["muB_MeV"],
                    "gc_effective_muS_MeV": gc_effective[energy]["muS_MeV"],
                    "gc_effective_muQ_MeV": gc_effective[energy]["muQ_MeV"],
                    "gc_effective_gammaS": gc_effective[energy]["gammaS"],
                    "gc_effective_R_fm": gc_effective[energy]["R_fm"],
                    **sce_effective[energy],
                }
            )

    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
