#!/usr/bin/env python3
"""Build a parameter summary table for the retained GC and SCE variants."""

from __future__ import annotations

import csv
from pathlib import Path

GC_LOW = Path("data/thermus_gc_from_sce_effective/parameter_grid.csv")
SCE_BLEND_ENDPOINTS = Path("data/thermus_sce_blended/endpoint_parameters.csv")
SCE_BLEND_WEIGHTS = Path("data/thermus_sce_blended/interpolation_grid.csv")
OUT = Path("data/model_parameter_summary_3to7p7GeV.csv")

TARGET_ENERGIES = [3.0, 3.2, 3.5, 3.9, 4.5, 5.2, 7.7]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def energy_map(rows: list[dict[str, str]]) -> dict[float, dict[str, str]]:
    return {float(r["energy_GeV"]): r for r in rows}


def gc_map() -> dict[float, dict[str, str]]:
    return {
        float(r["energy_GeV"]): {
            "T_MeV": r["T_MeV"],
            "muB_MeV": r["muB_MeV"],
            "muQ_MeV": r["muQ_MeV"],
            "gammaS": r["gammaS"],
            "Rc_fm": "",
            "R_fm": r["R_fm"],
            "note": r["note"],
        }
        for r in read_csv_rows(GC_LOW)
    }


def sce_blended_effective_map() -> dict[float, dict[str, str]]:
    endpoint_rows = energy_map(read_csv_rows(SCE_BLEND_ENDPOINTS))
    weights = energy_map(read_csv_rows(SCE_BLEND_WEIGHTS))

    low = endpoint_rows[3.0]
    high = endpoint_rows[7.7]
    out: dict[float, dict[str, str]] = {}
    for energy in TARGET_ENERGIES:
        w_high = float(weights[energy]["weight_7p7GeV"])
        w_low = 1.0 - w_high
        t = w_low * float(low["T_MeV"]) + w_high * float(high["T_MeV"])
        mub = w_low * float(low["muB_MeV"]) + w_high * float(high["muB_MeV"])
        muq = w_low * float(low["muQ_MeV"]) + w_high * float(high["muQ_MeV"])
        gammas = w_low * float(low["gammaS"]) + w_high * float(high["gammaS"])
        rc = w_low * float(low["Rc_fm"]) + w_high * float(high["Rc_fm"])
        r = w_low * float(low["R_fm"]) + w_high * float(high["R_fm"])
        out[energy] = {
            "T_MeV": f"{t:.6f}",
            "muB_MeV": f"{mub:.6f}",
            "muQ_MeV": f"{muq:.6f}",
            "gammaS": f"{gammas:.6f}",
            "Rc_fm": f"{rc:.6f}",
            "R_fm": f"{r:.6f}",
            "note": (
                "Effective endpoint-weighted parameter trace using the same accelerated 3.0/7.7 "
                "weights as the blended-yield model; blended SCE yields themselves are interpolated "
                "directly and are not generated from these intermediate parameters"
            ),
        }
    return out


def main() -> None:
    gc = gc_map()
    sce_blended = sce_blended_effective_map()

    with OUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "model",
                "T_MeV",
                "muB_MeV",
                "muQ_MeV",
                "gammaS",
                "Rc_fm",
                "R_fm",
                "note",
            ],
        )
        writer.writeheader()
        for energy in TARGET_ENERGIES:
            for model_name, mapping in [
                ("GC_effective_trace", gc),
                ("SCE_endpoint_interpolated_effective", sce_blended),
            ]:
                row = mapping[energy]
                writer.writerow(
                    {
                        "energy_GeV": f"{energy:g}",
                        "model": model_name,
                        "T_MeV": row["T_MeV"],
                        "muB_MeV": row["muB_MeV"],
                        "muQ_MeV": row["muQ_MeV"],
                        "gammaS": row["gammaS"],
                        "Rc_fm": row["Rc_fm"],
                        "R_fm": row["R_fm"],
                        "note": row["note"],
                    }
                )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
