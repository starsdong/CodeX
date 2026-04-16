#!/usr/bin/env python3
"""Build per-energy ThermalFIST fit inputs from BESII_Review CSV data."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

PARTICLE_TO_PDG = {
    "pi+": 211,
    "pi-": -211,
    "K+": 321,
    "K-": -321,
    "p": 2212,
    "pbar": -2212,
    "Lambda": 3122,
    "Lambda_bar": -3122,
    "Xi": 3312,
    "Xi_bar": -3312,
    "Omega": 3334,
    "Omega_bar": -3334,
    "phi": 333,
    "d": 1000010020,
    "t": 1000010030,
}


def to_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def combined_error(stat_err: float | None, sys_err: float | None) -> float | None:
    if stat_err is None and sys_err is None:
        return None
    if stat_err is None:
        return sys_err
    if sys_err is None:
        return stat_err
    return math.hypot(stat_err, sys_err)


def energy_tag(energy_gev: float) -> str:
    text = f"{energy_gev:g}"
    return text.replace(".", "p")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare ThermalFIST per-energy fit input files from CSV data.")
    parser.add_argument("--input", default="data/first_group_dn_dy_vs_energy.csv", help="Input CSV with experimental yields")
    parser.add_argument("--outdir", default="data/thermalfist_inputs", help="Output directory for *.dat files")
    parser.add_argument("--manifest", default="data/thermalfist_inputs/manifest.csv", help="Output manifest CSV")
    parser.add_argument("--min-points", type=int, default=4, help="Minimum number of points required per energy")
    parser.add_argument(
        "--lambda-to-p-br",
        type=float,
        default=0.639,
        help="Branching ratio used for Lambda->p+pi and anti-Lambda->anti-p+pi feeddown subtraction",
    )
    parser.add_argument(
        "--fit-particles",
        default="pi+,pi-,K+,K-,p,pbar,Lambda,Lambda_bar,Xi,Xi_bar",
        help="Comma-separated particle names to include in fit files",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    manifest_path = Path(args.manifest)

    outdir.mkdir(parents=True, exist_ok=True)
    fit_particles = {p.strip() for p in args.fit_particles.split(",") if p.strip()}

    by_energy: dict[float, list[tuple[int, float, float]]] = defaultdict(list)
    by_energy_particle: dict[float, dict[str, list[tuple[float, float]]]] = defaultdict(lambda: defaultdict(list))
    skipped = 0

    with input_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            obs = (row.get("observable") or "").lower()
            if "dn/dy" not in obs:
                continue

            particle = (row.get("particle") or "").strip()
            pdg = PARTICLE_TO_PDG.get(particle)
            if pdg is None:
                skipped += 1
                continue

            energy = to_float(row.get("energy_GeV"))
            value = to_float(row.get("value"))
            stat = to_float(row.get("stat_err"))
            sys = to_float(row.get("sys_err"))
            err = combined_error(stat, sys)

            if energy is None or value is None or err is None or err <= 0.0:
                skipped += 1
                continue

            by_energy_particle[energy][particle].append((value, err))
            if particle in fit_particles:
                by_energy[energy].append((pdg, value, err))

    # Proton feeddown correction from (anti-)Lambda decays:
    # p_corr = p_raw - BR * Lambda
    # pbar_corr = pbar_raw - BR * Lambda_bar
    # Error propagated in quadrature.
    corrected_rows = 0
    for energy, rows in by_energy.items():
        species = by_energy_particle.get(energy, {})
        lambda_pts = species.get("Lambda", [])
        lambda_bar_pts = species.get("Lambda_bar", [])

        if lambda_pts:
            lambda_val = sum(v for v, _ in lambda_pts) / len(lambda_pts)
            lambda_err = (sum(e * e for _, e in lambda_pts) / len(lambda_pts)) ** 0.5
            for i, (pdg, value, err) in enumerate(rows):
                if pdg == 2212:
                    new_value = value - args.lambda_to_p_br * lambda_val
                    new_err = math.hypot(err, args.lambda_to_p_br * lambda_err)
                    rows[i] = (pdg, new_value, new_err)
                    corrected_rows += 1

        if lambda_bar_pts:
            lambda_bar_val = sum(v for v, _ in lambda_bar_pts) / len(lambda_bar_pts)
            lambda_bar_err = (sum(e * e for _, e in lambda_bar_pts) / len(lambda_bar_pts)) ** 0.5
            for i, (pdg, value, err) in enumerate(rows):
                if pdg == -2212:
                    new_value = value - args.lambda_to_p_br * lambda_bar_val
                    new_err = math.hypot(err, args.lambda_to_p_br * lambda_bar_err)
                    rows[i] = (pdg, new_value, new_err)
                    corrected_rows += 1

    manifest_rows = []

    for energy in sorted(by_energy):
        rows = by_energy[energy]
        if len(rows) < args.min_points:
            continue

        outfile = outdir / f"sqrts_{energy_tag(energy)}GeV.dat"
        with outfile.open("w", encoding="utf-8") as out:
            out.write("# Auto-generated from BESII_Review data\n")
            out.write("#     is_fitted           pdg1           pdg2      feeddown1      feeddown2          value          error\n")
            for pdg, value, err in sorted(rows, key=lambda x: x[0]):
                out.write(f"{1:15d}{pdg:15d}{0:15d}{3:15d}{0:15d}{value:15.8g}{err:15.8g}\n")

        manifest_rows.append({
            "energy_GeV": f"{energy:g}",
            "fit_file": str(outfile),
            "n_points": str(len(rows)),
        })

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["energy_GeV", "fit_file", "n_points"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"wrote {len(manifest_rows)} fit files to {outdir}")
    print(f"wrote manifest: {manifest_path}")
    print(f"skipped rows: {skipped}")
    print(f"feeddown-corrected proton rows: {corrected_rows}")


if __name__ == "__main__":
    main()
