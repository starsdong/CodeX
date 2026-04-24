#!/usr/bin/env python3
"""Build per-energy THERMUS fit inputs from BESII data."""

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
    "Ks0": 310,
    "phi": 333,
    "p": 2212,
    "pbar": -2212,
    "Lambda": 3122,
    "Lambda_bar": -3122,
    "Xi": 3312,
    "Xi_bar": -3312,
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
    return f"{energy_gev:g}".replace(".", "p")


def parse_exclusions(spec: str) -> dict[float, set[str]]:
    mapping: dict[float, set[str]] = {}
    text = (spec or "").strip()
    if not text:
        return mapping
    for block in text.split(";"):
        item = block.strip()
        if not item or ":" not in item:
            continue
        e_txt, plist_txt = item.split(":", 1)
        try:
            energy = float(e_txt.strip())
        except ValueError:
            continue
        particles = {p.strip() for p in plist_txt.split(",") if p.strip()}
        if particles:
            mapping[energy] = particles
    return mapping


def parse_energies(spec: str) -> set[float] | None:
    text = (spec or "").strip()
    if not text:
        return None
    energies: set[float] = set()
    for item in text.split(","):
        token = item.strip()
        if not token:
            continue
        try:
            energies.add(float(token))
        except ValueError:
            continue
    return energies or None


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare THERMUS per-energy fit input files from CSV data.")
    parser.add_argument("--input", default="data/first_group_dn_dy_vs_energy.csv")
    parser.add_argument("--outdir", default="data/thermus_inputs")
    parser.add_argument("--manifest", default="data/thermus_inputs/manifest.csv")
    parser.add_argument("--min-points", type=int, default=4)
    parser.add_argument("--lambda-to-p-br", type=float, default=0.639)
    parser.add_argument(
        "--fit-particles",
        default="pi+,pi-,K+,K-,p,pbar,Lambda,Lambda_bar,Xi,Xi_bar,phi",
        help="Comma-separated particle names to include in fit files",
    )
    parser.add_argument(
        "--exclude-energy-particles",
        default="",
        help="Semicolon list like '62.4:Xi,Xi_bar;200:Omega' to exclude species only at selected energies",
    )
    parser.add_argument(
        "--energies",
        default="",
        help="Optional comma-separated energy list to keep, e.g. '3.0,7.7'",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    manifest_path = Path(args.manifest)
    outdir.mkdir(parents=True, exist_ok=True)

    fit_particles = {p.strip() for p in args.fit_particles.split(",") if p.strip()}
    exclusions = parse_exclusions(args.exclude_energy_particles)
    energy_filter = parse_energies(args.energies)
    by_energy_rows: dict[float, list[tuple[int, str, float, float]]] = defaultdict(list)
    by_energy_particle: dict[float, dict[str, tuple[float, float]]] = defaultdict(dict)
    skipped = 0

    with input_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "dn/dy" not in (row.get("observable") or "").lower():
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
            if energy_filter is not None and energy not in energy_filter:
                continue

            by_energy_particle[energy][particle] = (value, err)
            excluded_here = particle in exclusions.get(energy, set())
            if particle in fit_particles and not excluded_here:
                by_energy_rows[energy].append((pdg, particle, value, err))

    corrected_rows = 0
    for energy, rows in by_energy_rows.items():
        species = by_energy_particle.get(energy, {})

        lmb = species.get("Lambda")
        if lmb is not None:
            lmb_val, lmb_err = lmb
            for i, (pdg, particle, value, err) in enumerate(rows):
                if pdg == 2212:
                    rows[i] = (pdg, particle, value - args.lambda_to_p_br * lmb_val, math.hypot(err, args.lambda_to_p_br * lmb_err))
                    corrected_rows += 1

        almb = species.get("Lambda_bar")
        if almb is not None:
            almb_val, almb_err = almb
            for i, (pdg, particle, value, err) in enumerate(rows):
                if pdg == -2212:
                    rows[i] = (pdg, particle, value - args.lambda_to_p_br * almb_val, math.hypot(err, args.lambda_to_p_br * almb_err))
                    corrected_rows += 1

    manifest_rows = []
    for energy in sorted(by_energy_rows):
        rows = by_energy_rows[energy]
        if len(rows) < args.min_points:
            continue

        outfile = outdir / f"sqrts_{energy_tag(energy)}GeV.txt"
        with outfile.open("w", encoding="utf-8") as out:
            out.write("# id<TAB>descriptor<TAB>value<TAB>error\n")
            for pdg, particle, value, err in sorted(rows, key=lambda x: x[0]):
                out.write(f"{pdg}\t{particle}\t{value:.10g}\t{err:.10g}\n")

        manifest_rows.append(
            {
                "energy_GeV": f"{energy:g}",
                "fit_file": str(outfile),
                "n_points": str(len(rows)),
            }
        )

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
