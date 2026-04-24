#!/usr/bin/env python3
"""Build a literature thermal-fit table for central heavy-ion collisions."""

from __future__ import annotations

import csv
from pathlib import Path

import urllib.request
from urllib.parse import quote


OUT = Path("data/literature_thermal_fit_points_selected.csv")


def fetch_csv_lines(url: str) -> list[str]:
    text = urllib.request.urlopen(quote(url, safe=":/?=&")).read().decode("utf-8", "ignore")
    return text.splitlines()


def load_existing_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_sections(lines: list[str], xname: str) -> dict[str, list[tuple[float, float, float, float]]]:
    sections: dict[str, list[tuple[float, float, float, float]]] = {}
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if not ln or ln.startswith("#"):
            i += 1
            continue
        if ln.startswith(xname + ","):
            hdr = [h.strip() for h in ln.split(",")]
            label = hdr[1]
            i += 1
            rows: list[tuple[float, float, float, float]] = []
            while i < len(lines):
                l = lines[i].strip()
                if not l:
                    break
                if l.startswith("#"):
                    i += 1
                    continue
                if l.startswith(xname + ","):
                    break
                parts = [p.strip() for p in l.split(",")]
                if len(parts) >= 4:
                    try:
                        e = float(parts[0])
                    except ValueError:
                        i += 1
                        continue
                    ytxt = parts[1]
                    if ytxt not in {"", "-"}:
                        rows.append((e, float(ytxt), float(parts[2]), abs(float(parts[3]))))
                i += 1
            sections[label] = rows
            continue
        i += 1
    return sections


def main() -> None:
    fig11_up = "https://www.hepdata.net/download/table/ins871561/Figure 11 up/csv"
    fig11_down = "https://www.hepdata.net/download/table/ins871561/Figure 11 down/csv"

    rows: list[dict[str, str]] = []
    try:
        up = parse_sections(fetch_csv_lines(fig11_up), "Temp [GeV]")
        dn = parse_sections(fetch_csv_lines(fig11_down), "mu_B [GeV]")

        # Keep low-energy external points and STAR original-paper points from the same compilation;
        # exclude the generic "RHIC" section.
        for dataset in ["SIS", "AGS", "SPS", "this work"]:
            trows = up.get(dataset, [])
            mrows = dn.get(dataset, [])
            for tr, mr in zip(trows, mrows):
                et, tgev, te_plus, te_minus = tr
                em, mbgev, mb_plus, mb_minus = mr
                if abs(et - em) > 1e-9:
                    continue
                exp_group = "STAR (orig paper)" if dataset == "this work" else dataset
                rows.append(
                    {
                        "experiment_group": exp_group,
                        "collision_system": "A+A",
                        "energy_GeV": f"{et:g}",
                        "centrality": "central (see source refs)",
                        "fit_variant": "compiled point",
                        "T_MeV": f"{1000.0 * tgev:.6g}",
                        "T_err_plus_MeV": f"{1000.0 * te_plus:.6g}",
                        "T_err_minus_MeV": f"{1000.0 * te_minus:.6g}",
                        "muB_MeV": f"{1000.0 * mbgev:.6g}",
                        "muB_err_plus_MeV": f"{1000.0 * mb_plus:.6g}",
                        "muB_err_minus_MeV": f"{1000.0 * mb_minus:.6g}",
                        "muB_note": "",
                        "source": "STAR PRC83 Fig.11 literature compilation (Ref.[22])",
                        "source_doi_or_url": "10.17182/hepdata.96847.v2/t19,t20; https://www.hepdata.net/record/ins871561",
                    }
                )
    except Exception as exc:
        print(f"warning: falling back to existing {OUT} for external rows: {exc}")
        for row in load_existing_rows(OUT):
            if row["experiment_group"] in {"STAR (this work fit)", "ALICE"}:
                continue
            rows.append(row)

    # Add ALICE point from Andronic et al. Nature 561 (2018).
    rows.append(
        {
            "experiment_group": "ALICE",
            "collision_system": "Pb+Pb",
            "energy_GeV": "2760",
            "centrality": "central (see source refs)",
            "fit_variant": "SHM point from Nature 2018",
            "T_MeV": "156.5",
            "T_err_plus_MeV": "1.5",
            "T_err_minus_MeV": "1.5",
            "muB_MeV": "0.7",
            "muB_err_plus_MeV": "3.8",
            "muB_err_minus_MeV": "3.8",
            "muB_note": "",
            "source": "Andronic et al., Nature 561 (2018) 321-330",
            "source_doi_or_url": "10.1038/s41586-018-0491-6; https://arxiv.org/abs/1710.09425",
        }
    )

    # Add STAR points from our latest local THERMUS fit (same energies as STAR original-paper points).
    thermus_path = Path("data/thermus_fit_results.csv")
    with thermus_path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            e = float(row["energy_GeV"])
            if abs(e - 62.4) > 1e-6 and abs(e - 200.0) > 1e-6:
                continue
            rows.append(
                {
                    "experiment_group": "STAR (this work fit)",
                    "collision_system": "Au+Au",
                    "energy_GeV": f"{e:g}",
                    "centrality": "0-5% (fit input)",
                    "fit_variant": "THERMUS fit",
                    "T_MeV": row["T_MeV"],
                    "T_err_plus_MeV": row["T_err_MeV"],
                    "T_err_minus_MeV": row["T_err_MeV"],
                    "muB_MeV": row["muB_MeV"],
                    "muB_err_plus_MeV": row["muB_err_MeV"],
                    "muB_err_minus_MeV": row["muB_err_MeV"],
                    "muB_note": "",
                    "source": "Local THERMUS fit in this repository",
                    "source_doi_or_url": "data/thermus_fit_results.csv",
                }
            )

    rows.sort(key=lambda r: (float(r["energy_GeV"]), r["experiment_group"], r["fit_variant"]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        fields = [
            "experiment_group",
            "collision_system",
            "energy_GeV",
            "centrality",
            "fit_variant",
            "T_MeV",
            "T_err_plus_MeV",
            "T_err_minus_MeV",
            "muB_MeV",
            "muB_err_plus_MeV",
            "muB_err_minus_MeV",
            "muB_note",
            "source",
            "source_doi_or_url",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT} with {len(rows)} rows")


if __name__ == "__main__":
    main()
