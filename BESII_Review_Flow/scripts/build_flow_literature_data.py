#!/usr/bin/env python3
"""Build literature flow-data summary tables.

The raw inputs are downloaded publication/HEPData files kept under data/.
This script normalizes the compact excitation-function tables used by the
review workflow.  Older comparison points are vector-digitized from EPS
sources bundled with the corresponding arXiv submissions; those rows are
flagged in the extraction_method column.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def as_float(text: str) -> float:
    return float(text.strip())


def clean_particle(label: str) -> str:
    label = label.strip()
    replacements = {
        r"$\pi^{+}$": "pi+",
        r"$\pi^{-}$": "pi-",
        "anti-proton": "antiproton",
        "proton": "proton",
        "net-proton": "net-proton",
    }
    return replacements.get(label, label)


def parse_star_dv1_table(path: Path, table: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    particle: str | None = None
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#:"):
                continue
            if line.startswith(r"$\sqrt{s_{NN}}$"):
                match = re.search(r"for (.+?),stat", line)
                if not match:
                    raise ValueError(f"Could not parse particle from {line!r}")
                particle = clean_particle(match.group(1))
                continue
            if particle is None:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) != 6:
                continue
            sqrts, value, stat_p, stat_m, syst_p, syst_m = map(as_float, parts)
            rows.append(
                {
                    "observable": "dv1_dy_midrapidity",
                    "experiment": "STAR",
                    "collaboration": "STAR",
                    "system": "Au+Au",
                    "particle": particle,
                    "centrality": "10-40%",
                    "sqrt_s_NN_GeV": sqrts,
                    "beam_energy_A_GeV": "",
                    "value": value,
                    "stat_err_plus": abs(stat_p),
                    "stat_err_minus": abs(stat_m),
                    "syst_err_plus": abs(syst_p),
                    "syst_err_minus": abs(syst_m),
                    "source_citation": "L. Adamczyk et al. (STAR), Phys. Rev. Lett. 112, 162301 (2014)",
                    "source_doi": "10.1103/PhysRevLett.112.162301",
                    "data_source": "HEPData",
                    "data_source_doi": f"10.17182/hepdata.105867.v1/{table}",
                    "table_or_figure": table.upper(),
                    "extraction_method": "direct HEPData CSV",
                    "notes": "Cubic/linear midrapidity slope as published; HEPData errors are statistical and systematic.",
                }
            )
    return rows


def fixed_target_sqrts_from_kinetic(ekin_a_gev: float) -> float:
    mp = 0.9382720813
    e_total = ekin_a_gev + mp
    return math.sqrt(2 * mp * mp + 2 * mp * e_total)


def fixed_target_sqrts_from_momentum(p_a_gev_c: float) -> float:
    mp = 0.9382720813
    e_total = math.sqrt(p_a_gev_c * p_a_gev_c + mp * mp)
    return math.sqrt(2 * mp * mp + 2 * mp * e_total)


def star_fig3_y_to_dv1(y: float) -> float:
    return (y - 1133.0) * 0.05 / (1803.0 - 1133.0)


def add_digitized_dv1_rows(rows: list[dict[str, object]]) -> None:
    """Add E895 and NA49 proton points from STAR Fig. 3 EPS source.

    The STAR arXiv EPS places prior-experiment points in the same coordinate
    frame as the HEPData STAR values.  Coordinates and vertical error bars are
    read from fig3_new.eps.  The STAR caption states that prior-experiment cuts
    are comparable but not identical to STAR's 10-40% selection.
    """

    digitized = [
        {
            "experiment": "E895",
            "collaboration": "E895",
            "system": "Au+Au",
            "particle": "proton",
            "centrality": "comparable, not identical to STAR 10-40%",
            "sqrt_s_NN_GeV": fixed_target_sqrts_from_kinetic(6.0),
            "beam_energy_A_GeV": 6.0,
            "y": 2302.0,
            "err_px": 60.0,
            "citation": "H. Liu et al. (E895), Phys. Rev. Lett. 84, 5488 (2000)",
            "doi": "10.1103/PhysRevLett.84.5488",
            "source": "STAR Fig. 3 EPS comparison point",
        },
        {
            "experiment": "E895",
            "collaboration": "E895",
            "system": "Au+Au",
            "particle": "proton",
            "centrality": "comparable, not identical to STAR 10-40%",
            "sqrt_s_NN_GeV": fixed_target_sqrts_from_kinetic(8.0),
            "beam_energy_A_GeV": 8.0,
            "y": 1693.0,
            "err_px": 42.0,
            "citation": "H. Liu et al. (E895), Phys. Rev. Lett. 84, 5488 (2000)",
            "doi": "10.1103/PhysRevLett.84.5488",
            "source": "STAR Fig. 3 EPS comparison point",
        },
        {
            "experiment": "NA49",
            "collaboration": "NA49",
            "system": "Pb+Pb",
            "particle": "proton",
            "centrality": "comparable, not identical to STAR 10-40%",
            "sqrt_s_NN_GeV": fixed_target_sqrts_from_kinetic(40.0),
            "beam_energy_A_GeV": 40.0,
            "y": 1139.0,
            "err_px": 64.0,
            "citation": "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)",
            "doi": "10.1103/PhysRevC.68.034903",
            "source": "STAR Fig. 3 EPS comparison point",
        },
        {
            "experiment": "NA49",
            "collaboration": "NA49",
            "system": "Pb+Pb",
            "particle": "proton",
            "centrality": "comparable, not identical to STAR 10-40%",
            "sqrt_s_NN_GeV": fixed_target_sqrts_from_kinetic(158.0),
            "beam_energy_A_GeV": 158.0,
            "y": 1050.0,
            "err_px": 16.0,
            "citation": "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)",
            "doi": "10.1103/PhysRevC.68.034903",
            "source": "STAR Fig. 3 EPS comparison point",
        },
    ]

    err_scale = 0.05 / (1803.0 - 1133.0)
    for item in digitized:
        rows.append(
            {
                "observable": "dv1_dy_midrapidity",
                "experiment": item["experiment"],
                "collaboration": item["collaboration"],
                "system": item["system"],
                "particle": item["particle"],
                "centrality": item["centrality"],
                "sqrt_s_NN_GeV": item["sqrt_s_NN_GeV"],
                "beam_energy_A_GeV": item["beam_energy_A_GeV"],
                "value": star_fig3_y_to_dv1(item["y"]),
                "stat_err_plus": abs(item["err_px"] * err_scale),
                "stat_err_minus": abs(item["err_px"] * err_scale),
                "syst_err_plus": "",
                "syst_err_minus": "",
                "source_citation": item["citation"],
                "source_doi": item["doi"],
                "data_source": "arXiv EPS vector graphic",
                "data_source_doi": "10.48550/arXiv.1401.3043",
                "table_or_figure": "Fig. 3 of STAR PRL 112, 162301 EPS source",
                "extraction_method": "vector digitized from EPS coordinates",
                "notes": item["source"],
            }
        )


def na49_v2edep_y_to_v2(y: float) -> float:
    return 0.01 + (y - 440.0) * 0.05 / (1979.0 - 440.0)


def na49_v2edep_err(px: float) -> float:
    return abs(px * 0.05 / (1979.0 - 440.0))


def build_v2_cross_experiment_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    # Values from v2edep.eps in the NA49 arXiv source.  Energies are not
    # extracted from the shifted x positions; the physical sqrt(s_NN) values
    # are assigned from the cited experiments.
    items = [
        ("E877", "E877", "Au+Au", "charged hadrons", 4.72, "", 717, 62, "J. Barrette et al. (E877), Phys. Rev. C 56, 3254 (1997)", "10.1103/PhysRevC.56.3254", "standard"),
        ("NA49", "NA49", "Pb+Pb", "pi-", 8.76, 40.0, 861, 85, "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)", "10.1103/PhysRevC.68.034903", "standard/modified"),
        ("NA49", "NA49", "Pb+Pb", "pi-", 17.27, 158.0, 1130, "", "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)", "10.1103/PhysRevC.68.034903", "standard/modified"),
        ("NA49", "NA49", "Pb+Pb", "pi-", 8.76, 40.0, 804, 36, "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)", "10.1103/PhysRevC.68.034903", "cumulant"),
        ("NA49", "NA49", "Pb+Pb", "pi-", 17.27, 158.0, 1030, "", "C. Alt et al. (NA49), Phys. Rev. C 68, 034903 (2003)", "10.1103/PhysRevC.68.034903", "cumulant"),
        ("CERES", "CERES/NA45", "Pb+Au", "charged hadrons", 8.76, 40.0, 994, 77, "K. Filimonov et al. (CERES), arXiv:nucl-ex/0109017; J. Slivova, Nucl. Phys. A 715, 615c (2003)", "10.1016/S0375-9474(02)01545-2", "CERES"),
        ("CERES", "CERES/NA45", "Pb+Au", "charged hadrons", 17.27, 158.0, 1086, 77, "K. Filimonov et al. (CERES), arXiv:nucl-ex/0109017; J. Slivova, Nucl. Phys. A 715, 615c (2003)", "10.1016/S0375-9474(02)01545-2", "CERES"),
        ("STAR", "STAR", "Au+Au", "charged hadrons", 130.0, "", 1444, 80, "C. Adler et al. (STAR), Phys. Rev. C 66, 034904 (2002)", "10.1103/PhysRevC.66.034904", "STAR, low-pT corrected by NA49"),
        ("STAR", "STAR", "Au+Au", "charged hadrons", 200.0, "", 1604, 80, "L. Ray (STAR), Nucl. Phys. A 715, 45c (2003)", "10.1016/S0375-9474(02)01412-4", "STAR, low-pT corrected by NA49"),
        ("PHENIX", "PHENIX", "Au+Au", "charged hadrons", 200.0, "", 1795, 126, "S. Esumi et al. (PHENIX), Nucl. Phys. A 715, 599c (2003)", "10.1016/S0375-9474(02)01541-5", "PHENIX"),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", 130.0, "", 1610, 154, "S. Manly et al. (PHOBOS), Nucl. Phys. A 715, 611c (2003)", "10.1016/S0375-9474(02)01544-0", "PHOBOS"),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", 200.0, "", 1702, 154, "S. Manly et al. (PHOBOS), Nucl. Phys. A 715, 611c (2003)", "10.1016/S0375-9474(02)01544-0", "PHOBOS"),
    ]
    for experiment, collab, system, particle, sqrts, beam, y, err_px, citation, doi, method in items:
        rows.append(
            {
                "observable": "v2_near_midrapidity",
                "experiment": experiment,
                "collaboration": collab,
                "system": system,
                "particle": particle,
                "centrality": "mid-central, approx. 12-34% geometrical cross section",
                "sqrt_s_NN_GeV": sqrts,
                "beam_energy_A_GeV": beam,
                "value": na49_v2edep_y_to_v2(y),
                "stat_err_plus": na49_v2edep_err(err_px) if err_px != "" else "",
                "stat_err_minus": na49_v2edep_err(err_px) if err_px != "" else "",
                "syst_err_plus": "",
                "syst_err_minus": "",
                "source_citation": citation,
                "source_doi": doi,
                "data_source": "arXiv EPS vector graphic",
                "data_source_doi": "10.48550/arXiv.nucl-ex/0303001",
                "table_or_figure": "Fig. 24 / v2edep.eps of NA49 PRC 68, 034903",
                "extraction_method": "vector digitized from EPS coordinates",
                "notes": method + "; NA49 caption: near midrapidity, pions compared to charged-particle results.",
            }
        )
    return rows


def build_star_bes_eta0_rows() -> list[dict[str, object]]:
    base = [
        (7.7, 0.04515, 0.0005, 0.0023, "0.2<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub190 Fig. 8"),
        (11.5, 0.04750, 0.0002, 0.0024, "0.2<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub190 Fig. 8"),
        (19.6, 0.04995, 0.0001, 0.0025, "0.2<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub190 Fig. 8"),
        (27.0, 0.05105, 0.0001, 0.0026, "0.2<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub190 Fig. 8"),
        (39.0, 0.05275, 0.0000, 0.0026, "0.2<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub190 Fig. 8"),
        (62.4, 0.05556, 0.00012, "", "0.2<pT<2.0 GeV/c, |eta|~0.1 average", "STAR public data page pub190 Fig. 8, earlier 62.4 GeV data"),
        (200.0, 0.05810, 0.000256, "", "0.15<pT<2.0 GeV/c, |eta|~0.05 average", "STAR public data page pub108 Fig. 3"),
    ]
    rows = []
    for sqrts, value, stat, syst, notes, table_or_figure in base:
        citation = (
            "L. Adamczyk et al. (STAR), Phys. Rev. C 86, 054908 (2012)"
            if sqrts < 200
            else "B. I. Abelev et al. (STAR), Phys. Rev. C 77, 054901 (2008)"
        )
        doi = "10.1103/PhysRevC.86.054908" if sqrts < 200 else "10.1103/PhysRevC.77.054901"
        rows.append(
            {
                "observable": "v2_EP_eta0",
                "experiment": "STAR",
                "collaboration": "STAR",
                "system": "Au+Au",
                "particle": "inclusive charged hadrons",
                "centrality": "10-40%",
                "sqrt_s_NN_GeV": sqrts,
                "beam_energy_A_GeV": "",
                "value": value,
                "stat_err_plus": stat,
                "stat_err_minus": stat,
                "syst_err_plus": syst,
                "syst_err_minus": syst,
                "source_citation": citation,
                "source_doi": doi,
                "data_source": "STAR public data table",
                "data_source_doi": "",
                "table_or_figure": table_or_figure,
                "extraction_method": "midrapidity average of tabulated +/-eta bins",
                "notes": notes,
            }
        )
    return rows


def build_dv1_low_energy_rows() -> list[dict[str, object]]:
    """Additional low-energy directed-flow slope points from HEPData.

    These are kept separate from the original STAR PRL table because the
    low-energy entries come from later STAR/HADES fixed-target/SIS
    compilations with different centrality and kinematic choices.
    """

    rows: list[dict[str, object]] = []
    entries = [
        {
            "experiment": "HADES",
            "collaboration": "HADES",
            "system": "Au+Au",
            "particle": "proton",
            "centrality": "10-30%",
            "sqrts": 2.4,
            "value": 0.304,
            "stat": 0.001,
            "syst": 0.003,
            "citation": "J. Adamczewski-Musch et al. (HADES), Eur. Phys. J. A 59, 80 (2023)",
            "doi": "10.1140/epja/s10050-023-00936-6",
            "data_source": "STAR PLB 827 Fig. 5 HEPData compilation",
            "data_doi": "10.17182/hepdata.110656.v1/t25",
            "figure": "fig5_a $p$",
            "notes": "HADES/SIS18 proton point as included in the STAR 3 GeV excitation-function compilation.",
        },
        {
            "experiment": "STAR",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "proton",
            "centrality": "10-40%",
            "sqrts": 3.0,
            "value": 0.378858,
            "stat": 0.000239808,
            "syst": 0.00725157,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Lett. B 827, 137003 (2022)",
            "doi": "10.1016/j.physletb.2022.137003",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.110656.v1/t25",
            "figure": "fig5_a $p$",
            "notes": "STAR fixed-target 3 GeV proton directed-flow slope.",
        },
        {
            "experiment": "STAR FXT",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "Lambda",
            "centrality": "10-30%",
            "sqrts": 4.5,
            "value": 0.106,
            "stat": 0.011,
            "syst": 0.010,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)",
            "doi": "10.1103/PhysRevC.103.034908",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.95903.v1/t7",
            "figure": "Figure 11",
            "notes": "STAR fixed-target 4.5 GeV Lambda slope.",
        },
        {
            "experiment": "STAR FXT",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "proton",
            "centrality": "10-25%",
            "sqrts": 4.5,
            "value": 0.084,
            "stat": 0.002,
            "syst": 0.010,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)",
            "doi": "10.1103/PhysRevC.103.034908",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.95903.v1/t7",
            "figure": "Figure 11",
            "notes": "STAR fixed-target 4.5 GeV proton slope.",
        },
        {
            "experiment": "STAR FXT",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "K0S",
            "centrality": "10-30%",
            "sqrts": 4.5,
            "value": -0.030,
            "stat": 0.010,
            "syst": 0.017,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)",
            "doi": "10.1103/PhysRevC.103.034908",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.95903.v1/t7",
            "figure": "Figure 11",
            "notes": "STAR fixed-target 4.5 GeV K0S slope.",
        },
        {
            "experiment": "STAR FXT",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "pi-",
            "centrality": "10-30%",
            "sqrts": 4.5,
            "value": -0.005,
            "stat": 0.004,
            "syst": 0.003,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)",
            "doi": "10.1103/PhysRevC.103.034908",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.95903.v1/t7",
            "figure": "Figure 11",
            "notes": "STAR fixed-target 4.5 GeV negative-pion slope.",
        },
        {
            "experiment": "STAR FXT",
            "collaboration": "STAR",
            "system": "Au+Au",
            "particle": "pi+",
            "centrality": "10-30%",
            "sqrts": 4.5,
            "value": -0.024,
            "stat": 0.004,
            "syst": 0.003,
            "citation": "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)",
            "doi": "10.1103/PhysRevC.103.034908",
            "data_source": "HEPData",
            "data_doi": "10.17182/hepdata.95903.v1/t7",
            "figure": "Figure 11",
            "notes": "STAR fixed-target 4.5 GeV positive-pion slope.",
        },
    ]
    for entry in entries:
        rows.append(
            {
                "observable": "dv1_dy_midrapidity",
                "experiment": entry["experiment"],
                "collaboration": entry["collaboration"],
                "system": entry["system"],
                "particle": entry["particle"],
                "centrality": entry["centrality"],
                "sqrt_s_NN_GeV": entry["sqrts"],
                "beam_energy_A_GeV": "",
                "value": entry["value"],
                "stat_err_plus": entry["stat"],
                "stat_err_minus": entry["stat"],
                "syst_err_plus": entry["syst"],
                "syst_err_minus": entry["syst"],
                "source_citation": entry["citation"],
                "source_doi": entry["doi"],
                "data_source": entry["data_source"],
                "data_source_doi": entry["data_doi"],
                "table_or_figure": entry["figure"],
                "extraction_method": "direct HEPData CSV",
                "notes": entry["notes"],
            }
        )
    return rows


def build_v2_expanded_rows() -> list[dict[str, object]]:
    """Expanded v2 excitation function from STAR FXT Fig. 14 HEPData."""

    rows: list[dict[str, object]] = []
    source = "M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021)"
    source_doi = "10.1103/PhysRevC.103.034908"
    data_doi = "10.17182/hepdata.95903.v1/t9"

    entries = [
        ("STAR FXT", "STAR", "Au+Au", "proton", "0-30%", 4.6, 0.02831, 0.006, "STAR fixed-target proton point."),
        ("STAR FXT", "STAR", "Au+Au", "pion", "0-30%", 4.4, 0.02355, 0.00388, "STAR fixed-target pion point."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 1.907, 0.071, 0.008, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 1.922, 0.030, 0.005, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 1.936, -0.010, 0.004, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 1.984, -0.058, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.053, -0.080, 0.012, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.142, -0.062, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.227, -0.061, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.309, -0.057, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.389, -0.052, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("FOPI", "FOPI", "Au+Au", "charged/Z=1 hadrons", "0.25 < b0 < 0.45", 2.503, -0.041, 0.007, "FOPI SIS charged/Z=1 hadron excitation function."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 2.369, -0.065, 0.007, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 2.683, -0.047, 0.005, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 3.305, -0.002, 0.004, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 3.827, 0.016, 0.003, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 4.287, 0.019, 0.006, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("EOS/E895/E877", "EOS/E895/E877", "Au+Au", "proton", "5-7 fm", 4.857, 0.022, 0.005, "AGS/Bevalac proton comparison set compiled in STAR FXT Fig. 14."),
        ("CERES", "CERES/NA45", "Pb+Au", "charged hadrons", "", 8.84, 0.028, "", "CERES charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("CERES", "CERES/NA45", "Pb+Au", "charged hadrons", "", 12.36, 0.037, "", "CERES charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("CERES", "CERES/NA45", "Pb+Au", "charged hadrons", "", 17.26, 0.043, "", "CERES charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("NA49", "NA49", "Pb+Pb", "charged hadrons", "12-34%", 8.84, 0.024, 0.003, "NA49 charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("NA49", "NA49", "Pb+Pb", "charged hadrons", "12-34%", 17.26, 0.034, 0.001, "NA49 charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 7.7, 0.0357, 0.0009, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 11.5, 0.0402, 0.0005, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 14.5, 0.0412, 0.0005, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 19.6, 0.0412, 0.0003, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 27.0, 0.0412, 0.0002, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 39.0, 0.0431, 0.0002, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "proton", "10-40%", 62.4, 0.0527, 0.0002, "STAR proton comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 7.7, 0.0334, 0.0005, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 11.5, 0.0348, 0.0002, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 14.5, 0.0365, 0.0002, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 19.6, 0.0355, 0.0001, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 27.0, 0.0366, 0.0001, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 39.0, 0.0381, 0.0001, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("STAR", "STAR", "Au+Au", "pion", "10-40%", 62.4, 0.0387, 0.0001, "STAR pion comparison set compiled in STAR FXT Fig. 14."),
        ("PHENIX", "PHENIX", "Au+Au", "proton", "20-30%", 64.4, 0.045, 0.005, "PHENIX proton comparison set compiled in STAR FXT Fig. 14."),
        ("PHENIX", "PHENIX", "Au+Au", "proton", "20-30%", 200.0, 0.060, 0.008, "PHENIX proton comparison set compiled in STAR FXT Fig. 14."),
        ("PHENIX", "PHENIX", "Au+Au", "pion", "20-30%", 60.4, 0.040, 0.002, "PHENIX pion comparison set compiled in STAR FXT Fig. 14."),
        ("PHENIX", "PHENIX", "Au+Au", "pion", "20-30%", 200.0, 0.040, 0.002, "PHENIX pion comparison set compiled in STAR FXT Fig. 14."),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", "0-40%", 19.6, 0.033, 0.005, "PHOBOS charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", "0-40%", 58.4, 0.044, 0.003, "PHOBOS charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", "0-40%", 130.0, 0.048, 0.003, "PHOBOS charged-hadron comparison set compiled in STAR FXT Fig. 14."),
        ("PHOBOS", "PHOBOS", "Au+Au", "charged hadrons", "0-40%", 200.0, 0.052, 0.003, "PHOBOS charged-hadron comparison set compiled in STAR FXT Fig. 14."),
    ]
    for experiment, collab, system, particle, centrality, sqrts, value, stat, notes in entries:
        rows.append(
            {
                "observable": "v2_midrapidity_expanded",
                "experiment": experiment,
                "collaboration": collab,
                "system": system,
                "particle": particle,
                "centrality": centrality,
                "sqrt_s_NN_GeV": sqrts,
                "beam_energy_A_GeV": "",
                "value": value,
                "stat_err_plus": stat,
                "stat_err_minus": stat,
                "syst_err_plus": "",
                "syst_err_minus": "",
                "source_citation": source,
                "source_doi": source_doi,
                "data_source": "HEPData",
                "data_source_doi": data_doi,
                "table_or_figure": "Figure 14",
                "extraction_method": "direct HEPData CSV",
                "notes": notes,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "observable",
        "experiment",
        "collaboration",
        "system",
        "particle",
        "centrality",
        "sqrt_s_NN_GeV",
        "beam_energy_A_GeV",
        "value",
        "stat_err_plus",
        "stat_err_minus",
        "syst_err_plus",
        "syst_err_minus",
        "source_citation",
        "source_doi",
        "data_source",
        "data_source_doi",
        "table_or_figure",
        "extraction_method",
        "notes",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def star_dv1_table_path(name: str) -> Path:
    direct = DATA / "raw_star_prl112_162301" / name
    if direct.exists():
        return direct
    matches = sorted((DATA / "raw_star_prl112_162301").glob(f"**/{name}"))
    if not matches:
        raise FileNotFoundError(f"Could not find {name} under data/raw_star_prl112_162301")
    return matches[0]


def main() -> None:
    dv1_rows = []
    dv1_rows.extend(parse_star_dv1_table(star_dv1_table_path("Table5.csv"), "t5"))
    dv1_rows.extend(parse_star_dv1_table(star_dv1_table_path("Table6.csv"), "t6"))
    add_digitized_dv1_rows(dv1_rows)
    dv1_rows.sort(key=lambda r: (str(r["particle"]), str(r["experiment"]), float(r["sqrt_s_NN_GeV"])))
    write_csv(DATA / "flow_dv1dy_vs_energy.csv", dv1_rows)

    write_csv(DATA / "flow_dv1dy_low_energy_expanded.csv", build_dv1_low_energy_rows())
    write_csv(DATA / "flow_v2_cross_experiment_vs_energy.csv", build_v2_cross_experiment_rows())
    write_csv(DATA / "flow_v2_star_bes_eta0_vs_energy.csv", build_star_bes_eta0_rows())
    write_csv(DATA / "flow_v2_low_energy_expanded.csv", build_v2_expanded_rows())


if __name__ == "__main__":
    main()
