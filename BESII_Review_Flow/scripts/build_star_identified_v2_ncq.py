#!/usr/bin/env python3
"""Build STAR identified-particle v2 NCQ-scaling tables and figures."""

from __future__ import annotations

import csv
import html.parser
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / ".matplotlib-cache"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures" / "star_identified_v2_ncq"
RAW3 = DATA / "raw_star_3gev_hepdata_csv" / "HEPData-ins1897294-v1-csv"
RAW45 = DATA / "raw_star_4p5gev_hepdata_csv" / "HEPData-ins1809043-v1-csv"
RAWLOW = DATA / "raw_star_lowenergy_ncq_hepdata_csv" / "HEPData-ins2907591-v1-csv"
RAWBES = DATA / "raw_star_bes_identified_v2_hepdata_csv" / "HEPData-ins1210464-v1-csv"
STAR200 = DATA / "star_pub108_charged_strange_v2_200gev_data.html"
STAR241 = DATA / "star_pub241_pid_v2_cpp_0_80.txt"

MASSES = {
    "pi+": 0.13957039,
    "pi-": 0.13957039,
    "pion": 0.13957039,
    "K+": 0.493677,
    "K-": 0.493677,
    "K0S": 0.497611,
    "p": 0.93827208816,
    "pbar": 0.93827208816,
    "proton": 0.93827208816,
    "phi": 1.019461,
    "Lambda": 1.115683,
    "Lambdabar": 1.115683,
    "Xi-": 1.32171,
    "Xibar+": 1.32171,
    "Omega-": 1.67245,
    "Omegabar+": 1.67245,
}

NQ = {
    "pi+": 2,
    "pi-": 2,
    "pion": 2,
    "K+": 2,
    "K-": 2,
    "K0S": 2,
    "phi": 2,
    "p": 3,
    "pbar": 3,
    "proton": 3,
    "Lambda": 3,
    "Lambdabar": 3,
    "Xi-": 3,
    "Xibar+": 3,
    "Omega-": 3,
    "Omegabar+": 3,
}

PARTICLE_LABEL = {
    "pi+": r"$\pi^+$",
    "pi-": r"$\pi^-$",
    "pion": r"$\pi$",
    "K+": r"$K^+$",
    "K-": r"$K^-$",
    "K0S": r"$K^0_S$",
    "p": r"$p$",
    "pbar": r"$\bar{p}$",
    "proton": r"$p$",
    "phi": r"$\phi$",
    "Lambda": r"$\Lambda$",
    "Lambdabar": r"$\bar{\Lambda}$",
    "Xi-": r"$\Xi^-$",
    "Xibar+": r"$\bar{\Xi}^+$",
    "Omega-": r"$\Omega^-$",
    "Omegabar+": r"$\bar{\Omega}^+$",
}

STYLE = {
    "pi+": ("#d62728", "o"),
    "pi-": ("#ff9896", "o"),
    "pion": ("#d62728", "o"),
    "K+": ("#1f77b4", "s"),
    "K-": ("#6baed6", "s"),
    "K0S": ("#17becf", "D"),
    "p": ("#111111", "^"),
    "pbar": ("#777777", "^"),
    "proton": ("#111111", "^"),
    "phi": ("#9467bd", "P"),
    "Lambda": ("#2ca02c", "v"),
    "Lambdabar": ("#98df8a", "v"),
    "Xi-": ("#ff7f0e", "X"),
    "Xibar+": ("#f5a142", "X"),
    "Omega-": ("#8c564b", "*"),
    "Omegabar+": ("#c49c94", "*"),
}

OPEN_MARKERS = {"pi-", "K-", "pbar", "Lambdabar", "Xibar+", "Omegabar+"}


@dataclass(frozen=True)
class Point:
    energy: float
    particle: str
    centrality: str
    pt: float
    v2: float
    stat: float
    syst: float
    glob: float
    nq: int
    mass: float
    mt_minus_m0_over_nq: float
    v2_over_nq: float
    stat_over_nq: float
    syst_over_nq: float
    glob_over_nq: float
    source: str
    source_doi: str
    table: str
    table_doi: str
    extraction: str


def fnum(text: str | None) -> float | None:
    if text is None:
        return None
    text = text.strip().replace("\u2212", "-")
    if not text:
        return None
    text = text.replace("+/-", "")
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    return float(match.group(0)) if match else None


def clean_particle(raw: str) -> str | None:
    s = raw.replace(" ", "")
    s = s.replace("{", "").replace("}", "")
    s = s.replace("$", "")
    s = s.replace("\\mathrm", "")
    compact = s.replace("\\", "")
    if "Omega" in s or "Omega" in raw:
        return "Omega-"
    if "bar" in compact and "Xi" in s:
        return "Xibar+"
    if "Xi" in s:
        return "Xi-"
    if "bar" in compact and "Lambda" in s:
        return "Lambdabar"
    if "Lambda" in s:
        return "Lambda"
    if "barp" in compact or "overlinep" in compact:
        return "pbar"
    if compact in {"p", "p3GeV", "p54GeV"} or re.search(r"(?:^|/n_q)p(?:3|27|54|\\d|GeV)", compact):
        return "p"
    if "K^0" in raw or "K0" in s or "K_s" in s or "K^{0}_{s}" in raw or "Kshort" in raw:
        return "K0S"
    if "phi" in compact:
        return "phi"
    if "pi^+" in compact or "pi+" in compact:
        return "pi+"
    if "pi^-" in compact or "pi-" in compact:
        return "pi-"
    if "Pion" in raw or "pions" in raw:
        return "pion"
    if "K^+" in raw or "K{^+}" in raw or "K+" in s:
        return "K+"
    if "K^-" in raw or "K{^-}" in raw or "K-" in s:
        return "K-"
    if "Proton" in raw:
        return "proton"
    return None


def kinetic_from_pt(pt: float, mass: float, nq: int) -> float:
    return (math.sqrt(pt * pt + mass * mass) - mass) / nq


def pt_from_scaled_kinetic(x: float, mass: float, nq: int) -> float:
    mt = nq * x + mass
    return math.sqrt(max(mt * mt - mass * mass, 0.0))


def make_point(
    *,
    energy: float,
    particle: str,
    centrality: str,
    pt: float,
    v2: float,
    stat: float = 0.0,
    syst: float = 0.0,
    glob: float = 0.0,
    source: str,
    source_doi: str,
    table: str,
    table_doi: str,
    extraction: str,
) -> Point:
    nq = NQ[particle]
    mass = MASSES[particle]
    return Point(
        energy=energy,
        particle=particle,
        centrality=centrality,
        pt=pt,
        v2=v2,
        stat=abs(stat),
        syst=abs(syst),
        glob=abs(glob),
        nq=nq,
        mass=mass,
        mt_minus_m0_over_nq=kinetic_from_pt(pt, mass, nq),
        v2_over_nq=v2 / nq,
        stat_over_nq=abs(stat) / nq,
        syst_over_nq=abs(syst) / nq,
        glob_over_nq=abs(glob) / nq,
        source=source,
        source_doi=source_doi,
        table=table,
        table_doi=table_doi,
        extraction=extraction,
    )


def point_from_scaled(
    *,
    energy: float,
    particle: str,
    centrality: str,
    x: float,
    y: float,
    stat: float,
    syst: float,
    source: str,
    source_doi: str,
    table: str,
    table_doi: str,
) -> Point:
    nq = NQ[particle]
    mass = MASSES[particle]
    pt = pt_from_scaled_kinetic(x, mass, nq)
    return Point(
        energy=energy,
        particle=particle,
        centrality=centrality,
        pt=pt,
        v2=y * nq,
        stat=abs(stat) * nq,
        syst=abs(syst) * nq,
        glob=0.0,
        nq=nq,
        mass=mass,
        mt_minus_m0_over_nq=x,
        v2_over_nq=y,
        stat_over_nq=abs(stat),
        syst_over_nq=abs(syst),
        glob_over_nq=0.0,
        source=source,
        source_doi=source_doi,
        table=table,
        table_doi=table_doi,
        extraction="direct NCQ-scaled HEPData; pT and v2 reconstructed from mass and nq",
    )


def read_hepdata_meta_and_rows(path: Path) -> tuple[dict[str, str], list[str], list[list[str]]]:
    meta: dict[str, str] = {}
    header: list[str] = []
    rows: list[list[str]] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#:"):
                body = line[2:].strip()
                if ":" in body:
                    key, value = body.split(":", 1)
                    meta[key.strip()] = value.strip()
                if body.startswith("SQRT(S)"):
                    parts = next(csv.reader([body]))
                    if len(parts) > 1:
                        meta["energy"] = parts[1]
                elif body.startswith("Au+Au [GEV]"):
                    parts = next(csv.reader([body]))
                    if len(parts) > 1:
                        meta["energy"] = parts[1]
                elif body.startswith("Au+Au,"):
                    parts = next(csv.reader([body]))
                    if len(parts) > 1:
                        meta["particle"] = parts[1]
                elif body.startswith("$pt$,"):
                    parts = next(csv.reader([body]))
                    if len(parts) > 1:
                        meta["particle"] = parts[1]
                elif body.startswith("Centrality"):
                    parts = next(csv.reader([body]))
                    if len(parts) > 1:
                        meta["centrality"] = parts[1]
                continue
            if not line.strip():
                continue
            if not header:
                header = next(csv.reader([line]))
            else:
                rows.append(next(csv.reader([line])))
    return meta, header, rows


def parse_star_bes() -> list[Point]:
    points: list[Point] = []
    source = "STAR PRC 88, 014902 (2013), HEPData"
    source_doi = "10.17182/hepdata.102408"
    for path in sorted(RAWBES.glob("Table*.csv"), key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]):
        text = path.read_text(encoding="utf-8", errors="replace")
        if text.count("#: SQRT(S)") + text.count("#: Au+Au [GEV]") > 1:
            continue
        meta, header, rows = read_hepdata_meta_and_rows(path)
        if not header or len(header) < 2:
            continue
        if "p_T" not in header[0] and "p_T" not in header[0].replace("{", ""):
            continue
        if "v_2" not in header[1]:
            continue
        if "Delta" in header[1] or "\\Delta" in header[1]:
            continue
        if "Reduced" in header[0] or "m_T" in header[0]:
            continue
        if meta.get("centrality") != "0-80":
            continue

        energy = fnum(meta.get("energy") or meta.get("keyword cmenergies"))
        particle = clean_particle(meta.get("particle", ""))
        if energy is None or particle not in MASSES:
            continue

        col = {name: i for i, name in enumerate(header)}
        table_doi = meta.get("table_doi", "")
        table_name = meta.get("name", path.stem)
        for row in rows:
            if len(row) < 2:
                continue
            pt = fnum(row[0])
            v2 = fnum(row[1])
            if pt is None or v2 is None:
                continue
            points.append(
                make_point(
                    energy=energy,
                    particle=particle,
                    centrality="0-80%",
                    pt=pt,
                    v2=v2,
                    stat=fnum(row[col["stat +"]]) or 0.0 if "stat +" in col and col["stat +"] < len(row) else 0.0,
                    syst=fnum(row[col["sys +"]]) or 0.0 if "sys +" in col and col["sys +"] < len(row) else 0.0,
                    glob=fnum(row[col["glob +"]]) or 0.0 if "glob +" in col and col["glob +"] < len(row) else 0.0,
                    source=source,
                    source_doi=source_doi,
                    table=table_name,
                    table_doi=table_doi,
                    extraction="direct HEPData v2(pT); NCQ transform computed locally",
                )
            )
    return dedupe(points)


def parse_cpp_array(line: str) -> list[float] | None:
    match = re.search(r"=\s*\{([^}]*)\}", line)
    if not match:
        return None
    values: list[float] = []
    for item in match.group(1).split(","):
        value = fnum(item)
        if value is not None:
            values.append(value)
    return values


def flush_star241_block(
    points: list[Point],
    *,
    energy: float | None,
    raw_particle: str | None,
    arrays: dict[str, list[float]],
) -> None:
    particle_map = {
        "Phi": "phi",
        "XiM": "Xi-",
        "XiP": "Xibar+",
        "OmegaM": "Omega-",
        "OmegaP": "Omegabar+",
        "Lambda": "Lambda",
        "antiLambda": "Lambdabar",
        "K0S": "K0S",
        "Proton": "proton",
        "antiProton": "pbar",
        "PiM": "pi-",
        "PiP": "pi+",
        "KM": "K-",
        "KP": "K+",
    }
    particle = particle_map.get(raw_particle or "")
    pt = arrays.get("pt_bin_center", [])
    v2 = arrays.get("v2_values", [])
    stat = arrays.get("v2_stat_error", [])
    syst_low = arrays.get("v2_syst_low_error", [])
    syst_high = arrays.get("v2_syst_high_error", [])
    if energy is None or particle not in MASSES or not pt or len(pt) != len(v2):
        return

    source = "STAR PRC 93, 014907 (2016), STAR public data"
    source_doi = "10.1103/PhysRevC.93.014907"
    for i, (pt_value, v2_value) in enumerate(zip(pt, v2)):
        low = abs(syst_low[i]) if i < len(syst_low) else 0.0
        high = abs(syst_high[i]) if i < len(syst_high) else 0.0
        points.append(
            make_point(
                energy=energy,
                particle=particle,
                centrality="0-80%",
                pt=pt_value,
                v2=v2_value,
                stat=stat[i] if i < len(stat) else 0.0,
                syst=max(low, high),
                source=source,
                source_doi=source_doi,
                table="STAR public data 0-80% table",
                table_doi="",
                extraction="STAR public C++ v2(pT) table; NCQ transform computed locally",
            )
        )


def parse_star_pub241() -> list[Point]:
    points: list[Point] = []
    current_energy: float | None = None
    current_particle: str | None = None
    arrays: dict[str, list[float]] = {}

    for line in STAR241.read_text(encoding="utf-8", errors="replace").splitlines():
        energy_match = re.search(r"beam energy:\s*([0-9.]+)\s*GeV", line)
        if energy_match:
            flush_star241_block(points, energy=current_energy, raw_particle=current_particle, arrays=arrays)
            current_energy = float(energy_match.group(1))
            current_particle = None
            arrays = {}
            continue

        particle_match = re.search(r"Particle species:\s*([A-Za-z0-9]+)", line)
        if particle_match:
            flush_star241_block(points, energy=current_energy, raw_particle=current_particle, arrays=arrays)
            current_particle = particle_match.group(1)
            arrays = {}
            continue

        array_match = re.search(r"Double_t\s+([A-Za-z0-9_]+)\[\d+\]", line)
        if array_match:
            values = parse_cpp_array(line)
            if values is not None:
                arrays[array_match.group(1)] = values

    flush_star241_block(points, energy=current_energy, raw_particle=current_particle, arrays=arrays)
    return dedupe(points)


def low_energy_particle_from_path(path: Path) -> str | None:
    name = path.stem
    if "pi^{+}" in name:
        return "pi+"
    if "pi^{-}" in name:
        return "pi-"
    if "K^{+}" in name:
        return "K+"
    if "K^{-}" in name:
        return "K-"
    if "K^{0}_{S}" in name:
        return "K0S"
    if name.endswith("_p"):
        return "p"
    if name.endswith("_Lambda"):
        return "Lambda"
    return None


def parse_star_lowenergy_prl() -> list[Point]:
    points: list[Point] = []
    source = "STAR PRL 135, 072301 (2025), HEPData"
    source_doi = "10.17182/hepdata.159489"
    for path in sorted(RAWLOW.glob("Fig2_*GeV_*.csv")):
        meta, header, rows = read_hepdata_meta_and_rows(path)
        if not header or len(header) < 2:
            continue
        if "p_{T}" not in header[0]:
            continue
        if "v_{2}" not in header[1]:
            continue
        energy = fnum(meta.get("keyword cmenergies"))
        particle = low_energy_particle_from_path(path)
        if energy is None or particle not in MASSES:
            continue

        col = {name: i for i, name in enumerate(header)}
        table_doi = meta.get("table_doi", "")
        table_name = meta.get("name", path.stem)
        for row in rows:
            if len(row) < 2:
                continue
            pt = fnum(row[0])
            v2 = fnum(row[1])
            if pt is None or v2 is None:
                continue
            points.append(
                make_point(
                    energy=energy,
                    particle=particle,
                    centrality="10-40%",
                    pt=pt,
                    v2=v2,
                    stat=fnum(row[col["stat +"]]) or 0.0 if "stat +" in col and col["stat +"] < len(row) else 0.0,
                    syst=fnum(row[col["sys +"]]) or 0.0 if "sys +" in col and col["sys +"] < len(row) else 0.0,
                    source=source,
                    source_doi=source_doi,
                    table=table_name,
                    table_doi=table_doi,
                    extraction="direct HEPData v2(pT); NCQ transform computed locally",
                )
            )
    return dedupe(points)


def parse_scaled_star3_and_54(*, include_3gev: bool = True) -> list[Point]:
    points: list[Point] = []
    source = "STAR PLB 827, 137003 (2022), HEPData"
    source_doi = "10.17182/hepdata.110656"
    for path in sorted(RAW3.glob("fig4_*GeV*.csv")):
        name = path.name
        if "27GeV" in name:
            continue
        meta, header, rows = read_hepdata_meta_and_rows(path)
        if not header or len(header) < 2:
            continue
        match_energy = re.search(r"(\d+(?:\.\d+)?)\s*GeV", header[1])
        if not match_energy:
            match_energy = re.search(r"(\d+(?:\.\d+)?)GeV", path.stem)
        if not match_energy:
            continue
        energy = float(match_energy.group(1))
        if abs(energy - 54.0) < 0.05:
            energy = 54.4
        if energy not in {3.0, 54.4}:
            continue
        if energy == 3.0 and not include_3gev:
            continue
        particle = clean_particle(header[1])
        if particle not in MASSES:
            continue
        table_doi = meta.get("table_doi", "")
        table_name = meta.get("name", path.stem)
        for row in rows:
            if len(row) < 6:
                continue
            x = fnum(row[0])
            y = fnum(row[1])
            if x is None or y is None:
                continue
            points.append(
                point_from_scaled(
                    energy=energy,
                    particle=particle,
                    centrality="10-40%",
                    x=x,
                    y=y,
                    stat=fnum(row[2]) or 0.0,
                    syst=fnum(row[4]) or 0.0,
                    source=source,
                    source_doi=source_doi,
                    table=table_name,
                    table_doi=table_doi,
                )
            )
    return dedupe(points)


def parse_star_4p5() -> list[Point]:
    points: list[Point] = []
    path = RAW45 / "Figure13.csv"
    source = "STAR PRC 103, 034908 (2021), HEPData"
    source_doi = "10.17182/hepdata.95903"
    table_doi = ""
    table_name = "Figure 13"
    dataset = ""
    header: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                header = []
                continue
            if line.startswith("#: table_doi:"):
                table_doi = line.split(":", 2)[2].strip()
                continue
            if line.startswith("#: Dataset,"):
                dataset = line.split(",", 1)[1].strip()
                header = []
                continue
            if line.startswith("#:"):
                continue
            if not header:
                header = next(csv.reader([line]))
                continue
            if not dataset.startswith("STAR Au+Au"):
                continue
            row = next(csv.reader([line]))
            if len(row) < 6:
                continue
            x = fnum(row[0])
            y = fnum(row[1])
            if x is None or y is None:
                continue
            if "pions" in dataset:
                particle = "pion"
            elif re.search(r"\bp$", dataset):
                particle = "p"
            else:
                continue
            points.append(
                point_from_scaled(
                    energy=4.5,
                    particle=particle,
                    centrality="0-30%",
                    x=x,
                    y=y,
                    stat=fnum(row[2]) or 0.0,
                    syst=fnum(row[4]) or 0.0,
                    source=source,
                    source_doi=source_doi,
                    table=table_name,
                    table_doi=table_doi,
                )
            )
    return dedupe(points)


class TableParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_h3 = False
        self.current_heading = ""
        self._h3_parts: list[str] = []
        self.in_table = False
        self.in_cell = False
        self._cell_parts: list[str] = []
        self._row: list[str] = []
        self._rows: list[list[str]] = []
        self.tables: list[tuple[str, list[list[str]]]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "h3":
            self.in_h3 = True
            self._h3_parts = []
        elif tag == "table":
            self.in_table = True
            self._rows = []
        elif self.in_table and tag == "tr":
            self._row = []
        elif self.in_table and tag in {"td", "th"}:
            self.in_cell = True
            self._cell_parts = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "h3":
            self.in_h3 = False
            self.current_heading = " ".join("".join(self._h3_parts).split())
        elif tag == "table" and self.in_table:
            self.in_table = False
            self.tables.append((self.current_heading, self._rows))
        elif self.in_table and tag == "tr":
            if self._row:
                self._rows.append(self._row)
        elif self.in_table and tag in {"td", "th"}:
            self.in_cell = False
            self._row.append(" ".join("".join(self._cell_parts).split()))

    def handle_data(self, data: str) -> None:
        if self.in_h3:
            self._h3_parts.append(data)
        elif self.in_cell:
            self._cell_parts.append(data)


def parse_star200() -> list[Point]:
    parser = TableParser()
    parser.feed(STAR200.read_text(encoding="iso-8859-1", errors="replace"))
    points: list[Point] = []
    source = "STAR PRC 77, 054901 (2008), STAR public data"
    source_doi = "10.1103/PhysRevC.77.054901"

    target_specs = {
        "Kshort 0-80%": ("K0S", 1, 2, 3),
        "Lambda 0-80%": ("Lambda", 1, 2, 3),
        "Cascade 0-80%": ("Xi-", 1, 2, None),
        "Omega 0-80%": ("Omega-", 1, 2, None),
        "Proton 0-80%": ("proton", 1, 2, None),
        "Pion 0-80%": ("pion", 1, 2, None),
    }

    for heading, rows in parser.tables:
        if "Figure 7" not in heading or not rows:
            continue
        header = " ".join(rows[0])
        spec = None
        for key, value in target_specs.items():
            if key in header:
                spec = value
                break
        if spec is None:
            continue
        particle, value_col, stat_col, syst_col = spec
        for row in rows[1:]:
            if len(row) <= value_col:
                continue
            pt = fnum(row[0])
            v2 = fnum(row[value_col])
            if pt is None or v2 is None:
                continue
            stat = fnum(row[stat_col]) or 0.0 if len(row) > stat_col else 0.0
            syst = fnum(row[syst_col]) or 0.0 if syst_col is not None and len(row) > syst_col else 0.0
            points.append(
                make_point(
                    energy=200.0,
                    particle=particle,
                    centrality="0-80%",
                    pt=pt,
                    v2=v2,
                    stat=stat,
                    syst=syst,
                    source=source,
                    source_doi=source_doi,
                    table="STAR public data Figure 7",
                    table_doi="",
                    extraction="STAR public v2(pT) table; NCQ transform computed locally",
                )
            )
    return dedupe(points)


def dedupe(points: list[Point]) -> list[Point]:
    by_key: dict[tuple[float, str, str, float], Point] = {}
    for p in points:
        key = (round(p.energy, 3), p.particle, p.centrality, round(p.pt, 6))
        by_key.setdefault(key, p)
    return sorted(by_key.values(), key=lambda p: (p.energy, p.particle, p.pt))


def write_outputs(points: list[Point]) -> None:
    DATA.mkdir(exist_ok=True)
    raw_path = DATA / "star_identified_v2_pt.csv"
    ncq_path = DATA / "star_identified_v2_ncq.csv"
    fields = list(Point.__dataclass_fields__)
    with raw_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for p in points:
            writer.writerow(p.__dict__)

    ncq_fields = [
        "energy",
        "particle",
        "centrality",
        "nq",
        "mass",
        "pt",
        "v2",
        "mt_minus_m0_over_nq",
        "v2_over_nq",
        "stat_over_nq",
        "syst_over_nq",
        "glob_over_nq",
        "source",
        "source_doi",
        "table",
        "table_doi",
    ]
    with ncq_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ncq_fields)
        writer.writeheader()
        for p in points:
            writer.writerow({field: getattr(p, field) for field in ncq_fields})


def yerr_total(p: Point) -> float:
    return math.sqrt(p.stat_over_nq**2 + p.syst_over_nq**2 + p.glob_over_nq**2)


def set_value_ylim(ax: plt.Axes, points: list[Point]) -> None:
    yvals = [p.v2_over_nq for p in points]
    ymin = min(yvals)
    ymax = max(yvals)
    pad = max(0.006, 0.14 * (ymax - ymin if ymax > ymin else 0.02))
    ymin -= pad
    ymax += pad
    if ymin > 0.0:
        ymin = -0.004
    if ymax < 0.0:
        ymax = 0.004
    ax.set_ylim(ymin, ymax)


def energy_label(energy: float) -> str:
    return f"{energy:g}".replace(".", "p")


def plot_energy(points: list[Point], energy: float) -> None:
    subset = [p for p in points if abs(p.energy - energy) < 1e-6]
    if not subset:
        return
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    order = [
        "pi+",
        "pi-",
        "pion",
        "K+",
        "K-",
        "K0S",
        "p",
        "pbar",
        "proton",
        "phi",
        "Lambda",
        "Lambdabar",
        "Xi-",
        "Xibar+",
        "Omega-",
        "Omegabar+",
    ]
    seen_labels: set[str] = set()
    for particle in order:
        rows = sorted([p for p in subset if p.particle == particle], key=lambda p: p.mt_minus_m0_over_nq)
        if not rows:
            continue
        color, marker = STYLE[particle]
        label = PARTICLE_LABEL[particle]
        ax.errorbar(
            [p.mt_minus_m0_over_nq for p in rows],
            [p.v2_over_nq for p in rows],
            yerr=[yerr_total(p) for p in rows],
            marker=marker,
            linestyle="none",
            color=color,
            markerfacecolor="white" if particle in OPEN_MARKERS else color,
            markeredgecolor=color,
            elinewidth=0.9,
            capsize=2.0,
            markersize=5.3 if particle != "Omega-" else 7.0,
            label=label if label not in seen_labels else None,
        )
        seen_labels.add(label)

    centralities = sorted({p.centrality for p in subset})
    ax.axhline(0.0, color="#555555", linewidth=0.9)
    ax.grid(True, color="#dddddd", linewidth=0.8)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.set_xlabel(r"$(m_T-m_0)/n_q$ (GeV/$c^2$)")
    ax.set_ylabel(r"$v_2/n_q$")
    ax.set_title(rf"STAR Au+Au $\sqrt{{s_{{NN}}}}={energy:g}$ GeV")

    set_value_ylim(ax, subset)
    ax.set_xlim(left=0.0)
    ax.text(
        0.03,
        0.94,
        ", ".join(centralities),
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        color="#555555",
    )
    ax.legend(frameon=False, fontsize=8.2, ncol=2, loc="lower right")
    fig.tight_layout()
    stem = f"star_identified_v2_ncq_{energy_label(energy)}GeV"
    fig.savefig(FIGURES / f"{stem}.png", dpi=300)
    fig.savefig(FIGURES / f"{stem}.pdf")
    plt.close(fig)


def plot_all_grid(points: list[Point]) -> None:
    energies = sorted({p.energy for p in points})
    if not energies:
        return
    FIGURES.mkdir(parents=True, exist_ok=True)
    ncols = 3
    nrows = math.ceil(len(energies) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.2, 3.8 * nrows), sharex=False, sharey=False)
    flat_axes = list(axes.flat if hasattr(axes, "flat") else [axes])
    order = [
        "pi+",
        "pi-",
        "pion",
        "K+",
        "K-",
        "K0S",
        "p",
        "pbar",
        "proton",
        "phi",
        "Lambda",
        "Lambdabar",
        "Xi-",
        "Xibar+",
        "Omega-",
        "Omegabar+",
    ]
    for ax, energy in zip(flat_axes, energies):
        subset = [p for p in points if abs(p.energy - energy) < 1e-6]
        for particle in order:
            rows = sorted([p for p in subset if p.particle == particle], key=lambda p: p.mt_minus_m0_over_nq)
            if not rows:
                continue
            color, marker = STYLE[particle]
            ax.errorbar(
                [p.mt_minus_m0_over_nq for p in rows],
                [p.v2_over_nq for p in rows],
                yerr=[yerr_total(p) for p in rows],
                marker=marker,
                linestyle="none",
                color=color,
                markerfacecolor="white" if particle in OPEN_MARKERS else color,
                markeredgecolor=color,
                elinewidth=0.7,
                capsize=1.5,
                markersize=4.5,
                label=PARTICLE_LABEL[particle],
            )
        ax.axhline(0.0, color="#555555", linewidth=0.8)
        ax.grid(True, color="#dddddd", linewidth=0.7)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.set_title(rf"$\sqrt{{s_{{NN}}}}={energy:g}$ GeV", fontsize=11)
        ax.set_xlabel(r"$(m_T-m_0)/n_q$")
        ax.set_ylabel(r"$v_2/n_q$")
        set_value_ylim(ax, subset)
    for ax in flat_axes[len(energies):]:
        ax.axis("off")
    handles = []
    labels = []
    present = {p.particle for p in points}
    for particle in order:
        if particle not in present:
            continue
        color, marker = STYLE[particle]
        handle = plt.Line2D(
            [],
            [],
            marker=marker,
            linestyle="none",
            color=color,
            markerfacecolor="white" if particle in OPEN_MARKERS else color,
            markeredgecolor=color,
            markersize=5.0 if particle not in {"Omega-", "Omegabar+"} else 6.6,
        )
        handles.append(handle)
        labels.append(PARTICLE_LABEL[particle])
    fig.legend(handles, labels, frameon=False, ncol=8, loc="upper center", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURES / "star_identified_v2_ncq_all_energies.png", dpi=300)
    fig.savefig(FIGURES / "star_identified_v2_ncq_all_energies.pdf")
    plt.close(fig)


def main() -> None:
    points = []
    if RAWLOW.exists():
        points.extend(parse_star_lowenergy_prl())
        points.extend(parse_scaled_star3_and_54(include_3gev=False))
    else:
        points.extend(parse_scaled_star3_and_54())
        points.extend(parse_star_4p5())
    points.extend(parse_star_pub241() if STAR241.exists() else parse_star_bes())
    points.extend(parse_star200())
    points = dedupe(points)
    write_outputs(points)
    for energy in sorted({p.energy for p in points}):
        plot_energy(points, energy)
    plot_all_grid(points)
    print(f"Wrote {len(points)} points to {DATA / 'star_identified_v2_ncq.csv'}")
    print(f"Wrote figures to {FIGURES}")


if __name__ == "__main__":
    main()
