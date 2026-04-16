#!/usr/bin/env python3
"""Build LaTeX/PDF summary table for thermal-fit points and references."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


IN_CSV = Path("data/literature_thermal_fit_points_selected_avg.csv")
IN_RAW_CSV = Path("data/literature_thermal_fit_points_selected.csv")
OUT_TEX = Path("data/thermal_fit_points_summary.tex")
OUT_PDF = Path("data/thermal_fit_points_summary.pdf")


def fmt_num(v: str) -> str:
    try:
        return f"{int(round(float(v)))}"
    except ValueError:
        return v


def pm(v: str, e: str) -> str:
    return f"${fmt_num(v)} \\pm {fmt_num(e)}$"


def tex_escape(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("_", "\\_")
        .replace("#", "\\#")
        .replace("{", "\\{")
        .replace("}", "\\}")
    )


def main() -> None:
    rows: list[dict[str, str]] = []
    with IN_CSV.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    raw_rows: list[dict[str, str]] = []
    with IN_RAW_CSV.open(newline="", encoding="utf-8") as f:
        raw_rows = list(csv.DictReader(f))

    star_orig: dict[float, dict[str, str]] = {}
    star_fit: dict[float, dict[str, str]] = {}
    non_star: list[dict[str, str]] = []

    for r in rows:
        g = r["experiment_group"]
        e = float(r["energy_GeV"])
        if g == "STAR (this work fit)":
            star_fit[e] = r
        elif g.startswith("STAR (orig"):
            star_orig[e] = r
        else:
            non_star.append(r)

    out_rows: list[dict[str, object]] = []
    for r in non_star:
        g = r["experiment_group"]
        # AGS/SPS are handled below from non-averaged (raw) points.
        if g in {"AGS", "SPS"}:
            continue
        if g == "ALICE":
            refs = ["andronic_2018"]
            exp_name = "ALICE"
        elif g in {"SIS", "AGS"}:
            refs = ["cleymans_2006"]
            exp_name = "FOPI" if g == "SIS" else "E895/E866/E917"
        elif g == "SPS":
            refs = ["becattini_2006", "bravina_2002"]
            exp_name = "NA49"
        else:
            refs = ["cleymans_2006"]
            exp_name = g
        out_rows.append(
            {
                "energy": f"{float(r['energy_GeV']):g}",
                "t_orig": pm(r["T_MeV"], r["T_err_plus_MeV"]),
                "mb_orig": pm(r["muB_MeV"], r["muB_err_plus_MeV"]),
                "t_fit": "--",
                "mb_fit": "--",
                "exp": exp_name,
                "refs": refs,
            }
        )

    # Use all raw AGS/SPS points (no averaging), as requested.
    for r in raw_rows:
        g = r["experiment_group"]
        if g not in {"AGS", "SPS"}:
            continue
        exp_name = "E895/E866/E917" if g == "AGS" else "NA49"
        out_rows.append(
            {
                "energy": f"{float(r['energy_GeV']):g}",
                "t_orig": pm(r["T_MeV"], r["T_err_plus_MeV"]),
                "mb_orig": pm(r["muB_MeV"], r["muB_err_plus_MeV"]),
                "t_fit": "--",
                "mb_fit": "--",
                "exp": exp_name,
                "refs": ["cleymans_2006"],
            }
        )

    for e in sorted(set(star_orig.keys()) | set(star_fit.keys())):
        o = star_orig.get(e)
        f = star_fit.get(e)
        t_orig = pm(o["T_MeV"], o["T_err_plus_MeV"]) if o else "--"
        mb_orig = pm(o["muB_MeV"], o["muB_err_plus_MeV"]) if o else "--"
        t_fit = pm(f["T_MeV"], f["T_err_plus_MeV"]) if f else "--"
        mb_fit = pm(f["muB_MeV"], f["muB_err_plus_MeV"]) if f else "--"
        # Split R2 refs by source category as requested.
        if e <= 39.0:
            refs = ["adamczyk_2017_prc96", "adam_2020_prc102"]
        elif abs(e - 62.4) < 1e-9:
            refs = ["abelev_2009_prc79", "aggarwal_2011_prc83", "aggarwal_2011_prc83"]
        else:
            refs = ["abelev_2009_prc79", "abelev_2007_prl98", "aggarwal_2011_prc83"]

        out_rows.append(
            {
                "energy": f"{e:g}",
                "t_orig": t_orig,
                "mb_orig": mb_orig,
                "t_fit": t_fit,
                "mb_fit": mb_fit,
                "exp": "STAR",
                "refs": refs,
            }
        )

    out_rows.sort(key=lambda x: float(x["energy"]))

    # Add HADES Ar+KCl thermal-fit point at sqrt(sNN)=2.61 GeV.
    out_rows.append(
        {
            "energy": "2.61",
            "t_orig": "$70 \\pm 3$",
            "mb_orig": "$748 \\pm 8$",
            "t_fit": "--",
            "mb_fit": "--",
            "exp": "HADES (Ar+KCl)",
            "refs": ["hades_arkcl_2016"],
        }
    )

    # Add ALICE 5.02 TeV row using Andronic LHC freeze-out value.
    out_rows.append(
        {
            "energy": "5020",
            "t_orig": "$157 \\pm 2$",
            "mb_orig": "$1 \\pm 4$",
            "t_fit": "--",
            "mb_fit": "--",
            "exp": "ALICE",
            "refs": ["andronic_2018"],
        }
    )
    out_rows.sort(key=lambda x: float(x["energy"]))

    ref_text = {
        "andronic_2018": "A. Andronic, P. Braun-Munzinger, K. Redlich, and J. Stachel, Nature \\textbf{561}, 321 (2018), doi:10.1038/s41586-018-0491-6. \\textit{Comment: ALICE/LHC chemical freeze-out reference; 5.02 TeV row uses the LHC value from this work (no separate 5.02 best-fit table in this paper).}",
        "abelev_2009_prc79": "B. I. Abelev \\textit{et al.} (STAR Collaboration), Phys. Rev. C \\textbf{79}, 034909 (2009), doi:10.1103/PhysRevC.79.034909. \\textit{Comment: identified $\\pi^\\pm,K^\\pm,p,\\bar{p}$ yields used for STAR 62.4 and 200 GeV fits (Table VIII / arXiv:0808.2041).}",
        "abelev_2007_prl98": "B. I. Abelev \\textit{et al.} (STAR Collaboration), Phys. Rev. Lett. \\textbf{98}, 062301 (2007), doi:10.1103/PhysRevLett.98.062301. \\textit{Comment: $\\Lambda,\\bar{\\Lambda},\\Xi,\\bar{\\Xi}$ yields at 200 GeV used in this-work STAR fit inputs (Table I / arXiv:nucl-ex/0606014).}",
        "adamczyk_2017_prc96": "L. Adamczyk \\textit{et al.} (STAR Collaboration), Phys. Rev. C \\textbf{96}, 044904 (2017), doi:10.1103/PhysRevC.96.044904. \\textit{Comment: STAR original $(T,\\mu_B)$ points (GCEY, 0--5\\%) at 7.7--39 GeV; also $\\pi^\\pm,K^\\pm,p,\\bar{p}$ yield inputs via HEPData ins1510593.}",
        "adam_2020_prc102": "J. Adam \\textit{et al.} (STAR Collaboration), Phys. Rev. C \\textbf{102}, 034909 (2020), doi:10.1103/PhysRevC.102.034909. \\textit{Comment: 7.7--39 GeV $\\Lambda,\\bar{\\Lambda},\\Xi,\\bar{\\Xi}$ yield inputs via HEPData ins1738953.}",
        "aggarwal_2011_prc83": "M. M. Aggarwal \\textit{et al.} (STAR Collaboration), Phys. Rev. C \\textbf{83}, 024901 (2011); Erratum Phys. Rev. C \\textbf{107}, 049903 (2023), doi:10.1103/PhysRevC.83.024901. \\textit{Comment: STAR original 62.4 and 200 GeV points and 62.4 GeV hyperon-ratio inputs (HEPData ins871561).}",
        "cleymans_2006": "J. Cleymans, H. Oeschler, K. Redlich, and S. Wheaton, Phys. Rev. C \\textbf{73}, 034905 (2006), doi:10.1103/PhysRevC.73.034905. \\textit{Comment: Ref. 1; AGS/SPS entries are listed as individual (non-averaged) Table-I points in this table.}",
        "becattini_2006": "F. Becattini, J. Manninen, and M. Gazdzicki, Phys. Rev. C \\textbf{73}, 044905 (2006), doi:10.1103/PhysRevC.73.044905. \\textit{Comment: SPS thermal-fit set used in the non-conference averaged SPS points.}",
        "bravina_2002": "L. V. Bravina \\textit{et al.}, Phys. Rev. C \\textbf{66}, 014906 (2002), doi:10.1103/PhysRevC.66.014906. \\textit{Comment: SPS thermal-fit set used in the non-conference averaged SPS points.}",
        "hades_arkcl_2016": "G. Agakishiev \\textit{et al.} (HADES Collaboration), Eur. Phys. J. A \\textbf{52}, 178 (2016), doi:10.1140/epja/i2016-16178-x. \\textit{Comment: HADES Ar+KCl ($\\sqrt{s_{NN}}=2.61$ GeV) statistical-model reanalysis; table uses the commonly quoted Ar+KCl THERMUS fit point $T=70\\pm3$ MeV, $\\mu_B=748\\pm8$ MeV.}",
    }

    # Flat reference numbering in first-appearance order from table rows.
    ref_order: list[str] = []
    for r in out_rows:
        for k in r["refs"]:
            if k not in ref_order:
                ref_order.append(k)
    ref_num = {k: i + 1 for i, k in enumerate(ref_order)}

    body_lines = []
    for r in out_rows:
        row_refs: list[str] = []
        for k in r["refs"]:
            n = str(ref_num[k])
            if n not in row_refs:
                row_refs.append(n)
        body_lines.append(
            "  {energy} & {t_orig} & {mb_orig} & {t_fit} & {mb_fit} & {exp} & {ref} \\\\".format(
                energy=tex_escape(r["energy"]),
                t_orig=r["t_orig"],
                mb_orig=r["mb_orig"],
                t_fit=r["t_fit"],
                mb_fit=r["mb_fit"],
                exp=tex_escape(r["exp"]),
                ref=tex_escape(",".join(row_refs)),
            )
        )

    reduced_lines = []
    for r in out_rows:
        row_refs: list[str] = []
        for k in r["refs"]:
            n = str(ref_num[k])
            if n not in row_refs:
                row_refs.append(n)
        reduced_lines.append(
            "  {e} & {x} & {t} & {m} & {r} \\\\".format(
                e=tex_escape(r["energy"]),
                x=tex_escape(r["exp"]),
                t=(r["t_fit"] if (r["exp"] == "STAR" and r["t_fit"] != "--") else r["t_orig"]),
                m=(r["mb_fit"] if (r["exp"] == "STAR" and r["mb_fit"] != "--") else r["mb_orig"]),
                r=tex_escape(",".join(row_refs)),
            )
        )

    ref_items_latex = "\n".join(f"\\item [{i}] REFSLOT{i}" for i in range(1, len(ref_order) + 1))

    tex = r"""\documentclass[10pt]{article}
\usepackage[margin=0.8in,landscape]{geometry}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\setlength{\LTpre}{0pt}
\setlength{\LTpost}{0pt}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\begin{document}
\begin{center}
\Large Thermal-Fit Data Points
\end{center}
\vspace{0.1em}
\begin{center}
\normalsize
\textbf{Reduced Summary (STAR Uses This Work)}
\end{center}
\small
\begin{longtable}{c l c c c}
\toprule
Energy (GeV) & Experiment & T (MeV) & $\mu_B$ (MeV) & Ref \\
\midrule
\endfirsthead
\toprule
Energy (GeV) & Experiment & T (MeV) & $\mu_B$ (MeV) & Ref \\
\midrule
\endhead
""" + "\n".join(reduced_lines) + r"""
\bottomrule
\end{longtable}

\clearpage
\vspace{0.1em}
\begin{center}
\normalsize
\textbf{Full Summary}
\end{center}
\small
\begin{longtable}{c c c c c l l}
\toprule
Energy (GeV) & $T$ (original) & $\mu_B$ (original) & $T$ (this work, STAR) & $\mu_B$ (this work, STAR) & Experiment & Ref \\
\midrule
\endfirsthead
\toprule
Energy (GeV) & $T$ (original) & $\mu_B$ (original) & $T$ (this work, STAR) & $\mu_B$ (this work, STAR) & Experiment & Ref \\
\midrule
\endhead
""" + "\n".join(body_lines) + r"""
\bottomrule
\end{longtable}

\clearpage
\normalsize
\textbf{References}

\begin{itemize}
""" + ref_items_latex + r"""
\end{itemize}
\end{document}
"""

    # Replace placeholder reference slots with first-appearance numbered refs.
    for i, key in enumerate(ref_order, start=1):
        tex = tex.replace(f"REFSLOT{i}", ref_text[key], 1)

    OUT_TEX.write_text(tex, encoding="utf-8")

    # Build PDF from LaTeX.
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", OUT_TEX.name],
        cwd=str(OUT_TEX.parent),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    print(f"wrote {OUT_TEX}")
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
