#!/usr/bin/env python3
"""Add STAR 62.4 GeV 0-5% hyperon dN/dy points derived from Fig. 9 ratios."""

from __future__ import annotations

import csv
import math
from pathlib import Path

INPUT_CSV = Path("data/first_group_dn_dy_vs_energy.csv")
OUT_CROSSCHECK = Path("data/star_62p4_hyperon_consistency_check.csv")


def find_pi_minus_62p4(rows: list[dict[str, str]]) -> tuple[float, float]:
    candidates = [
        r
        for r in rows
        if r.get("particle") == "pi-"
        and r.get("energy_GeV") == "62.4"
        and "dn/dy" in (r.get("observable") or "").lower()
        and (r.get("centrality") or "").strip() == "0-5%"
    ]
    if not candidates:
        raise RuntimeError("Could not find 62.4 GeV 0-5% pi- row in data/first_group_dn_dy_vs_energy.csv")
    row = candidates[0]
    val = float(row["value"])
    stat = float(row["stat_err"]) if (row.get("stat_err") or "").strip() else 0.0
    sys = float(row["sys_err"]) if (row.get("sys_err") or "").strip() else 0.0
    return val, math.hypot(stat, sys)


def add_if_missing(rows: list[dict[str, str]], particle: str, value: float, err_total: float, note: str) -> bool:
    exists = any(
        r.get("particle") == particle
        and r.get("energy_GeV") == "62.4"
        and "dn/dy" in (r.get("observable") or "").lower()
        for r in rows
    )
    if exists:
        return False

    rows.append(
        {
            "paper_id": "ins871561",
            "particle": particle,
            "observable": "dN/dy",
            "energy_GeV": "62.4",
            "value": f"{value:.8g}",
            "stat_err": "",
            "sys_err": f"{err_total:.8g}",
            "centrality": "0-5%",
            "source": "ins871561/Figure9 (ratio-derived absolute yield)",
            "note": note,
        }
    )
    return True


def main() -> None:
    with INPUT_CSV.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    pi_val, pi_err = find_pi_minus_62p4(rows)

    # STAR Fig. 9 measured ratios for 62.4 GeV, 0-5% Au+Au:
    # Lambda_bar/Pi- = 0.0306 ± 0.0038
    # Lambda_bar/Lambda = 0.4656 ± 0.0172
    # Xi+/Pi- = 0.0048 ± 0.0007
    # Xi+/Xi- = 0.6539 ± 0.0444
    r_lb_pi, e_lb_pi = 0.0306, 0.0038
    r_lb_l, e_lb_l = 0.4656, 0.0172
    r_xb_pi, e_xb_pi = 0.0048, 0.0007
    r_xb_x, e_xb_x = 0.6539, 0.0444

    lbar = r_lb_pi * pi_val
    lbar_err = lbar * math.hypot(e_lb_pi / r_lb_pi, pi_err / pi_val)
    lamb = lbar / r_lb_l
    lamb_err = lamb * math.hypot(lbar_err / lbar, e_lb_l / r_lb_l)

    xibar = r_xb_pi * pi_val
    xibar_err = xibar * math.hypot(e_xb_pi / r_xb_pi, pi_err / pi_val)
    xi = xibar / r_xb_x
    xi_err = xi * math.hypot(xibar_err / xibar, e_xb_x / r_xb_x)

    added = 0
    added += int(
        add_if_missing(
            rows,
            "Lambda_bar",
            lbar,
            lbar_err,
            "Derived from STAR Fig.9 ratio Lambda_bar/Pi- (0-5%) and local pi- (0-5%); uncertainty propagated in quadrature.",
        )
    )
    added += int(
        add_if_missing(
            rows,
            "Lambda",
            lamb,
            lamb_err,
            "Derived via STAR Fig.9 ratios Lambda_bar/Pi- and Lambda_bar/Lambda (0-5%) with local pi- (0-5%); uncertainty propagated.",
        )
    )
    added += int(
        add_if_missing(
            rows,
            "Xi_bar",
            xibar,
            xibar_err,
            "Derived from STAR Fig.9 ratio Xi+/Pi- (0-5%) and local pi- (0-5%); uncertainty propagated in quadrature.",
        )
    )
    added += int(
        add_if_missing(
            rows,
            "Xi",
            xi,
            xi_err,
            "Derived via STAR Fig.9 ratios Xi+/Pi- and Xi+/Xi- (0-5%) with local pi- (0-5%); uncertainty propagated.",
        )
    )

    rows.sort(
        key=lambda r: (
            float(r.get("energy_GeV") or "0"),
            r.get("particle") or "",
            r.get("observable") or "",
        )
    )

    with INPUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Consistency cross-check against STAR Fig.6 energy-trend points
    # (not centrality-matched for Xi in that figure; included for transparency).
    fig6_reference = {
        "Lambda": 14.8,
        "Lambda_bar": 7.8,   # appears at x~63 GeV in Fig.6 up table
        "Xi": 1.63,          # Xi- in Fig.6 down
        "Xi_bar": 1.03,      # Xi+ appears at x~63 GeV in Fig.6 down
    }
    derived = {
        "Lambda": (lamb, lamb_err),
        "Lambda_bar": (lbar, lbar_err),
        "Xi": (xi, xi_err),
        "Xi_bar": (xibar, xibar_err),
    }
    with OUT_CROSSCHECK.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "energy_GeV",
                "particle",
                "centrality_used",
                "derived_value",
                "derived_total_err",
                "fig6_reference_value",
                "relative_diff_percent",
                "comment",
            ],
        )
        writer.writeheader()
        for particle in ["Lambda", "Lambda_bar", "Xi", "Xi_bar"]:
            val, err = derived[particle]
            ref = fig6_reference[particle]
            rel = 100.0 * (val - ref) / ref
            comment = "Fig.6 Xi values are 0-20% in caption; use with care for 0-5% fits." if particle.startswith("Xi") else "Fig.6 top-panel STAR points are quoted as 0-5%."
            writer.writerow(
                {
                    "energy_GeV": "62.4",
                    "particle": particle,
                    "centrality_used": "0-5% (ratio-derived)",
                    "derived_value": f"{val:.8g}",
                    "derived_total_err": f"{err:.8g}",
                    "fig6_reference_value": f"{ref:.8g}",
                    "relative_diff_percent": f"{rel:.4f}",
                    "comment": comment,
                }
            )

    print(f"Added {added} new 62.4 GeV hyperon rows to {INPUT_CSV}")
    print(f"Wrote consistency report to {OUT_CROSSCHECK}")


if __name__ == "__main__":
    main()
