#!/usr/bin/env python3
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt

SRC = Path('data/first_group_dn_dy_vs_energy.csv')
OUT_DIR = Path('data/particle_figures')
OUT_DIR.mkdir(parents=True, exist_ok=True)


def to_float(x):
    try:
        return float(x)
    except Exception:
        return None


def slug(s):
    s = s.strip().lower().replace('+', 'plus').replace('-', 'minus')
    s = re.sub(r'[^a-z0-9]+', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    return s or 'particle'


rows = list(csv.DictReader(SRC.open(encoding='utf-8')))
particles = sorted({r['particle'] for r in rows})

manifest = []
for particle in particles:
    pts = [r for r in rows if r['particle'] == particle]
    if not pts:
        continue

    # one curve per paper_id for this particle
    paper_ids = sorted({r['paper_id'] for r in pts}, key=lambda x: int(x))

    fig, ax = plt.subplots(figsize=(8, 5))
    for pid in paper_ids:
        p = [r for r in pts if r['paper_id'] == pid]
        p.sort(key=lambda r: to_float(r['energy_GeV']) or 0.0)
        xs = [to_float(r['energy_GeV']) for r in p]
        ys = [to_float(r['value']) for r in p]

        yerr = []
        for r in p:
            stat = to_float(r.get('stat_err')) or 0.0
            sys = to_float(r.get('sys_err')) or 0.0
            yerr.append(math.sqrt(stat * stat + sys * sys) if (stat or sys) else 0.0)

        obs = p[0].get('observable', 'dN/dy')
        ax.errorbar(xs, ys, yerr=yerr, marker='o', ms=4, lw=1.2, capsize=2, label=f'ID {pid} ({obs})')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\sqrt{s_{NN}}$ (GeV)')
    ax.set_ylabel(r'$dN/dy$')
    ax.set_title(f'{particle}: dN/dy vs sqrt(sNN)')
    ax.grid(True, which='both', ls='--', alpha=0.35)
    ax.legend(fontsize=8, loc='best')
    fig.tight_layout()

    out = OUT_DIR / f'dn_dy_vs_energy_{slug(particle)}.png'
    fig.savefig(out, dpi=240)
    plt.close(fig)

    manifest.append({'particle': particle, 'file': str(out), 'n_points': str(len(pts)), 'n_series': str(len(paper_ids))})

manifest_csv = OUT_DIR / 'manifest.csv'
with manifest_csv.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['particle', 'file', 'n_points', 'n_series'])
    w.writeheader()
    w.writerows(manifest)

print(f'wrote {len(manifest)} particle plots to {OUT_DIR}')
print(f'wrote manifest: {manifest_csv}')
