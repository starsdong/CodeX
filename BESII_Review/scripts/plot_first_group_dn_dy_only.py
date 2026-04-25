#!/usr/bin/env python3
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

in_csv = Path('data/first_group_yield_vs_energy.csv')
out_csv = Path('data/first_group_dn_dy_vs_energy.csv')
out_png = Path('data/first_group_dn_dy_vs_energy.png')

rows = list(csv.DictReader(in_csv.open(encoding='utf-8')))


def paper_sort_key(value: str):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)

keep = []
for r in rows:
    obs = (r.get('observable') or '').lower()
    if 'dn/dy' not in obs:
        continue
    if 'proxy' in obs:
        continue
    keep.append(r)

keep.sort(key=lambda r: (paper_sort_key(r['paper_id']), float(r['energy_GeV'])))

with out_csv.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(keep[0].keys()))
    w.writeheader()
    w.writerows(keep)

fig, ax = plt.subplots(figsize=(10, 6))
def series_key(r):
    # Connect low-energy ID207 positive particles with 200 GeV anchors from ID34.
    if r['particle'] in {'pi+', 'K+', 'p'} and r['paper_id'] in {'207', '34'}:
        return ('link_pos', r['particle'])
    return (r['paper_id'], r['particle'])


series_keys = sorted({series_key(r) for r in keep}, key=lambda x: (x[0], x[1]))
for k0, particle in series_keys:
    pts = [r for r in keep if series_key(r) == (k0, particle)]
    xs = [float(r['energy_GeV']) for r in pts]
    ys = [float(r['value']) for r in pts]
    yerr = []
    for r in pts:
        stat = float(r['stat_err']) if r.get('stat_err') else 0.0
        sys = float(r['sys_err']) if r.get('sys_err') else 0.0
        yerr.append(math.sqrt(stat * stat + sys * sys) if (stat or sys) else 0.0)
    if k0 == 'link_pos':
        label = f"{particle} (ID207 + 200 GeV anchor)"
    else:
        label = f"ID {k0}: {particle} ({pts[0]['observable']})"
    ax.errorbar(xs, ys, yerr=yerr, marker='o', ms=4, lw=1.3, capsize=2, label=label)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ (GeV)')
ax.set_ylabel(r'$dN/dy$')
ax.set_title(r'Central Au+Au: $dN/dy$ vs $\sqrt{s_{NN}}$ (Selected IDs)')
ax.grid(True, which='both', ls='--', alpha=0.35)
ax.legend(fontsize=8, loc='best')
fig.tight_layout()
fig.savefig(out_png, dpi=240)

print(f'wrote {out_csv}')
print(f'wrote {out_png}')
print('rows', len(keep), 'series', len(series_keys))
