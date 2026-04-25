#!/usr/bin/env python3
import csv
import json
import math
import re
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError


def to_float(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().replace(',', '')
    m = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', s)
    if not m:
        return None
    return float(m.group(1))


def parse_errors(yitem):
    stat = None
    sys = None
    for e in yitem.get('errors', []):
        v = to_float(e.get('symerror'))
        if v is None:
            continue
        label = (e.get('label') or '').lower()
        if 'stat' in label:
            stat = v
        elif 'sys' in label:
            sys = v
        elif stat is None:
            stat = v
        elif sys is None:
            sys = v
    return stat, sys


def fetch_json(url):
    last = None
    for i in range(10):
        try:
            with urlopen(url, timeout=45) as r:
                return json.load(r)
        except HTTPError as e:
            last = e
            if e.code == 429:
                time.sleep(4 + 3 * i)
                continue
            raise
    raise last

base = Path('data/first_group_dn_dy_vs_energy.csv')
rows = list(csv.DictReader(base.open(encoding='utf-8')))

# remove current ID237 rows; we'll rebuild all ID237 species consistently
rows = [r for r in rows if r['paper_id'] != '237']

species_237 = {
    'Ks0': 'K0S',
    'phi': 'Phi',
    'Lambda': 'Lambda',
    'Lambda_bar': 'AntiLambda',
    'Xi': 'Xi-',
    'Xi_bar': 'AntiXi+',
}
energies = [7.7, 11.5, 19.6, 27, 39]

for e in energies:
    for particle, tag in species_237.items():
        url = f'https://www.hepdata.net/download/table/ins1738953/{tag}%20dN/dy,%20Au+Au%20{e}%20GeV/json'
        t = fetch_json(url)
        best = None
        for r in t.get('values', []):
            npart = to_float((r.get('x') or [{}])[0].get('value'))
            if npart is None:
                continue
            yi = (r.get('y') or [{}])[0]
            if best is None or npart > best[0]:
                best = (npart, yi)
        if best is None:
            continue
        _, yi = best
        stat, sys = parse_errors(yi)
        rows.append({
            'paper_id': '237',
            'particle': particle,
            'observable': 'dN/dy',
            'energy_GeV': str(float(e)),
            'value': str(to_float(yi.get('value'))),
            'stat_err': '' if stat is None else str(stat),
            'sys_err': '' if sys is None else str(sys),
            'centrality': 'max Npart row (most central)',
            'source': f'ins1738953/{tag} dN/dy {e} GeV',
            'note': '',
        })
        time.sleep(1.2)

rows.sort(key=lambda r: (int(r['paper_id']), r['particle'], float(r['energy_GeV'])))

with base.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# Replot
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 6))
series_keys = sorted({(r['paper_id'], r['particle']) for r in rows}, key=lambda x: (int(x[0]), x[1]))
for pid, particle in series_keys:
    pts = [r for r in rows if r['paper_id'] == pid and r['particle'] == particle]
    xs = [float(r['energy_GeV']) for r in pts]
    ys = [float(r['value']) for r in pts]
    yerr = []
    for r in pts:
        stat = float(r['stat_err']) if r.get('stat_err') else 0.0
        sys = float(r['sys_err']) if r.get('sys_err') else 0.0
        yerr.append(math.sqrt(stat * stat + sys * sys) if (stat or sys) else 0.0)
    label = f"ID {pid}: {particle} ({pts[0]['observable']})"
    ax.errorbar(xs, ys, yerr=yerr, marker='o', ms=4, lw=1.2, capsize=2, label=label)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ (GeV)')
ax.set_ylabel(r'$dN/dy$')
ax.set_title(r'Central Au+Au: $dN/dy$ vs $\sqrt{s_{NN}}$ (Selected IDs)')
ax.grid(True, which='both', ls='--', alpha=0.35)
ax.legend(fontsize=8, loc='best')
fig.tight_layout()
fig.savefig('data/first_group_dn_dy_vs_energy.png', dpi=240)

print('updated', base)
print('updated data/first_group_dn_dy_vs_energy.png')
print('ID237 particles now:', sorted({r['particle'] for r in rows if r['paper_id']=='237'}))
print('ID237 rows:', sum(1 for r in rows if r['paper_id']=='237'))
