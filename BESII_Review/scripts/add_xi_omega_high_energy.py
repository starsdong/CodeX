#!/usr/bin/env python3
import csv
import json
import math
import re
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError

SRC = Path('data/first_group_yield_vs_energy.csv')


def to_float(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    m = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', str(x))
    return float(m.group(1)) if m else None


def fetch_json(url):
    last = None
    for i in range(10):
        try:
            with urlopen(url, timeout=45) as r:
                return json.load(r)
        except HTTPError as e:
            last = e
            if e.code == 429:
                time.sleep(2 + 2 * i)
                continue
            raise
    raise last


def parse_errors(yitem):
    stat = None
    sys = None
    for e in yitem.get('errors', []):
        v = to_float(e.get('symerror'))
        if v is None:
            continue
        lab = (e.get('label') or '').lower()
        if 'stat' in lab:
            stat = v
        elif 'sys' in lab:
            sys = v
        elif stat is None:
            stat = v
        elif sys is None:
            sys = v
    return stat, sys


rows = list(csv.DictReader(SRC.open(encoding='utf-8')))

# Drop old/replaced additions for this task.
rows = [r for r in rows if r.get('source') not in {
    'arXiv:nucl-ex/0307024/Table I',
    'arXiv:nucl-ex/0606014/Table I',
}]
rows = [r for r in rows if not (r.get('paper_id') == '237' and r.get('particle') in {'Omega', 'Omega_bar', 'Omega_sum'})]

# Low-energy Omega/Omega_bar from previous paper set (ID 237 HEPData, central=max Npart).
for e in [7.7, 11.5, 19.6, 27, 39]:
    vals = {}
    errs = {}
    for particle, tag in [('Omega', 'Omega-'), ('Omega_bar', 'AntiOmega+')]:
        url = f'https://www.hepdata.net/download/table/ins1738953/{tag}%20dN/dy,%20Au+Au%20{e}%20GeV/json'
        t = fetch_json(url)
        best = None
        for rr in t.get('values', []):
            npart = to_float((rr.get('x') or [{}])[0].get('value'))
            if npart is None:
                continue
            yi = (rr.get('y') or [{}])[0]
            if best is None or npart > best[0]:
                best = (npart, yi)
        if best is None:
            continue
        yi = best[1]
        stat, sys = parse_errors(yi)
        val = to_float(yi.get('value'))
        rows.append({
            'paper_id': '237',
            'particle': particle,
            'observable': 'dN/dy',
            'energy_GeV': str(float(e)),
            'value': str(val),
            'stat_err': '' if stat is None else str(stat),
            'sys_err': '' if sys is None else str(sys),
            'centrality': 'max Npart row (most central)',
            'source': f'ins1738953/{tag} dN/dy {e} GeV',
            'note': '',
        })
        vals[particle] = val
        errs[particle] = (stat or 0.0, sys or 0.0)

    if 'Omega' in vals and 'Omega_bar' in vals:
        stat = math.sqrt(errs['Omega'][0] ** 2 + errs['Omega_bar'][0] ** 2)
        sys = math.sqrt(errs['Omega'][1] ** 2 + errs['Omega_bar'][1] ** 2)
        rows.append({
            'paper_id': '237',
            'particle': 'Omega_sum',
            'observable': 'dN/dy',
            'energy_GeV': str(float(e)),
            'value': str(vals['Omega'] + vals['Omega_bar']),
            'stat_err': str(stat),
            'sys_err': str(sys),
            'centrality': 'max Npart row (most central)',
            'source': 'ins1738953/Omega + AntiOmega (summed)',
            'note': 'sum of Omega and anti-Omega',
        })

# High-energy from papers requested by user.
# nucl-ex/0307024 Table I, Au+Au 130 GeV, 0-10%.
# Quoted pT-independent systematic uncertainty is 10% (paper note).
for particle, val, stat in [
    ('Xi', 2.16, 0.09),
    ('Xi_bar', 1.81, 0.08),
    ('Omega_sum', 0.59, 0.14),
]:
    rows.append({
        'paper_id': '3024',
        'particle': particle,
        'observable': 'dN/dy',
        'energy_GeV': '130.0',
        'value': str(val),
        'stat_err': str(stat),
        'sys_err': str(0.10 * val),
        'centrality': '0-10%',
        'source': 'arXiv:nucl-ex/0307024/Table I',
        'note': 'includes 10% pT-independent systematic per paper',
    })

# nucl-ex/0606014 Table I, Au+Au 200 GeV, 0-5%.
for particle, val, stat, sys in [
    ('Lambda', 16.7, 0.2, 1.1),
    ('Lambda_bar', 12.7, 0.2, 0.9),
    ('Xi', 2.17, 0.06, 0.19),
    ('Xi_bar', 1.83, 0.05, 0.20),
    ('Omega_sum', 0.53, 0.04, 0.04),
]:
    rows.append({
        'paper_id': '6014',
        'particle': particle,
        'observable': 'dN/dy',
        'energy_GeV': '200.0',
        'value': str(val),
        'stat_err': str(stat),
        'sys_err': str(sys),
        'centrality': '0-5%',
        'source': 'arXiv:nucl-ex/0606014/Table I',
        'note': '',
    })

rows.sort(key=lambda r: (int(r['paper_id']), r['particle'], float(r['energy_GeV'])))
with SRC.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print('updated', SRC, 'rows', len(rows))
print('added/updated Xi/Xi_bar/Omega_sum high-energy points and Omega low-energy series')
