#!/usr/bin/env python3
import csv
import json
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


def fetch(url):
    last = None
    for i in range(10):
        try:
            with urlopen(url, timeout=50) as r:
                return json.load(r)
        except HTTPError as e:
            last = e
            if e.code == 429:
                time.sleep(2 + i * 2)
                continue
            raise
    raise last


def group_map(table, qkey):
    out = {}
    for q in table.get('qualifiers', {}).get(qkey, []):
        g = q.get('group')
        if g is None:
            continue
        out[int(g)] = (q.get('value') or '').strip().lower()
    return out

rows = list(csv.DictReader(SRC.open(encoding='utf-8')))
# rebuild IDs 207 and 237 from scratch
rows = [r for r in rows if r['paper_id'] not in {'207', '237'}]

# ID 207 exact species mapping
species_207 = {
    'pi-': {'pi-', 'pi−'},
    'pi+': {'pi+', 'pi＋'},
    'K-': {'ka-', 'k-', 'ka−', 'k−'},
    'K+': {'ka+', 'k+', 'ka＋', 'k＋'},
    'pbar': {'pbar', 'anti-proton', 'antiproton', 'p̄'},
    'p': {'proton'},
}
for e in [7.7, 11.5, 19.6, 27, 39]:
    t = fetch(f'https://www.hepdata.net/download/table/ins1510593/dN/dy,%20{e}%20GeV/json')
    gm = group_map(t, 'particle type')

    sel = {}
    for g, lbl in gm.items():
        clean = lbl.replace(' ', '')
        for sp, opts in species_207.items():
            if clean in {o.replace(' ', '') for o in opts}:
                sel[sp] = g

    for sp, gsel in sel.items():
        best = None
        for r in t.get('values', []):
            npart = to_float((r.get('x') or [{}])[0].get('value'))
            if npart is None:
                continue
            yis = [yi for yi in r.get('y', []) if int(yi.get('group', -1)) == gsel]
            if not yis:
                continue
            yi = yis[0]
            ynorm = to_float(yi.get('value'))
            if ynorm is None:
                continue
            if best is None or npart > best[0]:
                best = (npart, yi, ynorm)
        if best is None:
            continue
        npart, yi, ynorm = best
        stat_n, sys_n = parse_errors(yi)
        rows.append({
            'paper_id': '207',
            'particle': sp,
            'observable': 'dN/dy',
            'energy_GeV': str(float(e)),
            'value': str(ynorm * 0.5 * npart),
            'stat_err': '' if stat_n is None else str(stat_n * 0.5 * npart),
            'sys_err': '' if sys_n is None else str(sys_n * 0.5 * npart),
            'centrality': 'max Npart row (most central)',
            'source': f'ins1510593/dN/dy, {e} GeV',
            'note': 'derived from dN/dy/(0.5*Npart) using Npart',
        })

# ID 237 full strange set
species_237 = {
    'phi': 'Phi',
    'Lambda': 'Lambda',
    'Lambda_bar': 'AntiLambda',
    'Xi': 'Xi-',
    'Xi_bar': 'AntiXi+',
}
for e in [7.7, 11.5, 19.6, 27, 39]:
    for particle, tag in species_237.items():
        t = fetch(f'https://www.hepdata.net/download/table/ins1738953/{tag}%20dN/dy,%20Au+Au%20{e}%20GeV/json')
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
        time.sleep(0.7)

rows.sort(key=lambda r: (int(r['paper_id']), r['particle'], float(r['energy_GeV'])))
with SRC.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print('updated', SRC)
print('rows', len(rows))
