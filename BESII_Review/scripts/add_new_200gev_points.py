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
    m = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', str(x))
    return float(m.group(1)) if m else None


def fetch(url):
    last = None
    for i in range(10):
        try:
            with urlopen(url, timeout=40) as r:
                return json.load(r)
        except HTTPError as e:
            last = e
            if e.code == 429:
                time.sleep(2 + 2 * i)
                continue
            raise
    raise last


def combined_err(yitem):
    vals = []
    for e in yitem.get('errors', []):
        v = to_float(e.get('symerror'))
        if v is not None:
            vals.append(v)
    if not vals:
        return None, None
    if len(vals) == 1:
        return vals[0], None
    return vals[0], vals[1]

csv_path = Path('data/first_group_yield_vs_energy.csv')
rows = list(csv.DictReader(csv_path.open(encoding='utf-8')))
rows = [r for r in rows if r['paper_id'] not in {'34', '43', '137'}]

# ID 34: Figure 2-1, 0-5% Au+Au, y~0 point
f34 = fetch('https://www.hepdata.net/download/table/ins630160/Figure%202-1/json')
q = {int(x.get('group')): x.get('value', '') for x in f34.get('qualifiers', {}).get('$y$', []) if x.get('group') is not None}

best = None
for r in f34.get('values', []):
    x = to_float((r.get('x') or [{}])[0].get('value'))
    if x is None:
        continue
    if best is None or abs(x) < best[0]:
        best = (abs(x), x, r)

if best is not None:
    _, ynear, rr = best
    for yi in rr.get('y', []):
        g = int(yi.get('group', -1))
        label = q.get(g, '')
        if '0-5' not in label or 'au+au' not in label.lower():
            continue

        particle = None
        scale = 1.0
        ll = label.lower()
        if 'pi' in ll:
            particle = 'pi-'
            if '/5' in ll:
                scale = 5.0
        elif 'k^-' in ll or 'k-' in ll:
            particle = 'K-'
        elif '\\bar{p}' in label or 'bar{p}' in ll:
            particle = 'pbar'

        if particle is None:
            continue

        y = to_float(yi.get('value'))
        if y is None:
            continue
        stat, sys = combined_err(yi)
        if stat is not None:
            stat *= scale
        if sys is not None:
            sys *= scale

        note = 'from y~0 rapidity point'
        if scale != 1.0:
            note += '; scaled from plotted quantity (/5)'

        rows.append({
            'paper_id': '34',
            'particle': particle,
            'observable': 'dN/dy at y~0',
            'energy_GeV': '200.0',
            'value': str(y * scale),
            'stat_err': '' if stat is None else str(stat),
            'sys_err': '' if sys is None else str(sys),
            'centrality': '0-5%',
            'source': 'ins630160/Figure 2-1',
            'note': note,
        })

# ID 43: Table 1, 0-5% centrality, dN/dy group
f43 = fetch('https://www.hepdata.net/download/table/ins651461/Table%201/json')
for r in f43.get('values', []):
    x = (r.get('x') or [{}])[0]
    lo = to_float(x.get('low'))
    hi = to_float(x.get('high'))
    if lo == 0 and hi == 5:
        for yi in r.get('y', []):
            if int(yi.get('group', -1)) == 1:  # dN/dy
                stat, sys = combined_err(yi)
                y = to_float(yi.get('value'))
                rows.append({
                    'paper_id': '43',
                    'particle': 'phi',
                    'observable': 'dN/dy',
                    'energy_GeV': '200.0',
                    'value': str(y),
                    'stat_err': '' if stat is None else str(stat),
                    'sys_err': '' if sys is None else str(sys),
                    'centrality': '0-5%',
                    'source': 'ins651461/Table 1',
                    'note': '',
                })
        break

# Add positive-charge 200 GeV anchors by charge-symmetry assumption from ID34 negatives
# These selected tables do not provide direct pi+, K+, p dN/dy points.
neg200 = {r['particle']: r for r in rows if r['paper_id'] == '34' and float(r['energy_GeV']) == 200.0}
assumed_map = [('pi+', 'pi-'), ('K+', 'K-'), ('p', 'pbar')]
for pos, neg in assumed_map:
    src = neg200.get(neg)
    if not src:
        continue
    rows.append({
        'paper_id': '34',
        'particle': pos,
        'observable': 'dN/dy at y~0',
        'energy_GeV': '200.0',
        'value': src['value'],
        'stat_err': src.get('stat_err', ''),
        'sys_err': src.get('sys_err', ''),
        'centrality': src.get('centrality', '0-5%'),
        'source': src.get('source', 'ins630160/Figure 2-1'),
        'note': 'assumed charge symmetry at 200 GeV from corresponding negative species',
    })

rows.sort(key=lambda r: (int(r['paper_id']), r['particle'], float(r['energy_GeV'])))
with csv_path.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print('updated', csv_path)
print('added ID34/43 200 GeV dN/dy points (+ positive-charge symmetry anchors); ID137 not added (no explicit dN/dy table)')
