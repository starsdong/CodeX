#!/usr/bin/env python3
import csv
import json
import math
import re
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError

import matplotlib.pyplot as plt

OUT_DIR = Path('data')
OUT_DIR.mkdir(exist_ok=True)


def fetch_json(url: str):
    last_err = None
    for i in range(6):
        try:
            with urlopen(url, timeout=40) as r:
                return json.load(r)
        except HTTPError as e:
            last_err = e
            if e.code == 429:
                time.sleep(1.5 * (i + 1))
                continue
            raise
    if last_err is not None:
        raise last_err
    raise RuntimeError(f'Failed to fetch {url}')


def to_float(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if not s:
        return None
    s = s.replace(',', '')
    m = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', s)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def parse_errors(yitem):
    stat = None
    sys = None
    for e in yitem.get('errors', []):
        label = (e.get('label') or '').lower()
        v = to_float(e.get('symerror'))
        if v is None:
            continue
        if 'stat' in label:
            stat = v
        elif 'sys' in label:
            sys = v
        elif stat is None:
            stat = v
        elif sys is None:
            sys = v
    return stat, sys


def group_map(table, qkey):
    out = {}
    arr = table.get('qualifiers', {}).get(qkey, [])
    for q in arr:
        g = q.get('group')
        if g is None:
            continue
        out[int(g)] = q.get('value', '')
    return out


rows = []

# ID 188 proxy: k_s(P_T^s) at first point from energy tables (not integrated dN/dy)
for tnum in [80, 81, 82, 83, 84, 85]:
    url = f'https://www.hepdata.net/download/table/ins1378002/Table%20{tnum}/json'
    try:
        t = fetch_json(url)
    except Exception:
        continue
    q = t.get('qualifiers', {})
    e = to_float((q.get('SQRT(S)/NUCLEON') or [{}])[0].get('value', ''))
    c = (q.get('CENTRALITY') or [{}])[0].get('value', '')
    vals = t.get('values', [])
    if not vals:
        continue
    p0 = vals[0]
    x = to_float((p0.get('x') or [{}])[0].get('value'))
    yitem = (p0.get('y') or [{}])[0]
    y = to_float(yitem.get('value'))
    stat, sys = parse_errors(yitem)
    rows.append({
        'paper_id': '188',
        'particle': 'ks proxy (Omega+antiOmega)/(2phi)',
        'observable': 'k_s(P_T^s) at first bin (proxy)',
        'energy_GeV': e,
        'value': y,
        'stat_err': stat,
        'sys_err': sys,
        'centrality': c,
        'source': f'ins1378002/Table {tnum}',
        'note': f'proxy at P_T^s={x}; not integrated dN/dy',
    })

# ID 207: dN/dy for pi+/-, K+/-, p, pbar derived from dN/dy/(0.5*Npart)
species_aliases = {
    'pi+': ['pi+'],
    'pi-': ['pi-'],
    'K+': ['ka+', 'k+'],
    'K-': ['ka-', 'k-'],
    'p': ['proton', 'p'],
    'pbar': ['pbar', 'antiproton'],
}
for e in [7.7, 11.5, 19.6, 27, 39]:
    t = fetch_json(f'https://www.hepdata.net/download/table/ins1510593/dN/dy,%20{e}%20GeV/json')
    gmap = group_map(t, 'particle type')
    glabel = {g: (name or '').strip().lower() for g, name in gmap.items()}

    for sp, aliases in species_aliases.items():
        gsel = None
        for g, lname in glabel.items():
            if any(a in lname for a in aliases):
                gsel = g
                break
        if gsel is None:
            continue

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
        y = ynorm * 0.5 * npart
        stat = stat_n * 0.5 * npart if stat_n is not None else None
        sys = sys_n * 0.5 * npart if sys_n is not None else None
        rows.append({
            'paper_id': '207',
            'particle': sp,
            'observable': 'dN/dy',
            'energy_GeV': float(e),
            'value': y,
            'stat_err': stat,
            'sys_err': sys,
            'centrality': 'max Npart row (most central)',
            'source': f'ins1510593/dN/dy, {e} GeV',
            'note': 'derived from dN/dy/(0.5*Npart) using Npart',
        })

# ID 237: dN/dy at max Npart (central)
species_237 = {
    'phi': 'Phi',
    'Lambda': 'Lambda',
    'Lambda_bar': 'AntiLambda',
    'Xi': 'Xi-',
    'Xi_bar': 'AntiXi+',
}
for e in [7.7, 11.5, 19.6, 27, 39]:
    for particle, table_tag in species_237.items():
        url = f'https://www.hepdata.net/download/table/ins1738953/{table_tag}%20dN/dy,%20Au+Au%20{e}%20GeV/json'
        t = fetch_json(url)
        best = None
        for r in t.get('values', []):
            x = to_float((r.get('x') or [{}])[0].get('value'))
            if x is None:
                continue
            yi = (r.get('y') or [{}])[0]
            if best is None or x > best[0]:
                best = (x, yi)
        if best is None:
            continue
        yi = best[1]
        y = to_float(yi.get('value'))
        stat, sys = parse_errors(yi)
        rows.append({
            'paper_id': '237',
            'particle': particle,
            'observable': 'dN/dy',
            'energy_GeV': float(e),
            'value': y,
            'stat_err': stat,
            'sys_err': sys,
            'centrality': 'max Npart row (most central)',
            'source': f'ins1738953/{table_tag} dN/dy {e} GeV',
            'note': '',
        })

# ID 261: proton dN/dy at 14.5 GeV, max Npart
q261 = fetch_json('https://www.hepdata.net/download/table/ins1748776/dN/dy%20as%20a%20function%20of%20$N_{part}$%20(p)/json')
best = None
for r in q261.get('values', []):
    x = to_float((r.get('x') or [{}])[0].get('value'))
    if x is None:
        continue
    yi = (r.get('y') or [{}])[0]
    if best is None or x > best[0]:
        best = (x, yi)
yi = best[1]
rows.append({
    'paper_id': '261',
    'particle': 'p',
    'observable': 'dN/dy',
    'energy_GeV': 14.5,
    'value': to_float(yi.get('value')),
    'stat_err': parse_errors(yi)[0],
    'sys_err': parse_errors(yi)[1],
    'centrality': 'max Npart row (most central)',
    'source': 'ins1748776/dNdy_vs_Npart(p)',
    'note': '',
})

# ID 264: deuteron dN/dy derived from normalized dN/dy/(0.5*Npart), max Npart
for tnum in [88, 89, 90, 91, 92, 93, 94, 95]:
    t = fetch_json(f'https://www.hepdata.net/download/table/ins1727273/Table%20{tnum}%20from%20Fig%208/json')
    h = [x.get('name', '') for x in t.get('headers', [])]
    if len(h) < 2:
        continue
    e = to_float(h[1])
    best = None
    for r in t.get('values', []):
        x = to_float((r.get('x') or [{}])[0].get('value'))
        if x is None:
            continue
        yi = (r.get('y') or [{}])[0]
        if best is None or x > best[0]:
            best = (x, yi)
    if best is None:
        continue
    npart, yi = best
    stat_n, sys_n = parse_errors(yi)
    ynorm = to_float(yi.get('value'))
    y = ynorm * 0.5 * npart if ynorm is not None else None
    stat = stat_n * 0.5 * npart if stat_n is not None else None
    sys = sys_n * 0.5 * npart if sys_n is not None else None
    rows.append({
        'paper_id': '264',
        'particle': 'd',
        'observable': 'dN/dy',
        'energy_GeV': e,
        'value': y,
        'stat_err': stat,
        'sys_err': sys,
        'centrality': 'max Npart row (most central)',
        'source': f'ins1727273/Table {tnum} Fig8',
        'note': 'derived from dN/dy/(0.5*Npart) using Npart',
    })

# ID 310: phi dN/dy at y closest to 0, central 0-10 (group 0)
q310 = fetch_json('https://www.hepdata.net/download/table/ins1897327/Data%20from%20Figure%203(b)/json')
best = None
for r in q310.get('values', []):
    x = to_float((r.get('x') or [{}])[0].get('value'))
    if x is None:
        continue
    yis = [yi for yi in r.get('y', []) if int(yi.get('group', -1)) == 0]
    if not yis:
        continue
    yi = yis[0]
    score = abs(x)
    if best is None or score < best[0]:
        best = (score, yi, x)
_, yi, xnear = best
rows.append({
    'paper_id': '310',
    'particle': 'phi',
    'observable': 'dN/dy at y~0',
    'energy_GeV': 3.0,
    'value': to_float(yi.get('value')),
    'stat_err': parse_errors(yi)[0],
    'sys_err': parse_errors(yi)[1],
    'centrality': '0-10%',
    'source': 'ins1897327/Figure3(b)',
    'note': f'point at y={xnear}',
})

# ID 330: triton integral dN/dy, central 0-10 group
q330 = fetch_json('https://www.hepdata.net/download/table/ins2152917/Figure1,%20Triton%20integral%20yields%20dN/dy/json')
cent = group_map(q330, 'Collision energy (GeV)')
g0 = None
for g, v in cent.items():
    if '0-10' in v.replace('\\', ''):
        g0 = g
        break
if g0 is None:
    g0 = 0
for r in q330.get('values', []):
    e = to_float((r.get('x') or [{}])[0].get('value'))
    yis = [yi for yi in r.get('y', []) if int(yi.get('group', -1)) == g0]
    if not yis:
        continue
    yi = yis[0]
    stat, sys = parse_errors(yi)
    rows.append({
        'paper_id': '330',
        'particle': 't',
        'observable': 'integral dN/dy',
        'energy_GeV': e,
        'value': to_float(yi.get('value')),
        'stat_err': stat,
        'sys_err': sys,
        'centrality': '0-10%',
        'source': 'ins2152917/Figure1 triton integral dN/dy',
        'note': '',
    })

# ID 353: Lambda dN/dy at y closest to 0, central 0-10 (group 0)
q353 = fetch_json('https://www.hepdata.net/download/table/ins2807679/Table%201/json')
best = None
for r in q353.get('values', []):
    x = to_float((r.get('x') or [{}])[0].get('value'))
    if x is None:
        continue
    yis = [yi for yi in r.get('y', []) if int(yi.get('group', -1)) == 0]
    if not yis:
        continue
    yi = yis[0]
    score = abs(x)
    if best is None or score < best[0]:
        best = (score, yi, x)
_, yi, xnear = best
rows.append({
    'paper_id': '353',
    'particle': 'Lambda',
    'observable': 'dN/dy at y~0',
    'energy_GeV': 3.0,
    'value': to_float(yi.get('value')),
    'stat_err': parse_errors(yi)[0],
    'sys_err': parse_errors(yi)[1],
    'centrality': '0-10%',
    'source': 'ins2807679/Table1',
    'note': f'point at y={xnear}',
})

# Save data
rows = [r for r in rows if r.get('energy_GeV') is not None and r.get('value') is not None]
rows.sort(key=lambda r: (int(r['paper_id']), r['energy_GeV']))

out_csv = OUT_DIR / 'first_group_yield_vs_energy.csv'
with out_csv.open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=[
        'paper_id', 'particle', 'observable', 'energy_GeV', 'value', 'stat_err', 'sys_err',
        'centrality', 'source', 'note'
    ])
    w.writeheader()
    w.writerows(rows)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ids = sorted({r['paper_id'] for r in rows}, key=int)
for pid in ids:
    pts = [r for r in rows if r['paper_id'] == pid]
    xs = [r['energy_GeV'] for r in pts]
    ys = [r['value'] for r in pts]
    yerr = []
    for r in pts:
        s = r['stat_err'] or 0.0
        y = r['sys_err'] or 0.0
        yerr.append(math.sqrt(s * s + y * y) if (s or y) else 0.0)
    label = f"ID {pid}: {pts[0]['particle']} ({pts[0]['observable']})"
    ax.errorbar(xs, ys, yerr=yerr, marker='o', ms=4, lw=1.2, capsize=2, label=label)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\sqrt{s_{NN}}$ (GeV)')
ax.set_ylabel('Yield Observable')
ax.set_title('Central Au+Au Yield vs Energy (Selected IDs)')
ax.grid(True, which='both', ls='--', alpha=0.35)
ax.legend(fontsize=7, ncol=1, loc='best')
fig.tight_layout()
out_png = OUT_DIR / 'first_group_yield_vs_energy.png'
fig.savefig(out_png, dpi=220)
print(f'wrote {out_csv}')
print(f'wrote {out_png}')
print(f'rows: {len(rows)}; ids: {", ".join(ids)}')
