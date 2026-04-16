#!/usr/bin/env python3
import csv
import re
import subprocess
from pathlib import Path
from urllib.request import urlopen

PDF = Path('/tmp/0808.2041.pdf')
TXT = Path('/tmp/0808.2041.txt')
SRC = Path('data/first_group_yield_vs_energy.csv')


def ensure_text():
    if not PDF.exists():
        data = urlopen('https://arxiv.org/pdf/0808.2041.pdf', timeout=60).read()
        PDF.write_bytes(data)
    subprocess.run(['pdftotext', '-layout', str(PDF), str(TXT)], check=True)


def extract_table8_text(text: str) -> str:
    i = text.find('TABLE VIII:')
    j = text.find('TABLE IX:')
    if i < 0 or j < 0 or j <= i:
        raise RuntimeError('Could not locate TABLE VIII/TABLE IX in text')
    return text[i:j]


def parse_row_values(line: str):
    # Row format: centrality + 6 columns each as value ± error
    num = r'([0-9]+(?:\.[0-9]+)?)'
    pat = re.compile(
        r'^\s*([0-9]+-\s*[0-9]+%)\s+'
        + rf'{num}\s*±\s*{num}\s+'
        + rf'{num}\s*±\s*{num}\s+'
        + rf'{num}\s*±\s*{num}\s+'
        + rf'{num}\s*±\s*{num}\s+'
        + rf'{num}\s*±\s*{num}\s+'
        + rf'{num}\s*±\s*{num}\s+'
    )
    m = pat.match(line)
    if not m:
        return None
    c = m.group(1).replace(' ', '')
    vals = [float(x) for x in m.groups()[1:]]
    # columns: pi-, pi+, K-, K+, pbar, p ; each has (val, err)
    out = {
        'centrality': c,
        'pi-': (vals[0], vals[1]),
        'pi+': (vals[2], vals[3]),
        'K-': (vals[4], vals[5]),
        'K+': (vals[6], vals[7]),
        'pbar': (vals[8], vals[9]),
        'p': (vals[10], vals[11]),
    }
    return out


def extract_energy_rows(tbl: str):
    # Keep central rows only: 200->0-5%, 130->0-6%, 62.4->0-5%
    targets = {
        '200.0': '0-5%',
        '130.0': '0-6%',
        '62.4': '0-5%',
    }

    marker_pat = re.compile(r'\b(200|130|62\.4)\s*GeV\b')
    lines = tbl.splitlines()

    # Build block ranges keyed by marker occurrence order
    markers = []
    for i, ln in enumerate(lines):
        m = marker_pat.search(ln)
        if m:
            markers.append((i, m.group(1)))

    blocks = []
    for k, (idx, en) in enumerate(markers):
        end = markers[k + 1][0] if k + 1 < len(markers) else len(lines)
        blocks.append((en, lines[idx:end]))

    wanted = {}
    for en_raw, blines in blocks:
        en = '200.0' if en_raw == '200' else ('130.0' if en_raw == '130' else '62.4')
        if en not in targets or en in wanted:
            continue
        target_c = targets[en]
        for ln in blines:
            row = parse_row_values(ln)
            if not row:
                continue
            if row['centrality'] == target_c:
                wanted[en] = row
                break

    missing = [e for e in targets if e not in wanted]
    if missing:
        raise RuntimeError(f'Missing target rows for energies: {missing}')
    return wanted


def main():
    ensure_text()
    text = TXT.read_text(encoding='utf-8', errors='replace')
    t8 = extract_table8_text(text)
    rows_by_e = extract_energy_rows(t8)

    rows = list(csv.DictReader(SRC.open(encoding='utf-8')))
    # Remove all existing ID 34 entries; re-add from Table VIII.
    rows = [r for r in rows if r.get('paper_id') != '34']

    for en in ['62.4', '130.0', '200.0']:
        row = rows_by_e[en]
        cent = row['centrality']
        for particle in ['pi-', 'pi+', 'K-', 'K+', 'pbar', 'p']:
            v, e = row[particle]
            rows.append({
                'paper_id': '34',
                'particle': particle,
                'observable': 'dN/dy',
                'energy_GeV': en,
                'value': str(v),
                'stat_err': str(e),
                'sys_err': '',
                'centrality': cent,
                'source': 'arXiv:0808.2041/Table VIII',
                'note': 'table gives total uncertainty (stat+sys)',
            })

    rows.sort(key=lambda r: (int(r['paper_id']), r['particle'], float(r['energy_GeV'])))
    with SRC.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print('updated', SRC, 'rows', len(rows))
    for en in ['62.4', '130.0', '200.0']:
        print('added ID34 @', en, 'centrality', rows_by_e[en]['centrality'])


if __name__ == '__main__':
    main()
