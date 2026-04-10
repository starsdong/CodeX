#!/usr/bin/env python3
import csv
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

BASE = Path('/Users/starsdong/Work/CodeX/SQM2026')
XLSX_PATH = BASE / 'registrations_final.xlsx'

NS = {'a': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


US_KEYWORDS = [
    'ucla', 'university of houston', 'university of illinois', 'indiana university',
    'ohio state university', 'yale university', 'university of california', 'uc berkeley',
    'lawrence berkeley', 'lawrence livermore', 'los alamos', 'brookhaven', 'bnl', 'lbl', 'lbnl',
    'argonne', 'purdue university', 'vanderbilt university', 'stony brook', 'duke university',
    'duke ', 'east carolina university', 'massachusetts institute of technology',
    'prairie view', 'texas southern university', 'university of kansas', 'wayne state university',
    'university of texas', 'uic', 'facility for rare isotope beams', 'kent state university',
    'university of colorado boulder', 'kavli institute for theoretical physics', 'aamu',
    'columbia university', 'lehigh university'
]

EUROPE_KEYWORDS = [
    'cern', 'infn', 'universita', 'università', 'turin', 'torino', 'padova', 'trieste', 'cagliari',
    'catania', 'salerno', 'gsi', 'heidelberg', 'technical university of munich', 'university of bonn',
    'frankfurt institute', 'goethe university', 'paris-saclay', 'strasbourg', 'subatech', 'lyon',
    'niels bohr', 'copenhagen', 'jyväskylä', 'jyvaskyla', 'lund university', 'warsaw', 'prague',
    'ctu', 'austrian academy of sciences', 'polish academy of sciences', 'wigner', 'eotvos',
    'eötvös', 'hun-ren', 'university of ljubljana', 'kyiv', 'fnspe'
]

ASIA_KEYWORDS = [
    'fudan', 'tsinghua', 'central china normal', 'ccnu', 'shandong university', 'peking university',
    'chinese academy of sciences', 'china normal', 'ustc', 'university of science and technology of china',
    'sungkyunkwan', 'pusan national university', 'korea university', 'yonsei university', 'kangwon',
    'riken', 'university of tokyo', 'university of tsukuba', 'hiroshima university',
    'indian institute', 'iiser', 'iit', 'aligarh muslim university', 'niser', 'bose institute',
    'university of jammu', 'mrpd government college talwara', 'jazan university',
    'federal urdu university', 'agriculture peshawar', 'institue of particle physics, central china normal university',
    '中国科学院大学', 'academia sinica', 'hari singh gour vishwavidyalaya'
]

OTHER_KEYWORDS = [
    'brazil', 'sao paulo', 'universidade', 'eswatini', 'egyptian', 'universidad mayor', 'chile'
]

EUROPE_TLDS = {'fr', 'de', 'it', 'pl', 'dk', 'fi', 'hu', 'cz', 'si', 'ch', 'at', 'nl', 'be', 'se', 'eu'}
ASIA_TLDS = {'cn', 'kr', 'jp', 'in', 'pk', 'sa', 'tw'}
OTHER_TLDS = {'br', 'za', 'eg', 'cl'}


def col_letters(cell_ref: str) -> str:
    m = re.match(r'([A-Z]+)', cell_ref)
    return m.group(1) if m else ''


def load_rows(xlsx_path: Path):
    with zipfile.ZipFile(xlsx_path) as zf:
        sst = []
        root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
        for si in root.findall('a:si', NS):
            sst.append(''.join((t.text or '') for t in si.findall('.//a:t', NS)))

        ws = ET.fromstring(zf.read('xl/worksheets/sheet1.xml'))
        rows = ws.find('a:sheetData', NS).findall('a:row', NS)

        out = []
        for row in rows[1:]:
            rowmap = {}
            for c in row.findall('a:c', NS):
                col = col_letters(c.attrib.get('r', ''))
                t = c.attrib.get('t')
                v = c.find('a:v', NS)
                val = ''
                if v is not None and v.text is not None:
                    val = v.text
                    if t == 's':
                        val = sst[int(val)]
                rowmap[col] = val.strip()
            out.append({
                'ID': rowmap.get('A', ''),
                'Name': rowmap.get('B', ''),
                'Email': rowmap.get('C', ''),
                'Affiliation': rowmap.get('D', ''),
                'Country': rowmap.get('E', ''),
            })
    return out


def normalized_tld(email: str) -> str:
    if '@' not in email:
        return ''
    host = email.split('@', 1)[1].lower()
    parts = host.split('.')
    if len(parts) < 2:
        return ''
    if parts[-1] in {'uk', 'jp', 'kr', 'cn', 'br', 'in', 'fr', 'de', 'it', 'ch', 'at', 'pl', 'es', 'nl', 'hu', 'se', 'fi', 'dk', 'be', 'tw', 'au', 'za', 'eg', 'si', 'sa'} and len(parts) >= 3:
        return '.'.join(parts[-2:])
    return parts[-1]


def contains_any(text: str, needles) -> bool:
    return any(n in text for n in needles)


def classify_region(affiliation: str, email: str):
    txt = (affiliation or '').strip().lower()
    tld = normalized_tld(email)

    if contains_any(txt, US_KEYWORDS) or '(us)' in txt or ' usa' in txt or 'united states' in txt:
        return 'US', 'affiliation keyword'
    if contains_any(txt, EUROPE_KEYWORDS) or '(fr)' in txt or '(de)' in txt or '(it)' in txt or '(pl)' in txt or '(cz)' in txt or '(dk)' in txt or '(at)' in txt:
        return 'Europe', 'affiliation keyword'
    if contains_any(txt, ASIA_KEYWORDS) or '(cn)' in txt or '(kr)' in txt or '(jp)' in txt or '(in)' in txt:
        return 'Asia', 'affiliation keyword'
    if contains_any(txt, OTHER_KEYWORDS):
        return 'Other/Unknown', 'affiliation keyword'

    # email domain fallback
    if tld in {'edu', 'gov'}:
        return 'US', 'email tld fallback'
    if tld.endswith('.cn') or tld.endswith('.kr') or tld.endswith('.jp') or tld.endswith('.in') or tld in ASIA_TLDS:
        return 'Asia', 'email tld fallback'
    if tld in EUROPE_TLDS or any(tld.endswith('.' + x) for x in EUROPE_TLDS):
        return 'Europe', 'email tld fallback'
    if tld in OTHER_TLDS or any(tld.endswith('.' + x) for x in OTHER_TLDS):
        return 'Other/Unknown', 'email tld fallback'

    return 'Other/Unknown', 'unclassified'


def pie_path_d(cx, cy, r, start_angle, end_angle):
    x1 = cx + r * math.cos(start_angle)
    y1 = cy + r * math.sin(start_angle)
    x2 = cx + r * math.cos(end_angle)
    y2 = cy + r * math.sin(end_angle)
    large_arc = 1 if (end_angle - start_angle) > math.pi else 0
    return f'M {cx:.3f} {cy:.3f} L {x1:.3f} {y1:.3f} A {r:.3f} {r:.3f} 0 {large_arc} 1 {x2:.3f} {y2:.3f} Z'


def make_svg(title, counts, out_path):
    total = sum(counts.values())
    width, height = 900, 600
    cx, cy, r = 280, 320, 180
    palette = {
        'US': '#2563eb',
        'Europe': '#059669',
        'Asia': '#d97706',
        'Other/Unknown': '#6b7280',
    }

    labels = [k for k, v in counts.items() if v > 0]
    start = -math.pi / 2
    slices = []
    for k in labels:
        frac = counts[k] / total if total else 0
        end = start + 2 * math.pi * frac
        slices.append((k, counts[k], frac, start, end))
        start = end

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="40" y="60" font-family="Arial, sans-serif" font-size="30" font-weight="700" fill="#111827">{title}</text>',
        f'<text x="40" y="92" font-family="Arial, sans-serif" font-size="18" fill="#374151">Total registrations: {total}</text>',
    ]

    for k, c, frac, s, e in slices:
        lines.append(f'<path d="{pie_path_d(cx, cy, r, s, e)}" fill="{palette.get(k, "#999999")}" stroke="#ffffff" stroke-width="2"/>')

    lx, ly = 540, 180
    for i, (k, c, frac, _, _) in enumerate(slices):
        y = ly + i * 50
        lines.append(f'<rect x="{lx}" y="{y-16}" width="20" height="20" fill="{palette.get(k, "#999999")}"/>')
        lines.append(f'<text x="{lx+30}" y="{y}" font-family="Arial, sans-serif" font-size="20" fill="#111827">{k}: {c} ({frac*100:.1f}%)</text>')

    lines.append('</svg>')
    out_path.write_text('\n'.join(lines), encoding='utf-8')


def main():
    rows = load_rows(XLSX_PATH)
    region_counts = Counter()
    region3_counts = Counter()

    details_path = BASE / 'geo_classification_details.csv'
    with details_path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['ID', 'Name', 'Email', 'Affiliation', 'Region', 'Method'])
        for r in rows:
            region, method = classify_region(r['Affiliation'], r['Email'])
            region_counts[region] += 1
            if region in {'US', 'Europe', 'Asia'}:
                region3_counts[region] += 1
            w.writerow([r['ID'], r['Name'], r['Email'], r['Affiliation'], region, method])

    # ensure stable order
    ordered_all = Counter({k: region_counts.get(k, 0) for k in ['US', 'Europe', 'Asia', 'Other/Unknown']})
    ordered_3 = Counter({k: region3_counts.get(k, 0) for k in ['US', 'Europe', 'Asia']})

    summary_path = BASE / 'geo_region_summary.csv'
    with summary_path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['Region', 'Count'])
        for k in ['US', 'Europe', 'Asia', 'Other/Unknown']:
            w.writerow([k, ordered_all[k]])

    make_svg('Institution Geography Distribution (US / Europe / Asia / Other)', ordered_all, BASE / 'geo_pie_all.svg')
    make_svg('Institution Geography Distribution (US / Europe / Asia only)', ordered_3, BASE / 'geo_pie_us_europe_asia.svg')

    print('Wrote:')
    print(summary_path)
    print(BASE / 'geo_pie_all.svg')
    print(BASE / 'geo_pie_us_europe_asia.svg')
    print(details_path)
    print('\nCounts:', dict(ordered_all))


if __name__ == '__main__':
    main()
