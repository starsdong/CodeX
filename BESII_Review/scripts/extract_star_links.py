#!/usr/bin/env python3
import argparse
import csv
import html
import re
import time
from urllib import request, error

ROW_RE = re.compile(
    r'<span\s+class\s*=\s*"starpublications-table-entry">\s*(.*?)\s*</span>\s*</td>\s*<td>\s*:\s*</td>\s*<td>(.*?)</td>',
    re.IGNORECASE | re.DOTALL,
)
HREF_RE = re.compile(r'href\s*=\s*"([^"]+)"', re.IGNORECASE)
ARXIV_ID_RE = re.compile(r'\((\d{4}\.\d{4,5}(?:v\d+)?)\)')


def strip_tags(s: str) -> str:
    s = re.sub(r'<br\s*/?>', ' ', s, flags=re.IGNORECASE)
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def first_href(s: str) -> str:
    m = HREF_RE.search(s)
    return html.unescape(m.group(1)) if m else ''


def all_hrefs(s: str):
    return [html.unescape(x) for x in HREF_RE.findall(s)]


def fetch(url: str, timeout: float) -> str:
    req = request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Data extraction for publication metadata)"})
    with request.urlopen(req, timeout=timeout) as resp:
        enc = resp.headers.get_content_charset() or 'utf-8'
        return resp.read().decode(enc, errors='replace')


def parse_fields(page_html: str):
    fields = {}
    for key_html, val_html in ROW_RE.findall(page_html):
        key = strip_tags(key_html).rstrip(':').strip().lower()
        if key and key not in fields:
            fields[key] = val_html
    return fields


def enrich_row(row: dict, timeout: float):
    out = dict(row)
    out.update({
        'reference': '',
        'status': '',
        'journal_link': '',
        'journal_link_text': '',
        'hepdata_link': '',
        'arxiv_id': '',
        'arxiv_abs_link': '',
        'inspire_id': '',
        'fetch_error': '',
    })

    url = (row.get('star_url') or '').strip()
    if not url:
        out['fetch_error'] = 'missing star_url'
        return out

    try:
        page = fetch(url, timeout)
    except (error.URLError, TimeoutError, OSError) as e:
        out['fetch_error'] = str(e)
        return out

    fields = parse_fields(page)

    if 'reference' in fields:
        out['reference'] = strip_tags(fields['reference'])

    if 'status' in fields:
        out['status'] = strip_tags(fields['status'])

    if 'journal article' in fields:
        jhtml = fields['journal article']
        out['journal_link'] = first_href(jhtml)
        out['journal_link_text'] = strip_tags(jhtml)

    if 'hepdata data' in fields:
        hhtml = fields['hepdata data']
        out['hepdata_link'] = first_href(hhtml)

    if not out['hepdata_link'] and 'inspire id' in fields:
        links = all_hrefs(fields['inspire id'])
        for link in links:
            if 'hepdata.net' in link.lower():
                out['hepdata_link'] = link
                break
        out['inspire_id'] = strip_tags(fields['inspire id']).split('/')[0].strip()

    if 'e-print archives' in fields:
        ehtml = fields['e-print archives']
        etext = strip_tags(ehtml)
        m = ARXIV_ID_RE.search(etext)
        if m:
            out['arxiv_id'] = m.group(1)
        for link in all_hrefs(ehtml):
            if 'arxiv.org/abs/' in link:
                out['arxiv_abs_link'] = link
                break

    return out


def main():
    ap = argparse.ArgumentParser(description='Enrich STAR publication CSV with journal/HEPData links from STAR pages.')
    ap.add_argument('--input', default='data/publications_from_docx.csv')
    ap.add_argument('--output', default='data/publications_enriched.csv')
    ap.add_argument('--timeout', type=float, default=20.0)
    ap.add_argument('--sleep', type=float, default=0.15)
    args = ap.parse_args()

    with open(args.input, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))

    enriched = []
    for i, row in enumerate(rows, 1):
        enriched.append(enrich_row(row, args.timeout))
        if i % 10 == 0 or i == len(rows):
            print(f'processed {i}/{len(rows)}')
        time.sleep(args.sleep)

    fieldnames = [
        'index', 'title', 'star_url', 'journal_raw', 'journal', 'year',
        'reference', 'status', 'journal_link', 'journal_link_text',
        'hepdata_link', 'arxiv_id', 'arxiv_abs_link', 'inspire_id', 'fetch_error'
    ]
    with open(args.output, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(enriched)

    errors = sum(1 for r in enriched if r.get('fetch_error'))
    got_journal = sum(1 for r in enriched if r.get('journal_link'))
    got_hepdata = sum(1 for r in enriched if r.get('hepdata_link'))
    print(f'wrote {len(enriched)} rows to {args.output}')
    print(f'journal_link: {got_journal}, hepdata_link: {got_hepdata}, errors: {errors}')


if __name__ == '__main__':
    main()
