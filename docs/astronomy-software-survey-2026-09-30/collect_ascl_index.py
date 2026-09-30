"""Collect ASCL record titles and links, not copied abstracts or source code."""
from pathlib import Path
import hashlib
import html
import json
import re
import time
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
rows, pages = [], []
expected = None
page = 1
while True:
    url = f'https://ascl.net/code/all/page/{page}/limit/100/order/title/listmode/full/dir/asc'
    data = urllib.request.urlopen(url, timeout=30).read()
    body = data.decode('utf-8')
    total = re.search(r'Results\s+(\d+)-(\d+)\s+of\s+(\d+)\s+\((\d+)\s+ASCL,\s+(\d+)\s+submitted\)', body)
    assert total, f'Pagination metadata missing on page {page}'
    start, end, count, registered, submitted = map(int, total.groups())
    if expected is None: expected = (count, registered, submitted)
    assert expected == (count, registered, submitted), 'Registry changed during collection; rerun required'
    matches = re.findall(r'<div class="item">(.*?)(?=<div class="item">|<div class="ascl_pagination">)', body, re.S)
    n = 0
    for chunk in matches:
        title = re.search(r'<span class="title">\s*<a href="([^"]+)">(.*?)</a>', chunk, re.S)
        assert title, f'Title missing in page {page}'
        ascl_id = re.search(r'\[ascl:(\d{4}\.\d{3})\]', chunk)
        label = html.unescape(re.sub('<[^>]+>', '', title.group(2))).strip()
        rows.append({'name': label.split(':', 1)[0], 'title': label,
                     'ascl_id': ascl_id.group(1) if ascl_id else None,
                     'record_url': 'https://ascl.net' + title.group(1),
                     'registry_status': 'registered' if ascl_id else 'submitted',
                     'source_page': url})
        n += 1
    # ASCL's final page reports the nominal page end (4300), not the actual
    # last record (4260). Validate against the published total instead.
    assert n == min(end, count) - start + 1, (page, n, start, end)
    pages.append({'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                  'sha256': hashlib.sha256(data).hexdigest(), 'records': n, 'bytes': len(data)})
    if page % 5 == 0: print(f'{len(rows)}/{count} ASCL titles indexed', flush=True)
    if end >= count: break
    page += 1
    time.sleep(0.15)

assert len(rows) == expected[0]
assert len({r['record_url'] for r in rows}) == len(rows), 'Duplicate or shifted page boundaries'
assert sum(r['registry_status'] == 'registered' for r in rows) == expected[1]
assert sum(r['registry_status'] == 'submitted' for r in rows) == expected[2]
(ROOT / 'ascl-index.json').write_text(json.dumps(rows, indent=2, ensure_ascii=False) + '\n')
(ROOT / 'ascl-index-provenance.json').write_text(json.dumps({'expected': {'total': expected[0], 'registered': expected[1], 'submitted': expected[2]}, 'pages': pages, 'scope': 'All browse records at collection; titles and links only; no adoption or maintenance assessment'}, indent=2) + '\n')
print(json.dumps({'indexed':len(rows), 'registered':expected[1], 'submitted':expected[2], 'pages':len(pages)}))
