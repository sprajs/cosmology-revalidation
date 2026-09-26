"""Read-only refresh of primary metadata and newly needed source documents."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, urllib.request

OUT = Path('runs/mapping/sources-2026-09-20')
OUT.mkdir(parents=True, exist_ok=True)
URLS = {
    'behroozi2013.pdf': 'https://arxiv.org/pdf/1207.6105',
    'madau2014.pdf': 'https://arxiv.org/pdf/1403.0007',
    'wiseman2026-abs.html': 'https://arxiv.org/abs/2601.13785',
    'chung2026-abs.html': 'https://arxiv.org/abs/2605.21586',
    'murakami2026-abs.html': 'https://arxiv.org/abs/2604.16597',
    'park2026-abs.html': 'https://arxiv.org/abs/2605.12596',
    'son2025-abs.html': 'https://arxiv.org/abs/2510.13121',
    'titan-dr1.html': 'https://titan-snia.github.io/dr1.html',
    'wiseman-repo-head.json': 'https://api.github.com/repos/wisemanp/des_sn_hosts/commits?per_page=1',
    'wiseman-repo-releases.json': 'https://api.github.com/repos/wisemanp/des_sn_hosts/releases',
    'park-inspire.json': 'https://inspirehep.net/api/literature?q=arxiv:2605.12596',
}
for name, url in URLS.items():
    path = OUT / name
    if path.exists():
        continue
    row = {'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'scientific-local-reproduction/1.0'}), timeout=45) as r:
            raw = r.read()
            row.update(status=r.status, final_url=r.url, content_type=r.headers.get('Content-Type'))
        path.write_bytes(raw)
        row.update(path=str(path), size=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    except Exception as exc:
        row['error'] = repr(exc)
    with (OUT / 'acquisition.jsonl').open('a') as f:
        f.write(json.dumps(row) + '\n')
    print(name, row.get('status', row.get('error')), row.get('size'))
