"""Retrieve nominated public project pages and retain inspectable evidence.

HTTP success establishes retrieval only, not maintenance, adoption, licensing,
or scientific validation. No downloaded content is executed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import urllib.request
import urllib.error

ROOT = Path(__file__).resolve().parent

class Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.skip = 0
        self.in_title = False
        self.title = []
        self.parts = []
        self.description = ''
        self.article = False
        self.article_parts = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag in ('script', 'style'):
            self.skip += 1
        if tag == 'title': self.in_title = True
        if tag == 'article': self.article = True
        if tag == 'meta' and attrs.get('name', '').lower() == 'description':
            self.description = attrs.get('content', '')
        if tag in ('p', 'br', 'div', 'h1', 'h2', 'h3', 'li', 'tr'):
            self.parts.append('\n')
            if self.article: self.article_parts.append('\n')

    def handle_endtag(self, tag):
        if tag in ('script', 'style'): self.skip = max(0, self.skip - 1)
        if tag == 'title': self.in_title = False
        if tag == 'article': self.article = False

    def handle_data(self, data):
        if not self.skip:
            self.parts.append(data)
            if self.article: self.article_parts.append(data)
            if self.in_title: self.title.append(data)

def retrieve(record):
    result = {'id': record['id'], 'name': record['name'], 'url': record['source_url'],
              'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    try:
        req = urllib.request.Request(record['source_url'], headers={'User-Agent': 'AstronomySoftwareSurvey/1.0 (public research catalogue; Python urllib)'})
        with urllib.request.urlopen(req, timeout=16) as response:
            content = response.read(4_000_001)
            result.update(http_status=response.status, final_url=response.url,
                          content_type=response.headers.get('Content-Type', ''), bytes=len(content),
                          truncated=len(content) > 4_000_000)
        parser = Text()
        parser.feed(content.decode('utf-8', errors='replace'))
        text = ''.join(parser.article_parts if 'github.com/' in record['source_url'] and parser.article_parts else parser.parts)
        text = re.sub(r'[ \t\r\f\v]+', ' ', text)
        text = re.sub(r'\n\s*\n+', '\n\n', text).strip()
        path = ROOT / 'source-extracts' / (record['id'] + '.txt')
        path.parent.mkdir(exist_ok=True)
        path.write_text(text)
        result.update(title=''.join(parser.title).strip(), description=parser.description,
                      text_sha256=hashlib.sha256(text.encode()).hexdigest(),
                      extract_path=str(path.relative_to(ROOT)), text_chars=len(text),
                      preview=text[:2400], retrieval='retrieved')
    except Exception as exc:
        result.update(retrieval='unresolved', error=f'{type(exc).__name__}: {exc}')
        if isinstance(exc, urllib.error.HTTPError): result['http_status'] = exc.code
    return result

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--retry-unresolved', action='store_true')
    args = ap.parse_args()
    records = json.loads((ROOT / 'catalogue.json').read_text())
    prior = {}
    if (ROOT / 'source-checks.json').exists():
        prior = {r['id']: r for r in json.loads((ROOT / 'source-checks.json').read_text())}
    todo = [r for r in records if r['id'] not in prior or r['source_url'] != prior[r['id']]['url'] or
            (args.retry_unresolved and prior[r['id']]['retrieval'] == 'unresolved')]
    with ThreadPoolExecutor(max_workers=10) as pool:
        futures = {pool.submit(retrieve, r): r for r in todo}
        for n, future in enumerate(as_completed(futures), 1):
            result = future.result()
            prior[result['id']] = result
            if n % 40 == 0:
                print(f'{n}/{len(todo)} checked', flush=True)
                (ROOT / 'source-checks.json').write_text(json.dumps(list(prior.values()), indent=2) + '\n')
    final = [prior[r['id']] for r in records]
    (ROOT / 'source-checks.json').write_text(json.dumps(final, indent=2) + '\n')
    print(json.dumps({'total': len(final), 'retrieved': sum(r['retrieval'] == 'retrieved' for r in final),
                      'unresolved': [{'id': r['id'], 'name': r['name'], 'error': r.get('error')} for r in final if r['retrieval'] != 'retrieved']}, indent=2))

if __name__ == '__main__': main()
