"""Render the source-linked survey and local searchable catalogue (stdlib only)."""
from collections import Counter
from pathlib import Path
import html
import json
import re

ROOT = Path(__file__).resolve().parent
records = json.loads((ROOT / 'catalogue.json').read_text())
checks = {r['id']: r for r in json.loads((ROOT / 'source-checks.json').read_text())}
registry = json.loads((ROOT / 'ascl-index.json').read_text())
provenance = json.loads((ROOT / 'ascl-index-provenance.json').read_text())
assert len(records) == len(checks) == 528
assert len(registry) == len({r['record_url'] for r in registry}) == 4260
assert Counter(r['registry_status'] for r in registry) == {'registered': 4105, 'submitted': 155}

for r in records:
    c = checks[r['id']]
    assert c['url'] == r['source_url']
    if c.get('http_status') != 200:
        status = 'Retrieval unresolved'
    elif c.get('text_chars', 0) < 450:
        status = 'Short page, metadata or access shell'
    else:
        status = 'Source text retrieved'
    r['evidence_status'] = status
    r['retrieval'] = {key: c.get(key) for key in ('retrieved_utc', 'http_status', 'final_url', 'title', 'text_chars', 'extract_path', 'text_sha256', 'error') if c.get(key) is not None}
    r['maintenance_audit'] = 'Not systematically audited; see any explicit lineage notes'

(ROOT / 'catalogue.json').write_text(json.dumps(records, indent=2, ensure_ascii=False) + '\n')
counts = Counter(r['domain'] for r in records)
quality = Counter(r['evidence_status'] for r in records)

def md(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')

out = ['# Curated astronomy and physics software catalogue', '',
       '528 entries across 26 groups. Compiled 30 September 2026, Europe/London.', '',
       'This is a capability inventory, including legacy and unmaintained software. An entry can be a package, a family, a suite, a standard or a data service. Entries are not independent numerical engines. Language labels identify principal implementation or interface families; they are not source-percentage audits. A retrieved source is not evidence of present maintenance, buildability, adoption or numerical correctness.', '',
       'Read [the analysis](REPORT.md), use [the searchable catalogue](catalogue.html), or inspect [retrieval evidence and limitations](EVIDENCE.md). The separate ASCL index has 4,260 title-and-link records, including 155 submitted records. The two collections overlap.', '',
       '| Domain | Entries |', '|---|---:|']
out += [f'| {md(k)} | {v} |' for k, v in counts.items()]
for domain in counts:
    out += ['', '## ' + domain, '', '| Software / source | Language family | Capability | Lineage, integration and evidence notes |', '|---|---|---|---|']
    for r in records:
        if r['domain'] != domain:
            continue
        extra = f" {r['evidence_status']}."
        if r.get('original_source_url'):
            extra += f" [Original nominated site]({r['original_source_url']}); source link uses ASCL."
        out.append(f"| [{md(r['name'])}]({r['source_url']}) · {r['id']} | {md(r['language'])} | {md(r['capability'])} | {md(r['note'])}{extra} |")
(ROOT / 'CATALOGUE.md').write_text('\n'.join(out) + '\n')

evidence = ['# Evidence, coverage and reproducibility', '',
    'Collected 30 September 2026 in Europe/London. UTC timestamps in the logs begin on 29 September because the collection crossed the local midnight boundary.', '',
    '## What was checked', '',
    f'- Curated entries: **{len(records)}**, in **{len(counts)}** groups; every entry has a nominated source and a retrieval attempt.',
    f'- HTTP 200 responses: **{sum(c.get("http_status") == 200 for c in checks.values())}**. This is a transport-level count, not a count of verified scientific implementations.',
    f'- Extracted text of at least 450 characters: **{quality["Source text retrieved"]}**. This threshold is only a content triage aid; it is not a scientific or identity-verification score.',
    f'- Short pages, metadata, redirects, access challenges or JavaScript shells: **{quality["Short page, metadata or access shell"]}**.',
    f'- Unresolved nominated-source retrievals: **{quality["Retrieval unresolved"]}**. Retrieval failure does not imply that software is unavailable elsewhere.',
    f'- ASCL breadth index: **{len(registry):,}** records from **{len(provenance["pages"])}** browse pages: **4,105 registered** and **155 submitted**. Unique record URLs and registry totals agree.', '',
    'The review used project repositories, author pages, observatory documentation, software papers and registry records. Searches covered all 26 catalogue groups, including explicit legacy and alternative-language sweeps. Capability and language summaries were curated; repository contents were not cloned or exhaustively audited. Individual packages were not installed, tested or benchmarked. Release dates and activity metrics were not systematically collected, so this is not a maintenance ranking.', '',
    '## Unresolved direct retrievals', '',
    '| Entry | Nominated source | Retrieval result |', '|---|---|---|']
for r in records:
    if r['evidence_status'] == 'Retrieval unresolved':
        evidence.append(f"| {md(r['name'])} | [Source]({r['source_url']}) | {md(checks[r['id']].get('error', 'Unresolved'))} |")
evidence += ['', 'The Wilson–Devinney lineage also has [ASCL record 2004.004](https://ascl.net/2004.004). SSE/BSE are also referenced by coupling frameworks such as [AMUSE](https://github.com/amusecode/amuse); this does not resolve the old author-site retrieval.', '',
    '## Short, redirected or limited source pages', '',
    'These entries remain useful, but the automatic source-text snapshot is limited. Several have substantial documentation elsewhere, or were supported by additional web-search results. This table avoids treating a repository title or an HTTP 200 access challenge as full documentation.', '',
    '| Entry | Extracted characters | Landing-page title |', '|---|---:|---|']
for r in records:
    if r['evidence_status'] == 'Short page, metadata or access shell':
        c = checks[r['id']]
        evidence.append(f"| [{md(r['name'])}]({r['source_url']}) | {c.get('text_chars', 0)} | {md(c.get('title') or '(no usable title)')} |")
evidence += ['', '## Registry fallbacks', '',
    'Where an original site was inaccessible, a stable ASCL record was used for identity and scope. Both URLs are retained. A registry record does not guarantee that source code or a working build remains available.', '',
    '| Entry | Registry source | Original nominated URL |', '|---|---|---|']
for r in records:
    if r.get('original_source_url'):
        evidence.append(f"| {md(r['name'])} | [ASCL]({r['source_url']}) | [Original]({r['original_source_url']}) |")
evidence += ['', '## Identity and lineage findings', '',
    '- GLEE and WSLAP+ are documented using author papers. Their implementation language and public source-distribution status are explicitly unconfirmed here.',
    '- The proprietary-language identification is an inference: IDL is the strongest match. Mark Sullivan’s personal use was not established by a sufficiently direct primary source.',
    '- IRAF community maintenance, ESO MIDAS releases, the IDL-library archive migration and AIPS maintenance show why age alone cannot determine lifecycle status.',
    '- GRChombo/GRTeclyn, Tudat/tudatpy, COMPASS, SoFiA 2, Gradus and faer have successor, archive or migration information in their landing pages. Those relationships appear in the catalogue notes.',
    '- Name collisions require domain identity: ISIS optical subtraction, ISIS X-ray spectroscopy and USGS ISIS are different systems; CLASS cosmology and GILDAS CLASS are different; pulsar PINT is not the Python units package.',
    '- Wrappers, forks, suites and shared engines must not be counted as independent implementations. In particular, Python-facing software often depends on native C/C++/Fortran kernels.', '',
    '## Files and reproducibility', '',
    '- `catalogue.json`: curated entries, source links, language labels, capability notes and source-retrieval metadata.',
    '- `source-checks.json`: retrieval attempts, final URLs, UTC timestamps, response metadata and extracted-text checksums.',
    '- `source-extracts/`: local research snapshots of extracted source-page text, for traceability; these are not software source distributions or an independently publishable source collection.',
    '- `ascl-index.json`: registry titles and record links. Abstracts and source code are not reproduced in this index.',
    '- `ascl-index-provenance.json`: browse-page URLs, timestamps, page checksums, counts and collection scope.',
    '- `research-searches.json`: retained thematic search results. Search snippets are leads; primary project sources are the main catalogue references.',
    '- `REPORT.md`, `CATALOGUE.md`, `OVERLAP-MAP.md`, `overlap-map.json`, `EVIDENCE.md`, `catalogue.html`: reading and exploration views.', '',
    'Rebuild the curated catalogue and presentation from this directory:', '',
    '```sh', 'python build_catalogue.py', 'python verify_sources.py', 'python build_overlap_map.py', 'python render_catalogue.py', '```', '',
    '`verify_sources.py` reuses cached checks when the nominated URL is unchanged; use `--retry-unresolved` to retry failures. It does not refresh every cached page. `collect_ascl_index.py` performs a fresh network collection and can change the index. The rendered catalogue intentionally validates the snapshot counts; update its assertions deliberately when making a later survey edition.', '',
    '## Limits of the sweep', '',
    'The ASCL index was validated for pagination completeness, record-URL uniqueness and registry counts, not for the correctness of every title or scientific claim. The ASCL and curated lists overlap and their totals must not be added. Submitted registry entries are explicitly labelled. Registry search operates on titles and names, not full abstracts.', '',
    'Coverage is broad but cannot be a census of private collaboration code, unpublished scripts, every instrument recipe or every proprietary extension. General physics is represented by major frameworks and capability families rather than every computational-physics package. License compatibility, source availability, build reproduction, algorithm ancestry and scientific validation remain per-component follow-up work. The Rust module map is an architectural proposal, not an upstream consensus or an implemented product.', '']
(ROOT / 'EVIDENCE.md').write_text('\n'.join(evidence))

data = json.dumps({'curated': records, 'registry': registry}, ensure_ascii=False).replace('</', '<\\/')
template = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Astronomy & physics software atlas</title>
<style>
:root{color-scheme:light;--ink:#182d3d;--muted:#536673;--line:#d6e0e4;--paper:#f3f6f6;--accent:#086b66;--pale:#e3f2ee}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 system-ui,sans-serif}a{color:#086966;text-underline-offset:3px}header{background:#142f3d;color:#f4fafb;padding:42px max(24px,calc((100vw - 1280px)/2)) 34px}header a{color:#9de1d2}.eyebrow{font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:#a9c6d1;margin:0 0 9px}h1{font:600 clamp(28px,4vw,44px)/1.15 system-ui;margin:0 0 14px;letter-spacing:-.035em}header p{max-width:920px;margin:12px 0;color:#d8e6ea}.stats{display:flex;flex-wrap:wrap;gap:12px 36px;margin-top:26px}.stat strong{font-size:27px;color:#fff;display:block}.stat span{font-size:13px;color:#bfd3dc}nav{display:flex;flex-wrap:wrap;gap:20px;margin-top:25px;font-size:14px}main{max-width:1328px;margin:auto;padding:28px 24px 50px}.tabs{display:flex;gap:10px;margin-bottom:16px}.tabs button{border:1px solid var(--line);background:white;padding:11px 18px;border-radius:7px;font:600 15px system-ui;cursor:pointer;color:var(--ink)}.tabs button[aria-pressed=true]{background:var(--accent);border-color:var(--accent);color:white}.controls{padding:20px;background:white;border:1px solid var(--line);border-radius:9px;display:grid;grid-template-columns:minmax(260px,2.3fr) minmax(220px,1.5fr) minmax(155px,1fr) auto;gap:14px;align-items:end}label{display:block;font-size:12px;font-weight:650;color:var(--muted);margin-bottom:6px}input,select{border:1px solid #aebfc7;background:white;color:var(--ink);border-radius:5px;font:15px system-ui;padding:10px;width:100%;min-height:43px}input:focus,select:focus,button:focus-visible{outline:3px solid #87caba;outline-offset:2px}button.action{background:white;color:var(--ink);border:1px solid #aebfc7;padding:10px 14px;border-radius:5px;cursor:pointer;min-height:43px;font:14px system-ui}button:disabled{opacity:.4;cursor:default}.context{color:var(--muted);font-size:14px;max-width:1100px}.resultbar{display:flex;align-items:center;justify-content:space-between;gap:16px;margin:25px 0 12px;flex-wrap:wrap}#count{font-weight:650}.pager{display:flex;align-items:center;gap:12px;font-size:13px}#results{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.card{background:white;border:1px solid var(--line);border-radius:8px;padding:20px;min-width:0}.card h2{font-size:20px;line-height:1.3;margin:0 0 8px}.card h2 a{text-decoration:none}.card h2 a:hover{text-decoration:underline}.domain{font-size:11px;letter-spacing:.045em;text-transform:uppercase;color:var(--muted);margin-bottom:10px}.language{display:inline-block;background:var(--pale);padding:3px 8px;border-radius:4px;color:#225e58;font-size:12px;margin-bottom:6px}.card p{font-size:14px;margin:10px 0}.note{color:var(--muted)}.evidence{border-top:1px solid #e5ecee;padding-top:10px;margin-top:14px;font-size:11px;color:var(--muted);overflow-wrap:anywhere}.empty{grid-column:1/-1;background:white;border:1px dashed #aebfc7;padding:32px;border-radius:8px}footer{color:var(--muted);font-size:12px;margin-top:28px;max-width:1050px}.bottom{justify-content:flex-end}#registry-note{display:none}select:disabled{background:#edf1f2;color:#6e7a80}@media(max-width:1000px){.controls{grid-template-columns:1fr 1fr}.controls>div:first-child{grid-column:1/-1}}@media(max-width:670px){header{padding:30px 20px}main{padding:22px 16px}#results{grid-template-columns:1fr}.controls{grid-template-columns:1fr}.controls>div:first-child{grid-column:auto}.tabs{flex-wrap:wrap}.stats{gap:12px 25px}.card{padding:17px}}@media print{.controls,.tabs,.pager{display:none}header{background:white;color:black}header p,header a,.stat span,.stat strong{color:black}#results{display:block}.card{break-inside:avoid;margin-bottom:12px}}
</style></head><body>
<header><p class="eyebrow">Research inventory · 30 September 2026</p><h1>Astronomy & physics software atlas</h1><p>A capability map for a unified Rust suite — from observatory pipelines and legacy task systems to cosmology, lensing, relativity and general numerical physics.</p><div class="stats"><div class="stat"><strong>528</strong><span>curated entries</span></div><div class="stat"><strong>26</strong><span>capability groups</span></div><div class="stat"><strong>4,260</strong><span>ASCL discovery records</span></div></div><nav><a href="OVERLAP-MAP.md">Where packages overlap</a><a href="REPORT.md">Survey analysis</a><a href="CATALOGUE.md">Full reading list</a><a href="EVIDENCE.md">Evidence & limitations</a><a href="catalogue.json" download>Curated JSON</a><a href="ascl-index.json" download>ASCL JSON</a></nav></header>
<main><div class="tabs" aria-label="Catalogue collection"><button id="curated-tab" aria-pressed="true">Curated capability catalogue</button><button id="registry-tab" aria-pressed="false">Full ASCL discovery index</button></div>
<div class="controls"><div><label for="search">Search names, capabilities or notes</label><input type="search" id="search" placeholder="Try lensing, IDL, Einstein, Fortran, calibration…"></div><div><label for="domain">Capability group</label><select id="domain"><option value="">All 26 groups</option></select></div><div><label for="language">Language family</label><select id="language"><option value="">All languages</option></select></div><button class="action" id="reset">Reset filters</button></div>
<p class="context" id="curated-note">Includes packages, suites, related families, standards and essential data services. Legacy software remains in scope. Source retrieval does not establish current maintenance or scientific correctness.</p>
<p class="context" id="registry-note">4,105 registered records and 155 submitted records. This separate index searches titles and names only; individual records have not all received the curated review. It overlaps with the curated catalogue. Language and domain filters are unavailable for this index.</p>
<div class="resultbar"><div id="count" role="status" aria-live="polite"></div><div class="pager"><button class="action previous">Previous</button><span class="page-label"></span><button class="action next">Next</button></div></div><div id="results"></div><div class="resultbar bottom"><div class="pager"><button class="action previous">Previous</button><span class="page-label"></span><button class="action next">Next</button></div></div>
<footer>This survey records discovery and capability scope, including migration and legacy notes. It is not a universal census or an audit of licenses, adoption, buildability, numerical accuracy or performance. The proposed Rust architecture is the survey author's synthesis. The two collection totals must not be added as independent software counts. All search and filtering happens locally in this page.</footer></main>
<script id="inventory-data" type="application/json">__DATA__</script><script>
'use strict';
const data=JSON.parse(document.getElementById('inventory-data').textContent);
const $=id=>document.getElementById(id);
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let mode='curated',page=0;const pageSize=60;
const domains=[...new Set(data.curated.map(r=>r.domain))];
for(const d of domains){const o=document.createElement('option');o.value=d;o.textContent=d.slice(3);$('domain').append(o)}
const languages=['Python','Fortran','C++','C','Julia','IDL','Rust','R','Java','Perl','Wolfram/Mathematica','Maple','S-Lang','Yorick','SPP/CL','Other / mixed'];
function langMatch(r,l){const s=r.language;if(l==='Other / mixed')return /mixed|not confirmed|consult|standard|database|service|format|specification|tcl|scheme|javascript/i.test(s);if(l==='Wolfram/Mathematica')return /Wolfram|Mathematica/i.test(s);if(l==='SPP/CL')return /SPP|\bCL\b/.test(s);if(l==='C')return /(?:^|[; /])C(?:$|[; /])/.test(s);if(l==='R')return /(?:^|[; ])R(?:$|[; ])/.test(s);if(l==='Java')return /\bJava\b/.test(s);return s.toLowerCase().includes(l.toLowerCase())}
for(const l of languages){const o=document.createElement('option');o.value=l;o.textContent=l;$('language').append(o)}
function results(){const terms=$('search').value.toLowerCase().trim().split(/\s+/).filter(Boolean);return data[mode].filter(r=>{const s=(mode==='curated'?[r.name,r.language,r.capability,r.note,r.domain,r.id]:[r.name,r.title,r.ascl_id,r.registry_status]).join(' ').toLowerCase();return terms.every(t=>s.includes(t))&&(mode==='registry'||((!$('domain').value||r.domain===$('domain').value)&&(!$('language').value||langMatch(r,$('language').value))))})}
function card(r){if(mode==='registry')return `<article class="card"><div class="domain">ASCL · ${esc(r.registry_status)}</div><h2><a href="${esc(r.record_url)}" target="_blank" rel="noopener noreferrer">${esc(r.name)}</a></h2><p>${esc(r.title)}</p><div class="evidence">${esc(r.ascl_id||'Submitted record; registered ASCL identifier not assigned in listing')} · Title-and-link discovery record</div></article>`;return `<article class="card"><div class="domain">${esc(r.domain.slice(3))}</div><h2><a href="${esc(r.source_url)}" target="_blank" rel="noopener noreferrer">${esc(r.name)}</a></h2><span class="language">${esc(r.language)}</span><p>${esc(r.capability)}</p><p class="note">${esc(r.note)}</p><div class="evidence">${esc(r.id)} · ${esc(r.evidence_status)}<br>${esc(r.source_kind)}${r.original_source_url?` · <a href="${esc(r.original_source_url)}" target="_blank" rel="noopener noreferrer">Original nominated site</a>`:''}</div></article>`}
function render(){const rows=results(),pages=Math.max(1,Math.ceil(rows.length/pageSize));page=Math.min(page,pages-1);$('count').textContent=`${rows.length.toLocaleString()} of ${data[mode].length.toLocaleString()} ${mode==='curated'?'curated entries':'registry records'}`;$('results').innerHTML=rows.length?rows.slice(page*pageSize,(page+1)*pageSize).map(card).join(''):'<div class="empty">No entries match these filters. Try fewer search terms or reset the filters.</div>';document.querySelectorAll('.page-label').forEach(e=>e.textContent=`Page ${page+1} of ${pages}`);document.querySelectorAll('.previous').forEach(e=>e.disabled=page===0);document.querySelectorAll('.next').forEach(e=>e.disabled=page+1>=pages)}
function switchMode(m){mode=m;page=0;$('curated-tab').setAttribute('aria-pressed',m==='curated');$('registry-tab').setAttribute('aria-pressed',m==='registry');$('domain').disabled=m==='registry';$('language').disabled=m==='registry';$('curated-note').style.display=m==='curated'?'block':'none';$('registry-note').style.display=m==='registry'?'block':'none';$('search').placeholder=m==='curated'?'Try lensing, IDL, Einstein, Fortran, calibration…':'Search ASCL names and titles…';render()}
$('curated-tab').onclick=()=>switchMode('curated');$('registry-tab').onclick=()=>switchMode('registry');
for(const id of ['search','domain','language'])$(id).addEventListener(id==='search'?'input':'change',()=>{page=0;render()});
$('reset').onclick=()=>{$('search').value='';$('domain').value='';$('language').value='';page=0;render()};
document.querySelectorAll('.previous').forEach(e=>e.onclick=()=>{page--;render();$('count').scrollIntoView({block:'start'})});document.querySelectorAll('.next').forEach(e=>e.onclick=()=>{page++;render();$('count').scrollIntoView({block:'start'})});
render();
</script></body></html>'''
(ROOT / 'catalogue.html').write_text(template.replace('__DATA__', data))
print(json.dumps({'curated': len(records), 'domains': len(counts), 'ascl': len(registry), 'source_evidence': dict(quality)}, indent=2))
