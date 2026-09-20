#!/usr/bin/env python3
"""Build the human reading index; scientific findings are curated in docs/."""
import json
from pathlib import Path
from collect import ROOT, now

CORE={'2411.05299','2510.13121','2601.13785','2605.21586'}
FOCUSED={'2605.12596','2604.16597','2606.09650','2511.07517','2506.05471','2401.02929','2401.02945','2406.05046','2112.03863','2202.04077','2503.14738','1804.05850','2406.05051','1908.10375','2207.05583'}
EXTRA={
 '2411.05299':['chung2025-published.pdf','chung2025-correction-published.pdf'],
 '2510.13121':['son2025-published.pdf'],
 '2601.13785':['wiseman2026-published.pdf'],
 '2605.21586':['chung2026-published.pdf'],
 '2406.05051':['popovic2024-dust-published.pdf'],
}

def build():
    papers=json.loads((ROOT/'catalog/papers.json').read_text())
    records=[]
    lines=['# Collected literature','',f'Collection date: 2026-09-20. Generated inventory: {now()}.','',
      'Core papers were read for their claims, methods and data dependencies. Related papers received focused section-level reading; other references were collected and screened for relevance. This is not a completed methodological peer review of every paper.','',
      '| Paper | Saved PDFs | Reading level |','|---|---|---|']
    for p in papers:
        aid=p['arxiv_id'];m=p['metadata'];title=m['title'][0]
        paths=p['pdf_paths']+[f'papers/pdf/{n}' for n in EXTRA.get(aid,[]) if (ROOT/'papers/pdf'/n).exists()]
        level='Core: claims, methods and data-dependency reading' if aid in CORE else ('Focused: selected sections/abstracts and data-availability review; full review remains' if aid in FOCUSED else 'Supporting: metadata/abstract and relevance screening; detailed review remains')
        records.append(dict(arxiv_id=aid,title=title,reading_level=level,pdf_paths=paths,scientific_reproduction=False))
        refs=', '.join(f'[{Path(x).name}](../{x})' for x in paths) or 'No arXiv PDF; see linked published copy / acquisition log'
        lines.append(f'| [{aid}: {title}](https://arxiv.org/abs/{aid}) | {refs} | {level} |')
    lines += ['', '## Additional primary documents','',
      '- [Childress, Wolf & Zahid 2014, *Ages of Type Ia supernovae over cosmic time*](https://doi.org/10.1093/mnras/stu1892): [published PDF](../papers/pdf/childress2014-published.pdf). Focused reading of DTD/SFH role; no model executed.',
      '- [Correction to Chung Paper I](https://doi.org/10.1093/mnras/stag1210): [published correction](../papers/pdf/chung2025-correction-published.pdf). Correction text read; author claim of unchanged quantitative results not independently tested.',
      '', '## Version notes','',
      '- The Chung 2026 published reply adds a cumulative-redshift/mock test relative to arXiv v1. Use the published PDF for the latest result.',
      '- Park Paper III has a verified journal DOI (10.1093/mnras/stag935), but the local full text is arXiv v1; equivalence to the journal version has not been checked.',
      '- Full author lists, abstracts and arXiv metadata are in papers.json; references.bib is an arXiv-oriented export, supplemented by additional-references.bib.',
      '- An unversioned arXiv endpoint was sometimes required after versioned URLs returned HTTP errors. The resolved version comes from saved landing metadata; the acquisition log records the exact requested URL.',
      '- Rendered pages in papers/rendered/ were inspected for W26 Figure 1, Son Figure 2 and the published Chung Figure 1. No points were digitized.', '']
    (ROOT/'catalog/literature.md').write_text('\n'.join(lines))
    (ROOT/'catalog/reading-ledger.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'catalog/additional-references.bib').write_text('''@article{Childress2014,
  author = {Childress, Michael J. and Wolf, Christian and Zahid, H. Jabran},
  title = {Ages of Type Ia supernovae over cosmic time},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year = {2014}, volume = {445}, pages = {1898--1911}, doi = {10.1093/mnras/stu1892}
}

@article{Chung2026Correction,
  title = {Correction to: Strong progenitor age bias in supernova cosmology -- I},
  author = {Chung, Chul and Park, Seunghyun and Son, Junhyuk and Cho, Hyejeon and Lee, Young-Wook},
  journal = {Monthly Notices of the Royal Astronomical Society}, year = {2026},
  volume = {550}, number = {2}, doi = {10.1093/mnras/stag1210}
}
''')
    print('Indexed',len(records),'arXiv works plus additional primary documents')

if __name__=='__main__':build()
