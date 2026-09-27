#!/usr/bin/env python3
"""Inspect cited public source/data bundles without fetching multi-GB R19 chains."""
from pathlib import Path
import argparse, datetime, hashlib, json, re, subprocess, tarfile, urllib.request, zipfile
from run import ROOT, DEFAULT_ARCHIVE, save, digest

BASE = 'https://raw.githubusercontent.com/benjaminrose/MC-Age/92713be96a89da991fe53bffcc596a5c0942fc37/'
URLS = {
    'C25_arxiv_source':'https://arxiv.org/src/2411.05299v2',
    'C26_arxiv_source':'https://arxiv.org/src/2605.21586v1',
    'MC_Age_repo':'https://api.github.com/repos/benjaminrose/MC-Age',
    'MC_Age_tree':'https://api.github.com/repos/benjaminrose/MC-Age/git/trees/master?recursive=1',
    'MC_Age_readme.rst':BASE+'readme.rst',
    'MC_Age_resources_readme.rst':BASE+'resources/readme.rst',
    'MC_Age_data_readme.rst':BASE+'data/readme.rst',
    'MC_Age_Zenodo_3875482.json':'https://zenodo.org/api/records/3875482',
    'C25_journal':'https://academic.oup.com/mnras/article/538/4/3340/8098234',
    'C26_journal':'https://academic.oup.com/mnras/article/551/3/stag1513/8771022',
}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh',action='store_true');ap.add_argument('--archive',type=Path,default=DEFAULT_ARCHIVE);args=ap.parse_args()
    out=ROOT/'.work/age-reconciliation/public-input-check';out.mkdir(parents=True,exist_ok=True);records=[]
    for name,url in URLS.items():
        f=out/name;r={'name':name,'url':url,'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
        try:
            if not f.exists() or args.refresh:
                # curl retry tolerates interrupted large arXiv source responses.
                if 'arxiv_source' in name:
                    subprocess.run(['curl','--silent','--show-error','--fail','--location','--retry','2','--max-time','45','--output',str(f)+'.partial',url],check=True)
                    Path(str(f)+'.partial').replace(f)
                else:
                    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'cosmology-revalidation scholarly data availability audit'}),timeout=30) as response:f.write_bytes(response.read())
                r['access']='downloaded_this_execution'
            else:r['access']='cached_from_direct_check_2026-09-27'
            r.update(status=200,bytes=f.stat().st_size,sha256=digest(f))
            if 'arxiv_source' in name:
                with tarfile.open(f) as t:
                    r['members']=t.getnames();links=[]
                    for m in t.getmembers():
                        if m.name.endswith('.tex'):links+=re.findall(r'https?://[^\s}\]]+',t.extractfile(m).read().decode(errors='replace'))
                    r['nonbibliographic_source_links']=[z for z in links if 'ui.adsabs' not in z and '#1' not in z and '#2' not in z]
            elif name=='MC_Age_tree':
                d=json.loads(f.read_text());r.update(commit_sha=d['sha'],file_count=len(d['tree']),truncated=d.get('truncated'),resource_paths=[x['path'] for x in d['tree'] if x['path'].startswith('resources/')])
            elif name=='MC_Age_repo':
                d=json.loads(f.read_text());r['repository_metadata']={k:d.get(k) for k in ['pushed_at','updated_at','default_branch','description','homepage']}
            elif 'Zenodo' in name:
                d=json.loads(f.read_text());r.update(title=d['metadata']['title'],publication_date=d['metadata'].get('publication_date'),files=[{'key':z['key'],'size':z['size']} for z in d.get('files',[])],description=d['metadata'].get('description'))
            elif name.endswith('.rst'):r['links']=re.findall(r'https?://[^\s`]+',f.read_text())
        except Exception as e:r['error']=str(e)
        records.append(r)
    f=args.archive/'data/host_ages/chung2025/staf497_supplemental_file.zip'
    with zipfile.ZipFile(f) as z:records.append({'name':'C25_frozen_journal_supplement','sha256':digest(f),'members':z.namelist(),'source':'Archived official journal supplement; not live refetched.'})
    for name in ['chung2025-published.txt','chung2026-published.txt']:
        f=args.archive/'papers/text'/name;records.append({'name':name,'source':'Archived published text inspected for availability statements','sha256':digest(f)})
    save(ROOT/'studies/host_ages/results/age_reconciliation/public-input-check.json',{'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code_sha256':digest(Path(__file__)),'records':records,'findings':['Official C25 supplement contains only table1.dat and table2.dat; arXiv source contains figures/TeX, not posterior chains.','The cited MC-Age tree documentation links historical R19 chains via Zenodo3875482, not a C25 reprocessed posterior bundle.','C26v1 source contains manuscript/figures and no mock code/age-PDF arrays. Published C26 availability statement says no new data associated.','Journal live endpoints returned403; published text and archived supplement inspected. This bounded non-location is not proof of absence.'],'scope':'Direct cited repository, source archives, linked Zenodo, journal pages and archived official supplement. No author contact; no unlinked repository completeness claim.'})

if __name__=='__main__':main()
