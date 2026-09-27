#!/usr/bin/env python3
"""Recheck public DES/Dovekie assets without altering the historical archive."""
import argparse, json, hashlib, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[4]
DEFAULT=Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def get(url):
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'cosmology-revalidation-public-data-audit'}),timeout=40) as r:return r.read()
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--archive',type=Path,default=DEFAULT);args=pa.parse_args()
    old=args.archive/'catalog/lfs_hydration.json'
    objects=[o for o in json.loads(old.read_text())['objects'] if o['repository']=='des-science__DES-SN5YR' and o['status']=='unresolved']
    dest=ROOT/'data/age-recovery/survey-check';dest.mkdir(parents=True,exist_ok=True)
    def check(o):
        row={k:o[k] for k in ['path','url','sha256','bytes','revision']}
        try:
            b=get(o['url']);h=hashlib.sha256(b).hexdigest()
            row.update(http_status=200,received_bytes=len(b),received_sha256=h,verified=h==o['sha256'] and len(b)==o['bytes'])
            if row['verified']:(dest/o['sha256']).write_bytes(b)
        except urllib.error.HTTPError as e:row.update(http_status=e.code,verified=False)
        except Exception as e:row.update(error=type(e).__name__+': '+str(e),verified=False)
        return row
    with ThreadPoolExecutor(max_workers=6) as ex:rows=list(ex.map(check,objects))
    repositories=[]
    for slug in ['des-science/DES-SN5YR','bap37/Dovekie']:
        url=f'https://api.github.com/repos/{slug}/commits/main';b=get(url);d=json.loads(b)
        (dest/(slug.replace('/','__')+'-commit.json')).write_bytes(b)
        repositories.append(dict(repository=slug,url=url,commit=d['sha'],date=d['commit']['committer']['date'],response_sha256=hashlib.sha256(b).hexdigest()))
    report=dict(time_utc=datetime.now(timezone.utc).isoformat(),scope='Current public endpoint check, not proof these assets are unavailable elsewhere; recheck pinned previously unresolved objects and repository heads.',repositories=repositories,objects_checked=len(rows),verified=sum(r['verified'] for r in rows),http_status_counts={str(s):sum(r.get('http_status')==s for r in rows) for s in set(r.get('http_status') for r in rows)},objects=rows,source_inventory_sha256=sha(old),code_sha256=sha(Path(__file__)))
    p=ROOT/'studies/host_ages/results/age_recovery/survey-availability.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='objects'},indent=2))
if __name__=='__main__':main()
