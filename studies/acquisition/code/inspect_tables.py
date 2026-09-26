#!/usr/bin/env python3
"""Inventory schemas, dimensions and IDs only; no scientific estimates or fits."""
import csv
from collections import Counter
import gzip
import json
from pathlib import Path
import re
from collect import ROOT

def table(path):
    lines=[s for s in path.read_text().splitlines() if s.strip() and not s.startswith('#')]
    if lines[0].startswith('VARNAMES:'):
        cols=lines[0].split()[1:]
        rows=[dict(zip(cols,s.split()[1:])) for s in lines[1:]]
        fmt='SNANA whitespace table (VARNAMES/SN), despite .csv suffix'
    else:
        reader=csv.DictReader(lines, delimiter=',' if ',' in lines[0] else ' ', skipinitialspace=True)
        rows=list(reader);cols=reader.fieldnames;fmt='CSV' if ',' in lines[0] else 'whitespace'
    result=dict(rows=len(rows),columns=cols,format=fmt)
    if 'IDSURVEY' in cols:result['survey_counts']=dict(Counter(r['IDSURVEY'] for r in rows))
    if 'CID' in cols:result['unique_CID']=len({r['CID'] for r in rows})
    return result

def fits_headers(path):
    result=[]
    with gzip.open(path,'rb') as f:
        while True:
            h={};first=f.read(2880)
            if not first:break
            block=first
            while True:
                cards=[block[i:i+80].decode('ascii') for i in range(0,2880,80)]
                for c in cards:
                    if c.startswith('END '):break
                    if c[8:10]=='= ':h[c[:8].strip()]=c[10:].split('/')[0].strip().strip("'").strip()
                if any(c.startswith('END ') for c in cards):break
                block=f.read(2880)
            result.append({k:v for k,v in h.items() if k in ('XTENSION','NAXIS1','NAXIS2','TFIELDS','EXTNAME')})
            axes=int(h.get('NAXIS','0'));size=abs(int(h['BITPIX']))//8 if axes else 0
            for i in range(1,axes+1):size*=int(h['NAXIS'+str(i)])
            size+=int(h.get('PCOUNT','0'));size*=int(h.get('GCOUNT','1'))
            if size:f.seek((size+2879)//2880*2880,1)
    return result

if __name__=='__main__':
    output={}
    for p in (ROOT/'sources/repos').glob('des-science*/4_DISTANCES_COVMAT/*HD.csv'):
        output[str(p.relative_to(ROOT))]=table(p)
    for p in (ROOT/'data/host_ages/chung2025').glob('table*.dat'):
        output[str(p.relative_to(ROOT))]=dict(data_rows=sum(bool(re.match(r'\s*\d+\s*&',s)) for s in p.read_text().splitlines()),format='LaTeX-style ampersand table; see original paper captions for HR conventions')
    for p in (ROOT/'sources/repos').glob('des-science*/0_DATA/**/*HEAD.FITS.gz'):
        output[str(p.relative_to(ROOT))]=dict(fits_extensions=fits_headers(p))
    for p in (ROOT/'data/des-diffimg').glob('*HEAD.FITS.gz'):
        output[str(p.relative_to(ROOT))]=dict(fits_extensions=fits_headers(p))
    for p in (ROOT/'data/dust').glob('GSWLC-*.dat.gz'):
        with gzip.open(p,'rt') as f:
            first=next(f); n=1+sum(bool(line.strip()) for line in f)
        output[str(p.relative_to(ROOT))]=dict(rows=n,columns_count=len(first.split()),schema='data/dust/GSWLC-2-column-description.pdf',missing_value=-99)
    for p in (ROOT/'sources/repos').glob('PantheonPlus*/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat'):
        lines=p.read_text().splitlines()
        output[str(p.relative_to(ROOT))]=dict(rows=len(lines)-1,unique_CID=len({l.split()[0] for l in lines[1:]}),columns=lines[0].split())
    p=ROOT/'sources/repos/CobayaSampler__bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt'
    output[str(p.relative_to(ROOT))]=dict(data_rows=sum(bool(l.strip()) and not l.startswith('#') for l in p.read_text().splitlines()))
    (ROOT/'catalog/table_inventory.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
