#!/usr/bin/env python3
"""Predeclared eight-object archive and duplicate baseline engineering gates."""
import collections
import decimal
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'runs/research_2026_09_26/raisin_simulation_timing_pilot'
BASE = OUT / 'fits_authorgrid/baseline'
COPY = OUT / 'fits_authorgrid/baseline_copy'
ARCH = ROOT / 'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/nir.FITRES.gz'
LC = ROOT / 'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz'
CIDS = [str(i) for i in range(1,9)]


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parse_fit(path):
    rows={}; strings={}
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        names=None
        for l in f:
            if l.startswith('VARNAMES:'):names=l.split()[1:]
            if l.startswith('SN:'):
                v=l.split()[1:];r=dict(zip(names,v));
                if r['CID'] in CIDS: rows[r['CID']]=r;strings[r['CID']]=l.strip()
    return rows,strings


def parse_lc(path):
    d=collections.defaultdict(list)
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        names=None
        for l in f:
            if l.startswith('VARNAMES:'):names=l.split()[1:]
            if l.startswith('OBS:'):
                r=dict(zip(names,l.split()[1:]))
                if r['CID'] in CIDS and r['DATAFLAG']=='1':d[r['CID']].append(r)
    return d


def f32(x): return struct.unpack('<f',struct.pack('<f',float(x)))[0]


def f32_ulp(x):
    v=f32(x); bits=struct.unpack('<I',struct.pack('<f',v))[0]
    if bits==0: return 2**-149
    return abs(struct.unpack('<f',struct.pack('<I',bits+1))[0]-v)


def half_unit(token):
    x=decimal.Decimal(token)
    return float(decimal.Decimal(5).scaleb(x.as_tuple().exponent-1))


def gate():
    a,astring=parse_fit(ARCH);b,bstring=parse_fit(BASE/'fit.FITRES.TEXT');c,cstring=parse_fit(COPY/'fit.FITRES.TEXT')
    al,bl,cl=parse_lc(LC),parse_lc(BASE/'fit.LCPLOT.TEXT'),parse_lc(COPY/'fit.LCPLOT.TEXT')
    failures=[]; ledger=[]
    copy_raw=(BASE/'fit.LCPLOT.TEXT').read_bytes()==(COPY/'fit.LCPLOT.TEXT').read_bytes()
    if not copy_raw:failures.append('baseline-copy LCPLOT bytes differ')
    if bstring!=cstring:failures.append('baseline-copy FITRES SN rows differ')
    for cid in CIDS:
        R={'CID':cid,'baseline_present':cid in b,'copy_present':cid in c}
        if cid not in b or cid not in c:
            failures.append(f'{cid}: FITRES missing');ledger.append(R);continue
        x,y=a[cid],b[cid]
        R['ERRFLAG_FIT']=y['ERRFLAG_FIT'];R['NDOF_archive']=x['NDOF'];R['NDOF_baseline']=y['NDOF']
        R['DLMAG_archive']=x['DLMAG'];R['DLMAG_baseline']=y['DLMAG']
        R['abs_DLMAG_diff']=abs(float(y['DLMAG'])-float(x['DLMAG']))
        R['FITCHI2_archive']=x['FITCHI2'];R['FITCHI2_baseline']=y['FITCHI2']
        R['abs_dataQ_diff']=abs(float(y['FITCHI2'])-float(x['FITCHI2']))
        R['dataQ_cap']=max(.01,.001*abs(float(x['FITCHI2'])))
        R['peak_f32_input']=f32(x['PKMJDINI']);R['peak_baseline']=y['PKMJD']
        R['peak_abs_f32_diff']=abs(float(y['PKMJD'])-f32(x['PKMJDINI']))
        R['archived_data_rows']=len(al[cid]);R['baseline_data_rows']=len(bl[cid])
        R['accepted_mask_archive']=sorted((q['BAND'],q['MJD']) for q in al[cid])
        R['accepted_mask_baseline']=sorted((q['BAND'],q['MJD']) for q in bl[cid])
        if y['ERRFLAG_FIT']!='0':failures.append(f'{cid}: ERRFLAG {y["ERRFLAG_FIT"]}')
        if y['NDOF']!=x['NDOF']:failures.append(f'{cid}: NDOF mismatch')
        if R['accepted_mask_archive']!=R['accepted_mask_baseline']:failures.append(f'{cid}: accepted mask mismatch')
        if R['abs_DLMAG_diff']>.001:failures.append(f'{cid}: DLMAG cap')
        if R['abs_dataQ_diff']>R['dataQ_cap']:failures.append(f'{cid}: dataQ cap')
        if not math.isfinite(float(y['DLMAG'])) or not math.isfinite(float(y['FITCHI2'])):failures.append(f'{cid}: nonfinite')
        for name,target in [('STRETCH',1.),('AV',0.),('RV',1.518),('PKMJDERR',0.),('STRETCHERR',0.),('AVERR',0.)]:
            if float(y[name])!=target:failures.append(f'{cid}: fixed {name}')
        if R['peak_abs_f32_diff']>half_unit(y['PKMJD'])+f32_ulp(x['PKMJDINI']):failures.append(f'{cid}: fixed peak f32')
        R['data_token_checks']=[]
        orig_by_key=collections.defaultdict(list)
        for q in al[cid]:orig_by_key[(q['BAND'],q['MJD'])].append(q)
        for q in bl[cid]:
            key=(q['BAND'],q['MJD'])
            if not orig_by_key[key]:continue
            old=orig_by_key[key].pop(0)
            for col in ['FLUXCAL','FLUXCAL_ERR']:
                diff=abs(float(q[col])-float(old[col]))
                cap=half_unit(q[col])+half_unit(old[col])+f32_ulp(old[col])
                okay=diff<=cap
                R['data_token_checks'].append({'band':key[0],'mjd':key[1],'column':col,'old':old[col],'new':q[col],
                                               'diff':diff,'source_print_and_f32_cap':cap,'pass':okay})
                if not okay:failures.append(f'{cid}: outgoing {col} outside print/f32 bound')
        if not all(y['zHEL']==x['zHEL'] for _ in [0]):failures.append(f'{cid}: zHEL changed')
        ledger.append(R)
    result={'pass':not failures,'failures':failures,'rows':ledger,'Rcopy_fitres_SN_identical':bstring==cstring,
            'Rcopy_LCPLOT_bytes_identical':copy_raw,'source_baseline_runner_sha256':sha(OUT/'baseline_runner_authorgrid.py'),
            'gate_source_sha256':sha(__file__),'baseline_freeze_sha256':sha(OUT/'baseline-freeze-authorgrid.json'),
            'native_baseline':json.loads((BASE/'execution.json').read_text()),
            'native_copy':json.loads((COPY/'execution.json').read_text())}
    (OUT/'baseline-gate-authorgrid.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'pass':result['pass'],'failures':failures,'copy_SN':result['Rcopy_fitres_SN_identical'],
                      'copy_LCPLOT':copy_raw,'max_abs_DLMAG':max(x.get('abs_DLMAG_diff',0) for x in ledger),
                      'max_abs_dataQ':max(x.get('abs_dataQ_diff',0) for x in ledger)},indent=2))


if __name__=='__main__': gate()
