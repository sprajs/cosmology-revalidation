#!/usr/bin/env python3
"""Independent raw-log/DAT verification; never imports Astra runner or parser."""
import csv,ctypes,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
AS=ROOT/'runs/research_2026_09_26/astra_design/raisin_timing_assets'
PREC=AS/'rounding_2021'
ABS=AS/'instrumentation_2021/fits/absent'
PROTO=json.loads((HERE/'protocol.json').read_text())


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def parse_dat(path):
    hdr={};obs=[];var=None
    for line in path.read_text().splitlines():
        if line.startswith('VARLIST:'):var=line.split()[1:]
        elif line.startswith('OBS:'):
            values=line.split()[1:];assert var and len(values)==len(var),(path,line)
            obs.append(dict(zip(var,values)))
        elif ':' in line and not line.startswith('#'):
            k,v=line.split(':',1);hdr[k]=v.strip()
    assert var==['MJD','FLT','FIELD','FLUXCAL','FLUXCALERR']
    assert int(hdr['NOBS'])==len(obs)
    return hdr,obs

def parse_fitres(path):
    keys=None;out={}
    for line in path.read_text().splitlines():
        if line.startswith('VARNAMES:'):keys=line.split()[1:]
        elif line.startswith('SN:'):
            v=line.split()[1:];assert keys and len(v)==len(keys)
            r=dict(zip(keys,v));assert r['CID'] not in out;out[r['CID']]=r
    return out

def parse_log(path):
    out=defaultdict(lambda:defaultdict(lambda:{'rows':[],'wrows':{},'wdiag':{},'objective':None}))
    for line in path.read_text().splitlines():
        if not line.startswith('CSP_'):continue
        x=line.split();tag=x[0]
        if tag=='CSP_ROW:':
            assert len(x)==19,(path,line)
            cid,itr=x[1],int(x[2]);out[cid][itr]['rows'].append({'ifit':int(x[3]),'epoch':int(x[4]),'band':x[5],
                'mjd':float(x[6]),'trest':float(x[7]),'model':float(x[8]),'magerr':float(x[9]),
                'data':float(x[10]),'error':float(x[11]),'z':float(x[12]),'mwebv':float(x[13]),
                'lam':float(x[14]),'fudge':float(x[15]),'modelerr':float(x[16]),'rest1':float(x[17]),'rest2':float(x[18])})
        elif tag=='CSP_OBJECTIVE:':
            assert len(x)==15,(path,line)
            cid,itr=x[1],int(x[2]);v=list(map(float,x[5:]))
            assert out[cid][itr]['objective'] is None
            out[cid][itr]['objective']={'n':int(x[3]),'usefitcov':x[4], 'chi2':v[0], 'chi2ini':v[1],
                'sigma':v[2], 'D':v[3],'shape':v[4],'AV':v[5],'peak':v[6],
                'peak_ini':v[7],'search_peak':v[8],'mjd_off':v[9]}
        elif tag=='CSP_WROW:':
            cid,itr,i=x[1],int(x[2]),int(x[3]);assert i not in out[cid][itr]['wrows']
            out[cid][itr]['wrows'][i]=list(map(float,x[4:]))
        elif tag=='CSP_WDIAG:':
            assert len(x)==5
            cid,itr,i=x[1],int(x[2]),int(x[3]);assert i not in out[cid][itr]['wdiag']
            out[cid][itr]['wdiag'][i]=float(x[4])
    return out

def final_state(states,cid):
    assert cid in states and sorted(states[cid])==[1,2,3],(cid,sorted(states.get(cid,{})))
    return states[cid][3]

def matrix(st):
    o=st['objective'];n=o['n']
    if st['wrows']:
        assert sorted(st['wrows'])==list(range(1,n+1))
        W=np.array([st['wrows'][i] for i in range(1,n+1)],float)
    else:
        assert sorted(st['wdiag'])==list(range(1,n+1))
        W=np.diag([st['wdiag'][i] for i in range(1,n+1)])
    assert W.shape==(n,n)
    return W

def rowmultiset(rows):
    return Counter((r['mjd'],r['band'],r['data'],r['error']) for r in rows)

def source_multiset(obs):
    return Counter((float(r['MJD']),r['FLT'],float(np.float32(float(r['FLUXCAL']))),
                    float(np.float32(float(r['FLUXCALERR'])))) for r in obs)

def all_state_signature(st):
    # Exact Python float/text representation after normalizing the CID at parse time.
    return st

def compare_source(original,case,record):
    oldh,oldo=original;h,o=case
    issues=[]
    if h.get('SNID')!=record['new_CID']:issues.append('SNID')
    for k in set(oldh)|set(h):
        if k=='SNID':continue
        if k not in oldh or k not in h:issues.append('header_key_'+k);continue
        expected=oldh[k]
        changes=[c for c in record['changes'] if c['field']==k]
        if changes:
            assert len(changes)==1 and changes[0]['epoch'] is None
            oldparts=oldh[k].split();parts=h[k].split()
            if oldparts[1:]!=parts[1:]:issues.append('header_suffix_'+k)
            if float(parts[0])!=float(changes[0]['value']):issues.append('header_changed_'+k)
        elif h[k]!=expected:issues.append('header_undeclared_'+k)
    if len(o)!=len(oldo):issues.append('obs_count')
    else:
        for idx,(before,after) in enumerate(zip(oldo,o),1):
            for field in before:
                changes=[c for c in record['changes'] if c['field']==field and c['epoch']==idx]
                if changes:
                    if len(changes)!=1 or float(after[field])!=float(changes[0]['value']):issues.append(f'changed_{field}_{idx}')
                elif after[field]!=before[field]:issues.append(f'undeclared_{field}_{idx}')
        for c in record['changes']:
            if c['epoch'] is not None and (c['epoch']<1 or c['epoch']>len(o) or c['field'] not in oldo[c['epoch']-1]):issues.append('unknown_change')
    return issues

def check_case(stage,record,states,fit,baseline,basefit,original_data,archiveQ):
    cid=record['new_CID'];orig=record['original_CID'];issues=[]
    path=ROOT/record['file']
    if sha(path)!=record['sha256']:issues.append('input_hash')
    h,obs=parse_dat(path)
    issues+=compare_source(original_data[orig],(h,obs),record)
    if cid not in fit:issues.append('missing_FITRES')
    if cid not in states:issues.append('missing_native_log');return {'CID':cid,'original_CID':orig,'kind':record['kind'],'issues':issues}
    s=final_state(states,cid);obj=s['objective'];n=obj['n'];rows=sorted(s['rows'],key=lambda r:r['ifit'])
    if len(rows)!=n or [r['ifit'] for r in rows]!=list(range(1,n+1)):issues.append('native_row_count_index')
    if n!=len(obs) or rowmultiset(rows)!=source_multiset(obs):issues.append('source_physical_multiset')
    if cid in fit:
        fr=fit[cid]
        if fr['ERRFLAG_FIT']!='0':issues.append('ERRFLAG_FIT')
        if int(fr['NDOF'])!=n-1:issues.append('NDOF')
        for key,val in [('STRETCH',1.),('AV',0.),('RV',1.518)]:
            if float(fr[key])!=val:issues.append('fixed_'+key)
        for key in ['STRETCHERR','AVERR','RVERR','PKMJDERR']:
            if float(fr[key])!=0:issues.append('fixed_error_'+key)
        if float(fr['PKMJD'])!=float(fr['PKMJDINI']) or fr['PKMJD']!=basefit[orig]['PKMJD']:
            issues.append('fixed_peak')
    if obj is None:issues.append('missing_objective')
    if obj['n']!=n or obj['usefitcov']!='T':issues.append('final_cov_flag')
    if obj['shape']!=1 or obj['AV']!=0 or obj['peak']!=float(np.float32(float(h['PEAKMJD']))):issues.append('objective_fixed')
    W=matrix(s)
    if not np.all(np.isfinite(W)) or np.max(np.abs(W-W.T))>1e-10:issues.append('W_finite_symmetric')
    wmin=float(np.linalg.eigvalsh(W).min())
    if wmin<=0:issues.append('W_nonpositive')
    try:C=np.linalg.inv(W);cmin=float(np.linalg.eigvalsh(C).min())
    except np.linalg.LinAlgError:C=None;cmin=float('nan');issues.append('C_inverse')
    if not np.isfinite(cmin) or cmin<=0:issues.append('C_nonpositive')
    d=np.array([r['data'] for r in rows]);m=np.array([r['model'] for r in rows]);delta=d-m
    q=float(delta@W@delta)
    qgap=abs(q-obj['chi2'])
    if qgap>1e-8:issues.append('objective_closure')
    dstep=abs(obj['D']-states[cid][2]['objective']['D'])
    if dstep>0.001:issues.append('iteration_D')
    denom=float(m@W@m);numer=float(m@W@d)
    if denom<=0 or numer<=0:issues.append('nonpositive_amplitude');ampgap=float('nan')
    else:ampgap=abs(-2.5*math.log10(numer/denom))
    if not np.isfinite(ampgap) or ampgap>0.001:issues.append('fixedC_amplitude_gap')
    b=final_state(baseline,orig);bobj=b['objective']
    deltaD=obj['D']-bobj['D'];deltaQ=obj['chi2']-bobj['chi2']
    capQ=max(0.01,0.001*archiveQ[orig])
    if abs(deltaD)>0.001:issues.append('distance_cap')
    if abs(deltaQ)>capQ:issues.append('dataQ_cap')
    if record['kind']=='clone_control':
        if record['changes']:issues.append('clone_changes')
        if cid not in fit or {k:v for k,v in fit[cid].items() if k!='CID'}!={k:v for k,v in basefit[orig].items() if k!='CID'}:
            issues.append('clone_FITRES')
        if states[cid]!=baseline[orig]:issues.append('clone_fullstates')
    return {'CID':cid,'original_CID':orig,'kind':record['kind'],'issues':issues,
      'n':n,'q_from_rows':q,'q_native':obj['chi2'],'objective_gap':qgap,'W_min_eig':wmin,'C_min_eig':cmin,
      'finalD_step':dstep,'fixedC_amp_D_gap':ampgap,'D':obj['D'],'dataQ':obj['chi2'],
      'deltaD':deltaD,'deltaQ':deltaQ,'dataQ_cap':capQ}

def check_times():
    ledger=json.loads((PREC/'precision-ledger.json').read_text())
    coords=[c for c in ledger['coordinates'] if c['field']=='MJD']
    assert len(coords)==32
    libc=ctypes.CDLL(None);libc.snprintf.restype=ctypes.c_int
    result=[]
    for c in coords:
        token=c['token'];center=np.float32(float(token))
        scan=[center]
        x=center
        for _ in range(8):x=np.nextafter(x,np.float32(-np.inf));scan.append(x)
        x=center
        for _ in range(8):x=np.nextafter(x,np.float32(np.inf));scan.append(x)
        found=[]
        for x in scan:
            buf=ctypes.create_string_buffer(128)
            libc.snprintf(buf,128,b'%.3f',ctypes.c_double(float(x)))
            if buf.value.decode()==token:found.append(float(x))
        found=sorted(set(found));issues=[]
        if len(found)!=1 or c['float32_state_count']!=1:issues.append('state_count')
        if len(found)==1:
            s=np.float32(found[0]);prev=np.nextafter(s,np.float32(-np.inf));nxt=np.nextafter(s,np.float32(np.inf))
            low=math.nextafter((float(prev)+float(s))/2,math.inf)
            high=math.nextafter((float(s)+float(nxt))/2,-math.inf)
            if low!=c['input_low'] or high!=c['input_high']:issues.append('inverse_cell')
            if found[0]!=c['float32_low'] or found[0]!=c['float32_high']:issues.append('float32_state')
            width=(float(s)+float(nxt))/2-(float(prev)+float(s))/2
        else:low=high=width=None
        result.append({'CID':c['CID'],'epoch':c['epoch'],'token':token,'states':found,'inverse_low':low,'inverse_high':high,'cell_width_day':width,'issues':issues})
    return result

def main():
    base_states=parse_log(ABS/'fit.log');base_fit=parse_fitres(ABS/'fit.FITRES.TEXT')
    assert set(base_states)==set(base_fit)=={str(i) for i in range(1,9)}
    orig={cid:parse_dat(ABS/'data/DES_RAISIN_SIM'/f'{cid}.DAT') for cid in base_fit}
    archived=AS/'author_20211111/nir.FITRES.gz'
    assert sha(archived)=='c0e6de446c2d4b2179a6570c0a1f166b5384612101fbe5002938a73b8daa3b44'
    archiveQ={};keys=None
    with gzip.open(archived,'rt') as f:
        for line in f:
            if line.startswith('VARNAMES:'):keys=line.split()[1:]
            elif line.startswith('SN:'):
                values=line.split()[1:];assert len(values)==len(keys)
                r=dict(zip(keys,values));archiveQ[r['CID']]=float(r['FITCHI2'])
    assert set(base_fit).issubset(archiveQ)
    precision=json.loads((PREC/'precision-ledger.json').read_text())
    for cid in orig:assert sha(ABS/'data/DES_RAISIN_SIM'/f'{cid}.DAT')==precision['objects'][cid]['input_sha256']
    stages={}
    for stage,expected in [('endpoints',248),('corners',40)]:
        folder=PREC/'fits'/stage
        ledger=json.loads((folder/'case-ledger.json').read_text())
        assert len(ledger)==expected and len({r['new_CID'] for r in ledger})==expected
        states=parse_log(folder/'fit.log');fit=parse_fitres(folder/'fit.FITRES.TEXT')
        if set(states)!={r['new_CID'] for r in ledger} or set(fit)!={r['new_CID'] for r in ledger}:
            raise AssertionError((stage,'case-ID closure',len(states),len(fit)))
        cases=[]
        for rec in ledger:
            try:cases.append(check_case(stage,rec,states,fit,base_states,base_fit,orig,archiveQ))
            except Exception as ex:cases.append({'CID':rec['new_CID'],'original_CID':rec['original_CID'],'kind':rec['kind'],'issues':['exception:'+repr(ex)]})
        stages[stage]=cases
    times=check_times()
    (HERE/'case-results.json').write_text(json.dumps(stages,indent=2)+'\n')
    (HERE/'timestamp-cells.json').write_text(json.dumps(times,indent=2)+'\n')
    summary={'protocol_sha256':sha(HERE/'protocol.json'),'counts':{},'timestamp':{'count':len(times),'failures':[x for x in times if x['issues']]}}
    for stage,cases in stages.items():
        fails=[r for r in cases if r['issues']]
        summary['counts'][stage]={'cases':len(cases),'failures':len(fails),'failure_ids':[r['CID'] for r in fails],
          'issue_counts':dict(Counter(i for r in cases for i in r['issues'])),
          'max_objective_gap':max((r.get('objective_gap',float('nan')) for r in cases),default=None),
          'max_abs_deltaD':max((abs(r.get('deltaD',float('nan'))) for r in cases),default=None),
          'max_abs_deltaQ':max((abs(r.get('deltaQ',float('nan'))) for r in cases),default=None),
          'max_finalD_step':max((r.get('finalD_step',float('nan')) for r in cases),default=None),
          'max_amp_D_gap':max((r.get('fixedC_amp_D_gap',float('nan')) for r in cases),default=None)}
    (HERE/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
