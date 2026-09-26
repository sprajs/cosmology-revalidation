"""Independent source-line verification of the historical timing state controls."""
from pathlib import Path
import collections, hashlib, json
import numpy as np
from review_csp_native_states import parse

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/astra_design/raisin_timing_assets'
OUT=ROOT/'runs/research_2026_09_26/raisin_historical_pair_root_review/native-states'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def table(path):
    names=None; rows={}
    for line in path.read_text().splitlines():
        if line.startswith('VARNAMES:'):names=line.split()[1:]
        elif line.startswith('SN:'):
            x=line.split()[1:]; assert len(x)==len(names)
            row=dict(zip(names,x)); assert row['CID'] not in rows
            rows[row['CID']]=row
    return rows

def main():
    OUT.mkdir(exist_ok=False)
    inputs=[]; groups={}; entries={}; checks=[]; tables={}
    reference=table(BASE/'baseline_2021/fits/baseline/fit.FITRES.TEXT')
    for name in ['absent','zero','minus','plus']:
        work=BASE/'instrumentation_2021/fits'/name
        path=work/'fit.log';inputs.append(path)
        records,ent=parse(path);assert len(records)==24
        groups[name]={};entries[name]=ent; tables[name]=table(work/'fit.FITRES.TEXT')
        for r in records:
            cid=r['CID'];groups[name].setdefault(cid,[]).append(r)
            source=work/'data/DES_RAISIN_SIM'/f'{cid}.DAT';inputs.append(source)
            lines=source.read_text().splitlines(); obs=[x.split() for x in lines if x.startswith('OBS:')]
            peak=float(np.float32(float(next(x.split()[1] for x in lines if x.startswith('PEAKMJD:')))))
            expected=collections.Counter((x[2],float(x[1]),float(np.float32(x[4])),float(np.float32(x[5]))) for x in obs)
            actual=collections.Counter((x['band'],x['values'][0],x['values'][4],x['values'][5]) for x in r['rows'])
            assert actual==expected
            arr=np.array([x['values'] for x in r['rows']]); W=np.array(r['weights'])
            if not r['cov']: W=np.diag(W[:,0])
            L=np.linalg.cholesky(W); residual=arr[:,4]-arr[:,2]
            q=float(np.sum((L.T@residual)**2));o=r['objective'];gap=abs(q-o[0]+o[1]+o[2])
            assert gap<=1e-8 and o[4]==1 and o[5]==0 and o[6]==peak and o[9]==0
            amplitude=float(arr[:,2]@W@arr[:,4]/(arr[:,2]@W@arr[:,2]));assert amplitude>0
            optimum=float(-2.5*np.log10(amplitude))
            checks.append(dict(condition=name,CID=cid,iteration=r['iteration'],Q_gap=gap,fixed_C_D_gap=optimum))
        assert set(groups[name])==set(map(str,range(1,9)))
        for cid,rs in groups[name].items():
            assert [x['iteration'] for x in rs]==[1,2,3]
            assert abs(rs[-1]['objective'][3]-rs[-2]['objective'][3])<=.001
            assert abs(next(x['fixed_C_D_gap'] for x in checks if x['condition']==name and x['CID']==cid and x['iteration']==3))<=.001
            row=tables[name][cid]
            for k in ['STRETCH','STRETCHERR','AV','AVERR','RV','RVERR','PKMJD','PKMJDERR','PKMJDINI','NDOF','ERRFLAG_FIT']:
                assert row[k]==reference[cid][k],(name,cid,k)
    assert groups['absent']==groups['zero']
    for name in ['absent','zero']:
        assert tables[name]==reference
        actual=BASE/'instrumentation_2021/fits'/name/'fit.LCPLOT.TEXT'
        old=BASE/'baseline_2021/fits/baseline/fit.LCPLOT.TEXT'
        inputs += [actual,old];assert actual.read_bytes()==old.read_bytes()
    objects=[]
    for cid in reference:
        shifts=[]
        start=float(next(e for e in entries['absent'] if e[0]==cid and e[1]=='1')[5])
        for name,target in [('minus',-.2),('plus',.2)]:
            shift=float(next(e for e in entries[name] if e[0]==cid and e[1]=='1')[5])-start
            assert abs(shift-target)<1e-10;shifts.append(shift)
        ds=[groups[n][cid][-1]['objective'][3] for n in groups]
        spread=max(ds)-min(ds);assert spread<=.001
        objects.append(dict(CID=cid,actual_start_shifts=shifts,final_D_spread=spread))
    result=dict(pass_declared_state_gates=True,callback_count=len(checks),max_Q_closure=max(x['Q_gap'] for x in checks),
                max_final_D_spread=max(x['final_D_spread'] for x in objects),objects=objects,checks=checks,
                scope='Reproducibility of the fixed historical three-iteration estimator. Does not establish convergence of a parameter-dependent normalized likelihood or population calibration.')
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    inputs += [Path(__file__).with_name('review_csp_native_states.py'),BASE/'baseline_2021/fits/baseline/fit.FITRES.TEXT']
    (OUT/'manifest.json').write_text(json.dumps(dict(inputs={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(inputs))},outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks','objects']},indent=2))

if __name__=='__main__':main()
