"""Root independent verification of frozen all42 native continuation gates."""
from pathlib import Path
import hashlib, json
import numpy as np
from scipy.linalg import solve_triangular
from review_csp_native_states import parse

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26'
SRC=BASE/'raisin_profile_solver_review/csp-filter-interpretation/native-design/convergence-diagnostic'
OUT=BASE/'csp_native_filter_response/root-state-review/continuation'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def groups(log):
    rs,es=parse(log);out={}
    for r in rs:out.setdefault(r['CID'],[]).append(r)
    return out,es

def physical(r):
    return [(x['band'],x['values'][0],x['values'][4],x['values'][5]) for x in r['rows']]

def precision(r):
    w=np.array(r['weights']);return w if r['cov'] else np.diag(w[:,0])

def cchange(a,b):
    W0=precision(a);W1=precision(b);C0=np.linalg.inv(W0);C1=np.linalg.inv(W1)
    L=np.linalg.cholesky(C0);left=solve_triangular(L,C1-C0,lower=True)
    whiten=solve_triangular(L,left.T,lower=True).T
    return dict(C_relative=float(np.linalg.norm(C1-C0)/np.linalg.norm(C0)),
                W_relative=float(np.linalg.norm(W1-W0)/np.linalg.norm(W0)),
                C_whitened_operator_norm=float(np.max(abs(np.linalg.eigvalsh(whiten)))))

def main():
    OUT.mkdir(exist_ok=False)
    logs={n:SRC/'fits'/n/'fit.log' for n in ['iter09_default','iter12_default','iter12_minus','iter12_plus']}
    logs['original3']=BASE/'csp_native_filter_response/cohort-execution/full/fits/nominal/fit.log'
    fits={};entries={};gaps=[]
    for name,path in logs.items():
        fits[name],entries[name]=groups(path)
        assert len(fits[name])==42
        for cid,records in fits[name].items():
            for r in records:
                W=precision(r);L=np.linalg.cholesky(W)
                x=np.array([v['values'] for v in r['rows']]);res=x[:,4]-x[:,2]
                q=float(np.sum((L.T@res)**2));gap=abs(q-sum([r['objective'][0],-r['objective'][1],-r['objective'][2]]))
                gaps.append(gap);assert gap<=1e-8
    ids=list(fits['original3']);assert all(set(fits[n])==set(ids) for n in fits)
    prefix={}
    for old,new,n in [('original3','iter09_default',3),('iter09_default','iter12_default',9)]:
        same=all(fits[old][cid][:n]==fits[new][cid][:n] for cid in ids)
        prefix[f'{old}_to_{new}']=same;assert same
    objects=[]
    for cid in ids:
        orig=fits['original3'][cid][-1];nine=fits['iter09_default'][cid][-1];twelve=fits['iter12_default'][cid]
        assert len(twelve)==12 and len(fits['iter09_default'][cid])==9
        last=twelve[-1];x=np.array([r['values'] for r in last['rows']]);W=precision(last);m=x[:,2];y=x[:,4]
        amplitude=float(m@W@y/(m@W@m));assert amplitude>0
        optimum=-2.5*np.log10(amplitude)
        ds=[rs[cid][-1]['objective'][3] for name,rs in fits.items() if name.startswith('iter12')]
        start_shifts=[]
        baseline_entry=float(next(e for e in entries['iter12_default'] if e[0]==cid and e[1]=='1')[5])
        for name,offset in [('iter12_minus',-.2),('iter12_plus',.2)]:
            actual=float(next(e for e in entries[name] if e[0]==cid and e[1]=='1')[5])-baseline_entry
            assert abs(actual-offset)<1e-10
            start_shifts.append(actual)
        mask=all(physical(r)==physical(last) for name,rs in fits.items() for r in rs[cid] if r['iteration']>=2)
        updates=[cchange(twelve[i-1],twelve[i]) for i in [10,11]]
        dlast=[twelve[i]['objective'][3]-twelve[i-1]['objective'][3] for i in [10,11]]
        row=dict(CID=cid,D3=orig['objective'][3],D9=nine['objective'][3],D12=last['objective'][3],
                 D12_minus_D3=last['objective'][3]-orig['objective'][3],D12_minus_D9=last['objective'][3]-nine['objective'][3],
                 last_D_updates=dlast,fixed_C_optimum_D_gap=float(optimum),start_final_D_spread=max(ds)-min(ds),
                 actual_postinit_start_shifts=start_shifts,same_physical_rows_from_iter2=mask,last_C_updates=updates,
                 final_phase_supported=bool(np.all((x[:,1]>=-20)&(x[:,1]<=70))),
                 final_positive_mean=bool(np.all(m>0)))
        row['pass_numerical_convergence']=bool(mask and abs(row['D12_minus_D9'])<=.001 and
            max(abs(v) for v in dlast)<=.001 and abs(optimum)<=.001 and row['start_final_D_spread']<=.001 and
            max(c['C_whitened_operator_norm'] for c in updates)<=.001)
        objects.append(row)
    shifts=[r['D12_minus_D3'] for r in objects]
    result=dict(pass_numerical_convergence=all(r['pass_numerical_convergence'] for r in objects),
                scope='Nominal current-input numerical recipe only. Original archived-Q and initial empirical-support failures remain. No filter response or physical likelihood validation.',
                prefixes_exact=prefix,max_callback_Q_gap=max(gaps),callback_count=len(gaps),objects=objects,
                D12_minus_D3_mean=float(np.mean(shifts)),D12_minus_D3_range=[min(shifts),max(shifts)])
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'manifest.json').write_text(json.dumps(dict(inputs={str(p.relative_to(ROOT)):sha(p) for p in list(logs.values())+[Path(__file__).with_name('review_csp_native_states.py')]},outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='objects'}|{'failed':[r for r in objects if not r['pass_numerical_convergence']]},indent=2))

if __name__=='__main__':main()
