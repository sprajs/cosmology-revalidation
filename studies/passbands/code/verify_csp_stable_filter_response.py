"""Root independent stable filter-response audit; no collaborator parser imports."""
from pathlib import Path
from collections import Counter, defaultdict
import json, hashlib, csv
import numpy as np
from scipy.linalg import solve_triangular
from review_csp_native_states import parse

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26'
NATIVE=BASE/'raisin_profile_solver_review/csp-filter-interpretation/native-design'
SRC=NATIVE/'stable-filter-response'
OUT=BASE/'csp_native_filter_response/root-state-review/stable-response'
PREP=BASE/'csp_native_filter_response/cohort-preparation/full/data'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def groups(path):
    records,entries=parse(path);ans={}
    for r in records:ans.setdefault(r['CID'],[]).append(r)
    return ans,entries

def W(r):
    x=np.array(r['weights']);return x if r['cov'] else np.diag(x[:,0])

def physical(r):
    return Counter((x['values'][0],x['values'][4],x['values'][5],x['band']) for x in r['rows'])

def c_norm(a,b):
    Ca=np.linalg.inv(W(a));Cb=np.linalg.inv(W(b));L=np.linalg.cholesky(Ca)
    z=solve_triangular(L,Cb-Ca,lower=True)
    z=solve_triangular(L,z.T,lower=True).T
    return float(np.max(abs(np.linalg.eigvalsh(z))))

def text_rows(path,tag,cid=None):
    return [s for s in path.read_text().splitlines() if s.startswith('VARNAMES:') or
            s.startswith(tag) and (cid is None or s.split()[1]==cid)]

def main():
    OUT.mkdir(exist_ok=False)
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    proto=json.loads((SRC/'protocol.json').read_text());ids=proto['membership']
    # Derive the physical-row mapping independently from paired source lines.
    # Native data flux/error are R4, while CSP_ROW MJD remains R8.
    mappings={};input_hashes={};changed_by_cid=Counter()
    for cid in ids:
        paths=[PREP/arm/'CSPDR3_RAISIN'/f'CSPDR3_{cid}.PKMJD.DAT' for arm in ['nominal','known_Jdw_to_j']]
        source=[p.read_text().splitlines() for p in paths];assert len(source[0])==len(source[1])
        mapping=defaultdict(list)
        for x,y in zip(*source,strict=True):
            if x.startswith('OBS:'):
                a=x.split();b=y.split();assert a[:2]+a[3:]==b[:2]+b[3:]
                if a[2]!=b[2]:assert (a[2],b[2])==('J','j');changed_by_cid[cid]+=1
                key=(float(a[1]),float(np.float32(a[4])),float(np.float32(a[5])),a[2]);mapping[key].append(b[2])
            else:assert x==y
        mappings[cid]=mapping
        input_hashes.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
    assert sum(changed_by_cid.values())==188 and set(changed_by_cid)==set(proto['affected_32'])
    arms={};entries={};paths={}
    for n in ['changed09_default','changed12_default','changed12_minus','changed12_plus','nominal12_copy']:
        paths[n]=SRC/'fits'/n/'fit.log';arms[n],entries[n]=groups(paths[n]);assert set(arms[n])==set(ids)
    for n in ['iter12_default','iter12_minus','iter12_plus']:
        paths[n]=NATIVE/'convergence-diagnostic/fits'/n/'fit.log';arms[n],entries[n]=groups(paths[n]);assert set(arms[n])==set(ids)
    qgaps=[];mask_count=0;control_count=0
    for name,bycid in arms.items():
        for cid,states in bycid.items():
            assert len(states)==(9 if name=='changed09_default' else 12)
            for i,r in enumerate(states):
                x=np.array([v['values'] for v in r['rows']]);w=W(r);L=np.linalg.cholesky(w);res=x[:,4]-x[:,2]
                gap=abs(float(np.sum((L.T@res)**2))-(r['objective'][0]-r['objective'][1]-r['objective'][2]));qgaps.append(gap);assert gap<=1e-8
                assert np.all(np.isfinite(x)) and np.all(x[:,2]>0)
                if i>=1:assert np.all((x[:,1]>=-20)&(x[:,1]<=70))
                nominal=arms['iter12_default'][cid][i]
                if name.startswith('changed'):
                    expected=Counter()
                    for key,n in physical(nominal).items():
                        targets=mappings[cid][key];assert len(targets)>=n
                        if len(set(targets))==1:targets=[targets[0]]*n
                        else:assert len(targets)==n
                        expected.update((*key[:3],b) for b in targets)
                    assert expected==physical(r);mask_count+=1
                    if cid in proto['controls_10']:
                        ref='iter12_minus' if name.endswith('_minus') else 'iter12_plus' if name.endswith('_plus') else 'iter12_default'
                        assert r==arms[ref][cid][i];control_count+=1
                elif name=='nominal12_copy':assert r==nominal
    assert all(arms['changed09_default'][c]==arms['changed12_default'][c][:9] for c in ids)
    original_dir=NATIVE/'convergence-diagnostic/fits/iter12_default';copy_dir=SRC/'fits/nominal12_copy'
    assert text_rows(original_dir/'fit.FITRES.TEXT','SN:')==text_rows(copy_dir/'fit.FITRES.TEXT','SN:')
    assert (original_dir/'fit.LCPLOT.TEXT').read_bytes()==(copy_dir/'fit.LCPLOT.TEXT').read_bytes()
    results=[]
    for cid in ids:
        n=arms['iter12_default'][cid];c=arms['changed12_default'][cid];r=c[-1]
        x=np.array([s['values'] for s in r['rows']]);w=W(r);m=x[:,2];y=x[:,4]
        amplitude=float(m@w@y/(m@w@m));assert amplitude>0
        delta_opt=float(-2.5*np.log10(amplitude));assert abs(delta_opt)<=.001
        last_D=[c[i]['objective'][3]-c[i-1]['objective'][3] for i in [10,11]]
        last_C=[c_norm(c[i-1],c[i]) for i in [10,11]]
        ds=[arms[a][cid][-1]['objective'][3] for a in ['changed12_default','changed12_minus','changed12_plus']]
        difference_9=r['objective'][3]-arms['changed09_default'][cid][-1]['objective'][3]
        assert max(abs(v) for v in last_D)<=.001 and max(last_C)<=.001 and max(ds)-min(ds)<=.001 and abs(difference_9)<=.001
        e0=float(next(e for e in entries['changed12_default'] if e[0]==cid and e[1]=='1')[5])
        for arm,offset in [('changed12_minus',-.2),('changed12_plus',.2)]:
            e=float(next(e for e in entries[arm] if e[0]==cid and e[1]=='1')[5]);assert abs(e-e0-offset)<=1e-12
        if cid in proto['controls_10']:
            for filename,tag in [('fit.FITRES.TEXT','SN:'),('fit.LCPLOT.TEXT','OBS:')]:
                assert text_rows(original_dir/filename,tag,cid)==text_rows(SRC/'fits/changed12_default'/filename,tag,cid)
        results.append(dict(CID=cid,affected=cid in changed_by_cid,changed_input_rows=changed_by_cid[cid],
                            old_D12=n[-1]['objective'][3],new_D12=r['objective'][3],
                            delta_D12=r['objective'][3]-n[-1]['objective'][3],
                            delta_D3=c[2]['objective'][3]-n[2]['objective'][3],
                            last_D_changes=last_D,last_C_whitened_norms=last_C,
                            D12_minus_D9=difference_9,start_final_D_spread=max(ds)-min(ds),
                            final_fixed_C_D_gap=delta_opt))
    a32=[r['delta_D12'] for r in results if r['affected']];a42=[r['delta_D12'] for r in results]
    answer=dict(pass_root_numerical_and_physical_row_checks=True,scope='Current released-input stable native processing response; not corrected cosmology, exact WIRC calibration, or empirical validation of initial extrapolation.',
                objects=results,changed_input_rows=188,physical_callback_checks=mask_count,exact_control_callback_checks=control_count,
                max_callback_Q_gap=max(qgaps),mean_32=float(np.mean(a32)),mean_42=float(np.mean(a42)),
                holding_highz_fixed_high_minus_low=-float(np.mean(a42)),range_affected=[min(a32),max(a32)],
                secondary_mean_D3_42=float(np.mean([r['delta_D3'] for r in results])))
    (OUT/'result.json').write_text(json.dumps(answer,indent=2)+'\n')
    with (OUT/'responses.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['CID','affected','changed_input_rows','old_D12','new_D12','delta_D12','delta_D3']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in results)
    input_hashes.update({str(p.relative_to(ROOT)):sha(p) for p in list(paths.values())+[SRC/'protocol.json',Path(__file__).with_name('review_csp_native_states.py')]})
    (OUT/'manifest.json').write_text(json.dumps(dict(inputs=input_hashes,outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(json.dumps({k:v for k,v in answer.items() if k!='objects'},indent=2))

if __name__=='__main__':main()
