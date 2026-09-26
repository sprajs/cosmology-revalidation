"""Prespecified mean-versus-weight bookkeeping for the stable all42 response."""
from pathlib import Path
from collections import defaultdict, deque
import json, hashlib, csv
import numpy as np
from verify_csp_stable_filter_response import ROOT, BASE, NATIVE, SRC, PREP, groups, W

REVIEW=BASE/'csp_native_filter_response/root-state-review'
OUT=REVIEW/'response-decomposition'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    protocol=REVIEW/'response-decomposition-protocol.json'
    assert sha(protocol)=='ba0dc23d96b63dd3660594b88adc959fb1405d885f870685e1e3f6936fc1f959'
    result_path=REVIEW/'stable-response/result.json'
    verified=json.loads(result_path.read_text());assert verified['pass_root_numerical_and_physical_row_checks']
    OUT.mkdir(exist_ok=False)
    nominal_path=NATIVE/'convergence-diagnostic/fits/iter12_default/fit.log'
    changed_path=SRC/'fits/changed12_default/fit.log'
    nominal,_=groups(nominal_path);changed,_=groups(changed_path)
    rows=[]
    for saved in verified['objects']:
        cid=saved['CID'];old=nominal[cid][-1];new=changed[cid][-1]
        names=[PREP/arm/'CSPDR3_RAISIN'/f'CSPDR3_{cid}.PKMJD.DAT' for arm in ['nominal','known_Jdw_to_j']]
        mapping=defaultdict(list)
        for a,b in zip(*(p.read_text().splitlines() for p in names),strict=True):
            if a.startswith('OBS:'):
                a=a.split();b=b.split();key=(float(a[1]),float(np.float32(a[4])),float(np.float32(a[5])),a[2]);mapping[key].append(b[2])
        new_indices=defaultdict(deque)
        for i,r in enumerate(new['rows']):
            v=r['values'];new_indices[(v[0],v[4],v[5],r['band'])].append(i)
        order=[]
        for r in old['rows']:
            v=r['values'];key=(v[0],v[4],v[5],r['band']);target=mapping[key].pop(0)
            order.append(new_indices[(*key[:3],target)].popleft())
        assert not any(new_indices.values())
        x=np.array([r['values'] for r in old['rows']]);z=np.array([r['values'] for r in new['rows']])[order]
        assert np.array_equal(x[:,4:6],z[:,4:6])
        y=x[:,4];w0=W(old);w1=W(new)[np.ix_(order,order)]
        delta=new['objective'][3]-old['objective'][3]
        m1_at_D0=z[:,2]*10**(.4*delta)
        a0=float(m1_at_D0@w0@y/(m1_at_D0@w0@m1_at_D0));assert a0>0
        fixed_mean=-2.5*np.log10(a0)
        a1=float(z[:,2]@w1@y/(z[:,2]@w1@z[:,2]));assert a1>0
        opt=-2.5*np.log10(a1);assert abs(opt)<=.001
        rows.append(dict(CID=cid,affected=saved['affected'],total_delta_D=delta,
                         fixed_old_C_mean_response=float(fixed_mean),
                         remaining_weight_update_response=float(delta-fixed_mean),
                         final_C_optimum_gap=float(opt)))
    means={str(n):{k:float(np.mean([r[k] for r in rows if n==42 or r['affected']])) for k in
                      ['total_delta_D','fixed_old_C_mean_response','remaining_weight_update_response']}
           for n in [32,42]}
    answer=dict(scope='Order-dependent deterministic processing decomposition; no physical bias identification or corrected cosmology.',
                pass_arithmetic=True,means=means,objects=rows)
    (OUT/'result.json').write_text(json.dumps(answer,indent=2)+'\n')
    with (OUT/'decomposition.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    inputs=[protocol,result_path,nominal_path,changed_path,Path(__file__).with_name('verify_csp_stable_filter_response.py'),Path(__file__).with_name('review_csp_native_states.py')]
    (OUT/'manifest.json').write_text(json.dumps(dict(inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(json.dumps({k:v for k,v in answer.items() if k!='objects'},indent=2))

if __name__=='__main__':main()
