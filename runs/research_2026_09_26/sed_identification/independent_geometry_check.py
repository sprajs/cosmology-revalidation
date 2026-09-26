"""Independent QR/primal-KKT checks; no outcome-residual arrays."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.linalg import helmert
from numpy.polynomial.legendre import leggauss, legvander

HERE = Path(__file__).resolve().parent
arrays = np.load(HERE/'design-matrices.npz', allow_pickle=False)
coeffs = np.load(HERE/'amplitude-constrained-coefficients.npz', allow_pickle=False)
results = json.loads((HERE/'result.json').read_text())
pairs = [(n,0) for n in range(1,9)]+[(n,1) for n in range(9)]+[(n,2) for n in range(9)]
small = [pairs.index(p) for p in [(1,0),(2,0),(3,0),(1,1)]]
gauge = helmert(4,full=False).T
ledger = []
for name,ids,drift in [('small_shared',small,False),('small_drift',small,True),('broad_shared',list(range(26)),False),('broad_drift',list(range(26)),True)]:
    x,wx=leggauss(48); p,wp=leggauss(40)
    wv=legvander(x,8); pv=legvander(np.tanh((15+30*p)/20),2)
    cols=np.stack([(wv[:,n,None]*pv[None,:,q]).ravel() for n,q in pairs],axis=1)[:,ids]
    weight=(wx[:,None]*wp[None,:]/4).ravel()
    metric=cols.T@(weight[:,None]*cols)
    if drift: metric=np.block([[metric,np.zeros_like(metric)],[np.zeros_like(metric),metric/3]])
    for cohort in ['discovery','validation']:
        s=arrays[cohort+'_S'][:,ids]
        if drift: s=np.column_stack([s,s*arrays[cohort+'_v'][:,None]])
        o=arrays[cohort+'_O']
        qs,_=np.linalg.qr(s,mode='reduced');qo,_=np.linalg.qr(o,mode='reduced')
        rho=np.linalg.svd(qs.T@qo,compute_uv=False)
        archived=results['families'][name][cohort]
        rhoerr=float(np.max(abs(rho-archived['canonical_correlations'])))
        assert rhoerr<1e-8
        assert archived['numerical_rank']==s.shape[1]
        directions=gauge.T@np.array(archived['canonical_griz_dimming_mag']).T
        targets=np.column_stack([o@directions,arrays[cohort+'_fixed']])
        for limit in [.01,.02,.05,.10]:
            theta=coeffs[f'{name}_{cohort}_rms{limit}']
            saved=archived['hard_rms_amplitude_frontier'][str(limit)]
            penalty=np.array(saved['penalty_lagrange_multiplier'])
            rms=np.sqrt(np.einsum('ik,ij,jk->k',theta,metric,theta))
            kkt=s.T@(s@theta-targets)+(metric@theta)*penalty
            kktrel=float(np.linalg.norm(kkt)/max(np.linalg.norm(s.T@targets),1))
            assert kktrel<1e-8,(name,cohort,limit,kktrel)
            assert rms.max()<=limit+1e-8
            assert np.allclose(rms,saved['rms_mag'],atol=1e-10)
            rem=np.sum((targets-s@theta)**2,axis=0)/np.sum(targets**2,axis=0)
            assert np.allclose(rem,saved['remaining_information_fraction'],atol=1e-10)
            ledger.append(dict(family=name,cohort=cohort,rms_limit=limit,qr_canonical_max_error=rhoerr,primal_kkt_relative_error=kktrel,max_rms=float(rms.max())))
out={'checks':ledger,'max_qr_canonical_error':max(x['qr_canonical_max_error'] for x in ledger),'max_primal_kkt_relative_error':max(x['primal_kkt_relative_error'] for x in ledger),'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'design-matrices.npz',HERE/'amplitude-constrained-coefficients.npz',HERE/'result.json']}}
(HERE/'independent-geometry-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['checks','inputs_sha256']},indent=2))
