"""Independent native-export row/hash and amplitude-profile arithmetic check."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import linear_sum_assignment, minimize_scalar
O=Path(__file__).resolve().parent;R=O.parents[3];G=O/'native-profile-gate-v2'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    manifest=json.loads((G/'manifest.json').read_text())
    for f,h in manifest['all_files_sha256'].items():assert sha(G/f)==h,f
    protocol=json.loads((G/'protocol.json').read_text())
    for f,h in protocol['input_hashes'].items():assert sha(R/f)==h,f
    ref=np.load(G/'reference/native.npz',allow_pickle=False)
    nom=np.load(G/'fixed/native.npz',allow_pickle=False)
    raw=[]
    for l in (O/'data/RSR_B/DES16C1cim.snana.dat').read_text().splitlines():
      q=l.split()
      if q and q[0]=='VARLIST:':cols=q[1:]
      elif q and q[0]=='OBS:':raw.append(dict(zip(cols,q[1:])))
    errors=[];matches=[]
    for band in 'griz':
      ii=np.flatnonzero(ref['band']==band);rr=[r for r in raw if r['FLT']==band]
      dist=abs(ref['MJD'][ii,None]-np.array([float(r['MJD']) for r in rr])[None,:])
      aa,bb=linear_sum_assignment(dist)
      assert len(aa)==len(ii) and max(dist[aa,bb],default=0)<1e-10
      for a,b in zip(aa,bb):
        i=ii[a];r=rr[b]
        dif=[abs(float(ref['data_flux'][i])-float(r['FLUXCAL'])),abs(float(ref['data_fluxerr'][i])-float(r['FLUXCALERR']))]
        errors.append(dif);matches.append((band,i,b))
        assert dif==[0.,0.],(band,i,dif)
    assert len(matches)==len(ref['MJD'])
    C=ref['C'];W=ref['W'];h=nom['model_flux'];y=ref['data_flux']
    sym=float(np.max(abs(C-C.T)));inverse=float(np.max(abs(C@W-np.eye(len(y)))))
    L=np.linalg.cholesky((C+C.T)/2);hw=np.linalg.solve(L,h);yw=np.linalg.solve(L,y)
    a=float(hw@yw/(hw@hw));Dref=float(nom['parameters'][0]);dstar=float(Dref-2.5*np.log10(a))
    score=lambda d:float(np.sum((yw-10**(-.4*(d-Dref))*hw)**2))
    opt=minimize_scalar(score,bounds=(Dref-.3,Dref+.3),method='bounded',options={'xatol':1e-12})
    assert opt.success and abs(opt.fun-score(dstar))<1e-8
    native_total=float(ref['total_chi2']);native_prior=float(ref['prior_chi2'])
    rw=np.linalg.solve(L,y-ref['model_flux']);closure=float(rw@rw+native_prior-native_total)
    assert abs(closure)<1e-8
    result=dict(status='PASS: independent raw-row and whitened-arithmetic checks',N_epochs=len(y),
      hashes_checked=len(manifest['all_files_sha256'])+len(protocol['input_hashes']),raw_flux_error_max_abs=np.max(errors,axis=0).tolist(),
      covariance_max_asymmetry=sym,inverse_identity_max_abs=inverse,native_total_reconstruction_error=closure,
      amplitude_from_Cholesky=a,DLMAG_analytic=dstar,DLMAG_direct_scalar=opt.x,scalar_vs_analytic_minimum_objective_difference=float(opt.fun-score(dstar)),
      source_sha256=sha(__file__),proof_manifest_sha256=sha(G/'manifest.json'))
    (O/'native-profile-independent-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
