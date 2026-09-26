"""Independent SVD predictive inverse for frozen observer/SED probes."""
from pathlib import Path
import numpy as np,json,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[2];D=ROOT/'runs/research_2026_09_26/sed_shared_probe';des=json.loads((D/'design.json').read_text());ref=json.loads((D/'score.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();count=0
for p,h in {**des['input_sha256'],**ref['input_sha256']}.items():assert sha(ROOT/p)==h,p;count+=1
assert sha(D/'design.npz')==des['design_sha256'];assert sha(D/'design.json')==ref['design_manifest_sha256']
a=np.load(D/'design.npz',allow_pickle=False);post=np.load(P/'expanded12/posterior-models.npz',allow_pickle=False);data=np.load(P/'expanded12/validation-projected-modes.npz',allow_pickle=False);assert np.array_equal(a['CID'],data['CID'])
A=a['calibration_design'];V=a['posterior_covariance'];U,s,_=np.linalg.svd(A@np.linalg.cholesky(V),full_matrices=False)
def inv(v):return v-U@((s*s/(1+s*s))[:,None]*(U.T@v)) if v.ndim==2 else v-U@((s*s/(1+s*s))*(U.T@v))
y=data['residual']-A@post['systematics_only_discovery_posterior_mean'];vectors=np.column_stack([a['observer'],a['sed'],a['observer']-a['sed'],a['original_observer']]);gram=vectors.T@inv(vectors);err=float(np.max(abs(gram-np.array(des['shared_gram']))));results={}
for i,key in [(0,'observer'),(1,'sed'),(3,'original_observer')]:
 v=vectors[:,i];m=float(v@inv(y));info=float(v@inv(v));g=m-info/2;results[key]={'M':m,'I':info,'G':g};err=max(err,*[abs(results[key][k]-ref['probes'][key][k]) for k in ['M','I','G']])
assert err<1e-9
out={'status':'PASS: independent SVD inverse using prior-root calibration design','hashes_verified':count+2,'maximum_arithmetic_difference':err,'probes':results,'difference_shared_squared_norm':float(gram[2,2]),'difference_fraction_of_observer_information':float(gram[2,2]/gram[0,0]),'shared_vector_cosine':float(gram[0,1]/np.sqrt(gram[0,0]*gram[1,1])),'scope':'These are fixed additive probes in the same null-conditioned covariance. Both full offsets worsen that predictive score; this is not a comparison of fully conditioned/trained physical observer and SED hypotheses. Shared priors change relative distinguishability and cannot identify a correction.'}
(P/'sed-shared-independent.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
