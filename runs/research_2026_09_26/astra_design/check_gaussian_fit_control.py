"""Independent CSV arithmetic and scope check; no refits or score selection."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;D=ROOT/'runs/research_2026_09_26/gaussian_fit_control'
a=pd.read_csv(D/'draws.csv',dtype={'CID':str});result=json.loads((D/'result.json').read_text());meta=pd.read_csv(D/'cohort.csv',dtype={'CID':str});assert len(a)==768 and a.CID.nunique()==6 and not a.duplicated(['CID','replicate']).any();assert set(a.replicate)==set(range(128));assert a.success.all() and not a.boundary_hit.any()
allchecks={};rows=[]
for name,cols in {'paired_product_difference':('nonlinear_matched_product','linear_matched_product'),'paired_information_difference':('nonlinear_information','linear_information'),'paired_gain_difference':('nonlinear_gain','linear_gain')}.items():
 q=(a[cols[0]]-a[cols[1]]).to_numpy().reshape(6,128);assert np.array_equal(a.replicate.to_numpy().reshape(6,128),np.tile(np.arange(128),(6,1)))
 sums=q.sum(axis=0);mean=float(sum(sums)/128);se=float(np.sqrt(np.sum((sums-mean)**2)/(127*128)));stored=result['six_object_sum_statistics'][name]
 assert abs(mean-stored['mean'])<1e-12 and abs(se-stored['Monte_Carlo_SE'])<1e-12
 allchecks[name]={'mean':mean,'MCSE':se,'max_arithmetic_discrepancy':max(abs(mean-stored['mean']),abs(se-stored['Monte_Carlo_SE']))}
 for cid,v in zip(a.CID.unique(),q):rows.append(dict(CID=cid,metric=name,mean=float(v.mean()),MCSE=float(v.std(ddof=1)/np.sqrt(128))))
assert np.allclose(a.linear_gain,a.linear_matched_product-.5*a.linear_information,rtol=0,atol=1e-12)
assert np.allclose(a.nonlinear_gain,a.nonlinear_matched_product-.5*a.nonlinear_information,rtol=0,atol=1e-12)
second=a.second_start_success.dropna();assert len(second)==48 and second.all();assert a.second_start_objective_difference.abs().max()<1e-5
chi=[float(a[f'{k}_chi2'].sum()/a.dof.sum()) for k in ['linear','nonlinear']];assert np.allclose(chi,result['linear_and_nonlinear_pooled_chi2_per_dof'],rtol=0,atol=1e-12)
report={'status':'PASS: saved-draw arithmetic and paired control interpretation','fits':len(a),'second_start_checks':len(second),'second_start_max_objective_difference':float(a.second_start_objective_difference.abs().max()),'pooled_chi2_per_dof':chi,'paired_six_object_sums':allchecks,'max_independent_mean_vs_SNANA_fraction':float(meta.native_vs_exported_mean_max_fraction.max()),'generator_scope':'Both generator mean and fitted nonlinear mean use the same independent sncosmo SALT3 implementation; covariance is exact exported native SNANA. This isolates nonlinear fitting within that mean family and is not SNANA fitting/selection reproduction.','limits':'Six previously fixed discovery objects,128 draws each, frozen C and published masks; no model-C feedback, priors, clipping or classifier regeneration, no survey-wide bias bound.','sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'draws.csv',D/'result.json',D/'cohort.csv',D/'executed_source.py',Path(__file__)]}}
(P/'gaussian-control-independent.json').write_text(json.dumps(report,indent=2)+'\n');pd.DataFrame(rows).to_csv(P/'gaussian-control-independent-per-object.csv',index=False);print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))
