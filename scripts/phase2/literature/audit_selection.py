#!/usr/bin/env python3
from pathlib import Path
import json,re,hashlib,gzip,sys
import numpy as np,pandas as pd
from astropy.io import fits
from selection_assets import PopulationPDF,HostEfficiency,DetectionEfficiency,importance_diagnostics,ROOT,AUX
OUT=ROOT/'phase2/literature';OUT.mkdir(parents=True,exist_ok=True)
SIM=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS'
rows=[];files=[];frames=[]
def record(p):
 files.append({'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
for kind in ['SNIa_SIMULATIONS','SNnonIa_SIMULATIONS']:
 for directory in sorted((SIM/kind).iterdir()):
  if not directory.is_dir():continue
  readme=next(directory.glob('*.README'));txt=readme.read_text();record(readme)
  stats=re.search(r'- TOTAL\s+(\d+)\s+(\d+)\s+(\d+)',txt)
  dmp=next(directory.glob('*.DUMP'));record(dmp)
  count=sum(line.startswith('SN:') for line in dmp.open())
  selection=next((l.strip() for l in dmp.open() if 'SELECTION:' in l),'unknown')
  heads=list(directory.glob('*HEAD.FITS.gz'));head_count=0;masks=set()
  for p in heads:
   record(p)
   with fits.open(p) as h:
    data=h[1].data;head_count+=len(data);masks.update(map(int,np.unique(data['SIM_SEARCHEFF_MASK'])))
    if kind=='SNIa_SIMULATIONS':
     names=[n for n in data.names if n.startswith('SIM_') and data[n].dtype.kind in 'fi']+['SNID','HOSTGAL_LOGMASS','HOSTGAL_COLOR','HOSTGAL_LOGSFR','HOSTGAL_LOGsSFR','PTROBS_MIN']
     df=pd.DataFrame({n:np.asarray(data[n]).astype(str if n=='SNID' else float) for n in names});df['realization']=directory.name
     photfile=next(directory.glob('*PHOT.FITS.gz'))
     with fits.open(photfile) as ph:df['field']=[str(x).strip() for x in ph[1].data['FIELD'][data['PTROBS_MIN'].astype(int)-1]]
     frames.append(df)
  generated,written,spec=map(int,stats.groups())
  assert count==head_count==written,(directory,count,head_count,written)
  rows.append({'kind':kind,'realization':directory.name,'generated':generated,'written':written,'dump_rows':count,'head_rows':head_count,'selection':selection,'effmask_values':sorted(masks),'global_acceptance':written/generated})
pop=PopulationPDF();host=HostEfficiency();det=DetectionEfficiency()
truth=pd.concat(frames,ignore_index=True)
truth['EBV_truth']=truth.SIM_AV/truth.SIM_RV
truth['logp0_grid_candidate']=pop.log_joint(truth.SIM_SALT2c,truth.SIM_SALT2x1,truth.SIM_RV,truth.EBV_truth,truth.SIM_SALT2beta,truth['SIM_HOSTLIB(LOGMASS_TRUE)'],truth['SIM_HOSTLIB(ZTRUE)'])
truth['host_eff_grid']=np.array([host.probability(r,g,f,t) for r,g,f,t in zip(truth['SIM_HOSTLIB(r_obs_auto)'],truth['SIM_HOSTLIB(obs_gr_auto)'],truth.field,truth.SIM_PEAKMJD)])
truth.to_csv(OUT/'selected-ia-truth.csv.gz',index=False,compression={'method':'gzip','mtime':0})
summary={'status':'empirical_asset_audit_not_complete_selection_likelihood','simulations':rows,'totals':{},'truth_rows':len(truth),'truth_warning':'SIM_SALT2c and SIM_SALT2mB describe pre-host-dust SALT source. Neither equals a fitted noiseless SALT summary. No final fitted-SALT selection mask is released in these HEAD/DUMP files.','p0_status':'candidate_grid_density_only; old source file equality and full forward-map absolute continuity unproven','finite_p0':int(np.isfinite(truth.logp0_grid_candidate).sum()),'nonfinite_p0':int((~np.isfinite(truth.logp0_grid_candidate)).sum()),'coherent_scatter_values':np.unique(truth.SIM_MAGSMEAR_COH).tolist(),'mB_minus_mu_plus_alpha_x1_minus_beta_c_quantiles':np.quantile(truth.SIM_SALT2mB-truth.SIM_DLMU+truth.SIM_SALT2alpha*truth.SIM_SALT2x1-truth.SIM_SALT2beta*truth.SIM_SALT2c,[0,.5,1]).tolist(),'host_eff_quantiles':dict(zip(['min','p01','median','p99','max'],np.quantile(truth.host_eff_grid,[0,.01,.5,.99,1]).tolist())),'assets':files}
for kind in ['SNIa_SIMULATIONS','SNnonIa_SIMULATIONS']:
 rs=[r for r in rows if r['kind']==kind];g=sum(r['generated'] for r in rs);w=sum(r['written'] for r in rs)
 summary['totals'][kind]={'realizations':len(rs),'generated':g,'written':w,'global_acceptance':w/g}
summary['p0_grid']={k:{'names':v[0],'axes':[[float(a[0]),float(a[-1]),len(a)] for a in v[1]]} for k,v in pop.maps.items()}
# A population-only importance check, never a cosmological change or selection validation.
from scipy.stats import norm
c=truth.SIM_SALT2c.to_numpy();p0=pop.density('SALT2c',c)
summary['colour_reweighting_diagnostics']={}
for shift in [0,.02,.05,.1]:
 target=norm.pdf(c,loc=-.0736862+shift,scale=.05236237)
 with np.errstate(divide='ignore'):logw=np.log(target)-np.log(p0)
 # Analytic Gaussian vs rounded table mismatch is explicitly present at shift=0.
 summary['colour_reweighting_diagnostics'][str(shift)]=importance_diagnostics(logw)
for p in [pop.path,host.path,det.path,AUX/'models/searcheff/SEARCHEFF_PIPELINE_LOGIC.DAT']:record(p)
(OUT/'selection-audit.json').write_text(json.dumps(summary,indent=2)+'\n')
pd.DataFrame(rows).to_csv(OUT/'simulation-denominators.csv',index=False)
print(json.dumps({k:v for k,v in summary.items() if k not in ['assets','simulations']},indent=2))
