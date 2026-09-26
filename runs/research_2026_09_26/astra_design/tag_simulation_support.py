from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parents[2];O=P/'simulation_design';sys.path.insert(0,str(ROOT/'scripts/salt_dust_audit'))
from snana_extinction import SnanaExtinction
lib=ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so';ext=SnanaExtinction(lib);waves=np.arange(2800.,8001.,10.)
for rv in [.7,1.,3.1,6.]:assert np.allclose(ext(waves,rv,option=-99),ext(waves,rv,option=99,historical=True),rtol=0,atol=1e-10)
rows=[]
for arm in ['P21','G10']:
 m=pd.read_csv(O/f'{arm}-cohort.csv',dtype={'CID':str});pick=np.rint(np.linspace(0,len(m)-1,8)).astype(int);(O/f'{arm}-engineering8-cids.txt').write_text('\n'.join(m.iloc[pick].CID)+'\n')
 for q in m.itertuples():
  av=float(q.AV);rv=float(q.RV);d={'arm':arm,'CID':q.CID,'generated_attempt_index':q.generated_attempt_index,'AV':av,'RV':rv}
  if arm=='G10' and av==-9 and rv==-9:d.update(support_class='no_explicit_host_screen_sentinel',A8000_historical=np.nan,min_grid_A_historical=np.nan)
  elif not np.isfinite(av) or av<0:d.update(support_class='invalid_AV',A8000_historical=np.nan,min_grid_A_historical=np.nan)
  elif av==0:d.update(support_class='zero_AV',A8000_historical=0.,min_grid_A_historical=0.)
  elif not np.isfinite(rv) or not .01<=rv<=8:d.update(support_class='RV_outside_audit_range',A8000_historical=np.nan,min_grid_A_historical=np.nan)
  else:
   vals=ext(waves,rv,av/rv,option=99,historical=True);d.update(support_class='negative_passive_attenuation' if vals.min() < -1e-12 else 'nonnegative_declared_grid',A8000_historical=float(vals[-1]),min_grid_A_historical=float(vals.min()))
  rows.append(d)
d=pd.DataFrame(rows);d.to_csv(O/'truth-support-ledger.csv',index=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();paths=[P/'simulation-residual-protocol.md',Path(__file__),lib,ROOT/'scripts/salt_dust_audit/snana_extinction.py']+list(O.glob('*cohort.csv'))+list(O.glob('*cids.txt'))+[O/'truth-support-ledger.csv']
report={'freeze':'Before any pilot simulated epoch residual scoring; cohort IDs unchanged from initial metadata-only freeze.','support_counts':d.groupby(['arm','support_class']).size().to_dict().__str__(),'sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
(O/'handoff-manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(d.groupby(['arm','support_class']).size())
