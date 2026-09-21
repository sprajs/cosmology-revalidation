"""Freeze common object/field splits before opening real population-fit outcomes."""
from pathlib import Path
import hashlib,json,datetime
import pandas as pd
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3];R=ROOT/'sources/repos/des-science__DES-SN5YR@1.3';O=ROOT/'phase2/hierarchy';O.mkdir(parents=True,exist_ok=True)
meta=R/'4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv';d=pd.read_csv(meta,dtype={'CID':str})
phot=R/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_PHOT.FITS.gz';head=R/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
with fits.open(head) as f:h=f[1].data.copy()
with fits.open(phot) as f:p=f[1].data
field={str(sn).strip():str(p['FIELD'][lo-1]).strip() for sn,lo in zip(h['SNID'],h['PTROBS_MIN'])}
d=d[['CID','IDSURVEY','zHD','HOST_LOGMASS']].copy();d['field']=[field.get(cid,'LOWZ') if s==10 else 'LOWZ' for cid,s in zip(d.CID,d.IDSURVEY)]
d['fold']=[int(hashlib.sha256(('DES-PHASE2-20260920:'+str(s)+':'+cid).encode()).hexdigest()[:16],16)%5 for cid,s in zip(d.CID,d.IDSURVEY)]
d['depth']=np.where(d.field.isin(['C3','X3']),'DEEP',np.where(d.IDSURVEY==10,'SHALLOW','LOWZ'))
d['z_band']=pd.cut(d.zHD,[0,.1,.3,.5,.7,1.3],labels=['0-.1','.1-.3','.3-.5','.5-.7','.7-1.3']).astype(str)
d['primary_train']=d.fold.ne(0);d['primary_test']=d.fold.eq(0)
d['field_transfer_test']=d.field.isin(['C3','X1'])
d['redshift_transfer_train']=d.zHD.lt(.5);d['redshift_transfer_test']=d.zHD.ge(.5)
out=O/'folds.csv';d.to_csv(out,index=False)
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
(O/'folds-manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Preregistered object-level predictive partitions; shared across all correction models','inputs_sha256':{str(x.relative_to(ROOT)):sha(x) for x in [meta,head,phot,Path(__file__)]},'outputs_sha256':{str(out.relative_to(ROOT)):sha(out)},'configuration':{'primary':'hash modulo5, fold0test, folds1-4train','secondary':'all5foldcrossvalidation ifcomputationallyfeasible; C3+X1fieldholdout; z>=.5redshiftholdout','selection_notice':'This freezes official sample IDs, not a statement all fits succeeded or selection is ignorable. Any later exclusion requires explicit reason and identical model sample.'},'counts':d.groupby(['fold','depth']).size().to_dict().__str__()},indent=2)+'\n')
print(d.groupby(['fold','depth']).size())
