"""Freeze equal real-data IDs/folds before any conditional held-out scores."""
import datetime,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[3];records=[];cohort=None
for name in ['conditioned-multistart-best','recovered-mask']:
    folder=R/'phase2/hierarchy/data'/name;d=np.load(folder/'data.npz');rows=pd.read_csv(folder/'rows.csv',dtype={'CID':str});keep=(d['survey']==10)&(d['pIa']>.999)
    pairs={str(cid):int(fold) for cid,fold in zip(rows.loc[keep,'CID'],d['fold'][keep])}
    if cohort is None:cohort=pairs
    else:assert pairs==cohort
    records.extend({'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [folder/'data.npz',folder/'rows.csv'])
out=R/'phase2/hierarchy/conditional-cohort.json';assert not out.exists(),'Cohort already frozen'
out.write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'definition':'Original DES membership with published SNNV19>0.999 and available coherent fit in both preregistered arms; no prediction-residual cut.','cid_to_fold':dict(sorted(cohort.items())),'n_total':len(cohort),'n_test':sum(x==0 for x in cohort.values()),'known_multimodal_case':{'CID':'1307748','fold':cohort['1307748'],'primary':'retained','sensitivity':'omit from training for tripp/host/flexible only; never alter test cohort'},'inputs':records,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n');print(len(cohort))
