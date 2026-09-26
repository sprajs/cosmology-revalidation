"""Read the complete native reference export preceding RAISIN_READY, before scores."""
import numpy as np

def check(path,nom,cid,expected_parameters):
    flux={};inv={};objective=None;count=0
    for line in path.read_text().splitlines():
        a=line.split()
        if not a:continue
        if a[0]=='RAISIN_READY:':
            count=int(a[1]);break
        if a[0]=='PHASE2_FLUX:':
            assert a[1]==cid
            flux[int(a[2])]=[a[3]]+[float(x) for x in a[4:]]
        elif a[0]=='PHASE2_COVINV:':
            assert a[1]==cid
            inv[int(a[2])]=[float(x) for x in a[3:]]
        elif a[0]=='PHASE2_OBJECTIVE:':
            assert a[1]==cid
            objective=(int(a[2]),np.array([float(x) for x in a[3:]]))
    assert count and objective is not None and objective[0]==count==len(nom['MJD'])
    assert np.array_equal(objective[1][2:6],expected_parameters)
    ff=[flux[i] for i in range(1,count+1)];a=np.array([r[1:] for r in ff])
    values=dict(band=np.array([r[0] for r in ff]),MJD=a[:,0],data_flux=a[:,4],data_fluxerr=a[:,5],model_flux=a[:,2],W=np.array([inv[i] for i in range(1,count+1)]))
    for k,v in values.items():assert np.array_equal(v,nom[k]),k
    return dict(status='PASS',CID=cid,epochs=count,exact_parameters=objective[1][2:6].tolist(),exact_fields=list(values),before_first_query=True)
