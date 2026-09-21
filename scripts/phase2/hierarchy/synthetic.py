"""Record and generate complete/selected synthetic data before real comparisons."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from core import distance_kinematic

ROOT=Path(__file__).resolve().parents[3]
def main():
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--n',type=int,default=1500);p.add_argument('--seed',type=int,default=2026092102);p.add_argument('--family',choices=['dust','gaussian'],default='dust');p.add_argument('--selected',action='store_true');a=p.parse_args()
    out=ROOT/'phase2/hierarchy/synthetic'/a.name;out.mkdir(parents=True,exist_ok=True)
    truth={'M':-19.3,'alpha':.15,'beta':2.2 if a.family=='dust' else 3.0,'mx':0.,'mc':-.06 if a.family=='dust' else .01,'sx':1.,'sc':.05 if a.family=='dust' else .11,'sm':.1,'mx_h':-.35,'rb':3.4,'tau':.085,'q0':-.5,'q1':1.2}
    (out/'configuration.json').write_text(json.dumps({'purpose':'Synthetic recovery before DES model outcomes','arguments':vars(a),'truth':truth,'selection':'synthetic_probit' if a.selected else 'complete','limit':23.4,'width':.3},indent=2)+'\n')
    rng=np.random.default_rng(a.seed);z=rng.uniform(.02,.95,a.n);h=rng.integers(0,2,a.n);mu=np.array(distance_kinematic(z,z,truth['q0'],truth['q1']))
    err=np.column_stack([.035+.05*z/.95,.12+.25*z/.95,.012+.025*z/.95]);rho=np.array([[1,.15,-.12],[.15,1,.05],[-.12,.05,1]])
    cov=err[:,:,None]*err[:,None,:]*rho;chol=np.linalg.cholesky(cov);obs=np.zeros((a.n,3));pending=np.ones(a.n,dtype=bool);generated=0
    while pending.any():
        ids=np.flatnonzero(pending);nn=len(ids);generated+=nn
        x=rng.normal(truth['mx']+truth['mx_h']*h[ids],truth['sx']);ci=rng.normal(truth['mc'],truth['sc'],nn)
        dust=rng.exponential(truth['tau'],nn) if a.family=='dust' else np.zeros(nn)
        mag=mu[ids]+truth['M']-truth['alpha']*x+truth['beta']*ci+truth['rb']*dust+rng.normal(0,truth['sm'],nn)
        candidate=np.column_stack([mag,x,ci+dust])+np.einsum('nij,nj->ni',chol[ids],rng.normal(size=(nn,3)))
        accept=rng.uniform(size=nn)<ndtr((23.4-candidate[:,0])/.3) if a.selected else np.ones(nn,dtype=bool)
        obs[ids[accept]]=candidate[accept];pending[ids[accept]]=False
    np.savez_compressed(out/'data.npz',y=obs,cov=cov,z=z,zhel=z,host_prob=h.astype(float),distance=mu,limit=np.ones(a.n)*23.4,width=np.ones(a.n)*.3)
    (out/'generation.json').write_text(json.dumps({'proposals':generated,'selected':a.n,'note':'Redshift/host/covariance conditioned: rejection loop fills each predetermined object until selected; not a volumetric-rate simulation.'},indent=2)+'\n')
    print(out)
if __name__=='__main__':main()
