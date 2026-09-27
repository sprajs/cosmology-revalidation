"""Add the released SN design's local Fisher information to a proposal only."""
import json
from pathlib import Path
import numpy as np
import camb
from scipy.linalg import cho_factor,cho_solve
from proposal import REF,STEPS
from adapter import ROOT,WORK
from acquire import RESULTS,sha


def main():
    path=ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz'
    with np.load(path) as data:
        z=data['zHD'];zhel=data['zHEL'];cov=data['covariance']
    factor=cho_factor(cov,lower=True)
    u=cho_solve(factor,np.ones(len(z)));a=u.sum()
    def distance(p):
        cp=camb.set_params(**p,dark_energy_model='ppf',mnu=.06,
            num_massive_neutrinos=1,nnu=3.044,omk=0.)
        bg=camb.get_background(cp)
        dl=bg.angular_diameter_distance(z)*(1+z)*(1+zhel)
        return 5*np.log10(dl)+25
    deriv={};half_errors={}
    for key,h in STEPS.items():
        vals=[]
        for scale in [1.,.5]:
            lo=REF.copy();hi=REF.copy();lo[key]-=h*scale;hi[key]+=h*scale
            col=(distance(hi)-distance(lo))/(2*h*scale);col-=col.mean()
            vals.append(col)
        deriv[key]=vals[-1]
        half_errors[key]=float(np.max(abs(vals[1]-vals[0])))
    records={}
    for model in ['lcdm','cpl']:
        p=WORK/f'proposal-{model}.covmat'
        names=p.open().readline().lstrip('#').split();c=np.loadtxt(p)
        jac=np.column_stack([deriv.get(k,np.zeros(len(z))) for k in names])
        fj=cho_solve(factor,jac);fisher=jac.T@fj-np.outer(jac.T@u,jac.T@u)/a
        root=np.linalg.cholesky(c)
        middle=np.eye(len(c))+root.T@fisher@root
        new=root@np.linalg.solve(middle,root.T);new=(new+new.T)/2
        np.linalg.cholesky(new)
        q=WORK/f'proposal-{model}-dovekie.covmat';np.savetxt(q,new,header=' '.join(names))
        records[model]={'path':str(q.relative_to(ROOT)),'sha256':sha(q),
             'original_sha256':sha(p),'local_sd_before':{k:float(np.sqrt(c[names.index(k),names.index(k)])) for k in STEPS if k in names},
             'local_sd_after':{k:float(np.sqrt(new[names.index(k),names.index(k)])) for k in STEPS if k in names}}
    out={'status':'passed','reference':REF,'rows':len(z),'SN_input_sha256':sha(path),
         'step_halving_max_absolute_derivative_difference':half_errors,'models':records,
         'scope':'Proposal preconditioning only. Observed SN magnitudes are unused; redshifts and fixed covariance inform the local distance Jacobian, after projecting one free magnitude. No scientific likelihood/prior change.',
         'code_sha256':sha(__file__)}
    (RESULTS/'sn-proposal.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
