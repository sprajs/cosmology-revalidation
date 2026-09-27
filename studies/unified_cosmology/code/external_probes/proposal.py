"""Transform released Planck proposal covariance from100thetaMC to H0.

This is a sampling proposal only. Its source's older BAO sample is not added
to the scientific likelihood. The implicit Jacobian preserves acoustic-scale
correlations at one declared point and does not impose a prior.
"""
import json
import time
import numpy as np
import camb
from adapter import WORK, PACKAGES
from acquire import HERE, RESULTS, sha

REF={'H0':67.36,'ombh2':.02237,'omch2':.12,'w':-1.,'wa':0.}
STEPS={'H0':.005,'ombh2':2e-6,'omch2':1e-5,'w':1e-4,'wa':1e-4}
RENAME={'omegabh2':'ombh2','omegach2':'omch2','theta':'H0','calPlanck':'A_planck',
        'acib217':'A_cib_217','xi':'xi_sz_cib','asz143':'A_sz','aps100':'ps_A_100_100',
        'aps143':'ps_A_143_143','aps143217':'ps_A_143_217','aps217':'ps_A_217_217',
        'aksz':'ksz_norm','kgal100':'gal545_A_100','kgal143':'gal545_A_143',
        'kgal143217':'gal545_A_143_217','kgal217':'gal545_A_217',
        'galfTE100':'galf_TE_A_100','galfTE100143':'galf_TE_A_100_143',
        'galfTE100217':'galf_TE_A_100_217','galfTE143':'galf_TE_A_143',
        'galfTE143217':'galf_TE_A_143_217','galfTE217':'galf_TE_A_217',
        'cal0':'calib_100T','cal2':'calib_217T'}


def theta(point):
    p=camb.set_params(**point,dark_energy_model='ppf',mnu=.06,
                      num_massive_neutrinos=1,nnu=3.044,omk=0.)
    return 100*camb.get_background(p).cosmomc_theta()


def derivatives(scale=1.):
    d={}
    for k,h in STEPS.items():
        a=REF.copy();b=REF.copy();a[k]+=h*scale;b[k]-=h*scale
        d[k]=(theta(a)-theta(b))/(2*h*scale)
    return d


def main():
    first=derivatives();half=derivatives(.5)
    errors={k:abs(first[k]/half[k]-1) for k in first}
    assert max(errors.values())<2e-4,errors
    out={'reference':REF,'reference_100thetaMC':theta(REF),
         'derivatives':half,'step_halving_relative_difference':errors,
         'scope':'Proposal covariance only; no prior or scientific data added.', 'models':{}}
    folder=PACKAGES/'data/planck_supp_data_and_covmats/covmats'
    for model,name in [('lcdm','base_plikHM_TTTEEE_lowl_lowE_lensing.covmat'),
                       ('cpl','base_w_wa_plikHM_TTTEEE_lowl_lowE_BAO.covmat')]:
        p=folder/name
        old_names=p.open().readline().lstrip('#').split()
        names=[RENAME.get(n,n) for n in old_names]
        assert len(set(names))==len(names)
        cov=np.loadtxt(p);jac=np.eye(len(names));i=names.index('H0');jac[i,:]=0.
        jac[i,i]=1/half['H0']
        for key in ['ombh2','omch2','w','wa']:
            if key in names:jac[i,names.index(key)]=-half[key]/half['H0']
        transformed=jac@cov@jac.T;transformed=(transformed+transformed.T)/2
        assert np.linalg.eigvalsh(transformed).min()>0
        back=np.linalg.solve(jac,np.linalg.solve(jac,transformed).T).T
        relative=np.max(abs(back-cov)/np.sqrt(np.outer(np.diag(cov),np.diag(cov))))
        assert relative<1e-10
        q=WORK/f'proposal-{model}.covmat'
        np.savetxt(q,transformed,header=' '.join(names))
        out['models'][model]={'source':str(p.relative_to(WORK)),'source_sha256':sha(p),
            'proposal':str(q.relative_to(WORK)),'proposal_sha256':sha(q),'parameters':names,
            'H0_proposal_sd':float(np.sqrt(transformed[i,i])),
            'inverse_transform_max_correlation_scaled_error':float(relative),
            'minimum_eigenvalue':float(np.linalg.eigvalsh(transformed).min())}
    out['code_sha256']=sha(__file__)
    (RESULTS/'proposal.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
