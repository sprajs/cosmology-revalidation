"""Recover immutable released posterior headers; test declared prior shape.

This reads author samples, not a new independent posterior. A Git HEAD recorded
in a chain does not certify the absence of local source changes at execution.
"""
import configparser
import hashlib
import json
import urllib.request
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.stats import truncnorm
from acquire import ROOT, WORK, RESULTS, HERE, sha

RELEASE='c9a4fcafc4cbd19bd750dee47fc76194a45c181f'
CSL='a87c01b0795eb643a97bc03ca38f4f9eb94a01b6'
BASE=WORK/'author-chains'


def fetch(url,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(url,timeout=120) as r:path.write_bytes(r.read())
    return {'url':url,'path':str(path.relative_to(ROOT)),'sha256':sha(path),'bytes':path.stat().st_size}


def chain(model):
    tag='w0wacdm' if model=='cpl' else 'lcdm'
    name=f'dovekie_bao_cmbactspt_{tag}_nautilus.txt'
    folder='flcdm' if model=='lcdm' else tag
    path=BASE/name
    record=fetch(f'https://raw.githubusercontent.com/des-science/DES-SN5YR/{RELEASE}/5_COSMOLOGY/chains/{folder}/{name}',path)
    tree=json.loads((ROOT/'.work/unified-cosmology/survey-selection/release-tree.json').read_text())['tree']
    blob=next(x['sha'] for x in tree if x['path']==f'5_COSMOLOGY/chains/{folder}/{name}')
    assert hashlib.sha1(f'blob {path.stat().st_size}\0'.encode()+path.read_bytes()).hexdigest()==blob
    record['github_blob_sha1']=blob
    header=[]
    with path.open() as f:
        names=f.readline().lstrip('#').split()
        for line in f:
            if not line.startswith('#'):break
            header.append(line)
    header=''.join(header);configs={}
    for kind in ['PARAMS','VALUES','PRIORS']:
        txt=header.split(f'## START_OF_{kind}_INI\n')[1].split(f'## END_OF_{kind}_INI')[0]
        txt=''.join(x[3:] for x in txt.splitlines(keepends=True))
        (BASE/f'{model}-{kind.lower()}.ini').write_text(txt)
        cfg=configparser.ConfigParser(interpolation=None);cfg.read_string(txt);configs[kind]=cfg
    return names,np.loadtxt(path),configs,record


def main():
    out={'status':'released_chain_headers_recovered_runtime_not_certified','release_commit':RELEASE,'header_cwd_git_revision':CSL,'models':{},'code_sha256':sha(__file__)}
    for model in ['lcdm','cpl']:
        names,x,cfg,source=chain(model);n=len(names)-3
        ww=np.exp(x[:,-3]-logsumexp(x[:,-3]));mu=ww@x[:,:n]
        cov=(x[:,:n]-mu).T@(ww[:,None]*(x[:,:n]-mu))
        pred=np.zeros(len(x));gaussians={}
        for i,key in enumerate(names[:n]):
            sec,name=key.split('--');bounds=list(map(float,cfg['VALUES'][sec][name].split()));lo,_,hi=bounds
            p=cfg['PRIORS'].get(sec,name,fallback='uniform '+str(lo)+' '+str(hi)).split()
            if p[0]=='gaussian':
                mean,sd=map(float,p[1:]);pred+=truncnorm.logpdf(x[:,i],(lo-mean)/sd,(hi-mean)/sd,loc=mean,scale=sd)
                gaussians[key]={'mean':mean,'sigma':sd,'truncation':[lo,hi]}
            else:pred+=-np.log(hi-lo)
        diff=pred-x[:,-2]
        med=float(np.median(diff));shape=float(np.max(abs(diff-med)))
        assert shape<1e-9, (model,shape)
        out['models'][model]={'source':source,'rows':len(x),'importance_weight_ESS':float(1/(ww@ww)),
            'weighted_mean_sd':{k:{'mean':float(m),'sd':float(s)} for k,m,s in zip(names[:n],mu,np.sqrt(np.diag(cov)))},
            'max_posterior_row':int(np.argmax(x[:,-1])), 'active_modules':cfg['PARAMS']['pipeline']['modules'].split(),
            'active_module_options':{k:dict(cfg['PARAMS'][k]) for k in cfg['PARAMS']['pipeline']['modules'].split()},
            'values':{k:dict(cfg['VALUES'][k]) for k in cfg['VALUES'].sections()},'gaussian_priors':gaussians,
            'prior_reconstruction':{'max_shape_error':shape,'constant_offset':med,'all_rows_checked':len(x)},
            'scope':'Weighted moments of published Nautilus samples; not independently inferred measurements.'}
        # A proposal only: density transform is not used to reweight samples.
        ind={k.split('--')[-1]:i for i,k in enumerate(names[:n])};h=x[:,ind['h0']]
        import camb
        p=camb.CAMBparams();p.set_cosmology(H0=67.,ombh2=.022,omch2=.12,mnu=.06,nnu=3.046,num_massive_neutrinos=3)
        columns={'H0':100*h,'ombh2':x[:,ind['omega_b']]*h*h,
            'omch2':(x[:,ind['omega_m']]-x[:,ind['omega_b']])*h*h-p.omnuh2,
            'logA':np.log(1e10*x[:,ind['a_s']]),'ns':x[:,ind['n_s']],'tau':x[:,ind['tau']]}
        if model=='cpl':columns.update(w=x[:,ind['w']],wa=x[:,ind['wa']])
        columns.update({b:x[:,ind[a]] for a,b in [('a_planck','A_planck'),('p_act','P_act'),('tcal','Tcal'),('ecal','Ecal'),('a_fg','A_fg')]})
        arr=np.array(list(columns.values())).T;avg=ww@arr;c=(arr-avg).T@(ww[:,None]*(arr-avg));np.linalg.cholesky(c)
        target=WORK/f'proposal-author-{model}.covmat';np.savetxt(target,c,header=' '.join(columns),fmt='%.17e')
        ncosmo=8 if model=='cpl' else 6
        cosmofile=WORK/f'proposal-author-{model}-cosmology.covmat'
        np.savetxt(cosmofile,c[:ncosmo,:ncosmo],header=' '.join(list(columns)[:ncosmo]),fmt='%.17e')
        out['models'][model]['proposal_only']={'cosmology_path':str(cosmofile.relative_to(ROOT)),'cosmology_sha256':sha(cosmofile),'path':str(target.relative_to(ROOT)),'sha256':sha(target),
            'parameters':list(columns),'transformed_mean':dict(zip(columns,map(float,avg))),
            'transformed_max_posterior_point':dict(zip(columns,map(float,arr[np.argmax(x[:,-1])]))),
            'warning':'Initialization/proposal only; not a prior or an independent posterior; nuisance priors/configuration differ from our modern target.'}
    out['configuration_findings']=[
       'Active module list has no Commander or simall; PlanckPy has two Gaussian low-temperature bins and cropped high-l TTTEEE.',
       'External P_act prior is Gaussian sigma0.1, not sigma0.003; shape of every recorded prior is independently reproduced.',
       'Header fixes ACT A_act=1; public ACT wrapper reads this directly, not A_planck.',
       'Header retains SPT primary internal priors; historical YAML is audited separately before claiming which priors ran.',
       'Header selects SPT lens_only variant despite other CMB measurements in joint pipeline.',
       'Header sets fixed YHe0.245341, N_eff3.046 and3 massive neutrinos; our independent modern target uses BBN helium,3.044 and1 massive.',
       'Public PlanckPy interface at recorded cwd Git HEAD does not pass A_planck. Git HEAD alone cannot certify runtime source; calibration response requires retained-row closure.'
    ]
    (RESULTS/'author-configuration.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v['prior_reconstruction'] for k,v in out['models'].items()},indent=2))

if __name__=='__main__':main()
