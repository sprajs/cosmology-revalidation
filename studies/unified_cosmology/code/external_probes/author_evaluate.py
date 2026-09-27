"""Conditional public-source reconstruction of released posterior density.

Use recovered CSL CAMB parameter extraction directly through a minimal read-only
block API. This is not a CosmoSIS installation or a claim of runtime identity.
Missing original Dovekie source and unrecorded installed package versions remain
explicit. Compare density differences, allowing only an additive constant.
"""
import os
os.environ.setdefault('CLIPY_NOJAX','1')
import sys,types,importlib.util,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import cho_factor,cho_solve,block_diag
from scipy.special import logsumexp
import camb,candl,candl_data,spt_candl_data
from acquire import sha,RESULTS,ROOT,PACKAGES
from author_config import BASE,CSL,chain,fetch


class Block(dict):
    def __getitem__(self,key):return super().__getitem__(tuple(x.lower() for x in key))
    def __setitem__(self,key,val):super().__setitem__(tuple(x.lower() for x in key),val)
    def has_value(self,section,key):return (section.lower(),key.lower()) in self
    def get_double(self,s,k,default=None):return float(self[s,k]) if self.has_value(s,k) else default
    def get_int(self,s,k,default=None):return int(self[s,k]) if self.has_value(s,k) else default
    def get_bool(self,s,k,default=None):
        if not self.has_value(s,k):return default
        return str(self[s,k]).lower() in ('t','true','1')
    def get_string(self,s,k,default=None):return str(self[s,k]).strip('"') if self.has_value(s,k) else default


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def source_camb(cfg,row,names):
    # This supplies storage/getters only; all CAMB parameter setup is original CSL.
    for key in ['cosmosis','cosmosis.datablock','cosmosis.datablock.cosmosis_py']:
        if key not in sys.modules:sys.modules[key]=types.ModuleType(key)
    sys.modules['cosmosis.datablock'].names=types.SimpleNamespace(cosmological_parameters='cosmological_parameters',halo_model_parameters='halo_model_parameters')
    sys.modules['cosmosis.datablock'].option_section='options'
    sys.modules['cosmosis.datablock.cosmosis_py'].errors=types.SimpleNamespace()
    module=load(BASE/'CSL/boltzmann/camb/camb_interface.py','author_camb')
    opt=Block()
    for k,v in cfg['PARAMS']['camb_act'].items():
        try:v=float(v) if ('.' in v or 'e' in v.lower()) else int(v)
        except ValueError:
            v=v.strip(chr(34))
            if v.lower() in ('t','true','f','false'):v=v.lower() in ('t','true')
        opt['options',k]=v
    conf,more=module.setup(opt);b=Block()
    for sec in cfg['VALUES'].sections():
        for k,v in cfg['VALUES'][sec].items():
            vv=list(map(float,v.split()));b[sec,k]=vv[len(vv)//2]
    for k,v in zip(names,row):
        if '--' in k:b[tuple(k.split('--'))]=float(v)
    h=b['cosmological_parameters','h0'];sec='cosmological_parameters'
    # Exact consistency.py relation at recorded revision, not an approximate93eV conversion.
    nnu=b[sec,'nnu'];omnuh2=b[sec,'mnu']*(nnu/3)**.75/94.06410581217612
    b[sec,'ombh2']=b[sec,'omega_b']*h*h
    b[sec,'omch2']=(b[sec,'omega_m']-b[sec,'omega_b'])*h*h-omnuh2
    params=module.extract_camb_params(b,conf,more)
    return params,more,b


def main():
    names,x,cfg,record=chain('cpl');w=np.exp(x[:,-3]-logsumexp(x[:,-3]))
    # Frozen deterministic retained-row design: max posterior plus weighted w quintiles.
    wi=names.index('cosmological_parameters--w');order=np.argsort(x[:,wi]);cum=np.cumsum(w[order])
    rows=[int(np.argmax(x[:,-1]))]+[int(order[np.searchsorted(cum,q)]) for q in [.1,.3,.7,.9]]
    out={'status':'conditional_reconstruction_not_runtime_identity','rows_selected':rows,'selection':'MAP plus weighted w quantiles .1,.3,.7,.9; no outcome selection',
         'chain_sha256':record['sha256'],'code_sha256':sha(__file__),'rows':[],'limitations':[
        'Recorded cwd Git revision does not certify clean runtime source.',
        'Header Dovekie module is absent from public CSL at recorded revision; released DES module and normalized release data substitute, explicitly.',
        'Installed author CAMB, ACT, candl package versions are not recorded; current pinned versions are used.',
        'Alternative A_planck response is a diagnostic of unresolved source identity, not a recovered author configuration.']}
    planckmod=load(BASE/'modules/planck_lite_py.py','author_planck')
    planck=planckmod.PlanckLitePy(data_directory=BASE/'CSL/likelihood/planck_py/data',year=2018,spectra='TTTEEE',use_low_ell_bins=True,ell_max_tt=1000,ell_max_te=600,ell_max_ee=600)
    spt=candl.Like(spt_candl_data.SPT3G_D1_TnE_lite,feedback=False)
    lens=candl.LensLike(candl_data.SPT3G_2018_Lens_only,feedback=False);lens.priors=[]
    from act_dr6_cmbonly import ACTDR6CMBonly
    act=ACTDR6CMBonly(packages_path=str(PACKAGES))
    import act_dr6_lenslike as al
    lensdata=al.load_data(ddir=str(BASE/'act-lensing/v1.1'),variant='actplanck_baseline',lens_only=False,like_corrections=True,apply_hartlap=True,mock=False,nsims_act=792,nsims_planck=400,trim_lmax=2998)
    # Read literal DESI numbers from downloaded source without executing its CosmoSIS framework.
    import ast
    tree=ast.parse((BASE/'CSL/likelihood/bao/desi-dr2/desi_dr2.py').read_text())
    sets=ast.literal_eval(next(n.value for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='DESI_DATA_SETS' for t in n.targets)))
    selected=['BGS','LRG1','LRG2','LRG3+ELG1','ELG2','QSO','Lya'];means=[];covs=[]
    for key in selected:
        d=sets[key];means.extend(np.atleast_1d(d['mean']));s=np.atleast_1d(d['sigma']);cc=np.diag(s*s)
        if len(s)==2:cc[0,1]=cc[1,0]=d['corr']*s.prod()
        covs.append(cc)
    bao_cov=block_diag(*covs);bao_inv=np.linalg.inv(bao_cov)
    with np.load(ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz') as f:
        z,zhel,mu,cov=[f[k] for k in ['zHD','zHEL','MU','covariance']]
    fac=cho_factor(cov,lower=True);prec=cho_solve(fac,np.eye(len(z)));u=prec@np.ones(len(z));A=u.sum()
    out['native_assets']={'planck_bins':len(planck.data_vector),'spt_internal_prior_parameters':[p.par_names for p in spt.priors],
          'spt_lensing_internal_prior_count':len(lens.priors),'actplanck_lensing_version':'1.1','actplanck_trim_lmax':2998,'SN_count':len(z)}
    out['source_files']={str(p.relative_to(ROOT)):sha(p) for p in sorted(BASE.rglob('*')) if p.is_file() and p.suffix in ['.py','.yaml','.npz','.dat'] and '__pycache__' not in p.parts}
    for index in rows:
        start=time.monotonic();params,more,b=source_camb(cfg,x[index],names)
        r=camb.get_results(params);dl=r.get_total_cls(lmax=9000,CMB_unit='muK');pp=r.get_lens_potential_cls(lmax=9000)[:,0];ell=np.arange(len(dl))
        dls=dict(zip(['TT','EE','BB','TE'],dl.T));dls.update(ell=ell,pp=pp,kk=pp*np.pi/2)
        pt={k.split('--')[-1]:float(v) for k,v in zip(names,x[index])};cal={'Tcal':pt['tcal'],'Ecal':pt['ecal'],'A_fg':pt['a_fg'],'tau':pt['tau']}
        def candl_val(obj):
            sl=slice(obj.ell_min,obj.ell_max+1)
            return float(obj.log_like(cal|{'Dl':{k:v[sl] for k,v in dls.items()}}))
        zz=np.linspace(0,3.01,150);dm=r.comoving_radial_distance(zz);da=dm/(1+zz);hh=r.h_of_z(zz);rd=r.get_derived_params()['rdrag'];dv=(dm*dm*zz/hh)**(1/3)
        # CosmoSIS GaussianLikelihood default interpolation is cubic (source pinned separately).
        predict=[]
        for key in selected:
            d=sets[key];zt=d['z_eff']
            for arr in ([dv] if d['kind']=='d_v' else [dm,1/hh]):predict.append(float(interp1d(zz,arr/rd,kind='cubic')(zt)))
        res=np.array(predict)-means
        snpred=5*np.log10((1+z)*(1+zhel)*interp1d(zz,da,kind='cubic')(z))+25
        delta=snpred-mu;delta-=delta.mean();snchi=float(delta@prec@delta-(delta@u)**2/A)
        clfac=np.zeros(len(ell));clfac[2:]=2*np.pi/(ell[2:]*(ell[2:]+1));clkk=pp*np.pi/2
        clpp=np.zeros(len(ell));clpp[2:]=pp[2:]*2*np.pi/(ell[2:]*(ell[2:]+1))**2
        assert np.allclose(clkk,al.pp_to_kk(clpp,ell),rtol=1e-14,atol=1e-30)
        likes={'SN':-.5*snchi,'DESI':-.5*float(res@bao_inv@res),
            'PlanckPy':float(planck.loglike(dl[2:,0],dl[2:,3],dl[2:,1],calPlanck=1.)),
            'ACT':float(act.loglike({k.lower():dls[k] for k in ['TT','TE','EE']},A_act=1.,P_act=pt['p_act'])),
            'SPT':candl_val(spt),'SPT_lens_only':candl_val(lens),
            'ACTPlanck_lens':float(al.generic_lnlike(lensdata,ell,clkk,ell,dl[:,0]*clfac,dl[:,1]*clfac,dl[:,3]*clfac,dl[:,2]*clfac,2998))}
        total=sum(likes.values());original=float(x[index,-1]-x[index,-2])
        altplanck=float(planck.loglike(dl[2:,0],dl[2:,3],dl[2:,1],calPlanck=pt['a_planck']))
        altact=float(act.loglike({k.lower():dls[k] for k in ['TT','TE','EE']},A_act=pt['a_planck'],P_act=pt['p_act']))
        rec={'row':index,'point':pt,'components':likes,'sum_unnormalized':total,'recorded_post_minus_prior':original,
             'density_difference':total-original,'seconds':time.monotonic()-start,
             'alternative_Aplanck_in_Planck_only_difference':total-likes['PlanckPy']+altplanck-original,
             'alternative_Aplanck_in_Planck_and_ACT_difference':total-likes['PlanckPy']-likes['ACT']+altplanck+altact-original,
             'CAMB_max_l':params.max_l,'rdrag':rd,'YHe':params.YHe}
        out['rows'].append(rec);(RESULTS/'author-density-reconstruction.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rec),flush=True)
    out['density_shape_ranges']={k:float(np.ptp([r[k] for r in out['rows']])) for k in ['density_difference','alternative_Aplanck_in_Planck_only_difference','alternative_Aplanck_in_Planck_and_ACT_difference']}
    out['exact_author_density_established']=False
    (RESULTS/'author-density-reconstruction.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
