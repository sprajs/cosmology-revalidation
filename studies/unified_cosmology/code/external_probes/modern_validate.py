"""Validate actual modern CMB data, selected covariance and wrapper precision."""
import os
os.environ.setdefault('CLIPY_NOJAX','1')
import json
import time
from pathlib import Path
import numpy as np
import yaml
from scipy.linalg import solve_triangular
from cobaya.model import get_model
import candl
import spt_candl_data
import candl_data
from adapter import reference_point,PACKAGES,WORK
from modern_adapter import modern_info
from acquire import sha,HERE,RESULTS


def quad(cov,residual):
    q=solve_triangular(np.linalg.cholesky(cov),residual,lower=True)
    return float(q@q)


def official_spt_test(module,filename):
    path=Path(module.__file__).parent/'tests'/filename
    spec=yaml.safe_load(path.read_text())
    dataset=getattr(module,spec['data_set_file'].split('.')[1])
    like=(candl.LensLike if spec['lensing'] else candl.Like)(dataset,feedback=False)
    table=np.loadtxt(path.parent/spec['test_spectrum'])
    params=spec['param_values']|{'Dl':dict(zip(['ell','TT','TE','EE','BB','pp','kk'],table.T))}
    if 'test_chisq' in spec:
        value=float(like.chi_square(params));expected=spec['test_chisq']
    else:
        value=float(like.log_like(params));expected=spec['test_logl']
    assert abs((value-expected)/expected)<1e-3
    return {'expected':expected,'actual':value,'difference':value-expected,
            'official_relative_tolerance':1e-3,'test_file_sha256':sha(path),
            'test_spectrum_sha256':sha(path.parent/spec['test_spectrum'])}


def candl_pars(obj,dls,point):
    lo,hi=obj.ell_min,obj.ell_max+1
    dl={(k if k in ['pp','kk','ell'] else k.upper()):v[lo:hi]
        for k,v in dls.items()}
    dl['kk']=dl['pp']*np.pi/2
    return dict(point,Dl=dl)


def main():
    out={'status':'passed','official_spt_tests':{
        'SPT2025_TnE':official_spt_test(spt_candl_data,'SPT3G_D1_TnE_lite_test.yaml'),
        'SPT2023_lensing':official_spt_test(candl_data,'SPT3G_2018_Lens_and_CMB_test.yaml')}}
    t=time.monotonic();m=get_model(modern_info());out['init_seconds']=time.monotonic()-t
    point=reference_point(m);t=time.monotonic();result=m.logposterior(point)
    out['first_evaluation_seconds']=time.monotonic()-t
    assert np.isfinite(result.logpost)
    likes=dict(zip(m.likelihood,map(float,result.loglikes)))
    out['reference_loglikes']=likes;out['reference_point']={k:float(v) for k,v in point.items()}
    out['reference_is_not_MAP']=True
    dls=m.provider.get_Cl(ell_factor=True)
    act=m.likelihood['act_dr6_cmbonly.ACTDR6CMBonly'];prediction=np.zeros_like(act.data_vec)
    used=[]
    for meta in act.spec_meta:
        pol=meta['pol'];inds=meta['idx'];used.extend(inds)
        v=dls[pol][meta['window'].values]/point['A_planck']**2/point['P_act']**pol.count('e')
        prediction[inds]=meta['window'].weight.T@v
    used=np.array(used);res=act.data_vec[used]-prediction[used]
    chi2=quad(act.covmat[np.ix_(used,used)],res)
    assert abs(chi2+2*likes['act_dr6_cmbonly.ACTDR6CMBonly'])<1e-8
    actualfile=PACKAGES/'data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits'
    out['ACT']={'active_bins':len(used),'all_covariance_rows':len(act.covmat),
        'independent_chi2':chi2,'actual_data_sha256':sha(actualfile),
        'observed_data_source':'NASA LAMBDA public real-data tarball, not bundled simulated fallback',
        'ell_cuts':act.ell_cuts}
    pl=m.likelihood['act_dr6_cmbonly.PlanckActCut'];pred=[]
    for j,k in enumerate(['tt','te','ee']):
        for b in pl.used_bins[j]:
            sl=slice(pl.blmin[b],pl.blmax[b]+1)
            pred.append(dls[k][sl]@pl.weights[sl]/point['A_planck']**2)
    residual=pl.X_data-np.array(pred);selected=np.diag(pl.cov)<1e9
    selected_chi2=quad(pl.cov[np.ix_(selected,selected)],residual[selected])
    deleted_contribution=float(np.sum(residual[~selected]**2/np.diag(pl.cov)[~selected]))
    assert abs(selected_chi2+deleted_contribution+2*likes['act_dr6_cmbonly.PlanckActCut'])<1e-8
    out['cropped_Planck']={'selected_bins':int(selected.sum()),'unused_bins':int((~selected).sum()),
        'independent_selected_chi2':selected_chi2,'finite_large_variance_discarded_bin_chi2':deleted_contribution,
        'scope':'Native wrapper replaces discarded bins by independent variance1e10; contribution retained and quantified.'}
    out['SPT']={}
    for key in ['SPT2025_TnE','SPT2023_lensing']:
        obj=m.likelihood[key].candl_like;pars=candl_pars(obj,dls,point)
        pred=obj.get_model_specs(pars)
        if key=='SPT2025_TnE':pred=obj.bin_model_specs(pred)
        chi=quad(np.asarray(obj.covariance),np.asarray(obj._data_bandpowers-pred))
        double=float(obj.log_like(pars));native=likes[key]
        assert not obj.priors
        assert abs(chi+2*double)<1e-8
        assert abs(native-double)<=abs(np.spacing(np.float32(double)))
        out['SPT'][key]={'bins':len(obj._data_bandpowers),'independent_chi2':chi,
           'native_double_loglike':double,'cobaya_wrapper_loglike':native,
           'float32_wrapper_rounding':native-double,'cleared_internal_prior_count':len(obj.priors),
           'required_nuisance_parameters':obj.required_nuisance_parameters}
    lens=m.likelihood['act_dr6_lenslike.ACTDR6LensLike'];d=lens.data
    from act_dr6_lenslike.act_dr6_lenslike import generic_lnlike,pp_to_kk
    cls=m.provider.get_Cl(ell_factor=False);ell=cls['ell']
    direct,binned=generic_lnlike(d,ell,pp_to_kk(cls['pp'],ell),ell,cls['tt'],cls['ee'],cls['te'],cls['bb'],return_theory=True)
    hartlap=(400-len(binned)-2)/(400-1)
    chi=hartlap*quad(d['cov'],d['data_binned_clkk']-binned)
    assert abs(chi+2*direct)<1e-8 and abs(direct-likes['act_dr6_lenslike.ACTDR6LensLike'])<1e-8
    na=len(d['bcents_act']);corr=d['cov']/np.sqrt(np.outer(np.diag(d['cov']),np.diag(d['cov'])))
    out['joint_ACT_Planck_lensing']={'bins':len(binned),'ACT_bins':na,'Planck_bins':len(binned)-na,
        'hartlap_factor':hartlap,'independent_chi2':chi,
        'max_abs_cross_experiment_correlation':float(np.max(abs(corr[:na,na:]))),
        'standalone_Planck_lensing_present':False,'CMB_dependent_responses_enabled':bool(d['likelihood_corrections'])}
    out['benchmarks']={}
    for key,step in [('H0',.1),('A_planck',.001),('P_act',.001),('A_fg',.1)]:
        times=[]
        for n in range(3):
            p=dict(point);p[key]+=(n+1)*step
            t=time.monotonic();r=m.logposterior(p);times.append(time.monotonic()-t)
            assert np.isfinite(r.logpost)
        out['benchmarks'][key]={'seconds':times,'median_seconds':float(np.median(times))}
    m.logposterior(point);rep=m.logposterior(point)
    assert np.array_equal(rep.loglikes,result.loglikes)
    out['deterministic_repeat']=True
    # Scale-cut sensitivity at exactly the same spectra; no extra CAMB call.
    from act_dr6_cmbonly import ACTDR6CMBonly
    a8500=ACTDR6CMBonly({'ell_cuts':{p:[600,8500] for p in ['TT','TE','EE']}},packages_path=str(PACKAGES))
    out['ACT_8500_sensitivity']={'loglike':float(a8500.loglike(dls,1.,1.)),
                               'active_bins':sum(len(x['idx']) for x in a8500.spec_meta)}
    # Increase numerical accuracy, preserving all scientific inputs.
    more=modern_info();more['theory']['camb']['extra_args'].update(AccuracyBoost=2,lAccuracyBoost=2,lSampleBoost=2)
    t=time.monotonic();mm=get_model(more);rr=mm.logposterior(point)
    out['accuracy_doubling']={'seconds_including_initialization':time.monotonic()-t,
        'loglike_differences':dict(zip(m.likelihood,map(float,rr.loglikes-result.loglikes))),
        'total_loglike_difference':float(np.sum(rr.loglikes-result.loglikes))}
    assert np.isfinite(rr.logpost)
    out['code_sha256']={name:sha(HERE/name) for name in ['modern_adapter.py','modern_validate.py']}
    out['design_sha256']=sha(HERE/'modern-design.json')
    out['scope']='Likelihood acquisition/numerical validation only. Modern likelihood factorization assumes negligible omitted cross-experiment covariance. Public source/configuration ambiguities remain explicit; no author posterior reproduction claimed.'
    (RESULTS/'modern-validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
