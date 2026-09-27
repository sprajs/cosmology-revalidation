"""Analytic backgrounds, independent GLS and synthetic typed-cache controls."""
import argparse
import copy
from contextlib import ExitStack
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
from unittest.mock import patch
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
import native_accuracy_history as h
from native_accuracy_qualify import QualifiedNativeAccuracy2
from luminosity_sensitivity import IntegratedLuminosity, DESIGN_PATH
from exact_correction import summarize


def refused(fn):
    try:fn()
    except (AssertionError,ValueError,RuntimeError,FileExistsError,KeyError):return True
    raise AssertionError('Invalid synthetic history input accepted.')


def analytic(z,H0,om,ore=9e-5):
    x=1+np.asarray(z);f=om*x**3+ore*x**4+1-om-ore
    fp=3*om*x*x+4*ore*x**3;fpp=6*om*x+12*ore*x*x
    return {'H_km_s_Mpc':H0*np.sqrt(f),'q':.5*x*fp/f-1,'j':1-x*fp/f+.5*x*x*fpp/f}


class Parameters(SimpleNamespace):
    def copy(self):return copy.deepcopy(self)


class Background:
    def __init__(self,parameters):
        from camb import constants
        self.Params=parameters;self.ident=id(self);self.calls=[]
        self.age=constants.Mpc/1000/constants.Gyr*quad(lambda z:1/(float(self.hubble_parameter(z))*(1+z)),0,np.inf,epsabs=1e-12,epsrel=1e-11)[0]
    def hubble_parameter(self,z):return analytic(z,self.Params.H0,self.Params.omegam)['H_km_s_Mpc']
    def get_background_densities(self,a):
        from camb import constants
        a=np.asarray(a);return {'tot':3*a**4*(self.hubble_parameter(1/a-1)*1000/constants.c)**2}
    def get_derived_params(self):return {'age':self.age,'rdrag':147+.5*(self.Params.omegam-.3)}
    def physical_time(self,z):assert z==0;return self.age
    def angular_diameter_distance(self,z):
        self.calls.append(('distance',self.ident,np.asarray(z).tolist()))
        return np.array([299792.458*quad(lambda t:1/float(self.hubble_parameter(t)),0,float(v),epsabs=1e-11,epsrel=1e-11)[0]/(1+v) for v in z])


EXTRA={'AccuracyBoost':2,'lAccuracyBoost':2,'lSampleBoost':2,'lens_potential_accuracy':4,
       'lmax':40,'lens_margin':5,'halofit_version':'mead2016','omk':0.,'dark_energy_model':'ppf'}
ACTUAL={k:v for k,v in dict(EXTRA,lmax=41).items() if k!='halofit_version'}
META={'finalized_theory_extra_args':ACTUAL,'CAMB_Params_max_l':46,'provider_length':50,
      'ell_factor':False,'units':'FIRASmuK2','pp':'dimensionless C_phi_phi'}


def parameters(point):
    from camb import model
    return Parameters(H0=point['H0'],omegam=(point['ombh2']+point['omch2']+.00064)/(point['H0']/100)**2,
        max_l=46,WantCls=True,WantScalars=True,NonLinear=model.NonLinear_lens,DoLensing=True,SourceWindows=[],WantTransfer=False)


def point_at(x,evolution='none'):
    value={'H0':68+.2*x,'ombh2':.0224,'omch2':.118+.0005*x,'ns':.965,'tau':.055,'logA':3.04}
    if evolution=='linear':value['epsilon']=.04+.01*x
    return value


def native_for(background,point,sn,kernel):
    exact=analytic(np.array([0.,.5,1.]),background.Params.H0,background.Params.omegam)
    derived={name+label:float(exact[name][i]) for name in ['q','j'] for i,label in enumerate(['0','05','1'])}
    da=background.angular_diameter_distance(sn.z)
    prediction=5*np.log10(da*(1+sn.z)*(1+sn.zhel))+25
    sigma={'none':0.,'linear':0.,'smooth01':.1,'smooth03':.3}[kernel.evolution]
    # Independent full inverse-covariance GLS with explicit common intercept.
    covariance=sn.factor@sn.factor.T+sigma*sigma*sn.basis@sn.basis.T
    inverse=np.linalg.inv(covariance);ones=np.ones(sn.n);u=inverse@ones;A=float(ones@u)
    residual=prediction+point.get('epsilon',0.)*np.log1p(sn.z)/np.log(2)-sn.observed
    centered=residual-residual.mean();chi2=float(centered@inverse@centered-(u@centered)**2/A)
    loglike=-.5*(chi2+np.linalg.slogdet(covariance)[1]+np.log(A)+(sn.n-1)*np.log(2*np.pi))
    derived.update(omegam=background.Params.omegam,rdrag=background.get_derived_params()['rdrag'],sn_chi2=chi2)
    return {'point':point,'native2_derived':derived,'native2_loglikes':{'released_sn':float(loglike)},
            'value':{'H':background.hubble_parameter(h.rt.GRID).tolist(),'metadata':copy.deepcopy(META)}}


def gls_tests(rng):
    rows=[]
    for i in range(16):
        n=20+i;z=np.sort(rng.uniform(.01,1.15,n));zhel=z+.0005
        R=rng.normal(size=(n,n));C=(R@R.T/n+np.eye(n))*.015
        residual=rng.normal(size=n)*.2;sn=IntegratedLuminosity(z,zhel,-residual,C)
        for mode,sigma in [('smooth01',.1),('smooth03',.3)]:
            kernel=h.LuminosityKernel(sn,mode,h.luminosity.REDSHIFTS);value=kernel.from_prediction(np.zeros(n),{})
            D=np.column_stack([np.ones(n),sn.basis]);P=np.linalg.inv(C)
            cov=np.linalg.inv(D.T@P@D+np.diag([0.]+[1/sigma**2]*4));mean=-cov@D.T@P@residual
            S=C+sigma*sigma*sn.basis@sn.basis.T;Q=np.linalg.inv(S);u=Q@np.ones(n);A=u.sum()
            score=residual@Q@residual-(u@residual)**2/A
            loglike=-.5*(score+np.linalg.slogdet(S)[1]+np.log(A)+(n-1)*np.log(2*np.pi))
            shifted=kernel.from_prediction(np.full(n,53.),{})
            errors={'mean':float(np.max(abs(np.array(value['coefficient_mean_mag'])-mean[1:]))),
                    'covariance':float(np.max(abs(kernel.coefficient_covariance-cov[1:,1:]))),
                    'density':abs(value['SN_loglike']-loglike),'score':abs(value['SN_chi2']-score),
                    'offset':float(np.max(abs(np.array(value['mean_mag'])-shifted['mean_mag'])))}
            assert max(errors.values())<2e-10
            rows.append(errors)
    return {'cases':32,'maximum_errors':{k:max(row[k] for row in rows) for k in rows[0]}}


def main_run():
    import camb
    import cobaya.model
    rng=np.random.default_rng(274029);rootbase=h.ROOT/'.work/unified-cosmology/native-accuracy-history-validation'
    rootbase.mkdir(parents=True,exist_ok=True);root=Path(tempfile.mkdtemp(prefix='case-',dir=rootbase))
    physical=[]
    def forbidden(*a,**kw):physical.append('forbidden');raise RuntimeError('No physical calls in synthetic validator.')
    with ExitStack() as guard:
        for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),
                           (camb,'set_params'),(cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
            guard.enter_context(patch.object(owner,name,forbidden))
        design=json.loads(h.DESIGN.read_text());old=json.loads(h.expansion.DESIGN.read_text())
        gls=gls_tests(rng);mixtures=h.luminosity.validate()
        z=np.linspace(.01,1.1,36);zhel=z+.0003*np.sin(17*z)
        R=rng.normal(size=(36,36));C=(R@R.T/36+np.eye(36))*.008
        base=Background(parameters(point_at(0)))
        observed=5*np.log10(base.angular_diameter_distance(z)*(1+z)*(1+zhel))+25+rng.normal(0,.04,len(z))
        sn=IntegratedLuminosity(z,zhel,observed,C)
        checks=[];templates={};natives={};shared_calls=[]
        for mode in ['none','linear','smooth01','smooth03']:
            kernel=h.LuminosityKernel(sn,mode,design['luminosity_redshifts']);kernel.declared_extra=EXTRA
            for x in [-1.5,0.,1.5]:
                point=point_at(x,mode);bg=Background(parameters(point));native=native_for(bg,point,sn,kernel)
                value=h.consume_background(bg,point,native,kernel,design,old)
                assert not value['failed_numerical_gates']
                exact=analytic(np.array(old['redshift_grid']),bg.Params.H0,bg.Params.omegam)
                errors={k:float(np.max(abs(np.array(value['history'][k])-exact[k]))) for k in exact}
                assert errors['q']<1e-7 and errors['j']<2e-5
                assert abs(value['scalar']['age_Gyr']-bg.age)<1e-7
                checks.append({'mode':mode,'x':x,'errors':errors,'numerical':value['numerical_checks'],
                               'SN_density_error':value['SN_loglike_absolute_error'],'SN_score_error':value['SN_chi2_absolute_error']})
            # Run the real construction/share function with only both physical
            # entry points replaced by explicit analytic objects. Exactly one
            # fake background must feed expansion and SN distance prediction.
            point=point_at(.3,mode);bg=Background(parameters(point));native=native_for(bg,point,sn,kernel)
            def fake_params(**kwargs):
                assert kwargs['lmax']==41 and kwargs['halofit_version']=='mead2016'
                assert all(kwargs[k]==2 for k in ['AccuracyBoost','lAccuracyBoost','lSampleBoost'])
                return parameters(point)
            def fake_background(pars):
                assert pars.WantTransfer;shared_calls.append(mode);return bg
            with patch.object(camb,'set_params',fake_params),patch.object(camb,'get_background',fake_background):
                counter={'background_invocations':0};value=h.calculate(point,native,kernel,design,old,counter)
            assert counter['background_invocations']==1 and not value['failed_numerical_gates']
            assert value['thermal_adapter']=={'input_WantTransfer':False,'background_WantTransfer':True}
            assert value['background_extra_args']==dict(EXTRA,lmax=41)
        # Independently verify conditional-noise total variance and B-q covariance.
        kernel=h.LuminosityKernel(sn,'smooth01',design['luminosity_redshifts'])
        means=np.array([kernel.from_prediction(observed+d,{})['mean_mag'] for d in [np.sin(z),-.3*np.cos(z),.2*z]])
        w=np.array([.2,.3,.5]);q=np.array([-.6,-.3,-.5]);average=w@means
        covariance=kernel.covariance+sum(wi*np.outer(mu-average,mu-average) for wi,mu in zip(w,means))
        cross=sum(wi*(mu-average)*(qi-w@q) for wi,mu,qi in zip(w,means,q))
        assert np.max(abs(cross-(means-average).T@(w*(q-w@q))))<1e-15
        # Deterministic paired +/- coefficient deviations reproduce the complete
        # conditional covariance exactly; they add zero covariance with fixed q.
        L=np.linalg.cholesky(kernel.coefficient_covariance);atoms=[];qa=[];wa=[]
        for mu,qi,wi in zip(means,q,w):
            for j in range(4):
                for sign in [-1,1]:atoms.append(mu+sign*2*kernel.basis@L[:,j]);qa.append(qi);wa.append(wi/8)
        atoms=np.array(atoms);wa=np.array(wa);qa=np.array(qa);delta=atoms-wa@atoms
        total_error=float(np.max(abs((delta*wa[:,None]).T@delta-covariance)))
        cross_error=float(np.max(abs((delta*wa[:,None]).T@(qa-wa@qa)-cross)))
        assert total_error<1e-13 and cross_error<1e-13
        # A complete2000-slot cache path uses analytic templates, not2000 physical
        # backgrounds. Fresh typed parent boundary remains explicitly mocked.
        mode='smooth01';kernel=h.LuminosityKernel(sn,mode,design['luminosity_redshifts']);kernel.declared_extra=EXTRA
        xs=np.linspace(-1.5,1.5,20)
        for j,x in enumerate(xs):
            point=point_at(float(x),mode);bg=Background(parameters(point));native=native_for(bg,point,sn,kernel)
            value=h.consume_background(bg,point,native,kernel,design,old)
            assert not value['failed_numerical_gates'];templates[h.rt.identity(point)]=value;natives[j]=native
        source=root/'source';source.mkdir();data_path=source/'sn.npz'
        np.savez(data_path,zHD=z,zHEL=zhel,MU=observed,covariance=C)
        sequence=np.tile(np.arange(20),25);rng.shuffle(sequence);sequence=np.tile(sequence,4)
        sequence[167]=sequence[166];sequence[1229]=sequence[1228]
        records=[];inputs={h.rt.relative(data_path):h.rt.digest(data_path)}
        locations=[{'chain':f'chain.{i//500+1}.txt','row':i%500,'expanded_index':3*(i%500)} for i in range(2000)]
        locations[167]=copy.deepcopy(locations[166]);locations[1229]=copy.deepcopy(locations[1228])
        for i,j in enumerate(sequence):
            native=copy.deepcopy(natives[int(j)]);point=native['point'];q=-.5*((point['H0']-68)/.2)**2
            native.update(index=i,status='finite_native_accuracy2',numerical_target_identity='synthetic_native2',
                          proposal_logpost=q,log_weight=0.,record_path=h.rt.relative(source/f'{i:05d}.json'))
            path=h.ROOT/native['record_path'];h.rt.write_new(path,native);native['record_sha256']=h.rt.digest(path)
            inputs[native['record_path']]=native['record_sha256'];records.append(native)
        config={'theory':{'camb':{'extra_args':EXTRA}},'likelihood':{'released_sn':{'external':h.ReleasedDistances,'data_file':str(data_path),'smooth_sigma':.1}},'params':{'epsilon':0.}}
        pure=[{'status':'finite','point':r['point'],'derived':r['native2_derived'],'log_weight':0.} for r in records]
        groups=np.repeat(np.arange(4),500);posterior=summarize(pure,groups)
        assert posterior['status']=='passed_importance_weight_gates'
        target=QualifiedNativeAccuracy2(2,'synthetic_native2',config,{'evolution':mode,'model':'lcdm'},'synthetic_proposal',tuple(records),
           tuple(r['point'] for r in records),groups,tuple(locations),np.array([r['proposal_logpost'] for r in records]),
           np.zeros(2000),np.full(2000,1/2000),{'status':'qualified_native_accuracy2_importance_posterior','posterior':posterior['posterior'],
           'importance_diagnostics':posterior},inputs,{'evidence':{'versions':{}}},source/'plan.json',source/'summary.json')
        # Exact exp-normalized weights, including roundoff convention.
        target=replace(target,normalized_weights=np.exp(np.zeros(2000)-np.log(2000)))
        validation=source/'mock-validation.json';validation.write_text('{"explicit_synthetic_boundary":true}')
        guard.enter_context(patch.object(h,'qualify_run',lambda work,summary:copy.deepcopy(target)))
        guard.enter_context(patch.object(h,'validation_guard',lambda:{h.rt.relative(validation):h.rt.digest(validation)}))
        cache=root/'cache';h.prepare(source,source/'summary.json',cache)
        counter={'logical_backgrounds':0}
        def synthetic_calculate(point,native,kernel,design,expansion_design,count):
            counter['logical_backgrounds']+=1;count['background_invocations']+=1
            return copy.deepcopy(templates[h.rt.identity(point)])
        with patch.object(h,'calculate',synthetic_calculate):result=h.execute(cache,root/'result.json')
        assert result['status']=='qualified_native_accuracy2_histories' and counter['logical_backgrounds']==2000
        replay=h.report(cache,root/'result.json');assert replay==h.rt.plain(h.rt.read(root/'result.json'))
        means=np.array([templates[h.rt.identity(r['point'])]['luminosity']['mean_mag'] for r in records]);weights=target.normalized_weights
        expected=kernel.covariance+((means-weights@means)*weights[:,None]).T@(means-weights@means)
        covariance_error=float(np.max(abs(expected-np.array(result['luminosity']['covariance_mag2']))))
        assert covariance_error<1e-14
        negative=[]
        for label,bad in [('wrong_type',{}),('wrong_accuracy',replace(target,numerical_accuracy=1)),
                          ('stale_weights',replace(target,normalized_weights=np.zeros(2000)))]:
            with patch.object(h,'qualify_run',lambda w,s:bad):assert refused(lambda:h.prepare(source,source/'summary.json',root/label))
            assert not (root/label).exists();negative.append(label)
        bad=copy.deepcopy(target);bad.native_configuration['likelihood']['released_sn']['external']=object
        with patch.object(h,'qualify_run',lambda w,s:bad):assert refused(lambda:h.prepare(source,source/'summary.json',root/'anchored'))
        negative.append('alternate_SN_class')
        bad=copy.deepcopy(target);bad.records[0]['value']['metadata']['finalized_theory_extra_args']['AccuracyBoost']=1
        with patch.object(h,'qualify_run',lambda w,s:bad):assert refused(lambda:h.prepare(source,source/'summary.json',root/'accuracy'))
        negative.append('wrong_finalized_accuracy')
        claim=cache/'00000.claim.json';original=claim.read_bytes();claim.write_text('{}')
        assert refused(lambda:h.report(cache));claim.write_bytes(original);negative.append('changed_claim_bytes')
        path=cache/'00000.json';original=path.read_bytes();row=h.rt.read(path);row['history']['q'][0]+=.1
        row['payload_sha256']=h.rt.identity(row);path.write_text(json.dumps(row))
        assert refused(lambda:h.report(cache));path.write_bytes(original);negative.append('resealed_changed_history')
        assert refused(lambda:h.execute(cache,root/'second.json'));negative.append('retry_existing_cache')
        # Native SN/derived mismatch cannot pass numerical gates even though
        # the background and every cached value are finite.
        point=point_at(0);bg=Background(parameters(point));native=native_for(bg,point,sn,kernel)
        native['native2_loglikes']['released_sn']+=.01
        mismatch=h.consume_background(bg,point,native,kernel,design,old)
        assert 'SN_loglike_absolute_error' in mismatch['failed_numerical_gates'];negative.append('native2_SN_mismatch')
        native['native2_derived']['q0']+=.01
        mismatch=h.consume_background(bg,point,native,kernel,design,old)
        assert 'q_native_scaled_error' in mismatch['failed_numerical_gates'];negative.append('native2_derivative_mismatch')
        baseline=h.LuminosityKernel(sn,'none',design['luminosity_redshifts'])
        assert baseline.from_prediction(np.zeros(sn.n),{})['mean_mag']==[0.]*8
        assert np.count_nonzero(baseline.covariance)==0
    assert physical==[]
    return {'status':'passed_synthetic_native_accuracy2_history_validation','physical_calls':0,'source_sha256':h.sources(),
       'independent_GLS':gls,'analytic_background_cases':checks,'shared_background_calls_by_mode':shared_calls,
       'total_covariance_atom_error':total_error,'cross_q0_atom_error':cross_error,
       'Gaussian_mixture_kernel_validation':mixtures,'typed_cache_slots':2000,'logical_synthetic_backgrounds':2000,
       'full_cache_covariance_error':covariance_error,'refusals':negative,'synthetic_work':h.rt.relative(root),
       'limits':'Analytic backgrounds and explicitly mocked fresh typed parent only. Real CAMB/model calls were forbidden; no observational history or numerical convergence claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=h.VALIDATION);a=p.parse_args();result=main_run()
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'physical_calls':0,'refusals':len(result['refusals'])}))
