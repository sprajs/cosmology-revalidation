"""Contemporary CAMB background sensitivity of the released calibrated SN model."""
from __future__ import annotations
import argparse
import ast
import hashlib
import inspect
import json
from pathlib import Path
import time
import urllib.request

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg
from scipy.integrate import quad
import camb

from calibration_interface import ReleasedCalibration, ROOT

HERE=Path(__file__).resolve().parent
WORK=ROOT/'.work/unified-cosmology/calibration-background-review'
SOURCES=HERE/'calibration-background-sources.json'
DESIGN=HERE/'calibration-background-design.json'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def relative(path):return str(Path(path).resolve().relative_to(ROOT))
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def acquire(allow_download=False):
    manifest=json.loads(SOURCES.read_text());records=[]
    for row in manifest['files']:
        path=ROOT/row['path']
        if not path.exists():
            if not allow_download:raise FileNotFoundError(str(path)+'; use --acquire')
            path.parent.mkdir(parents=True,exist_ok=True)
            with urllib.request.urlopen(row['url'],timeout=120) as response:body=response.read()
            assert hashlib.sha256(body).hexdigest()==row['sha256'];path.write_bytes(body)
        assert sha(path)==row['sha256'];records.append(row)
    return records


def source_semantics():
    # Inspect literal source expressions without executing third-party pipeline code.
    relations={}
    for epoch in ['pre_chain','current']:
        path=WORK/(epoch+'-consistency.py');module=ast.parse(path.read_text())
        expr=next(n.value for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='COSMOLOGY_CONSISTENCY_RELATIONS' for t in n.targets))
        items=ast.literal_eval(expr)
        assert ('omega_m','omega_b+omega_c+omega_nu') in items
        assert ('omch2','ommh2-ombh2-omnuh2') in items
        relations[epoch]={'Omega_m':'omega_b+omega_c+omega_nu', 'omch2':'ommh2-ombh2-omnuh2',
                          'neutrino_conversions':[list(x) for x in items if x[0] in ('mnu','omnuh2') and ('93.14' in str(x) or '94.064' in str(x))]}
        text=(WORK/(epoch+'-camb_interface.py')).read_text()
        assert 'Parameter massless_nu is being ignored' in text
        parsed=ast.parse(text)
        optional=[n for n in ast.walk(parsed) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='get_optional_params' and len(n.args)>=3 and isinstance(n.args[2],ast.List) and any(isinstance(v,ast.Constant) and v.value=='num_massive_neutrinos' for v in n.args[2].elts)]
        assert len(optional)==1
        names=[v.value for v in optional[0].args[2].elts if isinstance(v,ast.Constant)]
        assert 'massive_nu' not in names and 'num_massive_neutrinos' in names
        relations[epoch]['optional_CAMB_cosmology_keys']=names
        relations[epoch]['massless_nu_ignored']=True
    old=(WORK/'camb-1.3.5-model.py').read_text();constants=(WORK/'camb-1.3.5-constants.py').read_text()
    assert 'num_massive_neutrinos=1' in old and 'default_nnu = 3.046' in constants
    return {'verified_relations':relations,'historical_CAMB_example':'1.3.5 defaults one massive and Neff3.046',
            'attested_author_runtime':False,
            'reason':'Official chain header omits CSL/CAMB source identity. Contemporary source snapshots do not prove which version the authors installed.'}


def parameters(case,Omega_m,H0,accuracy=1):
    h=H0/100.;obh2=.048*h*h
    if case=='modern_one_species':
        p=camb.set_params(H0=H0,ombh2=obh2,omch2=Omega_m*h*h-obh2,
                          mnu=.06,num_massive_neutrinos=1,nnu=3.044,standard_neutrino_neff=3.044,
                          TCMB=2.7255,omk=0.,w=-1.,wa=0.,dark_energy_model='ppf')
        p.omch2=Omega_m*h*h-obh2-p.omnuh2
    elif case=='declared_three_species':
        p=camb.set_params(H0=H0,ombh2=obh2,omch2=Omega_m*h*h-obh2-.00083,
                          mnu=.00083*93.14,num_massive_neutrinos=3,nnu=3.046,standard_neutrino_neff=3.046,
                          TCMB=2.7255,omk=0.,w=-1.,wa=0.,dark_energy_model='ppf')
        # Explicit fine-grained public CAMB fields implement the literal header
        # interpretation; automatic Neff splitting need not give .046 massless.
        p.omnuh2=.00083;p.nu_mass_eigenstates=1;p.num_nu_massive=3;p.num_nu_massless=.046
        p.nu_mass_degeneracies=[3.];p.nu_mass_numbers=[3];p.nu_mass_fractions=[1.]
    elif case=='contemporary_interface_candidate':
        p=camb.set_params(H0=H0,ombh2=obh2,omch2=Omega_m*h*h-obh2-.00083,
                          mnu=.00083*93.14,num_massive_neutrinos=1,nnu=3.046,standard_neutrino_neff=3.046,
                          TCMB=2.7255,omk=0.,w=-1.,wa=0.,dark_energy_model='ppf')
    else:raise ValueError(case)
    p.Accuracy.AccuracyBoost=accuracy
    # CAMB's static validate requires a requested transfer/Cl output even though
    # get_background itself supports neither. Validate physical parameters first,
    # then disable spectral requests; validation performs no theory calculation.
    assert p.validate()
    p.WantCls=False;p.WantTransfer=False;p.WantDerivedParameters=False
    return p


def state(p):
    return {'H0':p.H0,'ombh2':p.ombh2,'omch2':p.omch2,'omnuh2':p.omnuh2,
            'actual_Omega_m_density_sum':(p.ombh2+p.omch2+p.omnuh2)/(p.H0/100)**2,
            'N_eff':p.N_eff,'num_nu_massive':p.num_nu_massive,'num_nu_massless':p.num_nu_massless,
            'mass_eigenstates':p.nu_mass_eigenstates,'nu_mass_degeneracies':list(p.nu_mass_degeneracies[:p.nu_mass_eigenstates]),
            'nu_mass_numbers':list(p.nu_mass_numbers[:p.nu_mass_eigenstates]),'nu_mass_fractions':list(p.nu_mass_fractions[:p.nu_mass_eigenstates]),
            'TCMB':p.TCMB,'YHe':p.YHe,'AccuracyBoost':p.Accuracy.AccuracyBoost}


def reference_distance(z,Omega_m,H0):
    x,w=leggauss(96);zz=z[:,None]*(x+1)/2
    return 299792.458/H0*z*np.sum(w/np.sqrt(Omega_m*(1+zz)**3+1-Omega_m),axis=1)/2/(1+z)


def local_curvature(rows,producer):
    s=np.array([.03,1.5]);design=[];values=[]
    for r in rows:
        x,y=(np.array([r['Omega_m_requested'],r['H0']])-[.33245,73.55])/s
        design.append([1.,x,y,x*x,x*y,y*y]);values.append(r['delta_loglike'])
    D=np.array(design);v=np.array(values);fit,res,rank,sv=linalg.lstsq(D,v,lapack_driver='gelsd');assert rank==6
    gradient=fit[1:3]/s;hessian=np.array([[2*fit[3],fit[4]],[fit[4],2*fit[5]]])/np.outer(s,s)
    p=producer['posterior'];C=np.array([[p['Omega_m']['sd']**2,p['H0_Omega_m_covariance']],
                                       [p['H0_Omega_m_covariance'],p['H0_km_s_Mpc']['sd']**2]])
    precision=linalg.solve(C,np.eye(2),assume_a='pos')-hessian
    if np.linalg.eigvalsh(precision).min()<=0:return {'status':'unsupported_nonpositive_local_precision'}
    moved=linalg.solve(precision,np.eye(2),assume_a='pos');shift=moved@gradient
    fraction=float(np.sqrt(moved[1,1]/C[1,1])-1)
    return {'status':'local_quadratic_diagnostic_not_recomputed_posterior','max_fit_residual_loglike':float(abs(D@fit-v).max()),
            'gradient_Omega_H0':gradient.tolist(),'hessian_Omega_H0':hessian.tolist(),
            'local_H0_sigma_relative_change':fraction,'local_mean_shift_Omega_H0':shift.tolist(),
            'one_percent_width_followup_flag':abs(fraction)>=.01}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--acquire',action='store_true');ap.add_argument('--quadrature',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-background-review.json')
    args=ap.parse_args();start=time.monotonic();source_records=acquire(args.acquire);semantics=source_semantics()
    producer=json.loads(args.quadrature.read_text());assert producer['status']=='qualified_conditional_SN_only_LCDM_quadrature'
    design=json.loads(DESIGN.read_text());release=ReleasedCalibration()
    package=Path(inspect.getfile(camb)).parent
    dependencies=[Path(__file__),SOURCES,DESIGN,args.quadrature,HERE/'calibration_interface.py',package/'model.py',package/'constants.py',package/'camblib.so']
    bindings={relative(p):sha(p) for p in dependencies}
    bindings.update({r['path']:r['sha256'] for r in source_records})
    bindings.update({p:r['sha256'] for p,r in release.input_records.items()})
    execution={'design':design,'source_bindings':bindings,'CAMB_version':camb.__version__}
    execution_key=hashlib.sha256(json.dumps(execution,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:20]
    execution_path=WORK/f'execution-design-{execution_key}.json'
    if execution_path.exists():assert json.loads(execution_path.read_text())==execution
    else:write(execution_path,execution)
    def forbidden(*args,**kwargs):raise RuntimeError('Spectral evaluation forbidden in background audit')
    camb.get_results=forbidden;camb.get_transfer_functions=forbidden
    records=[];validation=[];background_calls=0
    for case in design['cases']:
        for om in design['grid']['Omega_m']:
            for h0 in design['grid']['H0']:
                p=parameters(case,om,h0);background=camb.get_background(p,no_thermo=True);background_calls+=1
                da=background.angular_diameter_distance(release.z_hd_noncalibrator)
                ref=reference_distance(release.z_hd_noncalibrator,om,h0)
                got=release.evaluate(da);original=release.evaluate(ref)
                delta=release.theory_mu(da)-release.theory_mu(ref)
                white=linalg.solve_triangular(release.chol,delta,lower=True)
                white-=release.one_white*(release.one_white@white)/release.flat_m_precision
                record={'case':case,'Omega_m_requested':om,'H0':h0,'CAMB_state':state(p),
                        'relative_distance_min':float(np.min(da/ref-1)),'relative_distance_max':float(np.max(da/ref-1)),
                        'max_abs_distance_modulus_change_mag':float(abs(delta).max()),
                        'delta_chi2':got['chi2']-original['chi2'],'delta_loglike':got['loglike']-original['loglike'],
                        'projected_model_shift_Mahalanobis_norm':float(np.linalg.norm(white)),
                        'conditional_M_shift_mag':got['M_conditional_mean']-original['M_conditional_mean']}
                records.append(record)
                if om==.33245 and h0==73.55:
                    high=camb.get_background(parameters(case,om,h0,accuracy=2),no_thermo=True);background_calls+=1
                    thermal=camb.get_background(p.copy(),no_thermo=False);background_calls+=1
                    da_high=high.angular_diameter_distance(release.z_hd_noncalibrator)
                    da_thermal=thermal.angular_diameter_distance(release.z_hd_noncalibrator)
                    checkz=np.array([.01,.1,.3,.5,1.,1.5,2.,float(release.z_hd_noncalibrator.max())])
                    direct=np.array([299792.458/(1+z)*quad(lambda zz:1/float(background.hubble_parameter(zz)),0,z,epsabs=1e-13,epsrel=1e-11)[0] for z in checkz])
                    native=background.angular_diameter_distance(checkz)
                    validation.append({'case':case,'independent_H_integral_max_relative_error':float(abs(direct/native-1).max()),
                                       'AccuracyBoost2_distance_relative_max_change':float(abs(da_high/da-1).max()),
                                       'AccuracyBoost2_loglike_change':release.evaluate(da_high)['loglike']-got['loglike'],
                                       'thermal_vs_no_thermo_distance_relative_max_change':float(abs(da_thermal/da-1).max()),
                                       'thermal_vs_no_thermo_loglike_change':release.evaluate(da_thermal)['loglike']-got['loglike']})
    assert background_calls==33
    assert max(r['independent_H_integral_max_relative_error'] for r in validation)<1e-7
    curvature={case:local_curvature([r for r in records if r['case']==case],producer) for case in design['cases']}
    for p,digest in bindings.items():assert sha(ROOT/p)==digest
    result={'status':'completed_bounded_background_sensitivity_no_new_posterior','CAMB_version':camb.__version__,
            'source_semantics':semantics,'grid_results':records,'numerical_validation':validation,'local_uncertainty_diagnostics':curvature,
            'background_calls':background_calls,'spectrum_calls':0,'sampling_steps':0,'seconds':time.monotonic()-start,
            'execution_design_path':relative(execution_path),'execution_design_sha256':sha(execution_path),'input_source_sha256':bindings,
            'limits':['All CAMB calculations use the current pinned installation; neither the header interpretation nor contemporary-interface candidate is an exact historical reproduction.',
                      'The fixed baryon fraction .048 is held common to isolate late-time neutrino/radiation effects. It is not the full modern sampled cosmological model.',
                      'Nine-point local curvature diagnostics do not establish global posterior widths or tail coverage. No corrected uncertainty is substituted into the qualified matter+Lambda result.',
                      'If a full sensitivity is needed, interpolate only the small background-distance difference over a bounded Omega_m/H0 grid and validate independent held-out backgrounds before deterministic posterior quadrature.']}
    write(args.output,result);print(json.dumps({'status':result['status'],'local_uncertainty_diagnostics':curvature,'numerical_validation':validation,'seconds':result['seconds']},indent=2))


if __name__=='__main__':main()
