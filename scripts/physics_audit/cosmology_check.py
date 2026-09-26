"""Independent physical identities, units and domains for local cosmology code.

Run: .venv/bin/python scripts/physics_audit/cosmology_check.py
Uses the separate hierarchy interpreter only for its optional JAX checks.
No data fitting, external likelihood execution, or historical output mutation.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import sys

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM, Flatw0waCDM
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/cosmology"))
import core
sys.path.insert(0, str(ROOT / "scripts/mapping"))
import csfh_dtd


def checks():
    result = {}
    z = np.array([.001, .01, .1, .3, .6, 1., 2., 2.5])
    cases = [[.3, -1., 0.], [.353, -.42, -1.75], [.24, -.7, .3]]
    error, qerror, conservation = [], [], []
    for theta in cases:
        ast = Flatw0waCDM(H0=70., Om0=theta[0], w0=theta[1], wa=theta[2], Tcmb0=0)
        independent = ast.luminosity_distance(z).value / (299792.458 / 70.)
        actual = 10**(core.mu(z, theta)[0] / 5)
        error.append(float(np.max(np.abs(5*np.log10(actual / independent)))))
        h = 1e-5
        q_numeric = -1 + (1+z)*(np.log(core.efunc(z+h, theta)[0])-np.log(core.efunc(z-h, theta)[0]))/(2*h)
        qerror.append(float(np.max(np.abs(q_numeric-core.qvalue(z, theta)[0]))))
        om, w0, wa = theta
        de = lambda u: (1+u)**(3*(1+w0+wa))*np.exp(-3*wa*u/(1+u))
        lhs = (np.log(de(z+h))-np.log(de(z-h)))/(2*h)
        conservation.append(float(np.max(np.abs(lhs-3*(1+w0+wa*z/(1+z))/(1+z)))))
    result['cpl_distance_max_abs_mag'] = max(error)
    result['cpl_q_derivative_max_abs'] = max(qerror)
    result['cpl_conservation_max_abs'] = max(conservation)
    assert max(error)<1e-10 and max(qerror)<2e-8 and max(conservation)<2e-8

    bin_errors = []
    for q in [-1., -.5, 0., 1e-12, .5, 2.]:
        expected = np.log1p(z) if abs(q)<1e-10 else -np.expm1(-q*np.log1p(z))/q
        bin_errors.append(float(np.max(np.abs(core.integral(z, np.full(5,q), 'qbins')[0]-expected))))
    result['qbins_constant_q_integral_max_abs'] = max(bin_errors)
    assert max(bin_errors)<1e-10
    q = np.array([-.6, -.2, .3, .7, .1])
    log_e = sum((1+v)*np.clip(np.log((1+z)/(1+lo)),0,np.log((1+hi)/(1+lo))) for v,lo,hi in zip(q,core.EDGES[:-1],core.EDGES[1:]))
    result['qbins_efunc_max_abs'] = float(np.max(np.abs(core.efunc(z,q,'qbins')[0]-np.exp(log_e))))
    assert result['qbins_efunc_max_abs']<1e-12
    domain = {}
    for name, fn in [('integral',core.integral),('efunc',core.efunc),('qvalue',core.qvalue)]:
        for bad in [-.01,2.50001,np.nan]:
            try:
                fn([bad],q,'qbins')
            except ValueError:
                domain[f'{name}:{bad}'] = 'ValueError'
            else:
                raise AssertionError((name,bad,'accepted invalid q-bin redshift'))
    result['qbins_outside_domain'] = domain

    ba = core.BAO()
    theta = [.353,-.42,-1.75]; hrd = 10000.
    ast = Flatw0waCDM(H0=70.,Om0=theta[0],w0=theta[1],wa=theta[2],Tcmb0=0)
    rd = hrd/70.
    dm = ast.comoving_transverse_distance(ba.z).value/rd
    dh = 299792.458/ast.H(ba.z).value/rd
    dv = (ba.z*dm*dm*dh)**(1/3)
    expected = np.where(ba.data.kind=='DM_over_rs',dm,np.where(ba.data.kind=='DH_over_rs',dh,dv))
    result['bao_prediction_max_abs_ratio'] = float(np.max(np.abs(ba.prediction(theta,[hrd])[0]-expected)))
    assert result['bao_prediction_max_abs_ratio']<1e-10

    ages = {}
    for name,c in csfh_dtd.COSMOS.items():
        clock,_ = csfh_dtd.clock(c)
        za = np.array([0.,.1,.5,1.,2.,10.])
        expected = Flatw0waCDM(**c,Tcmb0=0).age(za).value
        ages[name] = float(np.max(np.abs(clock(-np.log1p(za))-expected)))
    result['clock_vs_independent_astropy_max_abs_gyr'] = ages
    assert max(ages.values())<2e-6
    r0 = FlatLambdaCDM(H0=70.,Om0=.3,Tcmb0=0)
    rr = FlatLambdaCDM(H0=70.,Om0=.3,Tcmb0=2.7255,m_nu=0)
    result['radiation_approximation_example'] = {
        'definition':'Flat LCDM H0=70 Om0=.3; Tcmb=2.7255 K massless neutrinos versus Tcmb=0; not a fitted cosmology',
        'z':z.tolist(),
        'radiation_minus_no_radiation_mu_mag':(rr.distmod(z).value-r0.distmod(z).value).tolist(),
        'radiation_minus_no_radiation_age_gyr':(rr.age(z).value-r0.age(z).value).tolist(),
    }
    direct_b13 = .180/(10**(-.997*(z-1.243))+10**(.241*(z-1.243)))
    direct_md14 = .015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
    result['csfh_formula_max_abs'] = float(max(np.max(abs(csfh_dtd.csfh(z,'B13')-direct_b13)),np.max(abs(csfh_dtd.csfh(z,'MD14')-direct_md14))))
    assert result['csfh_formula_max_abs']<1e-14

    # Adaptive quadrature uses the same declared cosmic clock but independent
    # integration and direct SFH/DTD expressions, avoiding the delay-grid rule.
    c = csfh_dtd.COSMOS['son_cpl_H63p6']
    clock,z_of_age = csfh_dtd.clock(c)
    grid_result = csfh_dtd.calculate(c,'B13','C14_smooth')
    comparisons = []
    for redshift in [0.,1.]:
        age = float(clock(-np.log1p(redshift)))
        def weight(delay):
            formed_z = float(z_of_age(age-delay))
            if not np.isfinite(formed_z) or formed_z>1000:
                return 0.
            sf = .180/(10**(-.997*(formed_z-1.243))+10**(.241*(formed_z-1.243)))
            u = delay/.3
            return sf*u**20/(1+u**21)
        area = quad(weight,0.,age,points=[.3],epsabs=1e-10)[0]
        mean = quad(lambda delay:delay*weight(delay),0.,age,points=[.3],epsabs=1e-10)[0]/area
        cdf = lambda upper:quad(weight,0.,upper,points=[.3] if upper>.3 else None,epsabs=1e-10)[0]/area
        median = brentq(lambda upper:cdf(upper)-.5,1e-6,age-1e-6,xtol=1e-10)
        row = grid_result.loc[np.isclose(grid_result.z,redshift)].iloc[0]
        comparisons.append({'z':redshift,'adaptive_mean_gyr':mean,'adaptive_median_gyr':median,
                            'grid_minus_adaptive_mean_gyr':float(row.mean_delay_gyr-mean),
                            'grid_minus_adaptive_median_gyr':float(row.median_delay_gyr-median)})
    result['csfh_dtd_adaptive_quadrature']={'clock':'shared declared radiation-free CPL63.6 clock; independent delay quadrature', 'cases':comparisons}
    assert max(abs(v['grid_minus_adaptive_mean_gyr']) for v in comparisons)<1e-7
    assert max(abs(v['grid_minus_adaptive_median_gyr']) for v in comparisons)<2e-5

    # Uniform star formation with a tau^-1 DTD has exact mean and median.
    tmin,tmax=.04,13.
    norm=quad(lambda delay:1/delay,tmin,tmax,epsabs=1e-12)[0]
    mean=quad(lambda delay:1.,tmin,tmax)[0]/norm
    result['normalized_delay_control']={'mean_gyr':mean,'analytic_mean_gyr':(tmax-tmin)/np.log(tmax/tmin),'median_gyr':float(np.sqrt(tmin*tmax)),'cdf_at_median':float(np.log(np.sqrt(tmin*tmax)/tmin)/norm)}
    assert abs(result['normalized_delay_control']['cdf_at_median']-.5)<1e-12

    tab=pd.read_csv(core.ROOT/'sources/repos/CobayaSampler__sn_data/PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+')
    result['current_data_domains']={'pantheon_zHD_max':float(tab.loc[tab.zHD>.01,'zHD'].max())}
    sampled={}
    for path in sorted((ROOT/'phase2/hierarchy').glob('**/data*.npz')):
        with np.load(path,allow_pickle=False) as d:
            if 'z' in d: sampled[str(path.relative_to(ROOT))]={'n':len(d['z']),'min':float(d['z'].min()),'max':float(d['z'].max())}
    result['current_data_domains']['hierarchy_files']=sampled
    assert result['current_data_domains']['pantheon_zHD_max']<=2.5
    assert all(v['min']>0 and v['max']<=1.3 for v in sampled.values())
    return result


def jax_checks():
    exe=ROOT/'phase2/hierarchy/.venv/bin/python'
    program=r'''
import importlib.util, json, numpy as np, jax, jax.numpy as jnp
s=importlib.util.spec_from_file_location('hc','scripts/phase2/hierarchy/core.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
z=jnp.array([.01,.15,.35,.6,1.,1.3]); errors=[]
for q in [-1.,-.5,0.,1e-12,.5,2.]:
 chi=np.log1p(z) if abs(q)<1e-10 else -np.expm1(-q*np.log1p(z))/q
 expected=25+5*np.log10((1+z)*chi*m.C_KMS/70.)
 errors.append(float(np.max(np.abs(np.asarray(m.distance_qbins(z,z,jnp.full(4,q)))-expected))))
assert max(errors)<1e-10
bad=jnp.array([-.01,1.30001,np.nan]); out=jax.jit(m.distance_qbins)(bad,bad,jnp.zeros(4))
assert np.all(np.isnan(np.asarray(out)))
jac=np.asarray(jax.jacrev(lambda q:m.distance_qbins(z,z,q))(jnp.zeros(4)))
assert np.isfinite(jac).all()
print(json.dumps({'constant_q_max_abs_mag':max(errors),'jit_invalid_redshift_nan':bool(np.isnan(out).all()),'q0_gradients_finite':bool(np.isfinite(jac).all())}))
'''
    p=subprocess.run([str(exe),'-c',program],cwd=ROOT,text=True,capture_output=True,check=True)
    return json.loads(p.stdout)


def main():
    result={'interpretation':'Numerical implementation of stated physical equations and domains; not empirical validation of cosmology, SFH, DTD or a full external likelihood.',**checks(),'jax_qbins':jax_checks()}
    paths=[Path(__file__),ROOT/'scripts/cosmology/core.py',ROOT/'scripts/mapping/csfh_dtd.py',ROOT/'scripts/phase2/hierarchy/core.py']
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result['python']=sys.version.split()[0]
    out=ROOT/'runs/physics_audit/cosmology-check.json';out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['current_data_domains','radiation_approximation_example']},indent=2))


if __name__=='__main__': main()
