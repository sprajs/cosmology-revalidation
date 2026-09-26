#!/usr/bin/env python3
"""Physical invariants and regression checks for the observation chain.

Run with phase2/env-official/bin/python from the repository root.
No existing scientific results are overwritten.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import sys
from datetime import datetime, timezone
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts/phase2/independent_flux'))
from engine import Engine, read_inputs
OUT=ROOT/'runs/physics_audit'
OUT.mkdir(exist_ok=True, parents=True)

# A deliberately non-blackbody analytic SED; physical conservation must not
# depend on a thermal spectrum or a SALT fit. L_lambda units are erg/s/Angstrom.
hc=1.9864458571489286e-8 # erg Angstrom, exact SI constants converted
c=2.99792458e18 # Angstrom/s
lum=lambda w: 1e39*np.exp(-.5*((w-5500.)/1400.)**2)
emit_lo,emit_hi=500.,16000.
L=quad(lum,emit_lo,emit_hi,epsabs=0,epsrel=1e-12)[0]
records=[]
for z in [0.,.02,.4,1.1,2.]:
    dl=1e27
    f=lambda w:lum(w/(1+z))/(4*np.pi*dl**2*(1+z))
    lo,hi=emit_lo*(1+z),emit_hi*(1+z)
    bol=quad(f,lo,hi,epsabs=0,epsrel=1e-12)[0]
    expected=L/(4*np.pi*dl**2)
    # Photon integrals must be identical in frequency and wavelength units.
    trans=lambda w:np.exp(-.5*((w-7200.)/1000.)**2)
    nlam=quad(lambda w:w/hc*trans(w)*f(w),lo,hi,epsabs=0,epsrel=1e-11)[0]
    # dnu integration transformed to a dimensionless coordinate for quadrature.
    nu_lo,nu_hi=c/hi,c/lo
    nnu=quad(lambda x: (f(c/(x*1e14))*(c/(x*1e14))**2/c)
             /(hc/c*x*1e14)*trans(c/(x*1e14))*1e14,
             nu_lo/1e14,nu_hi/1e14,epsabs=0,epsrel=1e-11)[0]
    records.append({'z':z,'bolometric_relative_error':bol/expected-1,
                    'photon_frequency_wavelength_relative_error':nnu/nlam-1})
assert max(abs(r['bolometric_relative_error']) for r in records)<1e-11
assert max(abs(r['photon_frequency_wavelength_relative_error']) for r in records)<1e-10

# O08: identical AB bands with flat L_nu. Common constants cancel.
weight=lambda w: w*np.exp(-.5*((w-5500)/900)**2)
ref=quad(lambda w:weight(w)/w**2,1000,12000,epsabs=1e-14)[0]
kchecks=[]
for z in [.01,.3,1.2,2.]:
    shifted=quad(lambda w:weight(w)/(w/(1+z))**2,1000,12000,epsabs=1e-14)[0]
    measured=-2.5*np.log10(shifted/(1+z)/ref)
    exact=-2.5*np.log10(1+z)
    kchecks.append({'z':z,'K_mag':measured,'difference_mag':measured-exact})
assert max(abs(r['difference_mag']) for r in kchecks)<1e-12

engine=Engine();_,_,pars=read_inputs()
before=json.loads((OUT/'flux-before.json').read_text())
reference=[]
for cid in sorted(pars):
    q=np.load(ROOT/f'phase2/official/portable_pilot/objective_{cid}.npz')
    p=q['parameters_x0_x1_c_t0'].copy();p[0]=np.log(p[0])
    z=float(q['zHEL'][0]);g=engine.prepare(z,float(q['MWEBV'][0]))
    actual=engine.flux(p,q['band'],q['MJD'],z,g,interpolation='sncosmo')
    old=np.array(before[cid]['flux'])
    np.testing.assert_array_equal(actual,old)
    pdim=p.copy();pdim[0]-=.4*np.log(10)*.2
    dim=engine.flux(pdim,q['band'],q['MJD'],z,g,interpolation='sncosmo')
    np.testing.assert_allclose(dim,actual*10**(-.4*.2),rtol=3e-15,atol=0)
    pdistance=p.copy();pdistance[0]-=np.log(4.)
    farther=engine.flux(pdistance,q['band'],q['MJD'],z,g,interpolation='sncosmo')
    np.testing.assert_allclose(farther,actual/4,rtol=3e-15,atol=0)
    reference.append({'CID':cid,'epochs':len(actual),'unchanged_exactly':True,
                      'max_abs_delta_flux':float(np.max(abs(actual-old)))})
rejected=[]
for z in [1.,2.]:
    grid=engine.prepare(z,0.)
    assert not grid['g']['support_valid']
    try:engine.flux(np.array([0.,0.,0.,0.]),np.array(['g']),np.array([0.]),z,grid)
    except ValueError as e:rejected.append({'z':z,'band':'g','message':str(e)})
    else:raise AssertionError('Incomplete band was silently evaluated')
    # A supported filter remains usable although an unused filter is invalid.
    assert np.isfinite(engine.flux(np.array([0.,0.,0.,0.]),np.array(['z']),np.array([0.]),z,grid)).all()

inputs=[ROOT/'scripts/phase2/independent_flux/engine.py',Path(__file__),OUT/'flux-before.json']
inputs.extend(sorted((ROOT/'phase2/official/portable_pilot').glob('objective_*.npz')))
inputs.extend(p for p in (ROOT/'phase2/official/portable_pilot/assets').iterdir() if p.is_file())
inputs.extend(sorted((ROOT/'sources/repos/des-science__DES-SN5YR@1.3/2_LCFIT_MODEL/SALT3.DES5YR').glob('salt3_lc_*.dat.gz')))
result={'utc':datetime.now(timezone.utc).isoformat(),'status':'passed',
 'spectral_conservation':records,'flat_fnu_K_correction':kchecks,
 'reference_fluxes':reference,'total_reference_epochs':sum(r['epochs'] for r in reference),
 'gray_dimming_and_inverse_square_amplitude_checks':'passed',
 'unsupported_bands_rejected':rejected,'unused_unsupported_band_does_not_block_supported_band':True,
 'scope':'Mathematical conservation checks and existing in-domain flux regression; no survey calibration or cosmology refit.',
 'packages':{p:importlib.metadata.version(p) for p in ['numpy','scipy','astropy','sncosmo','extinction']},
 'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
(OUT/'observation-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'reference_objects':len(reference),'reference_epochs':result['total_reference_epochs'],
 'spectral_conservation_max_error':max(abs(r['bolometric_relative_error']) for r in records),
 'photon_integral_max_error':max(abs(r['photon_frequency_wavelength_relative_error']) for r in records),
 'invalid_bands_rejected':len(rejected)},indent=2))
