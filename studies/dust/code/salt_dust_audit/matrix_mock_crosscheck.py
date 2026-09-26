#!/usr/bin/env python3
"""Independent, bounded cross-review of released low-RV truth and red flux.

Use phase2/env-official/bin/python. Does not import the scripts under review.
"""
from pathlib import Path
import gzip
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import extinction

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/salt_dust_audit/matrix_mock_crossreview'


def old_8000(av, rv):
    # Independent Python transcription of historical GALextinct optical branch
    # from phase2/official/build/SNANA-2fe0f56/src/MWgaldust.c, lines 365-461.
    y = 1e4/8000 - 1.82
    a = np.polynomial.polynomial.polyval(y,[1,.104,-.609,.701,1.137,-1.718,-.827,1.647,-.505])
    b = np.polynomial.polynomial.polyval(y,[0,1.952,2.908,-3.989,-7.985,11.102,5.491,-10.805,3.347])
    ratio = np.polynomial.polynomial.polyval(8.,[
        .0855929205,1.91547833,-1.65101945,.750611119,-.200041118,.0330155576,
        -.00346344458,.000230741420,-.00000943018242,.000000214917977,-.00000000208276810])
    return av*(a+b/rv)*ratio


def exact_8000(av, rv):
    return np.array([extinction.fitzpatrick99(np.array([8000.]),float(a),float(r))[0] for a,r in zip(av,rv)])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    baseline = ROOT/'sources/repos/des-science__DES-SN5YR@1.3'
    base = baseline/'1_SIMULATIONS/SNIa_SIMULATIONS'
    reference = pd.read_csv(ROOT/'runs/salt_dust_audit/snana_mock_support/per_realization.csv').set_index('realization')
    paths = sorted(base.glob('*/*_HEAD.FITS.gz'))
    assert len(paths)==25
    results,examples,inputs = [],[],[]
    for p in [paths[0],paths[12],paths[24]]:
        inputs.append(p)
        with fits.open(p,memmap=False) as h:
            d = h[1].data
            rv = np.array(d['SIM_RV'],dtype=float)
            av = np.array(d['SIM_AV'],dtype=float)
            snid = np.char.strip(np.asarray(d['SNID']).astype(str))
        assert np.all(rv>0) and np.all(av>=0)
        old,exact = old_8000(av,rv),exact_8000(av,rv)
        counts = {'N_written_simulated':len(rv),'N_RV_lt1':int(sum(rv<1)),
                  'N_RV_lt2':int(sum(rv<2)), 'N_A8000_old_lt0':int(sum(old<0)),
                  'N_A8000_exact_lt0':int(sum(exact<0))}
        for name,n in counts.items():
            assert reference.loc[p.parent.name,name]==n,(p,name)
        results.append({'realization':p.parent.name,**counts})
        # One stored truth example per inspected realization, chosen by minimum RV.
        i = int(np.argmin(rv))
        examples.append({'realization':p.parent.name,'SNID':str(snid[i]),'SIM_RV':float(rv[i]),
                         'SIM_AV':float(av[i]),'A8000_historical_mag':float(old[i]),
                         'A8000_exact_mag':float(exact[i]),'historical_component_flux_multiplier':float(10**(-.4*old[i]))})
        assert old[i]<0 and exact[i]<0

    # Independent photon-weighted broadband integration for x1=c=0, phase zero.
    # Read the released M0 surface and filter directly; neither sncosmo nor the
    # source agent's response/model constructor is used in this calculation.
    template = baseline/'2_LCFIT_MODEL/SALT3.DES5YR/salt3_template_0.dat.gz'
    band = ROOT/'phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits.gz'
    inputs += [template,band]
    t = np.loadtxt(template)
    at_zero = t[t[:,0]==0]
    assert len(at_zero)>10
    sed = CubicSpline(at_zero[:,1],at_zero[:,2])
    with fits.open(band) as h:
        f = h['FilterTrans'].data
        wave_band = np.array(f['wavelength (A)'],dtype=float)
        transmission = np.array(f['DES-i'],dtype=float)
    nonzero = np.flatnonzero(transmission>0)
    wave = np.arange(wave_band[nonzero[0]-1],wave_band[nonzero[-1]+1]+.25,.25)
    trans = np.interp(wave,wave_band,transmission)
    rest = wave/1.1
    assert rest.min()>=at_zero[:,1].min() and rest.max()<=at_zero[:,1].max()
    integrand = wave*trans*sed(rest)
    attenuation = extinction.fitzpatrick99(rest,.1*.4,.4)
    baseflux = np.trapezoid(integrand,wave)
    dustflux = np.trapezoid(integrand*10**(-.4*attenuation),wave)
    delta = -2.5*np.log10(dustflux/baseflux)
    target = json.loads((ROOT/'runs/salt_dust_audit/f99_support/summary.json').read_text())['example']['extinction_mag']
    assert delta<-.02 and abs(delta-target)<1e-4
    example = {'z':.1,'rest_phase_day':0,'band':'DES-i','RV':.4,'EBV':.1,
               'independent_extinction_mag':float(delta),'source_agent_example_mag':float(target),
               'difference_mag':float(delta-target),'flux_ratio':float(dustflux/baseflux),
               'method':'Direct 0.25-A photon-weighted integration of released phase-zero M0 and filter; cubic SED interpolation, linear filter interpolation.'}
    inputs += [Path(__file__), Path(extinction.__file__),
               ROOT/'phase2/official/build/SNANA-2fe0f56/src/MWgaldust.c',
               ROOT/'phase2/official/build/SNANA-2fe0f56/src/genmag_SEDtools.c',
               ROOT/'phase2/official/build/SNANA-2fe0f56/src/genmag_SALT2.c',
               ROOT/'runs/salt_dust_audit/snana_mock_support/per_realization.csv',
               ROOT/'runs/salt_dust_audit/f99_support/summary.json']
    report = {'scope':'Independent bounded spot-check, not full 25-realization recount or end-to-end BBC/cosmology rerun.',
              'inspected_realizations':results,'stored_truth_examples':examples,'broadband_example':example,
              'RV_roots_at_8000':{'historical':float(brentq(lambda r:old_8000(1.,r),.3,1)),
                                  'exact':float(brentq(lambda r:exact_8000(np.array([1.]),np.array([r]))[0],.3,1))},
              'interpretation':'Negative A8000 means this extrapolated host-extinction component brightens a monochromatic SED. It is not a real observed negative extinction or a measured cosmological bias. The broadband demonstration is constructed, separate from the actual stored mock examples.',
              'inputs':[{'path':str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),
                         'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]}
    (OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))


if __name__=='__main__':
    main()
