#!/usr/bin/env python3
"""Pinned SNooPy VegaB synthetic-photometry reference check; no fit/KCOR."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import integrate, interpolate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
P = json.loads((HERE/'comparison-protocol.json').read_text())
SN = HERE/'snpy/snpy/filters'
OLD = HERE.parent/'snpy/snpy/filters/filters/LCO'
BD = ROOT/'phase2/official/inputs/SNDATA_ROOT/standards/bd_17d4708_stisnic_007.dat'
RAISIN = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.input'
checks = {
 'snpy_source': SN/'__init__.py',
 'snpy_standard_registry': SN/'standards/Vega/standards.dat',
 'vegaB_spectrum': SN/'standards/Vega/alpha_lyr_stis_005.ascii',
 'raisin_bd17_spectrum': BD,
 'raisin_kcor_input': RAISIN,
}
for k, path in checks.items():
    assert hashlib.sha256(path.read_bytes()).hexdigest() == P['input_hashes'][k], path
assert 'VegaB  alpha_lyr_stis_005.ascii' in checks['snpy_standard_registry'].read_text()
ref = np.loadtxt(BD)
vega = np.loadtxt(checks['vegaB_spectrum'])
assert np.all(np.diff(ref[:,0]) > 0) and np.all(np.diff(vega[:,0]) > 0)
oldacq = json.loads((HERE.parent/'snpy-acquisition.json').read_text())
newacq = json.loads((HERE/'curve-acquisition.json').read_text())
source_hash = {a['path'].split('/')[-1]:a['sha256'] for a in oldacq+newacq if a['path'].endswith('.dat')}
source_band = {'Y':'Swope', 'y':'Dupont', 'J':'Swope', 'j':'Swope', 'H':'Swope', 'h':'Dupont', 'WIRC-J':'Dupont'}
mirror_name = {'J':'Jrc1_SWO_TAM_scan_atm.dat', 'j':'Jrc2_SWO_TAM_scan_atm.dat'}
result = {'protocol_sha256':hashlib.sha256((HERE/'comparison-protocol.json').read_bytes()).hexdigest(), 'bands':{}}

def snpy_response(band, spec):
    swave, flux = spec.T
    x, s = band.T
    assert x[0] >= swave[0] and x[-1] <= swave[-1]
    i_min = np.flatnonzero(swave-x[0] > 0)[0]
    i_max = np.flatnonzero(swave-x[-1] > 0)[0]
    i_min = max(i_min-5, 0)
    i_max = min(i_max+5, len(swave)-1)
    wave = swave[i_min:i_max+1]
    f = flux[i_min:i_max+1]
    filt = interpolate.splev(wave, interpolate.splrep(x, s, k=1, s=0))
    filt = np.where((wave < x[0]) | (wave > x[-1]), 0, filt)
    # The 1/hc factor cancels in the BD17/Vega ratio.
    return float(integrate.simpson(filt*f*wave, x=wave))

def union_response(band, spec):
    x, s = band.T
    wave, f = spec.T
    u = np.unique(np.r_[x, wave[(wave > x[0]) & (wave < x[-1])]])
    return float(np.trapezoid(np.interp(u,x,s)*np.interp(u,wave,f)*u,u))

for item in P['filter_mapping']:
    label=item['result_label']; name=item['snpy_curve']; instr=source_band[label]
    d = (OLD if name in {'J_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'} else SN/'filters/LCO')/instr/name
    assert hashlib.sha256(d.read_bytes()).hexdigest() == source_hash[name], d
    mirror=ROOT/'phase2/official/inputs/SNDATA_ROOT/filters/CSP/CSP_TAMU_20180316'/(mirror_name.get(label,name))
    assert d.read_bytes()==mirror.read_bytes(), (label,'curve mirror mismatch')
    registry = (OLD/instr/'filters.dat').read_text()
    assert any(line.split()[:3]==[item['snpy_filter'].split('/')[1],name,'VegaB=0'] for line in registry.splitlines())
    if item['raisin_magref'] is not None:
        assert any(line.startswith('FILTER: CSP-'+label+' ') and item['raisin_magref'] in line for line in RAISIN.read_text().splitlines())
    band=np.loadtxt(d)
    v=snpy_response(band,vega); b=snpy_response(band,ref)
    uv=union_response(band,vega); ub=union_response(band,ref)
    mag=float(-2.5*np.log10(b/v))
    mag_union=float(-2.5*np.log10(ub/uv))
    reported=float(item['raisin_magref']) if item['raisin_magref'] is not None else None
    diff=mag-reported if reported is not None else None
    half_print=0.5*10**(-len(item['raisin_magref'].split('.')[1])) if reported is not None else None
    result['bands'][label]={
       'snpy_filter':item['snpy_filter'], 'snpy_curve':name,
       'vega_photon_integral_without_hc':v, 'bd17_photon_integral_without_hc':b,
       'bd17_mag_vegaB0_simpson':mag, 'bd17_mag_vegaB0_union_trapz':mag_union,
       'quadrature_delta_mag':mag_union-mag,
       'raisin_magref':reported, 'difference_mag':diff,
       'half_print_digit_mag':half_print,
       'matches_print_rounding':abs(diff)<=half_print if diff is not None else None,
       'source_sha256':source_hash[name]
    }
(HERE/'comparison-result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2))
