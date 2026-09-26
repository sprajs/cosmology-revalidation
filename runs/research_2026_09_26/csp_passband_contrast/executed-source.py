"""Frozen synthetic CSP instrument contrasts with the archived KCOR spectral SED.

No observed SN fluxes, fitted distances, dust inference or cosmology are scored.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
from astropy.io import fits
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/research_2026_09_26/csp_passband_contrast'
FILTERS = ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
KCOR = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
PAIRS = {'J_WIRC_minus_RC1': ('Jrc1_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'),
         'H_WIRC_minus_RetroCam': ('H_SWO_TAM_scan_atm.dat','H_DUP_TAM_scan_atm.dat'),
         'Y_WIRC_minus_RetroCam': ('Y_SWO_TAM_scan_atm.dat','Y_DUP_TAM_scan_atm.dat'),
         'J_RC2_minus_RC1': ('Jrc1_SWO_TAM_scan_atm.dat','Jrc2_SWO_TAM_scan_atm.dat')}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, x):
    p.write_text(json.dumps(x, indent=2, allow_nan=False)+'\n')


def prepare():
    OUT.mkdir(exist_ok=False)
    paths = [KCOR, FILTERS/'README', Path(__file__)]
    paths += [FILTERS/f for f in sorted({s for pair in PAIRS.values() for s in pair})]
    dump(OUT/'protocol.json', {
        'status': 'Frozen before synthetic integral outcomes; known metadata label merging motivated this diagnostic.',
        'phase_days': [-7,0,7,10,15,20,30], 'redshifts': [0,.01,.03,.05,.08],
        'pairs': PAIRS,
        'SED': 'Exact SN SED extension of released CSP KCOR; source input names Hsiao07. No new spectral model fitted.',
        'reference': 'Exact BD17 PrimarySED extension, zero differential assigned to the reference star. Static natural magnitude offsets excluded; phase/redshift differences cancel them.',
        'integral': 'Photon counts integral lambda*T(lambda)*f_lambda(lambda/(1+z))/(1+z). Signed tabulated transmission tails retained, zero outside each table. Piecewise linear SED/throughput, exact degree3 Gauss2 on union knots checked with Gauss4 and independent 0.5A trapezoid.',
        'gates': 'All support within SED/reference wavelengths; positive finite integrals; Gauss2/4 delta-mag difference <=1e-10; independent0.5A difference <=1e-4mag.',
        'interpretation': 'Conditional passband shape sensitivity only. Does not establish original transformation, natural zero points, actual SN SED, trained-model independence, selection or corrected distances.',
        'inputs_sha256': {str(p.relative_to(ROOT)):sha(p) for p in paths},
    })


def counts(wave, sed, band, order):
    lo, hi = band[0,0], band[-1,0]
    assert wave[0] <= lo and wave[-1] >= hi
    if order == 0:
        x = np.linspace(lo, hi, int(np.ceil((hi-lo)/.5))+1)
        return np.trapezoid(x*np.interp(x,band[:,0],band[:,1])*np.interp(x,wave,sed),x)
    knots = np.unique(np.r_[band[:,0], wave[(wave>lo)&(wave<hi)]])
    mid = (knots[1:]+knots[:-1])/2
    half = (knots[1:]-knots[:-1])/2
    nodes, weights = leggauss(order)
    x = mid[:,None]+half[:,None]*nodes
    value = x*np.interp(x,band[:,0],band[:,1])*np.interp(x,wave,sed)
    return float(np.sum(half*(value @ weights)))


def run():
    p = json.loads((OUT/'protocol.json').read_text())
    for name,digest in p['inputs_sha256'].items():
        assert sha(ROOT/name)==digest,name
    with fits.open(KCOR) as hdus:
        hdr = hdus['SN SED'].header
        wave = hdr['LMIN']+np.arange(hdr['NBL'])*hdr['LBIN']
        times = hdr['TMIN']+np.arange(hdr['NBT'])*hdr['TBIN']
        spectra = np.array(hdus['SN SED'].data.field(0),float).reshape(len(times),len(wave))
        ref_wave = np.array(hdus['PrimarySED'].data.field(0),float)
        ref_flux = np.array(hdus['PrimarySED'].data['BD17'],float)
    bands = {name:np.loadtxt(FILTERS/name) for pair in PAIRS.values() for name in pair}
    assert all(np.isfinite(x).all() and np.all(np.diff(x[:,0])>0) for x in bands.values())
    ref = {o:{name:counts(ref_wave,ref_flux,b,o) for name,b in bands.items()} for o in [2,4,0]}
    rows=[]
    for phase in p['phase_days']:
        ix=np.flatnonzero(times==phase);assert len(ix)==1
        spectrum=spectra[ix[0]]
        for z in p['redshifts']:
            sc={o:{name:counts(wave*(1+z),spectrum/(1+z),b,o) for name,b in bands.items()} for o in [2,4,0]}
            assert all(v>0 and np.isfinite(v) for group in sc.values() for v in group.values())
            for name,(a,b) in PAIRS.items():
                dm={o:float(-2.5*np.log10((sc[o][b]/ref[o][b])/(sc[o][a]/ref[o][a]))) for o in [2,4,0]}
                rows.append({'pair':name,'phase':phase,'z':z,'delta_mag_equal_BD17':dm[2],
                             'gauss2_4_gap':abs(dm[2]-dm[4]),'trapezoid_halfA_gap':abs(dm[2]-dm[0])})
    g=max(r['gauss2_4_gap'] for r in rows);t=max(r['trapezoid_halfA_gap'] for r in rows)
    assert g<=1e-10 and t<=1e-4,(g,t)
    for r in rows:
        atzero=next(v for v in rows if v['pair']==r['pair'] and v['phase']==0 and v['z']==r['z'])
        atz0=next(v for v in rows if v['pair']==r['pair'] and v['phase']==r['phase'] and v['z']==0)
        r['phase_contrast_vs_phase0']=r['delta_mag_equal_BD17']-atzero['delta_mag_equal_BD17']
        r['redshift_contrast_vs_z0']=r['delta_mag_equal_BD17']-atz0['delta_mag_equal_BD17']
    with (OUT/'grid.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summaries={}
    for name in PAIRS:
        group=[r for r in rows if r['pair']==name]
        summaries[name]={k:[min(r[k] for r in group),max(r[k] for r in group)] for k in ['delta_mag_equal_BD17','phase_contrast_vs_phase0','redshift_contrast_vs_z0']}
    result={'pass':True,'grid_points':len(rows),'max_gauss2_4_gap':g,'max_halfA_trapezoid_gap':t,
            'signed_negative_transmission_points':{k:int(np.sum(b[:,1]<0)) for k,b in bands.items()},
            'summaries':summaries,'protocol_sha256':sha(OUT/'protocol.json'),'scope':p['interpretation']}
    dump(OUT/'result.json',result)
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    dump(OUT/'manifest.json',{'files_sha256':{f.name:sha(f) for f in OUT.iterdir() if f.is_file() and f.name!='manifest.json'}})
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    {'prepare':prepare,'run':run}[sys.argv[1]]()
