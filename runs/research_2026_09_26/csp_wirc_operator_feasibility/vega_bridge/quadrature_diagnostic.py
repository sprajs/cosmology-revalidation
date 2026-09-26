#!/usr/bin/env python3
"""Post-result exact piecewise-linear quadrature; no source/version selection."""
import hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
P=json.loads((HERE/'comparison-protocol.json').read_text())
ORIGINAL=json.loads((HERE/'comparison-result.json').read_text())
S=HERE/'snpy/snpy/filters'
O=HERE.parent/'snpy/snpy/filters/filters/LCO'
vega=np.loadtxt(S/'standards/Vega/alpha_lyr_stis_005.ascii')
bd=np.loadtxt(ROOT/'phase2/official/inputs/SNDATA_ROOT/standards/bd_17d4708_stisnic_007.dat')
assert hashlib.sha256((S/'standards/Vega/alpha_lyr_stis_005.ascii').read_bytes()).hexdigest()==P['input_hashes']['vegaB_spectrum']

def exact_piecewise(band,spec):
    x,s=band.T;w,f=spec.T
    assert w[0]<=x[0] and w[-1]>=x[-1]
    u=np.unique(np.r_[x,w[(w>x[0])&(w<x[-1])]])
    left,right=u[:-1],u[1:]
    mid=(left+right)*.5;half=(right-left)*.5
    q=half/np.sqrt(3.)
    a,b=mid-q,mid+q
    fun=lambda z:np.interp(z,x,s)*np.interp(z,w,f)*z
    return float(np.sum(half*(fun(a)+fun(b))))
R={'protocol_sha256':hashlib.sha256((HERE/'quadrature-diagnostic-protocol.json').read_bytes()).hexdigest(),'bands':{}}
for item in P['filter_mapping']:
    label=item['result_label'];name=item['snpy_curve'];instr=item['snpy_filter'].split('/')[0]
    base=O if name in {'J_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'} else S/'filters/LCO'
    band=np.loadtxt(base/instr/name)
    vb=exact_piecewise(band,vega);bb=exact_piecewise(band,bd)
    mag=float(-2.5*np.log10(bb/vb))
    old=ORIGINAL['bands'][label]
    R['bands'][label]={'bd17_mag_vegaB0_exact_piecewise':mag,
       'minus_frozen_simpson':mag-old['bd17_mag_vegaB0_simpson'],
       'minus_union_trapz':mag-old['bd17_mag_vegaB0_union_trapz'],
       'minus_raisin_magref':mag-old['raisin_magref'] if old['raisin_magref'] is not None else None}
(HERE/'quadrature-diagnostic-result.json').write_text(json.dumps(R,indent=2,sort_keys=True)+'\n')
print(json.dumps(R,indent=2))
