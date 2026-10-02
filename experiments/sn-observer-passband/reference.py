#!/usr/bin/env python3
"""Independent high-precision/coordinate controls, never a production solver."""
import hashlib
import json
import math
from pathlib import Path
import sys
import mpmath as mp
import numpy as np

store=Path(sys.argv[1]);native=json.loads((store/'native.out').read_text());mp.mp.dps=80
checks=0;max_relative={};failures=[]
def compare(name,actual,expected,relative,absolute):
    global checks
    checks+=1;actual=mp.mpf(actual);expected=mp.mpf(expected);delta=abs(actual-expected)
    rel=delta/abs(expected) if expected else delta
    max_relative[name]=max(max_relative.get(name,0),float(rel))
    if delta>mp.mpf(absolute)+mp.mpf(relative)*abs(expected) or (expected>0 and actual<=0):failures.append(name)
for r in native['background']:
    z,o=mp.mpf(r['zHD']),mp.mpf(r['zHEL']);radial=mp.log1p(z);shape=(1+o)*radial
    compare('radial',r['radial'],radial,2e-12,1e-24)
    compare('shape',r['shape'],shape,2e-12,1e-24)
    compare('distance',r['DL_mpc'],shape*mp.mpf(299792.458)/70,2e-11,1e-9)
source_rows=[tuple(map(float,line.split())) for line in (store/'rows.tsv').read_text().splitlines()]
for i,(source_id,z,o) in enumerate(source_rows):
    a,b=native['background'][2*i:2*i+2]
    if (a['index'],a['zHD'],a['zHEL'],b['index'],b['zHD'],b['zHEL'])!=(2*i,z,o,2*i+1,z,z):failures.append('source_pair_order')
for i in range(0,642,2):
    a,b=native['background'][i:i+2]
    compare('observer_ratio',a['shape']/b['shape'],(1+mp.mpf(a['zHEL']))/(1+mp.mpf(a['zHD'])),2e-12,1e-24)
if any(native['background'][642][k]!=0 for k in ['radial','shape','DL_mpc']):failures.append('zero_limit')
low=native['background'][643]
compare('low_z_slope',low['shape']/low['zHD'],1+mp.mpf(low['zHD'])/2,2e-12,1e-24)
pairs=[tuple(map(float,line.split())) for line in (store/'passband-metre.tsv').read_text().splitlines()]
angstroms=[tuple(map(float,line.split())) for line in (store/'passband-angstrom.tsv').read_text().splitlines()]
if len(pairs)!=910 or len(angstroms)!=910:failures.append('passband_count')
for (a,ta),(b,tb) in zip(angstroms,pairs):
    checks+=1
    if float(mp.mpf(a)*mp.mpf('1e-10'))!=b or ta!=tb:failures.append('source_unit_or_transmission_identity')
x=np.array([a for a,b in pairs]);t=np.array([b for a,b in pairs]);h=mp.mpf('6.62607015e-34');c=mp.mpf(299792458)
width=mp.mpf(x[-1])-mp.mpf(x[0]);area=mp.mpf(0);weighted=mp.mpf(0)
for (a,ta),(b,tb) in zip(pairs,pairs[1:]):
    a,b,ta,tb=map(mp.mpf,[a,b,ta,tb]);d=b-a
    area+=d*(ta+tb)/2
    weighted+=d*((2*a+b)*ta+(a+2*b)*tb)/6
def frequency(n,z):
    nodes,weights=np.polynomial.legendre.leggauss(n)
    a,b=x[:-1],x[1:];lo,hi=float(c)/b,float(c)/a
    nu=(hi-lo)[:,None]/2*nodes+(hi+lo)[:,None]/2
    lam=float(c)/nu
    transmission=t[:-1,None]+(t[1:]-t[:-1])[:,None]*(lam-a[:,None])/(b-a)[:,None]
    f=1/(4*math.pi*1e40*(1+z));fnu=f*float(c)/nu**2
    dw=(hi-lo)[:,None]/2*weights
    return tuple(float(np.sum(values*dw)) for values in [fnu,fnu*transmission,fnu*transmission/(float(h)*nu)])
for r in native['photometry']:
    f=1/(4*mp.pi*mp.mpf(1e20)**2*(1+mp.mpf(r['z'])))
    expected=[f*width,f*area,f*weighted/(h*c)]
    fine=frequency(16,r['z']);coarse=frequency(8,r['z'])
    for k,e,ff,cc in zip(['flux','energy','photons'],expected,fine,coarse):
        compare('native_'+k,r[k],e,2e-12,1e-300)
        compare('frequency_'+k,ff,e,2e-13,1e-300)
        compare('frequency_refinement_'+k,cc,ff,2e-13,1e-300)
environment={}
for module in [mp,np]:
    path=Path(module.__file__);environment[module.__name__]={'version':module.__version__,'module_path':str(path),'module_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
environment['python']={'path':sys.executable,'version':sys.version,'sha256':hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest()}
result={'status':'accepted' if not failures else 'failed','checks':checks,'max_relative':max_relative,'failures':failures,'reference_identity':'80decimal coasting log1p and exact wavelength-linear moments; independent frequency GaussLegendre8/16 refinement','response_normalization':{'integral_T_dlambda_metre':float(area),'integral_lambda_T_dlambda_metre_squared':float(weighted),'raw_peak_preserved':max(t)},'environment':environment,'limits':'No source interpolation uncertainty, measured calibration distribution, photometric reduction, SN fit or posterior is established.'}
print(json.dumps(result));raise SystemExit(bool(failures))
