"""Independent bisection and simultaneous-interval arithmetic, no imported fit code."""
from pathlib import Path
import json,hashlib,math
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parents[2];D=ROOT/'runs/research_2026_09_26/bao_shape/curvature_certificate'
m=json.loads((D/'manifest.json').read_text());checks=0
for p,h in m['inputs_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h;checks+=1
for p,h in m['outputs_sha256'].items():assert hashlib.sha256((D/p).read_bytes()).hexdigest()==h;checks+=1
r=json.loads((D/'result.json').read_text());B=r['BGS'];H=r['highest_transverse'];t=math.log(1+B['z']);a=B['z']/(1+B['z']);records=[];maxerr=0.
for item in r['radial_results']:
 def calculation(k):
  bu=B['value']+k*B['sigma'];gl=item['g']-k*item['sigma_g'];dl=H['value']-k*H['sigma']
  assert min(bu,gl,dl)>0
  upperdm=math.sqrt(bu**3/(a*gl));theta=t*gl/dl;assert dl>upperdm and theta<math.pi/2
  lower=(a*gl*(dl*math.sin(theta))**2)**(1/3)
  return lower-bu,lower,upperdm
 if calculation(0)[0]<=0:k=0.
 else:
  lo=0.;hi=3.
  assert calculation(hi)[0]<0
  for _ in range(70):
   mid=(lo+hi)/2
   if calculation(mid)[0]>0:lo=mid
   else:hi=mid
  k=(lo+hi)/2
 val=calculation(k);tail=.5*math.erfc(k/math.sqrt(2));maxerr=max(maxerr,abs(k-item['certificate_k']),abs(min(1,8*tail)-item['six_radial_union_bound']))
 records.append(dict(z=item['z'],k=k,union8=min(1,8*tail),minimum_BGS_central=calculation(0)[1],upper_DM_at_boundary=val[2]))
assert maxerr<1e-10
out={'status':'PASS independent scalar bisection and eight one-sided marginal union bound','hashes_verified':checks,'maximum_numeric_difference':maxerr,'radial_results':records,'scope':'Finite six-radial family at fixed highest-z transverse datum; exploration outside this family is not covered. Arbitrary quoted cross-covariance allowed by union bound; Gaussian marginal calibration, common ruler, FLRW and first-antipode branch remain assumptions.'}
(P/'curvature-certificate-independent.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
