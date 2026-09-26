"""Synthetic algebra only; no RAISIN flux values or fitted quantities enter."""
import json, math
from pathlib import Path
P=Path(__file__).resolve().parent
q=.0075
epsilon=q/2.5

def signal(s,b,e):
    if b==0: return s**(1+e)
    return b**(1+e)*math.expm1((1+e)*math.log1p(s/b))

def k(s,b):
    if b==0: return math.log10(s)
    return (math.log(b)+(1+b/s)*math.log1p(s/b))/math.log(10)

def bias(s,b,sa,ba,e=epsilon):
    return -2.5*math.log10((signal(s,b,e)/s)/(signal(sa,ba,e)/sa))

checks=[]
max_deriv_error=0.
for s in [1e-6,.01,1.,100.,10000.]:
 for b in [0.,1.,10.]:
    sa=10000.;ba=1.
    eps=1e-5
    fd=(bias(s,b,sa,ba,eps)-bias(s,b,sa,ba,-eps))/(2*eps)
    expect=2.5*(k(sa,ba)-k(s,b))
    max_deriv_error=max(max_deriv_error,abs(fd-expect))
    checks.append(dict(source=s,background=b,calibrator_source=sa,calibrator_background=ba,exact_mag_bias=bias(s,b,sa,ba),first_order_mag_bias=q*(k(sa,ba)-k(s,b))))
assert max_deriv_error<1e-8
# Source-free background limit of incremental responsivity, evaluated stably.
s=1e-10;b=3.
limit=(1+epsilon)*b**epsilon
limit_error=abs(signal(s,b,epsilon)/s-limit)
assert limit_error<1e-10
# Same source/background regime cancels exactly after choosing its calibration pivot.
assert bias(1.,3.,1.,3.)==0.
# Multiplying response g by a common calibration factor cancels from the ratio.
# No-sky power law gives exact q per dex, fixing all signs and units.
assert abs(bias(1.,0.,10000.,0.)-4*q)<1e-14
result=dict(scope='Synthetic response-operator/algebra checks only; no measured RAISIN correction',q_mag_per_dex=q,epsilon=epsilon,max_first_derivative_error=max_deriv_error,faint_background_limit_error=limit_error,illustrative_four_dex=dict(delta_mag_unadjusted=4*q,correction_mag=-4*q,uncertainty_mag=4*.0006,flux_multiply=10**(.4*4*q),distance_multiply=10**(-4*q/5)),illustrative_004mag=dict(flux_multiply=10**(.4*.04),distance_multiply=10**(-.04/5)),cases=checks)
(P/'algebra-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
