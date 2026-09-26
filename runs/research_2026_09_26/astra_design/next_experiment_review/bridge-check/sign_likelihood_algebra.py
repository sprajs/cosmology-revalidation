"""Analytic/numerical sign-processing likelihood distinctions; no SN data used."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import norm
OUT=Path(__file__).resolve().parent
rows=[]
for f in [0.,.1,.5,1.,3.]:
    def folded(x):return norm.pdf(x-f)+norm.pdf(-x-f)
    def truncated(x):return norm.pdf(x-f)/ndtr(f)
    def sf(x):
        a=norm.pdf(x-f);b=norm.pdf(-x-f)
        return ((x-f)*a+(-x-f)*b)/max(a+b,1e-300)
    def st(x):return x-f-norm.pdf(f)/ndtr(f)
    r=dict(f_over_sigma=f,folded_norm=quad(folded,0,np.inf)[0],truncated_norm=quad(truncated,0,np.inf)[0],
           censored_norm=ndtr(-f)+quad(lambda x:norm.pdf(x-f),0,np.inf)[0],
           folded_expected_score=quad(lambda x:folded(x)*sf(x),0,15)[0],
           truncated_expected_score=quad(lambda x:truncated(x)*st(x),0,15)[0],
           folded_information=quad(lambda x:folded(x)*sf(x)**2,0,15)[0],
           truncated_information=quad(lambda x:truncated(x)*st(x)**2,0,15)[0],
           censored_information=quad(lambda x:norm.pdf(x-f)*(x-f)**2,0,15)[0]+norm.pdf(f)**2/ndtr(-f))
    rows.append(r)
for r in rows:
    assert max(abs(r[k]-1) for k in ['folded_norm','truncated_norm','censored_norm'])<1e-10
    assert max(abs(r[k]) for k in ['folded_expected_score','truncated_expected_score'])<1e-10
assert abs(rows[0]['folded_information'])<1e-15
assert abs(rows[0]['truncated_information']-(1-2/np.pi))<1e-12
assert abs(rows[0]['censored_information']-(.5+1/np.pi))<1e-12
out=dict(status='Synthetic likelihood normalization and expected-score PASS; no cause of released positive flux established.',
    units='sigma=1; information is multiplied by1/sigma^2 in dimensional form',
    zero_mean_positive_observation_mean=float(np.sqrt(2/np.pi)),results=rows,
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'sign-likelihood-algebra.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
