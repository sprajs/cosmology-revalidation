"""Independent analytically soluble check of raw importance-weight summaries."""
import hashlib
import json
from pathlib import Path
import numpy as np
from exact_correction import summarize,ROOT


def main():
    rng=np.random.default_rng(272642);n=40000;shift=np.array([.2,-.3])
    draws=rng.normal(size=(n,2))
    # Unit-covariance target N(shift,I), proposal N(0,I): exact log density ratio.
    logweights=draws@shift-.5*shift@shift
    records=[{'status':'finite','point':{'x':float(x),'y':float(y)},'derived':{},'log_weight':float(w)}
             for (x,y),w in zip(draws,logweights)]
    report=summarize(records,np.repeat(np.arange(4),n//4))
    mean_errors=np.array([report['posterior'][k]['mean'] for k in ['x','y']])-shift
    mcse=np.array([report['posterior'][k]['weighted_batch_mean_mcse'] for k in ['x','y']])
    sd_errors=np.array([report['posterior'][k]['sd'] for k in ['x','y']])-1
    expected_ess_fraction=float(np.exp(-shift@shift))
    assert np.all(abs(mean_errors)<4*mcse)
    assert np.max(abs(sd_errors))<.025
    assert abs(report['raw_weight_ess']/n-expected_ess_fraction)<.02
    assert report['status']=='passed_importance_weight_gates'
    output={'status':'passed','draws':n,'target_mean':shift.tolist(),
        'mean_errors':mean_errors.tolist(),'mean_mcse':mcse.tolist(),'sd_errors':sd_errors.tolist(),
        'observed_weight_ess_fraction':report['raw_weight_ess']/n,
        'analytic_asymptotic_weight_ess_fraction':expected_ess_fraction,
        'pareto_k':report['pareto_k'],
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('exact_correction.py')]},
        'scope':'Algebra and numerical summary check under known Gaussian density ratios, not validation of cosmological overlap or survey coverage.'}
    (ROOT/'studies/unified_cosmology/results/inference/correction-validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
