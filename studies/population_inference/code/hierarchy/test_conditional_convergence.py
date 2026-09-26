"""Exercise the pre-score gate without opening any held-out data."""
import numpy as np
from conditional import convergence_diagnostics

rng=np.random.default_rng(20260926)
samples={'parameter':rng.normal(size=(4,1000))}
extra={'accept_prob':np.full((4,1000),.9),
       'num_steps':np.full((4,1000),7),
       'diverging':np.zeros((4,1000),dtype=bool)}
assert convergence_diagnostics(samples,extra)['valid_for_scoring']
assert not convergence_diagnostics({'parameter':samples['parameter'][:2]},
                                   {k:v[:2] for k,v in extra.items()})['valid_for_scoring']
assert not convergence_diagnostics({'parameter':np.zeros((4,1000))},extra)['valid_for_scoring']
for key,value in [('accept_prob',np.zeros((4,1000))),
                  ('num_steps',np.full((4,1000),1023)),
                  ('diverging',np.ones((4,1000),dtype=bool))]:
    bad={**extra,key:value}
    assert not convergence_diagnostics(samples,bad)['valid_for_scoring']
print('Conditional convergence gate rejects immobility, insufficient chains, low acceptance, deep trees and divergences')
