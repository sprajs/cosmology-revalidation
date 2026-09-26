"""Recheck the explicitly bounded kernel after independent review."""
import json
from pathlib import Path
import numpy as np
from scipy.special import log_ndtr
import onefactor_sign_validation as base
import onefactor_sign_likelihood as kernel

original = base.OUT
base.OUT = original / 'v2'
protocol = json.loads((original / 'protocol.json').read_text())
protocol['amendment_sha256'] = base.SHA(base.OUT / 'amendment.json')
for path in [Path(kernel.__file__), Path(__file__)]:
    protocol['inputs_sha256'][str(path.relative_to(base.ROOT))] = base.SHA(path)
assert not (base.OUT / 'protocol.json').exists()
base.dump(base.OUT / 'protocol.json', protocol)
base.run()
scalar = kernel.log_orthant([-.2], [1.], [1000.], [1.])
exact = float(log_ndtr(-.2/np.sqrt(1+1000.**2)))
errors = {}
for name, call in [
    ('strong_multirow', lambda: kernel.log_orthant([0.,0.], [1.,1.], [1000.,1000.], [1.,1.])),
    ('remote_latent_mode', lambda: kernel.log_orthant([-100.,-100.], [1.,1.], [.4,.4], [1.,1.])),
    ('sparse_two_row_factorization', lambda: kernel.decompose_covariance(np.array([[1.,.1],[.1,1.]]))),
]:
    try:
        call()
    except ValueError as exc:
        errors[name] = str(exc)
    else:
        raise AssertionError(name+' must refuse unsupported domain')
assert scalar == exact and len(errors) == 3
base.dump(base.OUT / 'scope-check.json', {'scalar_logprob':scalar,'exact_scalar_logprob':exact,'explicit_refusals':errors})
(base.OUT / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
base.dump(base.OUT/'manifest.json',{'files_sha256':{str(f.relative_to(base.OUT)):base.SHA(f) for f in base.OUT.iterdir() if f.is_file() and f.name!='manifest.json'}})
print(json.dumps({'scalar_error':scalar-exact,'refusals':list(errors)}))
