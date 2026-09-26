"""Execute the predeclared classifier-cut secondary without changing64 records."""
from pathlib import Path
import json,hashlib
import residual64 as run
p=Path(__file__).resolve().parent
out=p/'residual43'
if out.exists():raise RuntimeError('Preserve completed run')
out.mkdir()
(out/'residual64-protocol.md').write_bytes((p/'residual43-protocol.md').read_bytes())
original=run.build
def subset():
    m,arms,info=original();keep=(m.PROB_SNNV19>.999).to_numpy()
    assert keep.sum()==43
    ids=set(m.loc[keep,'CID'])
    for arm,a in arms.items():
        for key in ['u','F','chi2']:a[key]=a[key][keep]
        a['templates']=[v for v,k in zip(a['templates'],keep) if k]
    return m.loc[keep].reset_index(drop=True),arms,[v for v in info if v['CID'] in ids]
run.build=subset;run.OUT=out
run.main()
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
(out/'wrapper-manifest.json').write_text(json.dumps({'wrapper_sha256':sha(__file__),
 'base_script_sha256':sha(p/'residual64.py'),'protocol_sha256':sha(p/'residual43-protocol.md'),
 'status':'Secondary after64 outcomes; exact published classifier cut, no residual exclusions.'},indent=2)+'\n')
