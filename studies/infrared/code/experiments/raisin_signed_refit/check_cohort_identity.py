"""Exact metric-reference versus mean-oracle row identity, independent certification gate."""
from pathlib import Path
import json,hashlib,numpy as np
O=Path(__file__).resolve().parent;C=O/'fixed-c-profile/cohort10'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
    index=json.loads((C/'active-attempts.json').read_text())['attempts'];records=[]
    for cid,rel in index.items():
      d=C/rel
      if not(d/'result.json').exists():continue
      ref=np.load(d/'reference/native.npz');nom=np.load(d/'fixed_reference/native.npz');end=np.load(d/'stream/native.npz')
      checks={k:bool(np.array_equal(ref[k],nom[k]) and np.array_equal(nom[k],end[k])) for k in ['band','MJD','data_flux','data_fluxerr']}
      assert all(checks.values()),(cid,checks)
      assert np.array_equal(nom['parameters'],end['parameters'])
      records.append(dict(CID=cid,attempt=rel,epochs=len(ref['MJD']),exact_data_rows_checks=checks,
        covariance_ref_nom_equality_NOT_REQUIRED=True,source_hashes={f:sha(d/f) for f in ['reference/native.npz','fixed_reference/native.npz','stream/native.npz']}))
    out=dict(status='PASS for every completed case',complete_count=len(records),records=records,source_sha256=sha(__file__))
    (C/'metric-oracle-identity-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],complete_count=len(records))))
if __name__=='__main__':main()
