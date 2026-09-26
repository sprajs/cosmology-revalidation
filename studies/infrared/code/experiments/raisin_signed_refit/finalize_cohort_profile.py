"""Final closed-file integrity inventory; preserves provisional engine manifests."""
from pathlib import Path
import json,hashlib
O=Path(__file__).resolve().parent;C=O/'fixed-c-profile/cohort10'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
  idx=json.loads((C/'active-attempts.json').read_text())['attempts'];ledger=[]
  for cid,rel in idx.items():
    d=C/rel;r=json.loads((d/'result.json').read_text());m=json.loads((d/'manifest.json').read_text());checks={};appended=False
    for f,h in m['files_sha256'].items():
      actual=sha(d/f)
      if actual==h:checks[f]='exact final hash';continue
      assert f=='runner.log',(cid,f)
      b=(d/f).read_bytes();suffix=(json.dumps(r,indent=2)+'\n').encode();assert b.endswith(suffix)
      assert hashlib.sha256(b[:-len(suffix)]).hexdigest()==h
      checks[f]='original hash matches exact prefix before final result JSON appended';appended=True
    ledger.append(dict(CID=cid,attempt=rel,provisional_manifest_all_entries_resolved=True,runner_final_result_append_verified=appended,
      provisional_manifest_sha256=sha(d/'manifest.json'),closed_log_sha256=sha(d/'runner.log'),entries=checks))
  (C/'provisional-manifest-closure.json').write_text(json.dumps(dict(records=ledger,source_sha256=sha(__file__),original_manifests_unchanged=True),indent=2)+'\n')
  sources=['cohort_profile_control.py','cohort_profile_engine.py','cohort_profile_resume.py','cohort_profile_engine_v2.py','cohort_profile_remaining_v3.py','cohort_profile_engine_v3.py','oracle_prefix_gate.py','native_mean_oracle.py','native_profile_gate_v2.py','cohort_profile_geometry.py','cohort_profile_geometry_v2.py','check_cohort_identity.py','report_cohort_profile.py','finalize_cohort_profile.py']
  m=dict(status='Final closed-artifact hashes, including original failed attempts and provisional manifests',files_sha256={str(f.relative_to(C)):sha(f) for f in sorted(C.rglob('*')) if f.is_file() and f.name!='cohort-manifest.json'},executed_source_sha256={str(O/f):sha(O/f) for f in sources})
  (C/'cohort-manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps(dict(final_cohort_files=len(m['files_sha256']),source_files=len(sources),cohort_manifest_sha256=sha(C/'cohort-manifest.json'))))
if __name__=='__main__':main()
