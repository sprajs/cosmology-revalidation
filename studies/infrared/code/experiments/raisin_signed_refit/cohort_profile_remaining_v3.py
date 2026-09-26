"""Remaining eight under unchanged scientific protocol plus strict initial stream gate."""
from pathlib import Path
import json,os,subprocess,time,sys
import paired_refit as p
O=p.OUT;C=O/'fixed-c-profile/cohort10';ENGINE=O/'cohort_profile_engine_v3.py'
def run():
    am=json.loads((C/'prefix-gate-amendment.json').read_text());assert p.sha(ENGINE)==am['engine_sha256'] and p.sha(Path(__file__))==am['controller_sha256']
    b=json.loads((C/'benchmark-v2-summary.json').read_text());assert b['native_mapping_and_stream_pass'] and b['forecast_total_seconds']<=900
    proto=json.loads((C/'protocol.json').read_text());index=json.loads((C/'active-attempts.json').read_text())['attempts'];records=[]
    for cid in proto['cohort'][2:]:
      executions=list(C.glob('*/execution.json'))+list(C.glob('*/v2/execution.json'))+list(C.glob('*/v3/execution.json'))
      active=sum(json.loads(f.read_text())['seconds'] for f in executions)
      if active>=900:records.append(dict(CID=cid,status='not run: total active resource cap'));continue
      d=C/index[cid];assert not(d/'runner.log').exists();env=os.environ.copy();env.update(RAISIN_PROFILE_CASE=str(d),RAISIN_PROFILE_CID=cid,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
      t=time.monotonic()
      with(d/'runner.log').open('x') as f:
        try:r=subprocess.run([sys.executable,str(ENGINE)],cwd=p.ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=max(1,900-active));rc=r.returncode
        except subprocess.TimeoutExpired:rc=-999
      r=dict(CID=cid,returncode=rc,seconds=time.monotonic()-t,case_protocol_sha256=p.sha(d/'protocol.json'),runner_log_sha256=p.sha(d/'runner.log'))
      p.save(d/'execution.json',r)
      if(d/'result.json').exists():
        q=json.loads((d/'result.json').read_text());r.update(status=q['status'],numerical_gate_pass=q['numerical_gate_pass'],delta_DLMAG=q['delta_DLMAG_B_minus_A'],native_vectors=q['total_oracle_calls'])
      else:r.update(status='native/stream gate failure; no profile result',numerical_gate_pass=False)
      records.append(r);print(json.dumps(r),flush=True)
      p.save(C/'remaining-v3-execution.json',dict(records=records,amendment_sha256=p.sha(C/'prefix-gate-amendment.json')))
if __name__=='__main__':run()
