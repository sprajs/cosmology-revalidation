"""Versioned retry after explicit per-object oracle-coordinate gate failure."""
from pathlib import Path
import argparse,json,os,subprocess,time,sys
import paired_refit as p
O=p.OUT;C=O/'fixed-c-profile/cohort10';ENGINE=O/'cohort_profile_engine_v2.py'

def prepare():
    proto=json.loads((C/'protocol.json').read_text());assert not(C/'oracle-coordinate-amendment.json').exists()
    amendment=dict(status='Benchmark adapter failure retained; C3 and remaining objects have no profile outcomes',
      failure='Generic wrapper omitted object-specific pars when constructing Oracle; inherited original pilot coordinates. C3 reference54rows versus stream57 triggered pre-profile assertion.',
      correction='Only replace Oracle(G/stream) by Oracle(G/stream,pars=pars). No model, data, masks, covariance, domain, score, or resolution changes.',
      benchmark_gate='Require native reference+mask gates and a saved result from a completed stream (runner exit0). Keep scientific edge/limit failures distinct.',
      original_protocol_sha256=p.sha(C/'protocol.json'),old_engine_sha256=p.sha(O/'cohort_profile_engine.py'),new_engine_sha256=p.sha(ENGINE),controller_sha256=p.sha(Path(__file__)),
      resource='Original900s total active runtime, including failed benchmark, retained. One worker;100000calls/object. Forecast reviewed before remaining8.',
      first_pilot_replay='C1cim new generic run exactly matches independently verified original; no rerun for one-line explicit-argument correction.')
    p.save(C/'oracle-coordinate-amendment.json',amendment)
    index={proto['cohort'][0]:str(Path(proto['cohort'][0]))}
    for cid in proto['cohort'][1:]:
      d=C/cid/'v2';d.mkdir();case=json.loads((C/cid/'protocol.json').read_text());case['engine_sha256']=p.sha(ENGINE)
      case['oracle_coordinate_amendment_sha256']=p.sha(C/'oracle-coordinate-amendment.json')
      old=str((O/'cohort_profile_engine.py').relative_to(p.ROOT));case['hashes'].pop(old)
      case['hashes'][str(ENGINE.relative_to(p.ROOT))]=p.sha(ENGINE)
      case['hashes'][str((C/'oracle-coordinate-amendment.json').relative_to(p.ROOT))]=p.sha(C/'oracle-coordinate-amendment.json')
      p.save(d/'protocol.json',case);index[cid]=str(Path(cid)/'v2')
    p.save(C/'active-attempts.json',dict(attempts=index,all_original_attempts_preserved=True,amendment_sha256=p.sha(C/'oracle-coordinate-amendment.json')))
    patch=subprocess.run(['diff','-u',str(O/'cohort_profile_engine.py'),str(ENGINE)],capture_output=True,text=True);(C/'oracle-coordinate-fix.patch').write_text(patch.stdout)
    print(p.sha(C/'oracle-coordinate-amendment.json'))

def run(stage):
    amendment=json.loads((C/'oracle-coordinate-amendment.json').read_text());assert p.sha(ENGINE)==amendment['new_engine_sha256'] and p.sha(Path(__file__))==amendment['controller_sha256']
    proto=json.loads((C/'protocol.json').read_text());index=json.loads((C/'active-attempts.json').read_text())['attempts']
    ids=[proto['cohort'][1]] if stage=='benchmark' else proto['cohort'][2:]
    if stage=='remaining':
      b=json.loads((C/'benchmark-v2-summary.json').read_text());assert b['native_mapping_and_stream_pass'] and b['forecast_total_seconds']<=900
    records=[]
    for cid in ids:
      executions=list(C.glob('*/execution.json'))+list(C.glob('*/v2/execution.json'))
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
    if stage=='benchmark':
      first=json.loads((C/proto['cohort'][0]/'execution.json').read_text());last=records[0]
      active=sum(json.loads(f.read_text())['seconds'] for f in list(C.glob('*/execution.json'))+list(C.glob('*/v2/execution.json')))
      forecast=active+8*max(first['seconds'],last['seconds']);passed=all(r['returncode']==0 and(C/index[r['CID']]/'result.json').exists() for r in records)
      p.save(C/'benchmark-v2-summary.json',dict(records=records,native_mapping_and_stream_pass=passed,forecast_total_seconds=forecast,total_active_seconds_so_far=active))
      print(json.dumps(dict(benchmark_native_mapping_and_stream_pass=passed,forecast_total_seconds=forecast)),flush=True)
    p.save(C/(stage+'-v2-execution.json'),dict(records=records,amendment_sha256=p.sha(C/'oracle-coordinate-amendment.json')))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','benchmark','remaining']);a=ap.parse_args();prepare() if a.action=='prepare' else run(a.action)
