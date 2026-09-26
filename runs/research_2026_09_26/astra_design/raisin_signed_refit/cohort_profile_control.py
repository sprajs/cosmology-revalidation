"""Freeze and run the complete source-defined ten-object restricted profile."""
from pathlib import Path
import argparse,csv,json,os,subprocess,time,sys
import numpy as np
import paired_refit as p
from check_results import table
O=p.OUT;C=O/'fixed-c-profile/cohort10';ENGINE=O/'cohort_profile_engine.py'

def prepare():
    C.mkdir(exist_ok=False)
    members=list(csv.DictReader((O/'cohort.csv').open()));members.sort(key=lambda r:r['CID'])
    assert len(members)==10
    A={r['CID']:r for r in table(O/'fits/full/fixed/A/fit.FITRES.TEXT')};B={r['CID']:r for r in table(O/'fits/full/fixed/B/fit.FITRES.TEXT')}
    pilot=json.loads((O/'fixed-c-profile/global-profile/protocol.json').read_text())
    protocol=dict(status='Frozen adaptive extension after independently verified first-object proof, before remaining-object profile outcomes',
      cohort=[r['CID'] for r in members],selection='All ten DES16 source-lineage objects, unchanged from paired-refit cohort. First pilot is known; no pristine-validation claim.',
      benchmark_pair=[r['CID'] for r in members[:2]],same_method=pilot,
      changes=['Per-object current/audit nominal FITRES closure and exact mean/C export before profiling.',
      'A to B exact native row multiset mapping preserves duplicate multiplicities; identical-key covariance swap invariance required.',
      'Same scientific domains, covariance anchors, native oracle, resolution, amplitude/clipping and exact float32-boundary gates as passed pilot.',
      'Every object retains all failure/edge/limit results. No new omissions or automatic domain expansion.',
      'Report all competing shape branches and descriptive deltaQ1/4/9 ridge ranges plus amplitude-inclusive finite-grid distance level sets; these are not confidence/posterior intervals.',
      'Select matched alternate B branch by proximity to the positive-arm minimum shape, independent of distance-shift sign; retain full branch list.'],
      expected_masks={r['CID']:{'A':int(float(A[r['CID']]['NDOF']))+3,'B':int(float(B[r['CID']]['NDOF']))+3} for r in members},
      resources=dict(workers=1,max_total_active_seconds=900,max_native_calls_per_object=100000,benchmark_before_remaining=True,
        forecast='Pilot50.7s grid plus about25s native initialization per object: about13min for10; review two-object measured forecast before other8.'),
      hashes={str(f.relative_to(p.ROOT)):p.sha(f) for f in [Path(__file__),ENGINE,O/'global_fixed_profile.py',O/'native_mean_oracle.py',O/'native_profile_gate_v2.py',O/'cohort.csv',O/'input-manifest.json',O/'fits/full/fixed/A/fit.FITRES.TEXT',O/'fits/full/fixed/B/fit.FITRES.TEXT',O/'fixed-c-profile/build-manifest.json',O/'fixed-c-profile/global-profile/result.json']})
    p.save(C/'protocol.json',protocol)
    for r in members:
      cid=r['CID'];d=C/cid;d.mkdir()
      case=dict(pilot);case.update(CID=cid,status='Frozen cohort case; same restricted profile and gates',engine_sha256=p.sha(ENGINE),cohort_protocol_sha256=p.sha(C/'protocol.json'),
        expected_A_epochs=protocol['expected_masks'][cid]['A'],expected_B_epochs=protocol['expected_masks'][cid]['B'],
        alternate_anchor_parameters=np.array([float(A[cid][k]) for k in ['DLMAG','STRETCH','AV','PKMJD']],float).astype(np.float32).astype(float).tolist(),
        hashes={str(f.relative_to(p.ROOT)):p.sha(f) for f in [ENGINE,O/'native_mean_oracle.py',O/'native_profile_gate_v2.py',C/'protocol.json',O/'data/RSR_A'/f'{cid}.snana.dat',O/'data/RSR_B'/f'{cid}.snana.dat',O/'fits/full/fixed/A/fit.FITRES.TEXT',O/'fits/full/fixed/B/fit.FITRES.TEXT',O/'fixed-c-profile/build-manifest.json']})
      p.save(d/'protocol.json',case)
    patch=subprocess.run(['diff','-u',str(O/'global_fixed_profile.py'),str(ENGINE)],capture_output=True,text=True)
    (C/'engine-versus-pilot.patch').write_text(patch.stdout)
    p.save(C/'input-manifest.json',dict(protocol_sha256=p.sha(C/'protocol.json'),files_sha256={str(f.relative_to(C)):p.sha(f) for f in C.rglob('*') if f.is_file() and f.name!='input-manifest.json'}))
    print(json.dumps(dict(protocol_sha256=p.sha(C/'protocol.json'),cohort=protocol['cohort'],expected_masks=protocol['expected_masks']),indent=2))

def execute(stage):
    protocol=json.loads((C/'protocol.json').read_text())
    for f,h in protocol['hashes'].items():assert p.sha(p.ROOT/f)==h,f
    ids=protocol['benchmark_pair'] if stage=='benchmark' else protocol['cohort'][2:]
    if stage=='remaining':
      benchmark=json.loads((C/'benchmark-summary.json').read_text());assert benchmark['native_mapping_all_pass']
      assert benchmark['forecast_total_seconds']<=900,'Resource forecast needs an explicit pre-remaining amendment'
    records=[]
    for cid in ids:
      existing=list(C.glob('*/execution.json'));active=sum(json.loads(f.read_text())['seconds'] for f in existing)
      if active>=900:records.append(dict(CID=cid,status='not run: total resource cap'));continue
      d=C/cid;assert not(d/'runner.log').exists()
      env=os.environ.copy();env.update(RAISIN_PROFILE_CASE=str(d),RAISIN_PROFILE_CID=cid,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
      t=time.monotonic()
      with(d/'runner.log').open('x') as log:
        try:r=subprocess.run([sys.executable,str(ENGINE)],cwd=p.ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,900-active));rc=r.returncode
        except subprocess.TimeoutExpired:rc=-999
      elapsed=time.monotonic()-t;rr=dict(CID=cid,returncode=rc,seconds=elapsed,runner_log_sha256=p.sha(d/'runner.log'),case_protocol_sha256=p.sha(d/'protocol.json'))
      p.save(d/'execution.json',rr)
      if(d/'result.json').exists():
        result=json.loads((d/'result.json').read_text());rr.update(status=result['status'],numerical_gate_pass=result['numerical_gate_pass'],delta_DLMAG=result['delta_DLMAG_B_minus_A'],native_vectors=result['total_oracle_calls'])
      else:rr.update(status='native/reference gate failure; no profile result',numerical_gate_pass=False)
      records.append(rr);print(json.dumps(rr),flush=True)
    if stage=='benchmark':
      ok=all((C/r['CID']/'mask-mapping.json').exists() and(C/r['CID']/'reference-closure.json').exists() for r in records)
      forecast=sum(r['seconds'] for r in records)+8*max(r['seconds'] for r in records)
      p.save(C/'benchmark-summary.json',dict(records=records,native_mapping_all_pass=ok,forecast_total_seconds=forecast,
        scope='Same algorithms and domains; failures remain cohort members. Mapping/runtime gate precedes remaining profile outcomes.'))
      print(json.dumps({'benchmark_native_mapping_all_pass':ok,'forecast_total_seconds':forecast}),flush=True)
    p.save(C/(stage+'-execution.json'),dict(records=records,protocol_sha256=p.sha(C/'protocol.json')))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','benchmark','remaining']);a=ap.parse_args();prepare() if a.action=='prepare' else execute(a.action)
