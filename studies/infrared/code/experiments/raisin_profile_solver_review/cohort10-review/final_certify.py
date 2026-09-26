"""Final closed hashes, aggregate fields and actual native amplitude replies."""
from pathlib import Path
import csv,hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parent
C=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/cohort10'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=json.loads((C/'cohort-manifest.json').read_text())
for f,h in manifest['files_sha256'].items():assert sha(C/f)==h,f
for f,h in manifest['executed_source_sha256'].items():assert sha(f)==h,f
cohort=json.loads((C/'protocol.json').read_text());index=json.loads((C/'active-attempts.json').read_text())['attempts']
summary=json.loads((C/'cohort-summary.json').read_text())['objects'];reported={r['CID']:r for r in summary};assert set(reported)==set(cohort['cohort'])
outputs={};raw_checks={};comparison=[]
for cid in cohort['cohort']:
    own=json.loads((OUT/(cid+'.json')).read_text());outputs[cid]=own;g=C/index[cid];rr=reported[cid]
    assert rr['attempt']==index[cid] and rr['rows_A']==own['epochs']['A'] and rr['rows_B']==own['epochs']['B'] and rr['added_negative']==own['negative_rows']
    assert rr['numerical_gate_pass']==own['numerical_gate_pass']
    for anchor in ['B','A']:assert abs(rr['delta_'+anchor+'_anchor']-own['paired'][anchor]['global_minimum_delta_DLMAG'])<1e-10
    assert abs(rr['matched_shape_delta']-own['paired']['B']['matched_mode_delta_DLMAG'])<1e-10
    for arm in ['A','B']:
        for level in [1,4,9]:
            ll=own['geometry']['Banchor_'+arm]['level_sets'][str(float(level))];pair=ll['amplitude_inclusive_distance_range']
            assert max(abs(rr[f'{arm}_level{level}_low']-pair[0]),abs(rr[f'{arm}_level{level}_high']-pair[1]))<1e-9
            assert rr[f'{arm}_level{level}_shape_edge']==ll['shape_box_touched'] and rr[f'{arm}_level{level}_AV_edge']==ll['AV_box_touched']
    # Read actual-D native replies, not merely the reported error column.
    with np.load(g/'native-profiles.npz') as z:
        coordinates=z['coordinates'];H=z['model_means'];par=z['reference_parameters']
    lookup={tuple(v):i for i,v in enumerate(coordinates)};n=len(H);m=H.shape[1]
    with (g/'competitive-amplitude-checks.csv').open() as f:checks=list(csv.DictReader(f))
    sequence=0;baseline=None;h=None;errors=[];pointwise=[];scale_report_errors=[]
    with (g/'stream/fit.log').open() as f:
        for line in f:
            if not line.startswith('RAISIN_MEAN:'):continue
            q=line.split(maxsplit=3);sequence+=1;assert int(q[1])==sequence and int(q[2])==m
            if 1<sequence<=n+1:continue  # Saved means were independently scored; raw full-vector pilot audit is retained.
            values=np.fromstring(q[3],sep=' ');p=values[:4];mean=values[4:]
            if sequence==1:baseline=mean.copy();assert np.array_equal(p,par);continue
            if sequence==n+2+2*len(checks):assert np.array_equal(p,par) and np.array_equal(mean,baseline);continue
            t=sequence-(n+2);r=checks[t//2];expected=par.copy();expected[1:3]=float(r['shape']),float(r['AV'])
            if t%2==0:
                assert np.array_equal(p,expected);h=mean.copy();assert np.array_equal(h,H[lookup[tuple(expected[1:3])]])
            else:
                expected[0]=float(r['DLMAG']);assert np.array_equal(p,expected)
                pred=h*10**(-.4*(expected[0]-par[0]));error=float(np.max(abs(mean-pred))/np.max(abs(pred)))
                errors.append(error);pointwise.append(float(np.max(abs((mean-pred)/pred))))
                scale_report_errors.append(abs(error-float(r['relative_error'])))
                assert error<2e-8 and scale_report_errors[-1]<1e-12
    result=json.loads((g/'result.json').read_text());assert sequence==n+2+2*len(checks)==result['total_oracle_calls']
    raw_checks[cid]=dict(raw_calls_counted=sequence,cached_mean_vectors=n,competitive_actual_amplitude_checks=len(checks),max_relative_error=max(errors),max_pointwise_relative_error=max(pointwise),max_reported_error_difference=max(scale_report_errors),baseline_replay_exact=True)
    print(json.dumps(dict(CID=cid,raw_calls=sequence,amplitude_checks=len(checks),max_error=max(errors))),flush=True)
zero={cid:o['zero_change_control'] for cid,o in outputs.items() if o['zero_change_control'] is not None}
assert set(zero)=={'DES16E2clk','DES16E2cqq'}
result=dict(status='PASS: complete ten-object saved-artifact verification',all_ten_complete=True,all_numerical_gates_verified=all(o['numerical_gate_pass'] for o in outputs.values()),final_files_verified=len(manifest['files_sha256']),executed_sources_verified=len(manifest['executed_source_sha256']),raw_native_call_count=sum(r['raw_calls_counted'] for r in raw_checks.values()),saved_mean_vectors=sum(r['cached_mean_vectors'] for r in raw_checks.values()),actual_amplitude_checks=sum(r['competitive_actual_amplitude_checks'] for r in raw_checks.values()),max_actual_amplitude_error=max(r['max_relative_error'] for r in raw_checks.values()),raw_checks=raw_checks,zero_change_controls=zero,max_Q_D_amplitude_error=np.max([m['Q_D_amplitude_max_error'] for o in outputs.values() for m in o['metric_errors'].values()],axis=0).tolist(),max_coarse_fine_Q=max(abs(m['coarse_minus_finest_Q']) for o in outputs.values() for m in o['metrics'].values()),max_coarse_fine_D=max(abs(m['coarse_minus_finest_DLMAG']) for o in outputs.values() for m in o['metrics'].values()),cohort_manifest_sha256=sha(C/'cohort-manifest.json'),summary_sha256=sha(C/'cohort-summary.json'),source_sha256=sha(Path(__file__)),interpretation='Conditional fixed-state numerical profile sensitivity; no unique-distance, selected-likelihood, population bias or cosmology correction claim.')
(OUT/'final-certification.json').write_text(json.dumps(result,indent=2)+'\n')
