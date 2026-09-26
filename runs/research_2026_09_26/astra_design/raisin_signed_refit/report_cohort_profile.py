"""Render all-ten conditional profile results without population-bias aggregation."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;C=O/'fixed-c-profile/cohort10'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
  protocol=json.loads((C/'protocol.json').read_text());index=json.loads((C/'active-attempts.json').read_text())['attempts'];rr=[]
  for cid in protocol['cohort']:
    d=C/index[cid];r=json.loads((d/'result.json').read_text());g=json.loads((d/'branch-level-geometry.json').read_text());row=dict(CID=cid,attempt=index[cid],rows_A=r['rows_A'],rows_B=r['rows_B'],added_negative=r['rows_B']-r['rows_A'],numerical_gate_pass=r['numerical_gate_pass'],computational_completion=r['computational_completion'],calls=r['total_oracle_calls'])
    for anchor in ['B','A']:
      row['delta_'+anchor+'_anchor']=r['delta_DLMAG_B_minus_A'][anchor]
    row.update(matched_shape_delta=g['paired']['B']['matched_mode_delta_DLMAG'],matched_shape_cost=g['paired']['B']['matched_B_mode_nearest_A_shape']['delta_Q'])
    for arm in ['A','B']:
      label='Banchor_'+arm;row[arm+'_minimum']=r['metrics'][label]['finest']['DLMAG']
      for lev in ['1.0','4.0','9.0']:
        ls=g['metrics'][label]['level_sets'][lev];stem=arm+'_level'+str(int(float(lev)))
        row[stem+'_low'],row[stem+'_high']=ls['amplitude_inclusive_distance_range'];row[stem+'_shape_edge']=ls['shape_box_touched'];row[stem+'_AV_edge']=ls['AV_box_touched']
    row['gate_flags']={key:{k:v for k,v in x.items() if k in ['coarse_fine_tolerance_pass','AV_edge_within9','shape_edge_within1']} for key,x in r['metrics'].items()}
    rr.append(row)
  with(C/'cohort-summary.csv').open('w') as f:
    keys=[k for k in rr[0] if k!='gate_flags'];w=csv.DictWriter(f,keys,extrasaction='ignore');w.writeheader();w.writerows(rr)
  (C/'cohort-summary.json').write_text(json.dumps(dict(objects=rr,source_sha256=sha(__file__),no_population_bias_estimate=True),indent=2)+'\n')
  fig,ax=plt.subplots(figsize=(10,7.2))
  for i,row in enumerate(rr):
    for arm,col,off in [('A','#416788',-.14),('B','#d06b38',.14)]:
      lo=row[arm+'_level1_low']-row['A_minimum'];hi=row[arm+'_level1_high']-row['A_minimum'];v=row[arm+'_minimum']-row['A_minimum']
      ax.plot([lo,hi],[i+off,i+off],color=col,lw=2);ax.plot(v,i+off,'o',color=col,ms=5)
    if not row['numerical_gate_pass']:ax.text(.98,i,'gate limited',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8,color='#993333')
  ax.set_yticks(range(10),[r['CID']+'  (+'+str(r['added_negative'])+' rows)' for r in rr]);ax.invert_yaxis();ax.axvline(0,color='gray',ls=':',lw=1)
  ax.set_xlabel('DLMAG relative to each positive-only minimum [mag]')
  ax.set_title('Fixed-covariance distance geometry after restoring measured negative epochs')
  ax.plot([],[],color='#416788',marker='o',label='A: author positive epochs');ax.plot([],[],color='#d06b38',marker='o',label='B: full signed epochs');ax.legend(loc='best',fontsize=9)
  ax.grid(axis='x',alpha=.2);fig.text(.02,.02,'Lines: amplitude-inclusive ΔQ ≤ 1 envelopes on the saved restricted grid. They are not confidence intervals.\nPoints: global minima of each arm’s own conditional quadratic; no evidence comparison between different row sets.',fontsize=9)
  fig.tight_layout(rect=[0,.07,1,1]);fig.savefig(C/'cohort-profile-geometry.png',dpi=170);fig.savefig(C/'cohort-profile-geometry.pdf');plt.close(fig)
  lines=['# Signed-epoch cohort profiles: conditional geometry, not a distance correction','',
    'The source-defined ten-object DES16 cohort has been evaluated with a common fixed-covariance rule and an exhaustive restricted native SNooPy shape/colour search. Restoring measured negative epochs can change a minimum, but broad or competing distance branches prevent treating the resulting point differences as measured bias corrections.',
    '',f"All {sum(x['computational_completion'] for x in rr)}/10 objects completed; {sum(x['numerical_gate_pass'] for x in rr)}/10 pass every declared numerical/support gate. Every source-defined object remains in the ledger. The two objects with identical accepted A/B rows provide exact-zero controls.",'',
    '![All-ten conditional distance geometry](cohort-profile-geometry.png)','',
    'The plotted intervals are finite-grid descriptive ΔQ≤1 sets, including analytic amplitude freedom. They are not confidence or posterior intervals. Each quadratic is centered on its own arm minimum; different row sets are never compared as a likelihood ratio.','',
    '| Object | A / B rows | Gate | ΔD, B anchor | ΔD, A anchor | Nearest-shape B branch ΔD (cost ΔQ) | B ΔQ≤1 range |','|---|---:|---|---:|---:|---:|---:|']
  for x in rr:
    lines.append(f"| {x['CID']} | {x['rows_A']} / {x['rows_B']} | {'pass' if x['numerical_gate_pass'] else 'limited'} | {x['delta_B_anchor']:+.6f} | {x['delta_A_anchor']:+.6f} | {x['matched_shape_delta']:+.6f} ({x['matched_shape_cost']:.3f}) | {x['B_level1_low']:.5f}–{x['B_level1_high']:.5f} |")
  lines += ['',
    'The branch-match rule was fixed in advance: choose the fine-grid B local mode/end point nearest in shape to A’s global minimum, irrespective of shift sign. It is a diagnostic branch correspondence, not a second preferred estimator. Gate-limited rows and all parameter-box flags remain in the machine-readable ledger. ΔQ=4 and9 envelopes and boundary flags are retained per object.',
    '', '## What is controlled', '',
    'A contains actual author-positive measurements; B adds the actual source-measured negative epochs, using the same author flux/error convention on shared rows. The same released header peak is fixed in both arms. This timing estimate comes from light-curve fitting; it is not an independent noiseless clock. The primary metric is the native B reference covariance, with A using the marginal covariance submatrix, independently inverted. The sensitivity uses full-B covariance evaluated at the pre-existing A/start1 coordinates. The reference covariance is deterministic; this does not claim its original iterative optimizer was globally optimal.',
    '',
    'Means come from the exact native observer-frame USRFUN with the released SNooPy grid/KCOR, host RV=1.518 and controlled modern MW law99. Shape covers every tabulated cell (approximately .7–1.3); empirical AV covers −1 to2. Negative AV is not passive dust. Positive amplitude is profiled analytically in the verified flat distance domain10–60; the native mean is replayed at competitive profiled distances. Coarse/fine grids, every local AV bracket, exact float32 algorithm boundaries, continuous refinements, clipping/amplitude guards and both covariance anchors are retained. The minimum gate requires AV edges beyond ΔQ9 and shape edges beyond ΔQ1; wider level sets can still touch shape limits. Coarse/fine tolerance tests certify minima, not complete convergence of every level-set endpoint.',
    '',
    'The oracle has object-specific parameters and accepted epochs. Initial native exports are checked for exact coordinates, rows, data/errors, mean and precision before stream queries for the remaining eight objects. The two benchmarks receive the same check retrospectively. A separate all-case certification requires exact reference↔fixed-reference↔final-stream band/MJD/data/error identity. Reference and fixed-reference covariances are not required to agree because their initialization differs by design.',
    '', '## Execution provenance', '',
    'The first generic C3cmy attempt omitted explicit object parameters and inherited the pilot defaults. Its accepted-row count mismatch stopped the run before profile scores. The failed attempt and original code remain. Version2 passes explicit parameters; version3 adds the initial native-state gate. All scientific domains, covariance rules, means and thresholds remain unchanged. The active-attempts index selects immutable result versions. No object is removed and no domain is expanded after seeing outcomes.',
    '', '## Scientific limits and next experiment', '',
    'This is a newly declared fixed-C conditional quadratic, not the historical iterative fitting likelihood and not the positive-selection-conditioned likelihood. Its normalizations are constant within an arm; comparing minima across different row vectors would not be a valid evidence test. Parameter-dependent covariance, timing uncertainty, SNooPy population support and cohort selection remain outside this calculation.',
    '',
    'DES16C1cim illustrates the central limitation: its large global-minimum shift switches between weakly separated shape/colour/distance branches. The signed high-shape branch costs only about0.116 quadratic units; adding negative rows at fixed A shape/AV moves distance in the opposite direction. A stable numerical minimum therefore cannot substitute for a proper distance likelihood retaining these branches.',
    '',
    'The next bias-recovery experiment must generate paired signed noise on a declared supported signal/cadence and separately vary frozen positive cadence versus realization-dependent sign cuts. It must retain failed/edge cases and use a numerically checked estimator or proper-prior likelihood. These ten conditional point differences are not a population bias estimate and do not justify a cosmology correction.',
    '', '## Artifacts', '',
    '`protocol.json`, both implementation amendments, `active-attempts.json`, `geometry-protocol.json`, `geometry-path-amendment.json`, and `metric-oracle-identity-protocol.json` preserve the declared method and timing. `cohort-summary.{json,csv}`, `branch-summary.json`, every active case’s `branch-level-geometry.json`, `result.json`, full `native-profiles.npz`, raw stream log and covariance exports preserve all branches. `cohort-manifest.json` provides final hashes. Independent review is recorded separately under `runs/research_2026_09_26/raisin_profile_solver_review/cohort10-review/`; consult its final verdict before calling this cohort independently verified.']
  report='\n'.join(lines)+'\n'
  resource=json.loads((C/'resource-numerical-summary.json').read_text())
  report=report.replace('It is a diagnostic branch correspondence, not a second preferred estimator.', 'It is a diagnostic branch correspondence, not a second preferred estimator. Because the matched branch uses a fine-grid point whereas A uses its refined minimum, even identical-row controls have tiny nonzero matched-branch offsets (up to .000367mag); their actual global A→B differences remain exactly zero.')
  report=report.replace('The failed attempt and original code remain.', 'The failed attempt and original code remain. Total active native/profile runtime, including this failed attempt, was '+f"{resource['active_runtime_seconds_including_failed_attempt']:.2f}"+' seconds, below the900-second cap.')
  report=report.replace('No object is removed and no domain is expanded after seeing outcomes.', 'No object is removed and no domain is expanded after seeing outcomes. Per-attempt engine manifests were written before the final JSON was printed to runner.log. The exact appended suffix reproduces every original log hash; original manifests are preserved and the final cohort manifest hashes the closed logs explicitly.')
  (C/'report.md').write_text(report)

if __name__=='__main__':main()
