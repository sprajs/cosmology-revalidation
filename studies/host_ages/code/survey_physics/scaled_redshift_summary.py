#!/usr/bin/env python3
"""Exploratory paired redshift-shape summaries, preserving joint dependence."""
import argparse
import json
from pathlib import Path
import numpy as np
from common import RESULTS, HERE, ROOT, sha
from scaled_bbc import SCALE

LABELS = ['frozen_nominal', 'age_nominal', 'age_literal', 'age_retained']
BINS = [0, 1, 2, 3]


def distribution(values):
    x = np.asarray(values, dtype=float)
    return {'supported_replicates': len(x),
            'sd_mag': float(x.std(ddof=1)) if len(x) > 1 else None,
            'percentile_95_mag': np.quantile(x,[.025,.975]).tolist() if len(x)>1 else None}


def vector(record, label):
    d = record.get('contrasts', {}).get(label+'_all_cases', {})
    if d.get('slope_mag_per_Gyr') is None:
        return None
    b = d['centered_redshift_means']
    return [b[i]['mean_mag'] for i in BINS]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--shards', type=int, default=16)
    args = parser.parse_args()
    records = {}
    dependencies = {}
    for variant in ['common-fixed', 'common-highmass-gauge']:
        point_path=RESULTS/f'scaled-response-k{args.shards:03d}-{variant}.json'
        point=json.loads(point_path.read_text())
        dependencies[str(point_path.relative_to(ROOT))]=sha(point_path)
        point_hilo={k: vector(point,k)[3]-vector(point,k)[0] for k in LABELS}
        modes={}
        for mode in ['joint','training_only']:
            summary_path=RESULTS/f'scaled-bootstrap-{mode}-k{args.shards:03d}-{variant}.json'
            summary=json.loads(summary_path.read_text())
            dependencies[str(summary_path.relative_to(ROOT))]=sha(summary_path)
            reps=[]
            for index,digest in summary['replicate_result_hashes'].items():
                p=SCALE/'bootstrap'/f'k{args.shards:03d}-{variant}'/mode/f'r{int(index):03d}'/'result.json'
                assert sha(p)==digest
                reps.append(json.loads(p.read_text()))
            highlow={k:[] for k in LABELS}
            changes={'retained_minus_frozen':[], 'retained_minus_nominal_refit':[]}
            complete=[]
            for r in reps:
                vec={k:vector(r,k) for k in LABELS}
                hl={k:v[3]-v[0] for k,v in vec.items() if v is not None and v[0] is not None and v[3] is not None}
                for k,v in hl.items():highlow[k].append(v)
                for key,before in [('retained_minus_frozen','frozen_nominal'),('retained_minus_nominal_refit','age_nominal')]:
                    if 'age_retained' in hl and before in hl:
                        changes[key].append(hl['age_retained']-hl[before])
                if all(v is not None and all(x is not None for x in v) for v in vec.values()):
                    complete.append([x for k in LABELS for x in vec[k]])
            modes[mode]={
                'replicates':len(reps), 'replicate_status_counts':summary['replicate_status_counts'],
                'high_minus_low':{k:distribution(v) for k,v in highlow.items()},
                'paired_changes':{k:distribution(v) for k,v in changes.items()},
                'joint_bin_covariance':{'ordering':[f'{k}:bin{i}' for k in LABELS for i in BINS],
                    'complete_replicates':len(complete),
                    'covariance_mag2':np.cov(np.array(complete),rowvar=False,ddof=1).tolist() if len(complete)>1 else None},
            }
        records[variant]={'point_high_minus_low_mag':point_hilo,
                         'point_retained_minus_frozen_mag':point_hilo['age_retained']-point_hilo['frozen_nominal'],
                         'point_retained_minus_nominal_refit_mag':point_hilo['age_retained']-point_hilo['age_nominal'],
                         'bootstrap':modes}
    out={'status':'computed', 'exploratory':True,
         'interpretation':'Post-pattern high-minus-low redshift contrast; globally centered paired distances, not within-bin age-slope survival. Conditional supported-replicate uncertainty, not coverage-validated or survey-systematic.',
         'design_sha256':sha(HERE/'scaled-redshift-summary-design.json'), 'code_sha256':sha(__file__),
         'dependencies_sha256':dependencies, 'variants':records}
    (RESULTS/f'scaled-redshift-summary-k{args.shards:03d}.json').write_text(json.dumps(out,indent=2)+'\n')
    for k,v in records.items():print(k,json.dumps({'points':v['point_high_minus_low_mag'],'joint':v['bootstrap']['joint']['high_minus_low'],'paired_change':v['bootstrap']['joint']['paired_changes']},indent=2))


if __name__=='__main__':main()
