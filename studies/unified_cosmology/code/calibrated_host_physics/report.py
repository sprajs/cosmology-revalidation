#!/usr/bin/env python3
"""Compact per-host results and paired sensitivity comparisons."""
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/unified-cosmology/calibrated-host-physics'
RESULT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source=WORK/'fit-records.json';records=json.loads(source.read_text())
    index={(r['targetid'],r['scenario']):r for r in records}
    ids=sorted({r['targetid'] for r in records})
    def bounds(tid,scenario):
        fit=index[(tid,scenario)]['fit']
        assert fit['status']=='compatible'
        return np.array([fit['lower']['age_Gyr'],fit['upper']['age_Gyr']])
    rows=[]
    for tid in ids:
        row=dict(targetid=tid,z=index[(tid,'primary')]['z'])
        for scenario in ['primary','clock_LCDM','original_free_emission','dust_tau4','mask800']:
            lo,hi=bounds(tid,scenario)
            row.update({scenario+'_lower_Gyr':float(lo),scenario+'_upper_Gyr':float(hi)})
        rows.append(row)
    table=RESULT/'host-age-feasibility.csv'
    with table.open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    primary=np.array([bounds(t,'primary') for t in ids])
    clock=np.array([bounds(t,'clock_LCDM') for t in ids])
    sensitivities={}
    for scenario in sorted({r['scenario'] for r in records}):
        values=np.array([bounds(t,scenario) for t in ids]);delta=values-primary
        sensitivities[scenario]=dict(median_width_Gyr=float(np.median(values[:,1]-values[:,0])),
            maximum_absolute_lower_endpoint_change_Gyr=float(abs(delta[:,0]).max()),
            maximum_absolute_upper_endpoint_change_Gyr=float(abs(delta[:,1]).max()),
            median_lower_endpoint_change_Gyr=float(np.median(delta[:,0])),
            median_upper_endpoint_change_Gyr=float(np.median(delta[:,1])))
    coverage=np.array([index[(t,'primary')]['retained_rest_width_A'] for t in ids])
    threshold_retries=sum(len(r['fit'][side]['retries']) for r in records for side in ['lower','upper'])
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='complete',hosts=len(ids),
        primary=dict(includes_age_le1_and_ge10=int(np.sum((primary[:,0]<=1)&(primary[:,1]>=10))),
            upper_near_library_ceiling_within_001Gyr=int(np.sum(primary[:,1]>=13.49)),
            lower_near_library_floor_within_001Gyr=int(np.sum(primary[:,0]<=.015)),
            interval_width_Gyr_percentiles=np.percentile(primary[:,1]-primary[:,0],[0,50,100]).tolist(),
            lower_age_Gyr_percentiles=np.percentile(primary[:,0],[0,50,100]).tolist()),
        clock=dict(upper_near_imposed_ceiling_within_001Gyr=int(sum(clock[j,1]>=index[(t,'clock_LCDM')]['age_ceiling_Gyr']-.01 for j,t in enumerate(ids))),
            median_width_Gyr=float(np.median(clock[:,1]-clock[:,0]))),
        paired_sensitivity=sensitivities,
        masked_retained_width_A_min=coverage.min(axis=0).tolist(),
        masked_retained_width_A_median=np.median(coverage,axis=0).tolist(),
        solver_rejected_attempts=threshold_retries,
        interpretation='The five compressed observed bands do not identify formed-mass stellar ages across this finite, flexible SFH/differential-dust family. This does not show that all55 true ages span the interval, nor that the full spectra or broader photometry lack additional information. Age-ceiling narrowing is conditional on the imposed cosmological clock.',
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),source,RESULT/'fits.json']},
        output_sha256={str(table.relative_to(ROOT)):sha(table)})
    (RESULT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
