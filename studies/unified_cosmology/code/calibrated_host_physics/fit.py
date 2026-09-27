#!/usr/bin/env python3
"""Execute all declared conditional stellar-mixture scenarios on observed hosts."""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone
import time
import numpy as np
from model import CODE,ROOT,WORK,INPUT,Library,problem,solve

RESULT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows):
    answer={}
    for scenario in sorted({r['scenario'] for r in rows}):
        selected=[r for r in rows if r['scenario']==scenario]
        compatible=[r for r in selected if r['fit']['status']=='compatible']
        lower=[r['fit']['lower']['age_Gyr'] for r in compatible]
        upper=[r['fit']['upper']['age_Gyr'] for r in compatible]
        width=np.array(upper)-np.array(lower)
        answer[scenario]=dict(hosts=len(selected),compatible=len(compatible),
            incompatible=sum(r['fit']['status']=='incompatible' for r in selected),
            numerical_failures=sum(r['fit']['status']=='numerical_failure' for r in selected),
            best_chi2_percentiles=np.percentile([r['fit']['best_chi2'] for r in selected],[0,50,100]).tolist(),
            lower_age_percentiles_Gyr=np.percentile(lower,[0,50,100]).tolist() if lower else None,
            upper_age_percentiles_Gyr=np.percentile(upper,[0,50,100]).tolist() if upper else None,
            width_percentiles_Gyr=np.percentile(width,[0,50,100]).tolist() if lower else None,
            includes_young_and_old=sum(l<=1. and u>=10. for l,u in zip(lower,upper)))
    return answer


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--refined',action='store_true')
    parser.add_argument('--scenarios',nargs='+')
    parser.add_argument('--targetids',nargs='+')
    args=parser.parse_args()
    plan=json.loads((CODE/'implementation-design.json').read_text())
    allids=sorted(p.stem for p in (INPUT/'spectra').glob('*.npz'))
    ids=args.targetids or (plan['refined_targetids'] if args.refined else allids)
    scenarios=args.scenarios or (['primary','original_free_emission'] if args.refined else plan['scenarios'])
    assert set(ids)<=set(allids) and set(scenarios)<=set(plan['scenarios'])
    lib=Library(refined=args.refined)
    suffix='-refined' if args.refined else ''
    if args.targetids or args.scenarios:suffix+='-subset'
    destination=WORK/('fits'+suffix);destination.mkdir(exist_ok=True,parents=True)
    dependencies=[CODE/'model.py',Path(__file__),CODE/'design.json',CODE/'implementation-design.json',
        RESULT/('stellar-library'+('-refined' if args.refined else '')+'.json'),lib.path]
    dependencies+=[INPUT/'spectra'/f'{tid}{ext}' for tid in ids for ext in ['.npz','.json']]
    dependency_hashes={str(p.relative_to(ROOT)):sha(p) for p in dependencies}
    started=time.monotonic();rows=[]
    for tid in ids:
        for scenario in scenarios:
            matrix,y,cov,ages,mass,parameters,info=problem(tid,lib,scenario)
            fit=solve(matrix,y,cov,ages,mass)
            case=destination/f'{tid}-{scenario}.npz'
            np.savez_compressed(case,matrix=matrix,y=y,cov=cov,ages=ages,mass=mass,
                best_coefficients=np.array(fit.pop('best_coefficients')),
                lower_coefficients=np.array(fit['lower'].pop('coefficients')) if fit['lower'] and 'coefficients' in fit['lower'] else [],
                upper_coefficients=np.array(fit['upper'].pop('coefficients')) if fit['upper'] and 'coefficients' in fit['upper'] else [])
            row=dict(info,fit=fit,fixture=str(case.relative_to(ROOT)),fixture_sha256=sha(case))
            rows.append(row)
        print(tid,len(rows),'cases complete',flush=True)
    output=WORK/('fit-records'+suffix+'.json');output.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    # Verify inputs stayed unchanged throughout the calculation.
    assert all(sha(ROOT/path)==digest for path,digest in dependency_hashes.items())
    summary=dict(created_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,
        status='complete' if all(r['fit']['status']!='numerical_failure' for r in rows) else 'complete_with_numerical_failures',
        distinct_hosts=len(ids),cases=len(rows),refined=args.refined,scenarios=summarize(rows),
        interpretation='Conditional fibre formed-mass age feasibility over a finite stellar/dust/response library. Not a Bayesian age posterior, empirical calibration error bound, progenitor age, supernova luminosity coefficient or cosmological constraint.',
        dependencies_sha256=dependency_hashes,output_sha256={str(output.relative_to(ROOT)):sha(output)},
        fixture_sha256={r['fixture']:r['fixture_sha256'] for r in rows})
    (RESULT/('fits'+suffix+'.json')).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if not k.endswith('sha256')},indent=2))


if __name__=='__main__':main()
