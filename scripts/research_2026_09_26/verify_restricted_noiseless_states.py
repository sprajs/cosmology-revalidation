"""Independent noiseless native state audit, using the root's raw log reader."""
from pathlib import Path
from collections import Counter
import argparse, hashlib, json, os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import numpy as np
from astropy.io import fits
from review_csp_native_states import parse

def science(path):
    return [s for s in path.read_text().splitlines() if s.startswith(('SN:','VARNAMES:'))]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('base',type=Path);ap.add_argument('generation');ap.add_argument('out',type=Path)
    ap.add_argument('--fits-directory',default='fits-restricted')
    ap.add_argument('--nir-fits-directory',help='Optional separate directory for a resumed NIR job; joint jobs remain unchanged.')
    ap.add_argument('--joint-only',action='store_true',help='Audit existing joint outputs without implying the NIR recovery gate has passed.')
    args=ap.parse_args();p=args.base;args.out.mkdir(exist_ok=False)
    source=p/args.generation/'noiseless'
    schema=json.loads((p/'instrumentation/schema-v2.json').read_text())['PROSP_EVENT']
    events=[dict(zip(schema,line.split()[1:])) for line in (source/'native.log').read_text().splitlines() if line.startswith('PROSP_EVENT ')]
    assert len(events)==1
    truth=events[0]; cid=truth['cid']; mu=float(truth['DLMU_true']);peak=float(truth['peak_true'])
    heads=list((source/'output/PTE').glob('*HEAD.FITS*'));phots=list((source/'output/PTE').glob('*PHOT.FITS*'))
    assert len(heads)==len(phots)==1
    with fits.open(heads[0]) as f:head=f[1].data.copy()
    with fits.open(phots[0]) as f:phot=f[1].data.copy()
    assert len(head)==1 and str(head[0]['SNID']).strip()==cid
    selected=phot[int(head[0]['PTROBS_MIN'])-1:int(head[0]['PTROBS_MAX'])]
    groups={};results=[];checks=[];inputs=[source/'native.log',heads[0],phots[0]]
    cases=[('noiseless_joint_active',117)]
    if not args.joint_only:cases.append(('noiseless_NIR_active',6))
    for name,n in cases:
        directory=args.nir_fits_directory if name=='noiseless_NIR_active' and args.nir_fits_directory else args.fits_directory
        work=p/directory/name; log=work/'native.log'; inputs+=[log,work/'fit.FITRES.TEXT']
        records,entries=parse(log);groups[name]=records
        assert len(records)==12 and {r['CID'] for r in records}=={cid}
        assert [r['iteration'] for r in records]==list(range(1,13))
        wanted=Counter((str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR']))
            for x in selected if n==117 or str(x['BAND']).strip() in ('J','H'))
        for r in records:
            actual=Counter((x['band'],x['values'][0],x['values'][4],x['values'][5]) for x in r['rows'])
            assert actual==wanted and r['n']==n
            x=np.array([s['values'] for s in r['rows']]); w=np.array(r['weights']);w=w if r['cov'] else np.diag(w[:,0])
            l=np.linalg.cholesky(w);q=float(np.sum((l.T@(x[:,4]-x[:,2]))**2))
            o=r['objective'];gap=abs(q-(o[0]-o[1]-o[2]));a=float(x[:,2]@w@x[:,4]/(x[:,2]@w@x[:,2])); assert a>0
            checks.append(dict(case=name,iteration=r['iteration'],Q_gap=gap,D=o[3],peak=o[6],
                phase_range=[float(x[:,1].min()),float(x[:,1].max())],fixed_C_D_gap=float(-2.5*np.log10(a))))
            assert gap<=1e-7 and o[4]==1 and o[5]==0
            assert np.isfinite(x).all() and -20<=x[:,1].min()<=x[:,1].max()<=70
            if name=='noiseless_NIR_active':assert o[6]==float(head[0]['PEAKMJD'])
        last=records[-1];o=last['objective'];x=np.array([s['values'] for s in last['rows']])
        residual=float(np.max(np.abs((x[:,4]-x[:,2])/x[:,5])))
        results.append(dict(case=name,truth_D=mu,D=o[3],D_minus_truth=o[3]-mu,
            peak_minus_truth=o[6]-peak,max_standardized_residual=residual,
            mean_recovery_pass=bool(abs(o[3]-mu)<=.001 and abs(o[6]-peak)<=.01 and residual<=.02)))
        if not name.endswith('_original'):
            summaries=[]
            for line in log.read_text().splitlines():
                v=line.split()
                if v and v[0].startswith('PROSP_SUPPORT'):
                    assert v[0]=='PROSP_SUPPORT' and len(v)==11
                    assert v[1]==cid and int(v[2])>0 and int(v[3])==int(v[4])==0
                    summaries.append(v)
            assert len(summaries)==1
    result=dict(exact_instrumentation_states=None,callbacks=len(checks),results=results,checks=checks,
        joint_only=args.joint_only,
        all_mean_recovery_gates_pass=all(r['mean_recovery_pass'] for r in results),
        scope='Restricted-estimator conditional noiseless recovery. Disabled identity is checked separately; no population, historical-execution or cosmological inference.')
    (args.out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs += [Path(__file__),Path(__file__).with_name('review_csp_native_states.py')]
    (args.out/'manifest.json').write_text(json.dumps({str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
