"""Independent arithmetic and raw-pixel audit of the frozen HST aperture ledger.

Consumes completed results only. Does not fit a noise scale or change selection.
Geometry/source-mask validation remains a separate pre-outcome gate.
"""
from pathlib import Path
import argparse, csv, hashlib, json, os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from astropy.io import fits

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rows(p):
    with p.open() as f:return list(csv.DictReader(f))

def close(a,b,label,relative=2e-7,absolute=1e-10):
    if not np.isfinite(a) or not np.isfinite(b) or abs(a-b)>absolute+relative*max(abs(a),abs(b)):
        raise AssertionError((label,a,b))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('execution',type=Path);ap.add_argument('out',type=Path)
    ap.add_argument('--summary',type=Path,help='Additive recovered primary summary, when original output serialization failed.')
    args=ap.parse_args();p=args.execution;args.out.mkdir(exist_ok=False)
    state=json.loads((p/'stage_a.json').read_text())
    main_path=args.summary or p/'stage_b_primary.json';result=json.loads(main_path.read_text())
    assert sha(p/'operators.npz')==state['operators_sha256']
    assert sha(p/'source_masks.npz')==state['source_masks_sha256']
    with np.load(p/'operators.npz',allow_pickle=False) as f:
        op={k:f[k] for k in ('flat','aperture','annulus','pam')}
    offsets={(r['candidate'],r['exposure']):(r['start'],r['end']) for r in state['operator_rows']}
    candidates={r['id']:r for r in state['rows']}
    exps=state['exposures'];index={r['filename']:i for i,r in enumerate(exps)}
    ledger=rows(p/'signed_exposure_ledger.csv');pairs=rows(p/'signed_pair_ledger.csv')
    common=rows(p/'common_template_ledger.csv');checks=[];read_paths=[];keyed={}
    for e in exps:
        raw=Path(e['path']);assert sha(raw)==e['sha256'];read_paths.append(raw)
        i=index[e['filename']]
        with fits.open(raw,memmap=True,do_not_scale_image_data=True) as f:
            for r in ledger:
                if r['filename']!=e['filename']:continue
                c=int(r['candidate']);lo,hi=offsets[c,i];flat=op['flat'][lo:hi]
                a,b,area=[op[k][lo:hi].astype(np.float64) for k in ('aperture','annulus','pam')]
                y,err=[f[k].data.ravel()[flat].astype(np.float64) for k in ('SCI','ERR')]
                dq=f['DQ'].data.ravel()[flat]
                valid=np.isfinite(y)&np.isfinite(err)&(err>0)&(dq==0)
                assert candidates[c].get(r['branch'])
                if r['branch']=='strict':assert valid.all()
                ca=float(np.sum(a[valid]*area[valid])/np.sum(a*area))
                cb=float(np.sum(b[valid])/np.sum(b))
                close(ca,float(r['aperture_coverage']),'coverage a',absolute=1e-7)
                close(cb,float(r['annulus_coverage']),'coverage b',absolute=1e-7)
                k=e['photflam']/exps[0]['photflam']
                a=a[valid]*area[valid];b=b[valid];y=y[valid];err=err[valid]
                sky=np.sum(b*y)/np.sum(b)
                flux=k*np.sum(a*(y-sky))
                w=k*(a-np.sum(a)*b/np.sum(b));variance=np.sum(np.square(w*err))
                close(flux,float(r['flux']),'raw flux',absolute=1e-9)
                close(variance,float(r['V']),'raw variance')
                close(np.sum(w),0,'constant image cancellation',absolute=1e-8)
                key=(r['branch'],c,e['visit'],e['ordinal']);assert key not in keyed
                keyed[key]=(float(r['flux']),float(r['V']))
                checks.append({'branch':r['branch'],'candidate':c,'filename':e['filename'],
                               'flux_gap':flux-float(r['flux']),'variance_gap':variance-float(r['V'])})
    assert len(checks)==len(ledger)
    for branch,entry in result['branches'].items():
        allowed=[c for c in candidates.values() if c.get(branch)]
        assert len(allowed)==entry['support_n']
        if len(allowed)<30:
            assert not any(r['branch']==branch for r in ledger)
            continue
        assert len([r for r in ledger if r['branch']==branch])==8*len(allowed)
        for visit in ('search','template'):
            pair=[r for r in pairs if r['branch']==branch and r['visit']==visit]
            assert {int(r['candidate']) for r in pair}=={c['id'] for c in allowed}
            for r in pair:
                c=int(r['candidate']);f2,v2=keyed[branch,c,visit,2];f4,v4=keyed[branch,c,visit,4]
                close(float(r['d']),f2-f4,'pair contrast');close(float(r['V']),v2+v4,'pair variance')
            d=np.array([float(r['d']) for r in pair]);v=np.array([float(r['V']) for r in pair]);z=d/np.sqrt(v)
            mu=np.sum(d/v)/np.sum(1/v)
            expected={'n':len(d),'sum_d':d.sum(),'sum_d2':d@d,'sum_V':v.sum(),
                'sum_d2_minus_sum_V':d@d-v.sum(),'sum_d2_over_sum_V':(d@d)/v.sum(),
                'signed_standardized_mean':z.mean(),'uncentered_mean_z2':np.mean(z*z),
                'weighted_intercept':mu,'centered_variance_ratio':np.sum((d-mu)**2/v)/(len(d)-1)}
            for name,value in expected.items():close(float(value),entry['visits'][visit][name],name,relative=1e-10)
    for r in common:
        branch,c=r['branch'],int(r['candidate'])
        x=np.array([keyed[branch,c,visit,j][0] for visit,j in [('search',2),('search',4),('template',2),('template',4)]])
        v=np.array([keyed[branch,c,visit,j][1] for visit,j in [('search',2),('search',4),('template',2),('template',4)]])
        cov=(v[2]+v[3])/4
        close(float(r['common_cov_offdiag']),cov,'shared reference covariance')
        close(float(r['separate_cov_offdiag']),0,'separate reference covariance')
        close(float(r['crossproduct_difference']),.5*(x[0]-x[1])*(x[3]-x[2])+.25*(x[2]-x[3])**2,'shared reference moment')
        close(float(r['same_template_difference_minus_search_pair']),0,'reference cancellation')
    out={'pass':True,'raw_pixel_exposure_rows':len(checks),'pair_rows':len(pairs),'shared_template_rows':len(common),
         'max_flux_gap':max((abs(r['flux_gap']) for r in checks),default=0),
         'max_variance_gap':max((abs(r['variance_gap']) for r in checks),default=0),
         'checks':checks,
         'scope':'Independent float64 raw-pixel, operator, pair, moment and shared-template arithmetic; no certification of physical noise independence, historical reduction or cosmological correction.'}
    (args.out/'result.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    paths=[Path(__file__),p/'stage_a.json',main_path,p/'operators.npz',p/'source_masks.npz',
           p/'signed_exposure_ledger.csv',p/'signed_pair_ledger.csv',p/'common_template_ledger.csv']+read_paths
    (args.out/'manifest.json').write_text(json.dumps({str(q.resolve()):sha(q) for q in paths},indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
