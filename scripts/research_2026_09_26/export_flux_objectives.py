"""Read-only export of output-only SNANA audit logs to a new destination."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser();p.add_argument('--log',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--expected-cids',type=Path)
    args=p.parse_args();text=args.log.read_text()
    assert 'ENDING PROGRAM GRACEFULLY.' in text
    if args.output.exists():raise RuntimeError('Preserve existing export destination')
    args.output.mkdir(parents=True)
    flux={};inverse={};context={};setup={};covset={}
    for line in text.splitlines():
        q=line.split()
        if not q:continue
        if q[0]=='PHASE2_FLUX:':flux.setdefault(q[1],{})[int(q[2])]=[q[3]]+[float(x) for x in q[4:]]
        elif q[0]=='PHASE2_COVINV:':inverse.setdefault(q[1],{})[int(q[2])]=[float(x) for x in q[3:]]
        elif q[0]=='PHASE2_OBJECTIVE:':context[q[1]]=(int(q[2]),[float(x) for x in q[3:]])
        elif q[0]=='PHASE2_SETUP:':setup[q[1]]=[float(x) for x in q[2:]]
        elif q[0]=='PHASE2_COVSET:':covset.setdefault(q[1],{})[int(q[2])]=[float(x) for x in q[3:]]
    if args.expected_cids:
        expected=set(args.expected_cids.read_text().split())
        assert set(context)==expected,(sorted(expected-set(context)),sorted(set(context)-expected))
    checks={};rows=[]
    for cid,(n,v) in context.items():
        f=[flux[cid][j] for j in range(1,n+1)];a=np.array([r[1:] for r in f])
        ci=np.array([inverse[cid][j] for j in range(1,n+1)])
        assert ci.shape==(n,n)
        C=np.linalg.inv(ci);res=a[:,4]-a[:,2];chi=float(res@ci@res)
        error=chi+v[1]-v[0];assert abs(error)<1e-7,(cid,error)
        assert np.linalg.eigvalsh(C).min()>0
        kw=dict(inverse_frozen_flux_covariance=ci,frozen_flux_covariance=C,
          parameters_x0_x1_c_t0=np.array(v[2:6]),chi2=np.array(v[0]),prior_chi2=np.array(v[1]),
          initial_search_peakmjd=np.array(v[6]),MJD=a[:,0],rest_phase=a[:,1],model_flux=a[:,2],
          model_magerr=a[:,3],data_flux=a[:,4],data_fluxerr=a[:,5],zHEL=a[:,6],MWEBV=a[:,7],
          band=np.array([r[0] for r in f]))
        if cid in setup:kw['setup']=np.array(setup[cid])
        if cid in covset:kw['covariance_components']=np.array([covset[cid][j] for j in range(1,n+1)])
        np.savez_compressed(args.output/f'objective_{cid}.npz',**kw)
        checks[cid]={'epochs':n,'objective_residual':error,'chi2':v[0],
          'prior_chi2':v[1],'min_covariance_eigenvalue':float(np.linalg.eigvalsh(C).min())}
        rows += [[cid,j,*r] for j,r in enumerate(f,1)]
    cols=['CID','fit_row_one_based','BAND','MJD','rest_phase','model_flux','model_magerr','data_flux','data_fluxerr','zHEL','MWEBV']
    pd.DataFrame(rows,columns=cols).to_csv(args.output/'epochs.csv',index=False,float_format='%.17g')
    sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
    report={'source_log':str(args.log.resolve()),'log_sha256':sha(args.log),'exporter_sha256':sha(__file__),
      'count':len(checks),'epochs':sum(v['epochs'] for v in checks.values()),'checks':checks,
      'scope':'Actual last accepted iteration dump of pinned SNANA output-only build; not historical runtime identification.',
      'outputs_sha256':{p.name:sha(p) for p in args.output.iterdir() if p.is_file()}}
    (args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'count':report['count'],'epochs':report['epochs'],'max_objective_error':max(abs(v['objective_residual']) for v in checks.values())},indent=2))

if __name__=='__main__':main()
