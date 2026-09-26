"""Trace paired forward experiments from generated attempts through flux fits."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np,pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'phase2/literature/simulations';OUT=ROOT/'phase2/hierarchy/forward';OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts/phase2/official'))
from audit_fits import read_fit,hessian
MODELS=['P21','BS21','G10','P21_dmplus010','P21_dmminus010','P21_dmz020','P21_rho000','P21_rho090','P21_noisetrue120']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--require-fits',action='store_true');args=ap.parse_args();reports=[];inputs=[];base=None
    for model in MODELS:
        version='PH2_pilot02_'+model;genpath=BASE/(version+'-generated.csv.gz');gen=pd.read_csv(genpath);inputs.append(genpath)
        manpath=BASE/'manifests'/(version+'.json');man=json.loads(manpath.read_text());assert man['exit_code']==0 and man['success_marker'];inputs.append(manpath)
        for obj in [man['input'],man['executable'],man['log'],*man['outputs']]:assert hashlib.sha256((ROOT/obj['path']).read_bytes()).hexdigest()==obj['sha256'],obj['path']
        if base is None:base=gen
        eq={c:bool(np.array_equal(base[c].to_numpy(),gen[c].to_numpy())) for c in ['GENZ','LIBID','GALID','PEAKMJD','RA','DEC']};assert all(eq.values())
        write=(gen.SIM_EFFMASK==5)&(gen.CUTMASK==4095);selected=gen.loc[write].copy();assert not selected.CID.duplicated().any()
        heads=fits.getdata(BASE/'outputs'/version/(version+'_HEAD.FITS'),1);assert set(np.asarray(heads['SNID']).astype(int))==set(selected.CID)
        stat=gen[['generated_attempt_index','CID','GENZ','LIBID','GALID','PEAKMJD','LOGMASS_TRUE','SALT2x1','SALT2c','AV','RV','MU']].copy();stat['written']=write;stat['fit_pass']=False;stat['basic_quality_pass']=False
        label='forward_'+model.lower();runrec=ROOT/f'phase2/official/diagnostics/forward-fit-{label}.json';log=ROOT/f'phase2/official/diagnostics/snana-{label.replace("_","-")}.log';fitpath=ROOT/f'phase2/official/results/snana_{label}.FITRES.TEXT'
        # Nominal was originally launched directly and has no helper manifest.
        ready=log.exists() and 'ENDING PROGRAM GRACEFULLY.' in log.read_text()
        record={'model':model,'version':version,'n_generated':len(gen),'n_written':int(write.sum()),'design_stream_pairing':eq,'fit_complete':ready,'n_duplicate_generated_cids':int(gen.CID.duplicated().sum())}
        if not ready:
            if args.require_fits:raise RuntimeError('Incomplete fitting arm: '+model)
        else:
            inputs.extend([fitpath,log]);f=read_fit(fitpath);f['CID']=f.CID.astype(int);raw=hessian(log)
            assert set(f.CID.astype(str))==set(raw);assert set(f.CID)<=set(selected.CID)
            quality=(f.x1.abs()<3)&(f.c.abs()<.3)&(f.x1ERR<1)&(f.PKMJDERR<2)&(f.cERR<1.5)&(f.FITPROB>.001)&(f.zHD>.025)&(f.zHD<1.2)
            f['basic_quality_pass']=quality
            pars=np.array([raw[str(cid)][0] for cid in f.CID]);cov=np.array([raw[str(cid)][1] for cid in f.CID]);jac=np.ones((len(f),4));jac[:,0]=-2.5/np.log(10)/pars[:,0];mcov=cov*jac[:,:,None]*jac[:,None,:];mpars=pars.copy();mpars[:,0]=-2.5*np.log10(pars[:,0])+10.635
            f['mB_double']=mpars[:,0];f['x1_double']=mpars[:,1];f['c_double']=mpars[:,2];f['t0_double']=mpars[:,3]
            f=f.merge(selected[['CID','generated_attempt_index','MU','LOGMASS_TRUE','SALT2x1','SALT2c','AV','RV']],on='CID',how='left',validate='one_to_one');assert f.generated_attempt_index.notna().all()
            stat['fit_pass']=stat.generated_attempt_index.isin(f.generated_attempt_index);stat['basic_quality_pass']=stat.generated_attempt_index.isin(f.loc[f.basic_quality_pass,'generated_attempt_index'])
            record.update(n_fit=len(f),n_basic_quality=int(quality.sum()),n_nonpositive_hessians=int((np.linalg.eigvalsh(mcov).min(axis=1)<=0).sum()))
            f.to_csv(OUT/(model+'-fitted.csv.gz'),index=False,compression={'method':'gzip','mtime':0});np.savez_compressed(OUT/(model+'-covariance.npz'),generated_attempt_index=f.generated_attempt_index.to_numpy(dtype=int),CID=f.CID.to_numpy(dtype=int),mag_parameters=mpars,mag_covariance=mcov)
        # Early failures reuse CID. Never join their status by CID alone.
        stat.to_csv(OUT/(model+'-attempts.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
        stat['zbin']=pd.cut(stat.GENZ,[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,2.],include_lowest=True)
        byz=stat.groupby('zbin',observed=True).agg(generated=('CID','size'),written=('written','sum'),fitted=('fit_pass','sum'),basic_quality=('basic_quality_pass','sum')).reset_index();byz.to_csv(OUT/(model+'-counts-by-z.csv'),index=False)
        reports.append(record)
    (OUT/'summary.json').write_text(json.dumps({'models':reports,'qualification':'Common design stream verified. Population/noise streams need not coincide. Generated-attempt index owns denominators. Basic fitted quality is not the BBC membership/classifier/systematic intersection. Counts alone cannot select a physical model because rates and parent population normalization are unknown.'},indent=2)+'\n')
    inputs.extend([Path(__file__),ROOT/'scripts/phase2/official/audit_fits.py']);sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    (OUT/'manifest.json').write_text(json.dumps({'inputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in inputs},'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in OUT.iterdir() if f.name!='manifest.json'}},indent=2)+'\n');print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
