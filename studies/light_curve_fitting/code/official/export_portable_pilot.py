#!/usr/bin/env python3
"""Package calibrated observations, fixed-parameter SNANA predictions and source contract.

This is a numerical comparison fixture, not an independent validation of SALT.
"""
from pathlib import Path
import gzip, hashlib, json, shutil
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'phase2/official'; DEST=OUT/'portable_pilot'
RELEASE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3'
SRC=OUT/'build/SNANA-current/src'

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',dtype={'CID':str})
    return d.drop(columns='VARNAMES:')
def main():
    DEST.mkdir(exist_ok=True)
    ids=(OUT/'inputs/pilot12.cid').read_text().split()
    # Input list has plain numeric CIDs; reject accidental non-ID tokens.
    ids=[x for x in ids if x.isdecimal()]
    headpath=RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    photpath=RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_PHOT.FITS.gz'
    with fits.open(headpath) as h, fits.open(photpath) as p:
        hd=h[1].data; ph=p[1].data
        hs=[]; ps=[]
        for row in hd:
            cid=str(row['SNID']).strip()
            if cid not in ids: continue
            rec={name:(row[name].item() if hasattr(row[name],'item') else row[name]) for name in hd.names}
            rec={k:(v.decode().strip() if isinstance(v,bytes) else v) for k,v in rec.items()}
            hs.append(rec)
            for j in range(int(row['PTROBS_MIN'])-1,int(row['PTROBS_MAX'])):
                r=ph[j]; d={name:(r[name].item() if hasattr(r[name],'item') else r[name]) for name in ph.names}
                d={k:(v.decode().strip() if isinstance(v,bytes) else v) for k,v in d.items()}
                d.update(CID=cid,source_phot_row_one_based=j+1);ps.append(d)
    pd.DataFrame(hs).to_csv(DEST/'heads.csv',index=False,float_format='%.17g')
    pd.DataFrame(ps).to_csv(DEST/'calibrated_observations.csv',index=False,float_format='%.17g')
    f=table(OUT/'results/snana_mask32_pilot.FITRES.TEXT')
    f.to_csv(DEST/'snana_parameters.csv',index=False,float_format='%.17g')
    from audit_fits import hessian
    precise=hessian(OUT/'diagnostics/snana-mask32-pilot.log')
    (DEST/'double_parameters_and_hessian.json').write_text(json.dumps({cid:{'parameters_x0_x1_c_t0':v[0].tolist(),'hessian_covariance':v[1].tolist()} for cid,v in precise.items()},indent=2)+'\n')
    cols=['CID','MJD','observer_phase','FLUXCAL','FLUXCALERR','DATA_MODEL','BAND','CHI2','IGNOREME']
    lc=pd.read_csv(OUT/'results/snana_mask32_pilot.LCPLOT.TEXT',sep=r'\s+',names=cols,dtype={'CID':str})
    lc.to_csv(DEST/'snana_predictions_and_mask.csv',index=False,float_format='%.17g')
    with gzip.open(RELEASE/'0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz','rt') as fh:
        public=pd.read_csv(fh,sep=r'\s+',comment='#',names=cols,dtype={'CID':str})
    public[public.CID.isin(ids)].to_csv(DEST/'published_predictions_and_mask.csv',index=False,float_format='%.17g')
    assets=DEST/'assets'; assets.mkdir(exist_ok=True)
    model=RELEASE/'2_LCFIT_MODEL/SALT3.DES5YR'
    for name in ['SALT3.INFO','salt3_template_0.dat.gz','salt3_template_1.dat.gz','salt3_color_correction.dat.gz','salt3_color_dispersion.dat.gz','salt3_lc_model_variance_0.dat.gz','salt3_lc_model_variance_1.dat.gz','salt3_lc_model_covariance_01.dat.gz','salt3_lc_variance_0.dat.gz','salt3_lc_variance_1.dat.gz','salt3_lc_covariance_01.dat.gz']:
        shutil.copy2(model/name,assets/name)
    calibration=OUT/'inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits.gz'
    shutil.copy2(calibration,assets/calibration.name)
    with fits.open(calibration) as cal:
        response=cal['FilterTrans'].data
        pd.DataFrame({n:response[n] for n in response.names}).to_csv(assets/'filter_transmission.csv',index=False,float_format='%.17g')
    cfg=OUT/'inputs/snana_mask32_pilot.nml';shutil.copy2(cfg,DEST/'runtime.nml')
    for name in ['genmag_SALT2.c','genmag_SALT2.h','genmag_SEDtools.c','genmag_extrap.c','MWgaldust.c','snlc_fit.F90','sntools_output_text.c']:
        shutil.copy2(SRC/name,assets/name)
    contract={
      'schema_version':1,'purpose':'Independent fixed-parameter flux/likelihood crosscheck; matching SNANA does not validate astrophysical assumptions.',
      'sample':{'n':len(ids),'CID':ids,'selection':'Frozen pilot: evenly spaced redshift ranks in original DES1635; no outcome selection.'},
      'reference_model':'SALT3.DES5YR; original release tag1.3',
      'numerical_precision':{'observations':'Original FITS binary precision, serialized17g; source row IDs retained.','parameters':'double_parameters_and_hessian.json contains true REAL8 fit parameters; snana_parameters.csv ordinary SNANA FITRES values can be float32.','predictions':'LCPLOT FLUXCAL and FLUXCALERR are printed 5 significant digits; MJD is float32 printed4dec for local and3dec for public. This fixture cannot test sub-1e-4 relative flux accuracy; begin within1e-3 relative for nonzero flux, then refine output.','DATA_MODEL':'1 accepted data; -1 rejected data; 0 model curve. Model curve fixed at the fit parameters. No independent refit yet.'},
      'units':{'flux':'27.5-2.5log10(FLUXCAL) magnitude','wavelength':'Angstrom','time':'MJD days observer frame','phase':'(MJD-t0)/(1+zHEL)','amplitude':'mx=-2.5log10(x0); mB=mx+10.635; internal SALT MAG_OFFSET=.27 must not be counted twice'},
      'known_limits':['Runtime NML absolute paths need relocation; all needed model and filter assets bundled.','SFD full sky maps are not bundled; exact effective MWEBV values from fit table must be used for fixed-parameter check.','LCPLOT model grid is an output comparison, not an exact per-epoch double-precision forward-model API.','Conditioning on published accepted mask is a separate arm; this fixture preserves both public and own mask32 masks.','Published model parameters are rounded; own true double parameters supplied.','No independent sncosmo prediction has yet been claimed.'],
      'source_files':{str(p.relative_to(ROOT)):digest(p) for p in [headpath,photpath,cfg,calibration,OUT/'build/SNANA-current/bin/snlc_fit.exe']},
      'file_sha256':{str(p.relative_to(DEST)):digest(p) for p in sorted(DEST.rglob('*')) if p.is_file() and p.name!='contract.json'}
    }
    (DEST/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    print(json.dumps({'objects':len(hs),'raw_observations':len(ps),'local_lcplot_rows':len(lc),'public_pilot_rows':int(public.CID.isin(ids).sum()),'directory':str(DEST)},indent=2))

if __name__=='__main__':main()
