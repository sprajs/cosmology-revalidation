"""Map-specified observer-frame foreground alternative through SALT response.

Uses existing fixed design, covariance and derivatives. This is neither a new
detector reduction nor selection/BBC/training closure, and CSFD is not assumed true.
"""
from pathlib import Path
import json
import hashlib
from datetime import datetime,timezone
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'runs/research_2026_09_26/deps'))
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy_healpix import HEALPix
import astropy_healpix
from scipy.linalg import solve_triangular

SOURCE=Path(__file__).read_bytes()
PLAN=(ROOT/'docs/research-2026-09-26/plan.md').read_bytes()


def main():
    out=ROOT/'runs/research_2026_09_26/foreground_response'
    if out.exists():raise FileExistsError(out)
    out.mkdir();(out/'executed_source.py').write_bytes(SOURCE);(out/'plan_snapshot.md').write_bytes(PLAN)
    mp=ROOT/'runs/salt_dust_audit/flux_response/matrices.npz'
    meta=ROOT/'runs/salt_dust_audit/flux_response/selected_objects.csv'
    hpfile=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    foldp=ROOT/'phase2/hierarchy/data/conditioned-multistart-best/rows.csv'
    maps=[ROOT/'data/dust/csfd-v2'/f'{n}.fits' for n in ['sfd_ebv','csfd_ebv','mask']]
    d=pd.read_csv(meta,dtype={'CID':str});assert set(d.IDSURVEY)=={10}
    fields=pd.read_csv(foldp,dtype={'CID':str});fields=fields[fields.IDSURVEY==10].set_index('CID')
    with fits.open(hpfile) as h:
        head=h[1].data
        lookup={str(s).strip():i for i,s in enumerate(head['SNID'])}
        indices=[lookup[c] for c in d.CID]
        d['RA']=head['RA'][indices];d['DEC']=head['DEC'][indices]
        d['HEAD_MWEBV']=head['MWEBV'][indices]
    d['field']=[fields.loc[c,'field'] for c in d.CID]
    sky=SkyCoord(d.RA.to_numpy()*u.deg,d.DEC.to_numpy()*u.deg,frame='icrs').galactic
    hp=HEALPix(nside=2048,order='ring',frame='galactic')
    pixels,weights=hp.bilinear_interpolation_weights(sky.l,sky.b)
    nearest=hp.lonlat_to_healpix(sky.l,sky.b)
    assert np.max(abs(weights.sum(axis=0)-1))<1e-12
    for name,p in zip(['sfd','csfd','mask'],maps):
        with fits.open(p,memmap=True) as f:
            assert f[1].header['NSIDE']==2048 and f[1].header['ORDERING'].strip()=='RING'
            a=f[1].data.field(0).reshape(-1)
            if name=='mask':
                d['mask_nearest']=a[nearest].astype(int)
                d['reliable']=np.all((a[pixels].astype(int)&4)!=0,axis=0)
                d['correction_footprint']=np.all((a[pixels].astype(int)&1)!=0,axis=0)
            else:d[name+'_ebv']=np.sum(a[pixels]*weights,axis=0)
    d['delta_ebv_cleaning']=.86*(d.csfd_ebv-d.sfd_ebv)
    d['delta_ebv_resampling']=.86*d.sfd_ebv-d.HEAD_MWEBV
    d['delta_ebv_direct_replacement']=.86*d.csfd_ebv-d.HEAD_MWEBV
    arrays=np.load(mp);rows=[]
    for row in d.itertuples():
        cid=row.CID;prefix=cid+'__'
        J=arrays[prefix+'jacobian_flux'];G=arrays[prefix+'nuisance_flux'][:,4]
        r=arrays[prefix+'flux_observed']-arrays[prefix+'flux_model']
        for weighting in ['measurement_only','measurement_plus_model']:
            C=arrays[prefix+weighting+'_covariance'];L=np.linalg.cholesky(C)
            jw=solve_triangular(L,J,lower=True);gw=solve_triangular(L,G,lower=True);rw=solve_triangular(L,r,lower=True)
            Q=np.linalg.qr(jw,mode='reduced')[0]
            gp=gw-Q@(Q.T@gw);rp=rw-Q@(Q.T@rw)
            resp=arrays[prefix+weighting+'_response'][:,4]
            standardizer=np.array([1,.16087,-3.11780,0])
            for mode in ['cleaning','resampling','direct_replacement']:
                e=getattr(row,'delta_ebv_'+mode);v=gp*e
                # Changing fixed E in the fitting model induces -response*deltaE.
                shift=-float(standardizer@resp)*e
                b=float(rp@v);info=float(v@v)
                rows.append({'CID':cid,'field':row.field,'zHEL':row.zHEL,'weighting':weighting,'mode':mode,
                    'reliable':row.reliable,'delta_EBV':e,'pre_BBC_delta_standardized_mag':shift,
                    'projected_template_information':info,'data_template_inner_product':b,
                    'fixed_map_loglike_gain':b-.5*info,'baseline_projected_chi2':float(rp@rp),
                    'epochs':len(r)})
    r=pd.DataFrame(rows);d.to_csv(out/'map_samples.csv',index=False);r.to_csv(out/'object_responses.csv',index=False)
    summaries=[]
    for (weighting,mode),part in r.groupby(['weighting','mode']):
        for subset in ['all','reliable']:
            p=part if subset=='all' else part[part.reliable]
            if p.empty:continue
            info=float(p.projected_template_information.sum());b=float(p.data_template_inner_product.sum())
            shift=p.pre_BBC_delta_standardized_mag.to_numpy()
            summaries.append({'weighting':weighting,'mode':mode,'subset':subset,'objects':len(p),
                'fields':int(p.field.nunique()),'epochs':int(p.epochs.sum()),
                'fixed_map_loglike_gain':float(p.fixed_map_loglike_gain.sum()),
                'fixed_map_residual_snr_under_conditional_noise':float(np.sqrt(info)),
                'unrestricted_template_amplitude_mle':b/info if info>0 else None,
                'conditional_template_amplitude_sd':1/np.sqrt(info) if info>0 else None,
                'pre_BBC_standardized_shift_mean_mag':float(shift.mean()),
                'pre_BBC_standardized_shift_rms_mag':float(np.sqrt(np.mean(shift**2))),
                'pre_BBC_standardized_shift_min_max_mag':[float(shift.min()),float(shift.max())]})
    field=r.groupby(['weighting','mode','field']).agg(objects=('CID','count'),
        mean_z=('zHEL','mean'),mean_shift=('pre_BBC_delta_standardized_mag','mean'),
        score_gain=('fixed_map_loglike_gain','sum'),information=('projected_template_information','sum'))
    field.to_csv(out/'field_summary.csv')
    report={'scope':'Externally specified foreground-map alternatives in a fixed 64-object SALT tangent calculation; no cosmology fit or empirical bias estimate.',
        'coordinate_source':'SN HEAD RA/DEC, not host coordinates',
        'map_scale':.86,'astropy_healpix':astropy_healpix.__version__,
        'refit_sign':'delta_theta=-[J whitened pseudoinverse G] delta_EBV for a fixed model foreground change at fixed observed flux',
        'three_modes':'cleaning=.86(CSFD-SFD); resampling=.86*SFD_HEALPIX-HEAD; direct_replacement=sum. Primary is cleaning, retaining original map sampling.',
        'calibration_limits':'Derivatives retain original SALT training/calibration, model covariance and accepted mask. No retraining, Poisson/detection redraw, classifier selection or BBC regeneration. A map comparison cannot certify CSFD truth or calibrate its uncertainty.',
        'population_limits':'Deterministic redshift-ranked sample; object/field overlap and selection persist. Do not extrapolate a 64-object range into a survey-wide physical bound.',
        'summaries':summaries}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    inputs=[mp,meta,hpfile,foldp,*maps]
    outputs=[out/f for f in ['map_samples.csv','object_responses.csv','field_summary.csv','summary.json']]
    (out/'manifest.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),
        'source_sha256':hashlib.sha256(SOURCE).hexdigest(),'plan_sha256':hashlib.sha256(PLAN).hexdigest(),
        'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
