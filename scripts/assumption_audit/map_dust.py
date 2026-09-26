"""Sample NASA/Chiang CSFD-v2 maps at released Pantheon+ coordinates.

This audits foreground-map structure, not a new standardized SN correction.
Requires astropy-healpix==1.1.2 in addition to the pinned project environment.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import astropy_healpix
from astropy_healpix import HEALPix

ROOT = Path(__file__).resolve().parents[2]


def digest(p):
    with open(p,'rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def stats(x):
    return {'mean':float(np.mean(x)), 'sd':float(np.std(x)),
            'q05_q50_q95':np.quantile(x,[.05,.5,.95]).tolist(), 'min':float(min(x)), 'max':float(max(x))}


def main():
    base=ROOT/'data/dust/csfd-v2'
    p=ROOT/'sources/repos/CobayaSampler__sn_data/PantheonPlus/Pantheon+SH0ES.dat'
    f=pd.read_csv(p,sep=r'\s+')
    sky=SkyCoord(f.RA.to_numpy()*u.deg,f.DEC.to_numpy()*u.deg,frame='icrs').galactic
    hp=HEALPix(nside=2048,order='ring',frame='galactic')
    pix=hp.lonlat_to_healpix(sky.l,sky.b)
    ind,weights=hp.bilinear_interpolation_weights(sky.l,sky.b)
    assert np.max(abs(weights.sum(axis=0)-1))<1e-12
    cols={}
    for name in ['sfd_ebv','csfd_ebv','mask']:
        with fits.open(base/f'{name}.fits',memmap=True) as hdus:
            h=hdus[1].header
            assert h['NSIDE']==2048 and h['ORDERING'].strip()=='RING'
            a=hdus[1].data.field(0).reshape(-1)
            assert len(a)==12*2048**2
            cols[name+'_nearest']=np.asarray(a[pix]).astype(float)
            if name!='mask':
                cols[name]=np.sum(a[ind]*weights,axis=0)
            else:
                cols['all_interpolation_pixels_cosmology']=np.all((a[ind].astype(int)&4)!=0,axis=0)
    for k,v in cols.items():f[k]=v
    f['galactic_l']=sky.l.deg;f['galactic_b']=sky.b.deg
    # Both maps retain the original SFD normalization. Pantheon uses 0.86.
    f['delta_EBV_scaled']=.86*(f.csfd_ebv-f.sfd_ebv)
    f['nearest_delta_EBV_scaled']=.86*(f.csfd_ebv_nearest-f.sfd_ebv_nearest)
    f['cosmology_mask']=(f.mask_nearest.astype(int)&4)!=0
    f['lss_correction_mask']=(f.mask_nearest.astype(int)&1)!=0
    f['analysis_selected']=f.zHD>.01
    out=ROOT/'runs/assumption_audit';out.mkdir(parents=True,exist_ok=True)
    dest=out/'foreground-map-samples.csv'
    save=['CID','IDSURVEY','zHD','RA','DEC','MWEBV','galactic_l','galactic_b',*cols.keys(),
          'delta_EBV_scaled','nearest_delta_EBV_scaled','cosmology_mask','lss_correction_mask','analysis_selected']
    f[save].to_csv(dest,index=False)
    result={'map_version':'Chiang CSFD v2, NASA LAMBDA, 2023-08-02',
            'statistic':'E(B-V) map differences only; not standardized magnitude residuals',
            'map_scale':.86,'astropy_healpix':astropy_healpix.__version__, 'samples':{}}
    for name,mask in [('all_selected',f.analysis_selected),
                      ('reliable_selected',f.analysis_selected & f.cosmology_mask & f.all_interpolation_pixels_cosmology)]:
        d=f[mask];v=d.delta_EBV_scaled.to_numpy()
        X=np.column_stack([np.ones(len(d)),d.MWEBV])
        fit=X@np.linalg.lstsq(X,v,rcond=None)[0]
        result['samples'][name]={
            'rows':len(d),'literal_unique_CIDs':int(d.CID.nunique()),
            'in_correction_footprint':int(d.lss_correction_mask.sum()),
            'in_cosmology_mask':int(d.cosmology_mask.sum()),
            'scaled_delta_EBV_mag':stats(v),
            'RMS_after_fitting_constant_and_global_MWEBV_scale':float(np.sqrt(np.mean((v-fit)**2))),
            'fraction_centered_variance_not_explained_by_global_scale':float(np.sum((v-fit)**2)/np.sum((v-v.mean())**2)),
            'bilinear_minus_nearest_delta_RMS':float(np.sqrt(np.mean((v-d.nearest_delta_EBV_scaled)**2))),
            'z_bins':[]}
        for lo,hi in zip([.01,.05,.1,.3,.6,1.],[.05,.1,.3,.6,1.,3.]):
            g=d[(d.zHD>=lo)&(d.zHD<hi)]
            if len(g):
                result['samples'][name]['z_bins'].append({'z_lo':lo,'z_hi':hi,'rows':len(g),
                    'literal_unique_CIDs':int(g.CID.nunique()),'delta_EBV_mag':stats(g.delta_EBV_scaled)})
    result['caveats']=[
        'The published mask is respected; even the reliable area retains reconstruction uncertainty and some one-halo residuals.',
        'Rows are not independent sky samples. Duplicate observations and common fields are retained; no naive standard errors or significance are assigned.',
        'A global scale cannot span this sampled foreground pattern, but that alone does not quantify missing standardized-distance covariance.',
        'Do not multiply delta_EBV by a constant R_B and subtract it from m_b_corr: bandpasses, redshift, light-curve colour, training and selection respond jointly.',
        'No final photometry, SN fit, host age or BAO distance was changed.',
    ]
    target=out/'foreground-map-summary.json';target.write_text(json.dumps(result,indent=2)+'\n')
    inputs=[p,Path(__file__),*[base/f'{n}.fits' for n in ['sfd_ebv','csfd_ebv','mask']],base/'readme.txt']
    (out/'foreground-map-manifest.json').write_text(json.dumps({
        'inputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in inputs},
        'outputs_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [target,dest]},
        'source_urls':{n:'https://lambda.gsfc.nasa.gov/data/foregrounds/CSFD/'+n for n in
                       ['sfd_ebv.fits','csfd_ebv.fits','mask.fits','readme.txt']},
        'environment':'project .venv plus astropy-healpix==1.1.2, installed separately',
    },indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
