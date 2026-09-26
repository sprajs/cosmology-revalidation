#!/usr/bin/env python3
"""Prepare original DES inputs and an explicitly versioned SNANA baseline."""
from pathlib import Path
import hashlib, json, re
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'phase2/official'
RELEASE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3'

def main():
    d=pd.read_csv(RELEASE/'4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv',dtype={'CID':str})
    h=.03754*(expit((d.HOST_LOGMASS.to_numpy()-10)/.001)-.5)
    reconstructed=-2.5*np.log10(d.x0)+.16087*d.x1-3.1178*d.c+h-d.biasCor_mu+29.95821
    ledger=d[['CID','IDSURVEY','zHD','MU']].copy()
    ledger['mx']=-2.5*np.log10(d.x0);ledger['stretch']=.16087*d.x1;ledger['color']=-3.1178*d.c
    ledger['host']=h;ledger['negative_bias']=-d.biasCor_mu;ledger['negative_M0']=29.95821
    ledger['reconstructed_mu']=reconstructed;ledger['residual']=reconstructed-d.MU
    ledger.to_csv(OUT/'results/distance-algebra.csv',index=False,float_format='%.17g')
    host={}
    for group in ['DES','Foundation','LOWZ']:
        head=fits.getdata(RELEASE/f'0_DATA/DES-SN5YR_{group}/DES-SN5YR_{group}_HEAD.FITS.gz',1)
        for row in head:host[str(row['SNID']).strip()]=float(row['HOSTGAL_LOGMASS'])
    hmass=d.CID.map(host)
    assert hmass.notna().all()
    head_step=.03754*(expit((hmass-10)/.001)-.5)
    ledger['head_logmass']=hmass;ledger['head_host_step']=head_step
    ledger['head_mass_reconstructed_mu']=reconstructed-h+head_step
    ledger['head_mass_residual']=ledger.head_mass_reconstructed_mu-d.MU
    ledger.to_csv(OUT/'results/distance-algebra.csv',index=False,float_format='%.17g')
    away=np.abs(d.HOST_LOGMASS.to_numpy()-10)>.02
    diag={'n':len(d),'parameters':{'alpha':.16087,'beta':3.1178,'gamma':.03754,'M0':-29.95821,'host_transition_width':.001},
          'max_abs_mag':float(np.max(np.abs(ledger.residual))), 'rms_mag':float(np.sqrt(np.mean(ledger.residual**2))),
          'away_rounded_host_boundary_n':int(away.sum()),'away_max_abs_mag':float(np.max(np.abs(ledger.residual[away]))),
          'away_rms_mag':float(np.sqrt(np.mean(ledger.residual[away]**2))),
          'head_host_mass_all_rows_max_abs_mag':float(np.abs(ledger.head_mass_residual).max()),
          'head_host_mass_all_rows_rms_mag':float(np.sqrt(np.mean(ledger.head_mass_residual**2))),
          'warning':'Uses published nuisance and bias-correction outputs; does not independently regenerate BBC.'}
    (OUT/'diagnostics/distance-algebra.json').write_text(json.dumps(diag,indent=2)+'\n')
    des=d[d.IDSURVEY==10].sort_values(['zHD','CID']).reset_index(drop=True)
    sample=des.iloc[np.rint(np.linspace(0,len(des)-1,12)).astype(int)]
    sample.to_csv(OUT/'inputs/pilot12.csv',index=False)
    groups=[('pilot12',sample,'desSMP_5yr'),('all_des',des,'desSMP_5yr'),
            ('lowz',d[(d.IDSURVEY!=10)&(d.IDSURVEY!=150)],'lowz_des5yr'),
            ('foundation',d[d.IDSURVEY==150],'found_des5yr')]
    for label,frame,base in groups:
        (OUT/f'inputs/{label}.cid').write_text('\n'.join(frame.CID)+'\n')
        src=(RELEASE/f'7_PIPPIN_FILES/base_files/lcfit/lcfit_{base}.nml').read_text().split('#END_YAML')[1]
        src=src.replace("'$SNDATA_ROOT/lcmerge/DES-SN5YR'",repr(str(RELEASE/'0_DATA')))
        src=src.replace("'$SNDATA_ROOT/models/SALT3/SALT3.DES5YR'",repr(str(RELEASE/'2_LCFIT_MODEL/SALT3.DES5YR')))
        src=re.sub(r"ROOTFILE_OUT\s*=.*\n",'',src)
        src=re.sub("'DES5YR_(?:DES|LOWZ|FOUND)_CHANGEME'",repr(str(OUT/f'results/snana_{label}')),src)
        src=src.replace('&SNLCINP',f"&SNLCINP\n SNCID_LIST_FILE = '{OUT}/inputs/{label}.cid'\n")
        (OUT/f'inputs/snana_{label}.nml').write_text(src)
    inputs=[RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz',RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_PHOT.FITS.gz',RELEASE/'4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv',RELEASE/'7_PIPPIN_FILES/base_files/lcfit/lcfit_desSMP_5yr.nml',OUT/'PLAN.md',OUT/'inputs/SNDATA_ROOT_2024-07-03.tar.gz']
    inputs+=list((RELEASE/'2_LCFIT_MODEL/SALT3.DES5YR').glob('*.gz'))
    inputs+=list((OUT/'inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR').glob('*'))
    manifest=[]
    for p in inputs:
        if p.is_file():
            hs=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
            manifest.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hs})
    (OUT/'inputs/source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(diag,indent=2))

if __name__=='__main__':main()
