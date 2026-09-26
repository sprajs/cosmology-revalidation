"""Header-only instrument and DQ-layout ledger for the four pilot RAWs."""
import json
from pathlib import Path
from astropy.io import fits

B=Path(__file__).resolve().parent
ROOTS=('idbx43p7q','idp247tnq','icxoi1bcq','icxoi4hgq')
KEYS=('ROOTNAME','EXPFLAG','OBSMODE','SCLAMP','DETECTOR','FILTER','APERTURE','SAMP_SEQ','NSAMP','SAMPZERO','SUBARRAY','SUBTYPE','CCDAMP','CCDGAIN','ATODGNA','ATODGNB','ATODGNC','ATODGND','READNSEA','READNSEB','READNSEC','READNSED','ZOFFCORR','BLEVCORR','NLINCORR','DQICORR','DARKCORR','FLATCORR','CRCORR','UNITCORR','CCDTAB','BPIXTAB','CRREJTAB','DARKFILE','NLINFILE','CRDS_CTX','CRDS_VER','CAL_VER','SAA_DARK','SAACRMAP')
out={'status':'complete','records':[],'scope':'FITS header cards only'}
for root in ROOTS:
    with fits.open(B/'files'/(root+'_raw.fits'),memmap=True,do_not_scale_image_data=True) as h:
        row={'root':root,'primary':{k:h[0].header.get(k) for k in KEYS},'dq_layout':[],'time_layout':[]}
        for x in h:
            header=x.header
            if header.get('EXTNAME')=='DQ':row['dq_layout'].append({'extver':header.get('EXTVER'),'naxis':header.get('NAXIS'),'naxis1':header.get('NAXIS1'),'naxis2':header.get('NAXIS2'),'pixvalue':header.get('PIXVALUE'),'bitpix':header.get('BITPIX')})
            if header.get('EXTNAME')=='TIME':row['time_layout'].append({'extver':header.get('EXTVER'),'naxis':header.get('NAXIS'),'pixvalue':header.get('PIXVALUE'),'bunit':header.get('BUNIT')})
        out['records'].append(row)
(B/'state-ledger.json').write_text(json.dumps(out,indent=2)+'\n')
