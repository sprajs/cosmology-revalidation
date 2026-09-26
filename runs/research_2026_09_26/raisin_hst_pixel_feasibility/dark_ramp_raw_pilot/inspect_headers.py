"""Inspect FITS headers, structure and checksum status, never image values."""
import hashlib,json,warnings
from pathlib import Path
from astropy.io import fits

B=Path(__file__).resolve().parent
P=json.loads((B/'verified-products.json').read_text())
KEYS=('EXTNAME','EXTVER','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','SAMPTIME','SAMPSEQ','SAMP_SEQ','NSAMP','SAMPZERO','DELTATIM','ROUTTIME','TIME','EXPSTART','EXPEND','EXPTIME','DETECTOR','SUBARRAY','SUBTYPE','CCDGAIN','ATODGNA','ATODGNB','ATODGNC','ATODGND','READNSEA','READNSEB','READNSEC','READNSED','CAL_VER','DARKFILE','CCDTAB','BPIXTAB','NLINFILE','ZSIGFILE','PFLTFILE','DFLTFILE','LFLTFILE','UNITCORR','NLINCORR','ZOFFCORR','BLEVCORR','DQICORR','DARKCORR','FLATCORR','CHECKSUM','DATASUM')
out={'status':'partial','records':[]}
for src in P['records']:
    path=B/'files'/src['filename']
    row={'filename':src['filename'],'size':path.stat().st_size,'expected_size':src['size'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'hdus':[]}
    assert row['size']==row['expected_size']
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        with fits.open(path,mode='readonly',memmap=True,do_not_scale_image_data=True,checksum=False,lazy_load_hdus=False) as hdul:
            hdul.verify('exception')
            for index,hdu in enumerate(hdul):
                h=hdu.header
                hdur={'index':index,'class':type(hdu).__name__,'header':{k:h.get(k) for k in KEYS if k in h},'header_card_count':len(h.cards),'fileinfo':{k:v for k,v in (hdu.fileinfo() or {}).items() if k in ('hdrLoc','datLoc','datSpan')}}
                # FITS checksum verification is integrity-only; no pixel statistic is computed.
                hdur['verify_checksum']=hdu.verify_checksum() if 'CHECKSUM' in h else None
                hdur['verify_datasum']=hdu.verify_datasum() if 'DATASUM' in h else None
                row['hdus'].append(hdur)
        row['warnings']=[str(w.message) for w in ws]
    out['records'].append(row)
    (B/'header-inspection.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
out['status']='complete'
(B/'header-inspection.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
