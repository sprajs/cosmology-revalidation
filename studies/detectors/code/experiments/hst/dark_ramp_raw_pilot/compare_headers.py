"""Compare RAW metadata only: first eight chronological reads, geometry, refs."""
import hashlib,json
from pathlib import Path
from astropy.io import fits

B=Path(__file__).resolve().parent
PAIRS=(('search','idbx43p7q','icxoi1bcq'),('template','idp247tnq','icxoi4hgq'))
PKEYS=('DETECTOR','SAMP_SEQ','NSAMP','SUBARRAY','SUBTYPE','SAMPZERO','CCDTAB','BPIXTAB','NLINFILE','DARKFILE','ZSIGFILE','ATODGNA','ATODGNB','ATODGNC','ATODGND','READNSEA','READNSEB','READNSEC','READNSED','ZOFFCORR','BLEVCORR','NLINCORR','DARKCORR','FLATCORR','DQICORR','UNITCORR','CAL_VER')
GKEYS=('BITPIX','NAXIS','NAXIS1','NAXIS2','LTV1','LTV2','LTM1_1','LTM2_2','BUNIT','BSCALE','BZERO')
def read(root):
    path=B/'files'/(root+'_raw.fits')
    with fits.open(path,memmap=True,do_not_scale_image_data=True,lazy_load_hdus=False) as h:
        primary={k:h[0].header.get(k) for k in PKEYS}
        groups=[]
        for sci in [x for x in h if x.header.get('EXTNAME')=='SCI']:
            hh=sci.header;extver=hh['EXTVER']
            samp=h['SAMP',extver].header;tim=h['TIME',extver].header
            groups.append({'extver':extver,'sampnum':hh.get('SAMPNUM'),'samptime':hh.get('SAMPTIME'),'deltatim':hh.get('DELTATIM'),'routtime':hh.get('ROUTTIME'),'samp_pixvalue':samp.get('PIXVALUE'),'time_pixvalue':tim.get('PIXVALUE'),'time_naxis':tim.get('NAXIS'),'samp_naxis':samp.get('NAXIS'),'science_geometry':{k:hh.get(k) for k in GKEYS}})
        groups.sort(key=lambda x:x['samptime'])
        return {'primary':primary,'groups':groups,'hdu_count':len(h),'file_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
out={'status':'complete','pairs':[],'scope':'FITS header metadata only, no SCI/ERR/DQ/TIME arrays read'}
for label,dark,science in PAIRS:
    a=read(dark);b=read(science)
    ar=a['groups'][:8];br=b['groups']
    fields=('sampnum','samptime','deltatim','samp_pixvalue','time_pixvalue','time_naxis','samp_naxis','science_geometry')
    matched={k:[x[k]==y[k] for x,y in zip(ar,br)] for k in fields}
    primary_equal={k:a['primary'][k]==b['primary'][k] for k in PKEYS}
    out['pairs'].append({'visit':label,'dark_root':dark,'science_root':science,'dark':a,'science':b,'first8_field_exact':matched,'first8_all_exact':all(all(z) for z in matched.values()),'primary_equal':primary_equal,'primary_mismatches':{k:[a['primary'][k],b['primary'][k]] for k in PKEYS if not primary_equal[k]}})
(B/'header-comparison.json').write_text(json.dumps(out,indent=2)+'\n')
