#!/usr/bin/env python3
"""Generic instrument response curves, source-pinned and explicitly approximate."""
from pathlib import Path
import json,hashlib,datetime,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
from astropy.io.votable import parse_single_table
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    arrays={};rows=[]
    for band in ['u','g','r','i','z','J','H','Ks']:
        p=W/'filters'/f'{band}.xml';tree=ET.parse(p)
        params={e.attrib['name']:e.attrib.get('value','') for e in tree.iter() if e.tag.split('}')[-1]=='PARAM'}
        t=parse_single_table(p).to_table();wave=np.asarray(t['Wavelength'],float);response=np.asarray(t['Transmission'],float)
        assert params['WavelengthUnit']=='Angstrom' and params['DetectorType']=='1'
        assert np.isfinite(wave).all() and np.isfinite(response).all() and (np.diff(wave)>0).all() and (response>=0).all()
        pivot=float(np.sqrt(np.trapezoid(response*wave,wave)/np.trapezoid(response/wave,wave)))
        assert abs(pivot/float(params['WavelengthPivot'])-1)<2e-3
        band=band.lower();arrays[band+'_wavelength_A']=wave;arrays[band+'_throughput']=response
        csv=W/'filters'/f'{band}-response.csv';pd.DataFrame({'wavelength_A':wave,'throughput':response}).to_csv(csv,index=False)
        rows.append({'band':band,'filter_id':params['filterID'],'url':'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='+params['filterID'],'description':params['Description'],'components':params.get('components'), 'detector':'photon counter','profile_reference':params.get('ProfileReference'),'pivot_A_direct':pivot,'pivot_A_SVO':float(params['WavelengthPivot']),'samples':len(wave),'input_path':str(p.relative_to(ROOT)),'input_sha256':sha(p),'output_path':str(csv.relative_to(ROOT)),'output_sha256':sha(csv)})
    np.savez(W/'filters/des-deep-responses.npz',**arrays)
    q=pd.read_csv(W/'des-deep-crosswalk.csv');coeff={}
    for b in ['U','G','R','I','Z','J','H','KS']:
        f=q[f'deep_BDF_FLUX_CALIB_{b}'];fd=q[f'deep_BDF_FLUX_DERED_CALIB_{b}'];ebv=q.deep_EBV_SFD98;ok=(f>0)&(fd>0)&(ebv>0)
        r=2.5*np.log10(fd[ok]/f[ok])/ebv[ok];coeff[b.lower()]={'n':int(ok.sum()),'median_A_over_EBV_SFD98':float(r.median()),'minmax': [float(r.min()),float(r.max())]}
    result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'band_identity':'Hartley+2022 sections2.1/2.2: DECam ugriz and VISTA VIRCam VIDEO JHKs for C3,E2,X3; not CFHT u or SDSS passbands.','approximation':'Public SVO full-system reference curves; not the exact exposure/CCD-weighted throughput of the deep coadds. Account for field calibration, throughput and atmosphere sensitivity. The SVO Vega calibration metadata is NOT applied to this AB-calibrated catalogue.','photometry_calibration':'DES Y3 source documentation and exact flux-mag identity give AB ZP30 for stored fluxes: divide by1000 to obtain AB nanomaggies (ZP22.5); all8 bands are AB even though NIR images were initially Vega-calibrated.','foreground':'FLUX_CALIB contains field-to-field calibration but no foreground removal; FLUX_DERED_CALIB additionally removes SFD98-based foreground. Do not deredden twice. Paper sec3.4 uses F99 Rv3.1 optical coefficients and Gonzalez-Fernandez2018 NIR coefficients; empirical per-band coefficients recovered below.','foreground_coefficients_from_released_flux_ratio':coeff,'filters':rows,'npz_path':'.work/host-transport/filters/des-deep-responses.npz','npz_sha256':sha(W/'filters/des-deep-responses.npz'),'code_sha256':sha(Path(__file__)),'paper_sha256':sha(W/'des/deep-paper.pdf')}
    (O/'filter-interface.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'coefficients':coeff,'bands':[r['filter_id'] for r in rows]},indent=2))
if __name__=='__main__':main()
