"""Post-outcome numerical source check of already inspected native CSP K tables.

Not an independent preregistration: exploratory native discrepancies were known.
"""
from pathlib import Path
import csv, hashlib, json
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_passband_contrast'
OUT=BASE/'native-control'
KCOR=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
INPUT=KCOR.with_suffix('.input')
FILT=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=False)
    rows=list(csv.DictReader((BASE/'grid.csv').open()))
    native=[];filters={};zps={}
    for line in INPUT.read_text().splitlines():
        a=line.split()
        if a and a[0]=='FILTER:':filters[a[1]]=a[2];zps[a[1]]=float(a[3])
    with fits.open(KCOR) as h:
        k=h['KCOR'];names=k.columns.names;ncol=k.header['TFIELDS'];nrow=k.header['NAXIS2']
        assert all(k.header[f'TFORM{i}']=='1E' for i in range(1,ncol+1))
        assert k.header['NAXIS1']==4*ncol
        # Duplicate K_YJ names preclude Astropy structured .data. Decode declared
        # contiguous IEEE big-endian float32 fields; do not rename or edit input.
        with KCOR.open('rb') as f:
            f.seek(k.fileinfo()['datLoc']);table=np.frombuffer(f.read(4*ncol*nrow),dtype='>f4').reshape(nrow,ncol)
        selected=table[(table[:,1]==0)&(table[:,2]==0)]
        controls={'J_RC2_minus_RC1':('CSP-J','CSP-j','K_Jj'),
                  'H_WIRC_minus_RetroCam':('CSP-H','CSP-h','K_Hh'),
                  'Y_WIRC_minus_RetroCam':('CSP-Y','CSP-y','K_Yy')}
        for pair,(a,b,name) in controls.items():
            assert names.count(name)==1
            for r in rows:
                if r['pair']!=pair or float(r['z'])!=0:continue
                match=selected[selected[:,0]==float(r['phase'])];assert len(match)==1
                value=float(match[0,names.index(name)])
                ours=float(r['delta_mag_equal_BD17'])+zps[b]-zps[a]
                native.append(dict(pair=pair,phase=r['phase'],native_K=value,synthetic_with_primary_offset=ours,difference=ours-value))
        tab=h['FilterTrans'].data;wave=np.array(tab.field(0),float);comparison={}
        for name in ['CSP-Y','CSP-y','CSP-J','CSP-j','CSP-H','CSP-h']:
            data=np.loadtxt(FILT/filters[name]);expected=np.interp(wave,data[:,0],data[:,1],left=0,right=0)
            observed=np.asarray(tab[name],float)
            comparison[name]=dict(max_abs_transmission_gap=float(np.max(abs(expected-observed))),rms_gap=float(np.sqrt(np.mean((expected-observed)**2))),source=filters[name])
        zp={str(r['Filter Name']).strip():float(r['Primary Mag']) for r in h['ZPoff'].data}
    with (OUT/'native-grid.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=native[0]);w.writeheader();w.writerows(native)
    result=dict(design='Post-outcome source/numerical check; discrepancies known before script, no claim of preregistration or native integral exactness.',
                native_table_decoder='Verified all30 contiguous1E fields, repeatedK_YJ bypassed by positional byte decode; original file unchanged.',
                max_mag_gap={pair:max(abs(r['difference']) for r in native if r['pair']==pair) for pair in controls},
                filter_table_interpolation=comparison,printed_primary_magnitudes=zp,
                scope='At z=AVwarp=0, compare synthetic contrasts plus explicit BD17 magnitude offset to shipped K-correction. Native wavelength discretization differs. No native WIRC J entry; no photometric correction inferred.')
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    paths=[KCOR,INPUT,BASE/'grid.csv',Path(__file__)]+[FILT/filters[n] for n in comparison]
    (OUT/'manifest.json').write_text(json.dumps({'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},'outputs':{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
