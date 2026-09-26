"""Offline descriptive quality/reference audit; no science pixels."""
import csv,hashlib,json
from collections import Counter
from pathlib import Path
from astropy.io import fits

B=Path(__file__).resolve().parent
Q=json.loads((B/'quality-header-results.json').read_text())
assert Q['status']=='complete' and len(Q['records'])==14
KEYS=('EXPFLAG','QUALITY','QUALCOM1','QUALCOM2','EXPSTART','EXPEND','EXPTIME','SAMP_SEQ','NSAMP','SUBARRAY','SUBTYPE','FILTER','SAA_TIME','SAA_EXIT','CCDAMP','CCDGAIN','ATODGNA','READNSEA','CCDTAB','OSCNTAB','BPIXTAB','NLINFILE','DARKFILE','CRDS_CTX','CAL_VER','ZOFFCORR','BLEVCORR','NLINCORR','DQICORR','CRCORR','DARKCORR','FLATCORR','UNITCORR')
rows=[{'root':r['root'],'visit':r['visit'],'header_sha256':r['header_sha256'],'source':r['source'],**{k:r['keys'].get(k) for k in KEYS}} for r in Q['records']]
with (B/'candidate-quality-ledger.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def table_info(filename):
    path=B/'reference_files'/filename
    with fits.open(path,memmap=True) as h:
        h.verify('exception')
        return h[0].header,{c.name:c.unit for c in h[1].columns},h[1].data
ccdh,units,ccd=table_info('t2c16200i_ccd.fits')
osc_h,osc_units,osc=table_info('q911321mi_osc.fits')
match=[(i,r) for i,r in enumerate(ccd,1) if r['CCDAMP'].strip()=='ABCD' and float(r['CCDGAIN'])==2.5]
assert len(match)==1
i,r=match[0]
ccd_row={'row':i,'CCDAMP':r['CCDAMP'].strip(),'CCDCHIP':int(r['CCDCHIP']),'CCDGAIN':float(r['CCDGAIN']),'CCDOFFST':[int(r['CCDOFST'+a]) for a in 'ABCD'],'BINAXIS':[int(r['BINAXIS1']),int(r['BINAXIS2'])],'ATODGAIN_electrons_per_DN':[float(r['ATODGN'+a]) for a in 'ABCD'],'READNSE_electrons':[float(r['READNSE'+a]) for a in 'ABCD'],'PEDIGREE':r['PEDIGREE'].strip(),'DESCRIP':r['DESCRIP'].strip()}
omatch=[(j,s) for j,s in enumerate(osc,1) if s['CCDAMP'].strip()=='ABCD' and int(s['NX'])==1024 and int(s['NY'])==1024]
assert len(omatch)==1
j,s=omatch[0]
osc_row={'row':j,'CCDAMP':s['CCDAMP'].strip(),'CCDCHIP':int(s['CCDCHIP']),'BIN':[int(s['BINX']),int(s['BINY'])],'NXY':[int(s['NX']),int(s['NY'])],'TRIM':[int(s['TRIMX1']),int(s['TRIMX2']),int(s['TRIMY1']),int(s['TRIMY2'])],'BIASSECTA':[int(s['BIASSECTA1']),int(s['BIASSECTA2'])],'BIASSECTB':[int(s['BIASSECTB1']),int(s['BIASSECTB2'])]}
out={'status':'complete descriptive metadata/reference audit, no cohort exclusions','candidate_count':len(rows),'expflag_counts':dict(Counter(r['EXPFLAG'] for r in rows)),'quality_blank_count':sum(r['QUALITY']=='' for r in rows),'qualcom1_blank_count':sum(r['QUALCOM1']=='' for r in rows),'qualcom2_blank_count':sum(r['QUALCOM2']=='' for r in rows),'all_same_sampling':len({(r['SAMP_SEQ'],r['NSAMP'],r['EXPTIME'],r['SUBARRAY']) for r in rows})==1,'all_same_ccdtab':len({r['CCDTAB'] for r in rows})==1,'all_same_oscntab':len({r['OSCNTAB'] for r in rows})==1,'ccdtab_observed_amp_gain_unique_row':ccd_row,'oscntab_fullframe_unique_row':osc_row,'pilot_template_log_evidence':{'file':'support_files/idp247tnq_log.txt','lines':[25,26,29,30,31],'tdf_down_readouts':9,'total_readouts':16,'flag':'INDETERMINATE','readout_time_log_seconds':2.91181,'exptime_log_seconds':702.938171},'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (B/'reference_files').glob('*.fits')}}
(B/'result.json').write_text(json.dumps(out,indent=2)+'\n')
