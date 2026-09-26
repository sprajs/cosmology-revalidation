from pathlib import Path
import csv,hashlib,json,re,collections
import numpy as np
from astropy.io import fits
P=Path(__file__).parent;root=Path.cwd();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=P/'metadata_protocol.json';table=P/'spectra/Lu2023_TableA1_extension.csv'
with table.open(newline='') as f: raw=list(csv.reader(f))
header=raw[0]; rows=[]
for a in raw[1:]:
    r={k:a[i] for i,k in enumerate(header) if header.count(k)==1}
    r['SNRH_first']=a[16];r['SNRH_second']=a[17]
    rows.append(r)
# Paper Table1 exact official names only, no other sections.
m20=root/'runs/research_2026_09_26/bayesn_distance_identification/Mandel-M20-paper.txt'
s=m20.read_text(); t=s[s.index('Table 1. Table of supernovae'):s.index('4     IMPLEMENTATION')]
m20names=set(re.findall(r'^\s+(SN\d{4}[A-Za-z]+)\s',t,re.M));assert len(m20names)==79
norm=lambda x:re.sub(r'^sn','',re.sub(r'\s+','',x).lower())
dr3=root/'runs/research_2026_09_26/astra_design/raisin_phase_identification/dr3_extension/objects.csv'
listed=root/'runs/research_2026_09_26/astra_design/raisin_phase_identification/objects.csv'
pub=root/'runs/research_2026_09_26/raisin_differential/frozen-membership.csv'
sets={}
for n,p in [('DR3_134',dr3),('RAISIN_CSP_list76',listed),('RAISIN_selected79',pub)]:
    with p.open() as f:rr=list(csv.DictReader(f))
    sets[n]={norm(x['CID']) for x in rr}
sets['M20_Table1_79']={norm(x) for x in m20names}
ledger=[];allheaders={};objs=collections.defaultdict(list)
for r in rows:
    base=r['filename']; fp=P/'spectra/observed_spectra/spec_fits'/f'{base}.fits';tp=P/'spectra/observed_spectra/spec_txt'/f'{base}.txt'
    with fits.open(fp) as f:
        d=np.asarray(f[0].data); hh=dict(f[0].header)
    txt=np.loadtxt(tp)
    assert d.shape[0]==3
    idx=np.searchsorted(d[0],txt[:,0]);assert np.all(idx<d.shape[1])
    assert np.array_equal(txt,d[:,idx].T,equal_nan=True)
    omitted=np.setdiff1d(np.arange(d.shape[1]),idx)
    w=txt[:,0];err=txt[:,2]
    # Never form a spectral color, signal statistic or template residual.
    assert np.all(np.isfinite(w)) and np.all(np.diff(w)>0)
    z=float(r['zhel']);phase=float(r['epoch']);tm=float(r['Tmax(MJD)']);mjd=float(r['MJD'])
    entry={**r,'wmin_um':float(w.min()),'wmax_um':float(w.max()),'rest_wmin_um':float(w.min()/(1+z)),'rest_wmax_um':float(w.max()/(1+z)),'pixels':len(w),'FITS_pixels':d.shape[1],'FITS_omitted_from_TXT':len(omitted),'omitted_FITS_all_nonfinite_flux':bool(np.all(~np.isfinite(d[1,omitted]))),'YUNITS':hh.get('YUNITS'),'XUNITS':hh.get('XUNITS'),'ERRSCALE':hh.get('ERRSCALE'),'DAIRMASS':hh.get('DAIRMASS'),'nonfinite_error_pixels':int((~np.isfinite(err)).sum()),'nonpositive_finite_error_pixels':int(((err<=0)&np.isfinite(err)).sum()),'FITS_TXT_equal_including_NaN':True,'phase_recomputed_from_rounded_table':(mjd-tm)/(1+z),'phase_delta_from_rounded_metadata':phase-(mjd-tm)/(1+z)}
    ledger.append(entry);objs[r['name']].append(r)
    if not allheaders:allheaders={'filename':base,'header':hh}
with (P/'spectrum-schema-ledger.csv').open('w',newline='') as f:
    dw=csv.DictWriter(f,fieldnames=list(ledger[0]));dw.writeheader();dw.writerows(ledger)
objledger=[]
for name,rr in sorted(objs.items()):
    info=rr[0]; samekeys=['zhel','Tmax(MJD)','e(Tmax)','sBV','e(sBV)','EBV_MW','spec_counts']
    assert all(len({r[k] for r in rr})==1 for k in samekeys)
    assert len(rr)==int(info['spec_counts'])
    row={'name':name,'spectra':len(rr),**{k:info[k] for k in samekeys},'M20_Table1_name_match':norm(name) in sets['M20_Table1_79'],'DR3_134_name_match':norm(name) in sets['DR3_134'],'RAISIN_CSP_list76_name_match':norm(name) in sets['RAISIN_CSP_list76'],'RAISIN_selected79_name_match':norm(name) in sets['RAISIN_selected79']}
    base=[r for r in rr if r['host_contamination']=='No' and .75<=float(r['sBV'])<=1.18 and float(r['e(Tmax)'])<=.5]
    for band in ['Y','J','H']:
        col='SNR'+band if band!='H' else 'SNRH_first'
        for suffix,snr in [('SNR10',10),('noSNR',-np.inf)]:
            q=[r for r in base if float(r[col])>=snr]
            early=[r for r in q if -7<=float(r['epoch'])<=7]
            late=[r for r in q if 10<=float(r['epoch'])<=20]
            later=[r for r in q if 20<=float(r['epoch'])<=35]
            row[f'{band}_{suffix}_early']=len(early);row[f'{band}_{suffix}_late']=len(late);row[f'{band}_{suffix}_later']=len(later)
            row[f'{band}_{suffix}_paired']=bool(early and late)
            row[f'{band}_{suffix}_later_paired']=bool(early and later)
    objledger.append(row)
with (P/'object-support-and-overlap.csv').open('w',newline='') as f:
    dw=csv.DictWriter(f,fieldnames=list(objledger[0]));dw.writeheader();dw.writerows(objledger)
summary={'protocol_sha256':sha(protocol),'scope':'metadata/schema only; no spectral integral or model residual evaluated','table_header':header,'duplicate_header_names':{k:header.count(k) for k in set(header) if header.count(k)>1},'duplicate_H_columns_identical':all(r['SNRH_first']==r['SNRH_second'] for r in rows),'spectra':len(rows),'objects':len(objs),'host_contamination_counts':dict(collections.Counter(r['host_contamination'] for r in rows)),'ranges':{k:[min(float(r[k]) for r in rows),max(float(r[k]) for r in rows)] for k in ['zhel','epoch','sBV','e(Tmax)','MJD']},'wavelength_ranges_um':{k:[min(r[k] for r in ledger),max(r[k] for r in ledger)] for k in ['wmin_um','wmax_um','rest_wmin_um','rest_wmax_um']},'pixels':sorted(set(r['pixels'] for r in ledger)),'FITS_extra_pixels':sum(r['FITS_omitted_from_TXT'] for r in ledger),'all_omitted_FITS_nonfinite_flux':all(r['omitted_FITS_all_nonfinite_flux'] for r in ledger),'ERRSCALE_range':[min(r['ERRSCALE'] for r in ledger if r['ERRSCALE'] is not None),max(r['ERRSCALE'] for r in ledger if r['ERRSCALE'] is not None)],'nonfinite_error_spectra':sum(r['nonfinite_error_pixels']>0 for r in ledger),'nonfinite_error_pixels':sum(r['nonfinite_error_pixels'] for r in ledger),'nonpositive_finite_error_pixels':sum(r['nonpositive_finite_error_pixels'] for r in ledger),'max_abs_phase_rounding_delta':max(abs(r['phase_delta_from_rounded_metadata']) for r in ledger),'name_overlap':{n:[name for name in objs if norm(name) in ss] for n,ss in sets.items()},'paired_support':{},'header_example':allheaders}
for band in ['Y','J','H']:
    summary['paired_support'][band]={}
    for suffix in ['SNR10','noSNR']:
        for l in ['', 'later_']:
            names=[r['name'] for r in objledger if r[f'{band}_{suffix}_{l}paired']]
            summary['paired_support'][band][suffix+'_'+l+'names']=names
            summary['paired_support'][band][suffix+'_'+l+'count']=len(names)
summary['inputs']={str(p):sha(p) for p in [protocol,table,m20,dr3,listed,pub]}
(P/'metadata-result.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['header_example','inputs','table_header']},indent=2))
