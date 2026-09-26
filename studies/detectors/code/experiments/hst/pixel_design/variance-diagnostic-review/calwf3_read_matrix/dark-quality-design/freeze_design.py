from pathlib import Path
import json,csv,hashlib,re
import numpy as np
R=Path.cwd();O=Path(__file__).resolve().parent;H=R/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
rows=list(csv.DictReader((H/'quality_reference_audit/candidate-quality-ledger.csv').open()))
logs=[]
for r in rows:
 f=O/'logs'/(r['root']+'_log.txt')
 if r['root']=='idp247tnq':f=H/'quality_reference_audit/support_files/idp247tnq_log.txt'
 s=f.read_text().split('CALBEG')[0]
 matches=re.findall(r'TDF down for (\d+) out of (\d+) readouts',s)
 pk=re.findall(r'(\d+) Science packets, (\d+) Engineering \(ULX\) packets, (\d+) Bad pkts',s)
 assert len(pk)==1
 assert len(matches)==(r['EXPFLAG']=='INDETERMINATE')
 down=int(matches[0][0]) if matches else None
 logs.append({'root':r['root'],'visit':r['visit'],'EXPFLAG':r['EXPFLAG'],'tdf_down':down,'total_reads':16,'prefix8_tdf_down_min':max(0,down-8) if down is not None else None,'prefix8_tdf_down_max':min(8,down) if down is not None else None,'science_packets':int(pk[0][0]),'engineering_packets':int(pk[0][1]),'bad_packets':int(pk[0][2]),'log_path':str(f.relative_to(R)),'log_sha256':sha(f)})
with (O/'log-quality-ledger.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=logs[0].keys());w.writeheader();w.writerows(logs)
pairs=[]
for v in ['search','template']:
 a=sorted([r for r in rows if r['visit']==v],key=lambda x:(float(x['EXPSTART']),x['root']))
 for i in range(0,len(a),2):pairs.append({'pair':f'{v}_{i//2+1}','visit':v,'earlier':a[i]['root'],'later':a[i+1]['root'],'kind':'primary_all14'})
pairs.append({'pair':'normal_search','visit':'search','earlier':'idbx41onq','later':'idbx43p7q','kind':'prespecified_sensitivity_one_pair'})
products={};source={}
for f in sorted((H/'dark_ramp_metadata/queries').glob('*products.json')):
 for r in json.loads(f.read_text()).get('data',[]):
  name=r.get('productFilename','');root=name.removesuffix('_raw.fits')
  if name.endswith('_raw.fits') and root in {x['root'] for x in rows}:
   x={'root':root,'filename':name,'uri':r['dataURI'],'size':int(r['size']),'existing_pilot':root in ['idbx43p7q','idp247tnq']}
   if root in products:assert products[root]==x
   products[root]=x;source[str(f.relative_to(R))]=sha(f)
assert len(products)==14 and sum(not x['existing_pilot'] for x in products.values())==12
comparison=json.loads((H/'dark_ramp_raw_pilot/header-comparison.json').read_text())
t=np.array([x['samptime'] for x in comparison['pairs'][0]['science']['groups'][1:]],float)
powers=[0,.4,1,3,6,10];hs=[]
for p in powers:
 u=np.abs((np.arange(7)-3)/3)**p;tb=np.sum(u*t)/sum(u);Q=np.sum(u*(t-tb)**2);h=u*(t-tb)/Q
 assert abs(sum(h))<1e-15 and abs(h@t-1)<1e-13
 hs.append(h)
np.savez(O/'fixed-operators.npz',times_seconds=t,powers=np.array(powers),h=np.array(hs))
protocol={'status':'frozen design before any dark SCI outcome; acquisition and executor require root release','scope':'Exploratory total digitized RAW paired second moments, never isolated electronic-noise covariance or science ERR correction','cohort':rows,'quality_logs':logs,'pairs':pairs,'acquisition':{'records':sorted(products.values(),key=lambda x:x['root']),'remaining_RAW_count':12,'remaining_expected_bytes':sum(x['size'] for x in products.values() if not x['existing_pilot']),'cap_bytes':600000000,'cap_seconds':180,'one_attempt_per_product':True,'no_substitution':True,'phases':'Acquire/verify missing NORMAL search member first; engineering uses only synthetic arrays plus its metadata/read decoding before cohort outcomes; remaining12 acquisition is independent of all pixel effects'},'time_operator':{'t_seconds_nominal':t.tolist(),'powers':powers,'primary_power':0,'definition':'u_i=abs((i-3)/3)^p; h_i=u_i*(t_i-tbar_u)/sum[u*(t-tbar_u)^2]; no data-selected weights','times_scope':'Recorded nominal schedule; does not assert correct physical exposure times under TDF flags','fixed_operator_npz_sha256':sha(O/'fixed-operators.npz')},'masks':{'geometry':'raw1024x1024 zero-based x,y5..1018 inclusive','primary':'Union of nonzero bad-pixel flags mapped by exact BPIXTAB source rules from3562029fi_bpx.fits and3562028ni_bpx.fits; same fixed detector mask for all14','reference_hash_gate':'Acquire official exact named references and freeze SHA/table-to-image mapping before SCI outcomes; no inferred scalar/image substitution for BPIXTAB table','raw_DQ':'retain in primary; report original flags/counts/contributions, no CR or residual-derived removal','geometry_only':'separate total-moment diagnostic before reference mask, no alternate primary selection','no':'no science-derived masks, no clipping, no finite-difference CR removal, no hot-pixel threshold on these ramps'},'estimands':{'difference':'later minus earlier, promoted float64 before subtraction, first7 nonzero chronological reads','read_matrix':'Gamma_pair_region=mean_pixels(delta_y delta_y^T)/2, uncentered; mean delta_y retained separately','slope':'mean_pixels[(h^T delta_y)^2]/2 and h^T Gamma h identity, units DN^2/s_nominal^2','gain':'secondary fixed per-amplifier CCDTAB A/D conversion, source orientation gate mandatory; never nominal CCDGAIN2.5 or final mean-gain FLT semantics','aggregation':'equal-pixel within fixed region; per-pair/visit ledgers; no pooled visits as exchangeable and no detection p','spatial':'fixed masked aperture-annulus contractions and spatial block summaries, retaining covariance rather than adding pixel variances'},'spatial_operator':{'centres':'256 detector centres(x,y)=(32+64i,32+64j),i,j0..15; zero-based pixel centres','radii_pixels':[.4/.12825,1.2/.12825,2/.12825],'fractional_area':'fixed geometric circle/square pixel overlaps, validated without SCI to absolute total-area and weight-sum1e-8; source or independent quadrature comparison','mask_rule':'m=fixed BPIXTAB union-good; w=m*a - sum(m*a)/sum(m*b)*m*b, no aperture coverage renormalization','support':'fixed-mask aperture area>=90% unmasked geometry and annulus>=75%; geometry/support determined beforeSCI, all rejected sites counted; no fallback mask choice','units':'raw DN contrast; not calibrated point-source flux or exact previous sky-plane aperture','blocks':'fixed16x16 detector blocks of64 pixels, clipped to active region; report block/pair values, no independent-pixel errors'},'engineering':{'reference_pair':['idbx41onq','idbx43p7q'],'before_all14_scoring':['unsigned FITS decoding/constant HDU expansion verified','chronological first8, EXTVER mapping and metadata identity','fixed-mask source mapping and coordinate/gain gates','synthetic test pass and independent implementation arithmetic','two independent float64 contraction routes agree1e-10 relative or1e-12 absolute'], 'outcome_policy':'engineering observed means/variance do not select objects/operators; any structural/numerical failure is preserved and fixed additively before scoring, not a new cohort'},'uncertainty':'descriptive conditional spatial blocks and3/4 primary disjoint pairs; one NORMAL pair is engineering sensitivity, no temporal population CI or white-read sigma claim','failures':'retain all14 plan, failed files/pairs/flags; missing required inputs stop complete-cohort result, no replacement; no automatic fallback to NORMAL-only or relaxed reference mask','resources_analysis':{'one_worker':True,'RAM_live_arrays_bytes':268435456,'new_output_bytes':100000000,'engineering_seconds':60,'full_analysis_seconds':180},'outputs':['acquisition.json plus input/source hashes','HDU/read-order/quality/reference/gain mapping ledger','fixed mask and support summary','read-matrix.npz with per-pair/quadrant/block Gamma and mean vectors','pair-quadrant-power.csv forDN and gain-secondary','aperture-pair-power.csv with signed contrast,secondmoment,support','raw-DQ/extreme-digital-value contribution ledger without deletions','normal-search-sensitivity.csv','arithmetic-gates.json','report.md separating total raw variability from calibrated/electronic covariance'], 'source_input_hashes':source|{str(f.relative_to(R)):sha(f) for f in [H/'quality_reference_audit/candidate-quality-ledger.csv',H/'quality_reference_audit/result.json',H/'dark_ramp_raw_pilot/header-comparison.json',H/'calwf3_variance_source/source-3.7.3/cridcalc.c',H/'calwf3_variance_source/source-3.7.3/noiscalc.c']}}
(O/'experiment-protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps({'sha256':sha(O/'experiment-protocol.json'),'primary_pairs':len(pairs)-1,'remaining_bytes':protocol['acquisition']['remaining_expected_bytes'],'tdf_counts':{str(k):sum(r['tdf_down']==k for r in logs) for k in [None,5,9,16]}}))
