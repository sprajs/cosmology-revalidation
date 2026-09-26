"""Post-outcome source-model contrasts of already frozen sufficient statistics.
Never opens a RAW FITS or changes an extraction/mask/projection.
"""
from pathlib import Path
import hashlib,json,csv
import numpy as np
O=Path(__file__).resolve().parent;H=O.parents[4];D=H/'dark_ramp_execution';P=O.parent/'dark-quality-design/fixed-operators.npz'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
paths=[D/'read-matrix.npz',D/'pair-quadrant-power.csv',D/'decomposition.csv',D/'raw-dq-digital-contributions.csv',D/'aperture-pair-power.csv',P,H/'calwf3_variance_source/source-3.7.3/cridcalc.c',H/'calwf3_native_replay/replay-result.json']
proto={'status':'Post-outcome interpretation, not a new observation test. Only saved sufficient statistics and declared operators; no RAW SCI reads, new masks, optimized weights, significance or rescaling.', 'reference':'Source linfit uses hardcoded21electronCDS and wf3.mean_gain. Compare reported read-only term and a hypothetical independent-single-read model sigma=CDS/sqrt2, with no Poisson/dark/CR/spatial component. Native mean_gain=2.35 is distinct from amplifier gains.', 'amplifier_gains':{'B':2.37,'C':2.31,'A':2.34,'D':2.38}, 'source_mean_gain':2.35,'inputs':{str(p):sha(p) for p in paths}}
assert not (O/'protocol.json').exists();(O/'protocol.json').write_text(json.dumps(proto,indent=2)+'\n')
with np.load(P,allow_pickle=False) as x:t=x['times_seconds'];hs=x['h'];powers=x['powers']
with np.load(D/'read-matrix.npz',allow_pickle=False) as x:G=x['quadrant_gamma'];means=x['quadrant_mean'];quad=list(x['quadrant_names'].astype(str))
qcsv=list(csv.DictReader((D/'pair-quadrant-power.csv').open()));pairs=list(dict.fromkeys(r['pair'] for r in qcsv));out=[];refs=[]
for i,(power,h) in enumerate(zip(powers,hs)):
 w=np.abs(np.arange(7)-3.)/3.;w=w**power;S=w.sum();Sx=w@t;Sxx=w@(t*t);den=S*Sxx-Sx*Sx
 h0=w*(S*t-Sx)/den;assert np.max(abs(h-h0))<1e-15
 vrep=21**2/(2.35**2)*S/den;vind=(21/2.35)**2/2*(h@h);pois=h@np.minimum.outer(t,t)@h
 refs.append(dict(power=float(power),source_mean_gain_read_var_DN2_per_s2=float(vrep),white_single_read_var_DN2_per_s2=float(vind),reported_over_white=float(vrep/vind),cumulative_Poisson_coefficient_s_inverse=float(pois),source_Poisson_coefficient_s_inverse=float(1/(t[-1]-t[0]))))
 for j,pair in enumerate(pairs):
  for k,q in enumerate(quad):
   v=float(h@G[j,0,k]@h);mean=float(h@means[j,0,k]);source=next(r for r in qcsv if r['pair']==pair and r['quadrant']==q and r['scope']=='masked' and float(r['power'])==power);assert abs(v-float(source['dn2_per_nominal_s2']))<1e-12
   out.append(dict(pair=pair,kind=source['kind'],quadrant=q,power=float(power),total_DN2_per_s2=v,reported_source_read_only_DN2_per_s2=float(vrep),total_over_source_read_only=v/vrep,total_over_white_single_read=v/vind,half_squared_mean=mean*mean/2,coherent_fraction=mean*mean/2/v,nominal_amp_electron2_per_s2=v*proto['amplifier_gains'][q]**2))
with (O/'source-reference-contrasts.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
primary=[r for r in out if r['power']==0 and r['kind']=='primary_all14'];normal=[r for r in out if r['power']==0 and r['pair']=='normal_search'];dec=list(csv.DictReader((D/'decomposition.csv').open()));digital=list(csv.DictReader((D/'raw-dq-digital-contributions.csv').open()));ap=list(csv.DictReader((D/'aperture-pair-power.csv').open()));ap0=[r for r in ap if float(r['power'])==0];apgroup=[]
for pair in pairs:
 a=[r for r in ap0 if r['pair']==pair];v=np.array([float(r['half_squared_dn2_per_nominal_s2']) for r in a]);apgroup.append(dict(pair=pair,N_fixed_supported_apertures=len(a),mean_half_squared_DN2_per_s2=float(v.mean()),median_half_squared_DN2_per_s2=float(np.median(v)),max_half_squared_DN2_per_s2=float(v.max()),note='Descriptive existing fixed-aperture contractions; no new pixel covariance or effective independent N inferred.'))
result={'reference_models':refs,'primary_p0':{'raw_total_range':[min(r['total_DN2_per_s2'] for r in primary),max(r['total_DN2_per_s2'] for r in primary)],'ratio_to_source_reported_read_only_range':[min(r['total_over_source_read_only'] for r in primary),max(r['total_over_source_read_only'] for r in primary)],'coherent_fraction_max':max(r['coherent_fraction'] for r in primary)},'normal_search_p0':normal,'power10_over_power0_range':[min(r['total_DN2_per_s2']/next(s['total_DN2_per_s2'] for s in out if s['pair']==r['pair'] and s['quadrant']==r['quadrant'] and s['power']==0) for r in out if r['power']==10),max(r['total_DN2_per_s2']/next(s['total_DN2_per_s2'] for s in out if s['pair']==r['pair'] and s['quadrant']==r['quadrant'] and s['power']==0) for r in out if r['power']==10)],'source_white_model_power10_over_power0':refs[-1]['white_single_read_var_DN2_per_s2']/refs[0]['white_single_read_var_DN2_per_s2'],'raw_DQ_endpoint_flagged_counts':sum(int(r['selected_pixels']) for r in digital),'CR_count':'NOT IDENTIFIED: RAW DQ zero is not a CR detection/inventory; no CR rejection was run in this experiment.','aperture_summary':apgroup,'no_attribution_or_ERR_rescaling':True}
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('normal_search_p0','aperture_summary')},indent=2))
