"""Descriptive branch and amplitude-inclusive level geometry; no CI labels."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
O=Path(__file__).resolve().parent;C=O/'fixed-c-profile/cohort10'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def analyze(d):
    result=json.loads((d/'result.json').read_text());arr=np.load(d/'native-profiles.npz',allow_pickle=False)
    records=list(csv.DictReader((d/'shape-profiles.csv').open()));coords=arr['coordinates'];means=arr['model_means'];vals=arr['metric_results']
    anchorCs=[arr['covariance_B'],arr['covariance_A_anchor']];metrics={};indices=[arr['A_indices'],np.arange(len(arr['data_flux']))]
    Dref=arr['reference_parameters'][0];knots=np.array(json.loads((d/'protocol.json').read_text())['domain']['shape_knots'])
    for k,label in enumerate(['Banchor_A','Banchor_B','Aanchor_A','Aanchor_B']):
      if label not in result['metrics']:continue
      r=[{key:(value if key in ['metric','stage'] else float(value)) for key,value in row.items()} for row in records if row['stage']=='fine' and row['metric']==label]
      if not r:continue
      minima=[]
      for j in range(len(r)):
        if (j==0 or r[j]['Q']<=r[j-1]['Q']) and(j==len(r)-1 or r[j]['Q']<=r[j+1]['Q']):minima.append(r[j])
      minimum=result['metrics'][label]['finest']['Q'];ix=indices[k%2];Cv=anchorCs[k//2][np.ix_(ix,ix)];W=np.linalg.inv(Cv)
      h=means[:,ix];q=np.einsum('ni,ij,nj->n',h,W,h,optimize=True)
      profQ=vals[:,k,0];ahat=vals[:,k,2];levelsets={}
      for threshold in [1.,4.,9.]:
        valid=np.isfinite(profQ)&np.isfinite(ahat)&(q>0)&(ahat>0)&(profQ<=minimum+threshold)
        if not np.any(valid):continue
        da=np.sqrt(np.maximum(0,minimum+threshold-profQ[valid])/q[valid]);alo=ahat[valid]-da;ahi=ahat[valid]+da
        # Keep the declared distance domain and native magnitude-validity domain explicit.
        alo_domain=10**(-.4*(60-Dref));ahi_domain=10**(-.4*(10-Dref))
        clipped_distance=bool(np.any(alo<=alo_domain) or np.any(ahi>=ahi_domain))
        native_alo=np.max(1e-5/means[valid],axis=1);native_ahi=np.min(1e9/means[valid],axis=1)
        clipped_native=bool(np.any(alo<=native_alo) or np.any(ahi>=native_ahi))
        alo=np.maximum(np.maximum(alo,alo_domain),native_alo);ahi=np.minimum(np.minimum(ahi,ahi_domain),native_ahi)
        Dlo=Dref-2.5*np.log10(ahi);Dhi=Dref-2.5*np.log10(alo)
        selected=coords[valid];ridge=[x for x in r if x['Q']<=minimum+threshold]
        levelsets[str(threshold)]=dict(amplitude_inclusive_distance_range=[float(Dlo.min()),float(Dhi.max())],
          shape_profile_ridge_distance_range=[min(x['DLMAG'] for x in ridge),max(x['DLMAG'] for x in ridge)],
          finite_native_grid_coordinates=int(valid.sum()),shape_box_touched=bool(np.any(abs(selected[:,0]-knots[0])<5e-7) or np.any(abs(selected[:,0]-knots[-1])<5e-7)),
          AV_box_touched=bool(np.any(abs(selected[:,1]+1)<1e-9) or np.any(abs(selected[:,1]-2)<1e-9)),distance_box_touched=clipped_distance,native_magnitude_box_touched=clipped_native,
          interpretation='Union over saved native coordinates plus analytic positive-amplitude freedom; descriptive restricted-objective level set, not confidence/posterior interval.')
      metrics[label]=dict(minimum_Q=minimum,cached_minimum_Q=float(np.nanmin(profQ)),local_modes=[dict(x,delta_Q=x['Q']-minimum) for x in minima],level_sets=levelsets)
    paired={}
    for anchor in ['B','A']:
      a=anchor+'anchor_A';b=anchor+'anchor_B'
      if a not in metrics or b not in metrics:continue
      ap=result['metrics'][a]['finest'];candidates=metrics[b]['local_modes']
      bp=min(candidates,key=lambda r:abs(r['shape']-ap['shape']))
      paired[anchor]=dict(global_minimum_delta_DLMAG=result['delta_DLMAG_B_minus_A'][anchor],
        matched_B_mode_nearest_A_shape=bp,matched_mode_delta_DLMAG=bp['DLMAG']-ap['DLMAG'],matching_rule='Nearest shape to A minimum among all fine-grid local modes/endpoints, independent of shift sign.')
    out=dict(CID=d.name,numerical_gate_pass=result['numerical_gate_pass'],metrics=metrics,paired=paired,
      source_sha256=sha(__file__),input_hashes={f:sha(d/f) for f in ['result.json','native-profiles.npz','shape-profiles.csv','protocol.json']})
    (d/'branch-level-geometry.json').write_text(json.dumps(out,indent=2)+'\n');return out

def main():
    protocol=json.loads((C/'protocol.json').read_text());out={};rows=[]
    for cid in protocol['cohort']:
      d=C/cid
      if not(d/'result.json').exists():rows.append(dict(CID=cid,status='no profile result'));continue
      r=analyze(d);out[cid]=r;main=r['paired'].get('B',{})
      rows.append(dict(CID=cid,status='profile saved',numerical_gate_pass=r['numerical_gate_pass'],
        global_minimum_delta=main.get('global_minimum_delta_DLMAG'),matched_mode_delta=main.get('matched_mode_delta_DLMAG'),
        matched_mode_delta_Q=main.get('matched_B_mode_nearest_A_shape',{}).get('delta_Q')))
    (C/'branch-summary.json').write_text(json.dumps(dict(cases=rows,results=out,source_sha256=sha(__file__)),indent=2)+'\n')
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
