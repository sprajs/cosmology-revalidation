"""Compare all numerical starts and profiles without choosing a scientific shift sign."""
import json,csv
from pathlib import Path
import numpy as np
import paired_refit as p
from check_results import table,gate
O=p.OUT

def mask(d,cid):
    return sorted((r['BAND'],r['MJD'],r['FLUXCAL'],r['FLUXCAL_ERR']) for r in table(d/'fit.LCPLOT.TEXT','OBS:') if r['CID']==cid and r['DATAFLAG']=='1')

def main():
    summaries=[];details=[]
    for timing in ['free','fixed']:
        for arm in ['R','A','B','H']:
            ds=[O/'fits/full'/timing/arm]+[O/'fits/multistart'/str(s)/timing/arm for s in [.85,1.15]]
            rr=[{r['CID']:r for r in gate(d,arm)} for d in ds]
            for cid in rr[0]:
                rows=[x[cid] for x in rr];masks=[mask(d,cid) for d in ds]
                same=all(x==masks[0] for x in masks[1:]);chi=np.array([r['FITCHI2'] for r in rows]);best=int(np.argmin(chi))
                out=dict(CID=cid,timing=timing,arm=arm,same_accepted_mask=same,best_start_by_objective=[1.,.85,1.15][best],objective_improvement=float(chi[0]-chi.min()),objective_range=float(np.ptp(chi)),
                    DLMAG_range=float(np.ptp([r['DLMAG'] for r in rows])),AV_range=float(np.ptp([r['AV'] for r in rows])),STRETCH_range=float(np.ptp([r['STRETCH'] for r in rows])),PKMJD_range=float(np.ptp([r['PKMJD'] for r in rows])),
                    all_ERRFLAG_zero=all(r['ERRFLAG_FIT']==0 for r in rows),all_positive_covariance=all(r['covariance_min_eigenvalue']>0 for r in rows),all_inside_grid=all(r['within_stretch_grid'] for r in rows),
                    comparable_objective=same,conditionally_preferred_DLMAG=rows[best]['DLMAG'] if same else '',conditionally_preferred_AV=rows[best]['AV'] if same else '')
                summaries.append(out)
                for start,r in zip([1.,.85,1.15],rows):details.append(dict(start=start,**r))
    p.csvsave(O/'multistart-comparison.csv',summaries);p.csvsave(O/'multistart-fit-ledger.csv',details)
    profiles=[];pr=json.loads((O/'stability-protocol.json').read_text())
    for j in pr['jobs']:
        if j['kind']!='profile':continue
        d=(p.ROOT/j['nml']).parent
        if not(d/'fit.FITRES.TEXT').exists() and not(d/'fit').exists():continue
        r=table(d/'fit.FITRES.TEXT')[0];cid=r['CID'];original=O/'fits/full/free/B'
        profiles.append(dict(CID=cid,offset_cells=j['offset_cells'],fixed_shape=j['fixed_shape'],same_mask_as_original=mask(d,cid)==mask(original,cid),**{k:float(r[k]) for k in ['STRETCH','AV','DLMAG','PKMJD','NDOF','FITCHI2','ERRFLAG_FIT']}))
    if profiles:p.csvsave(O/'shape-profile.csv',profiles)
    # Range of paired effects across same-start solutions, not an uncertainty interval.
    pairs=[]
    lut={(r['start'],r['timing'],r['arm'],r['CID']):r for r in details}
    for timing in ['free','fixed']:
        for cid in sorted({r['CID'] for r in details}):
            for before,after,label in [('A','B','B_minus_A'),('R','H','H_minus_R')]:
                values=[lut[s,timing,after,cid]['DLMAG']-lut[s,timing,before,cid]['DLMAG'] for s in [1.,.85,1.15]]
                pairs.append(dict(CID=cid,timing=timing,comparison=label,initial_delta=values[0],start085_delta=values[1],start115_delta=values[2],minimum_delta=min(values),maximum_delta=max(values),range_delta=max(values)-min(values)))
    p.csvsave(O/'paired-start-sensitivity.csv',pairs)
    out=dict(status='Numerical stability audit, no assertion of global convergence',N_native_fits=len(details),N_same_mask=sum(r['same_accepted_mask'] for r in summaries),N_total_cases=len(summaries),
        max_parameter_range={k:max(r[k+'_range'] for r in summaries) for k in ['DLMAG','AV','STRETCH','PKMJD']},
        objective_improvement_same_mask_max=max(r['objective_improvement'] for r in summaries if r['same_accepted_mask']),
        N_DLMAG_range_over_0p001=sum(r['DLMAG_range']>.001 for r in summaries),N_DLMAG_range_over_0p01=sum(r['DLMAG_range']>.01 for r in summaries),
        N_profile_points=len(profiles),all_ERRFLAG_zero=all(r['all_ERRFLAG_zero'] for r in summaries),all_positive_covariance=all(r['all_positive_covariance'] for r in summaries),
        most_unstable=sorted(summaries,key=lambda r:r['DLMAG_range'],reverse=True)[:12],protocol_sha256=p.sha(O/'stability-protocol.json'),source_sha256=p.sha(Path(__file__)))
    p.save(O/'stability-result.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()
