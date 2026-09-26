"""Fixed-cohort DES/PS1 stellar extinction branch sensitivity, not a release refit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'sources/repos/bap37__Dovekie'


def design(x,degree):
    return np.column_stack([(x-.6)**j for j in range(degree+1)])


def main():
    out=ROOT/'runs/research_2026_09_26/calibration_stars'
    out.mkdir(exist_ok=False)
    source=Path(__file__).read_bytes();(out/'executed_source.py').write_bytes(source)
    plan=(ROOT/'docs/research-2026-09-26/plan.md').read_bytes();(out/'protocol.md').write_bytes(plan)
    obsfile=BASE/'output_observed_apermags+AV/DES_observed.csv'
    d=pd.read_csv(obsfile);d['row_id']=np.arange(len(d))
    assert not d[['RA','DEC']].duplicated().any()
    cells=[f'{int(np.floor(r/2))}:{int(np.floor((c+90)/2))}' for r,c in zip(d.RA,d.DEC)]
    d['cell']=cells;d['test']=[int(hashlib.sha256(c.encode()).hexdigest()[:8],16)%2==1 for c in cells]
    catalogs={};cuts={};variables={}
    for name,factor in [('stored_extinction',1),('zero_extinction',0)]:
        m={s+'-'+b:d[s+'-'+b].to_numpy()-factor*d[s+'-'+b+'_AV'].to_numpy() for s in ['DES','PS1'] for b in 'griz'}
        x=m['PS1-g']-m['PS1-i'];cut=(x>.25)&(x<1)
        for b,limit in zip('griz',[14.3,14.4,14.6,14.1]):cut&=m['PS1-'+b]>limit
        cuts[name]=cut;variables[name]=(x,m)
    common=cuts['stored_extinction']&cuts['zero_extinction']
    synfiles=[BASE/f'output_synthetic_magsaper/synth_{s}_shift_0.000.txt' for s in ['DES','PS1']]
    sd,sp=[pd.read_csv(p,sep=r'\s+') for p in synfiles]
    assert not sd.standard.duplicated().any() and not sp.standard.duplicated().any()
    s=sd.merge(sp,on=['standard','standard_catagory'],validate='one_to_one',suffixes=('_DES','_PS1'))
    sx=s['PS1-g']-s['PS1-i']
    scut=(sx>.25)&(sx<1.6)&((s['PS1-g']-s['PS1-r'])>.25)
    synthrows=[]
    for cat in ['stis_ngsl_v2','calspec23']:
        mask=scut&(s.standard_catagory==cat)
        for b in 'griz':
            xx=sx[mask].to_numpy();yy=(s['DES-'+b]-s['PS1-'+b])[mask].to_numpy()
            # Source one-pass clipping fits before its cut; retain all source-cut
            # standards for this explicitly independent descriptive sensitivity.
            coef=np.linalg.lstsq(design(xx,1),yy,rcond=None)[0]
            catalogs[cat,b]=coef
            synthrows.append(dict(library=cat,band=b,n=int(mask.sum()),intercept_pivot06=coef[0],slope=coef[1],
                rms=float(np.sqrt(np.mean((yy-design(xx,1)@coef)**2)))))
    rows=[];pred=[];offsetrows=[]
    for arm in ['fixed_intersection','own_cuts']:
        for mode,(x,m) in variables.items():
            for b in 'griz':
                y=m['DES-'+b]-m['PS1-'+b]
                mask=(common if arm=='fixed_intersection' else cuts[mode])&(abs(d['DES-'+b]-d['PS1-'+b])<1)&np.isfinite(x)&np.isfinite(y)
                for degree in [1,2]:
                    X=design(x,degree);coef=np.linalg.lstsq(X[mask],y[mask],rcond=None)[0]
                    train=mask&~d.test.to_numpy();test=mask&d.test.to_numpy()
                    ct=np.linalg.lstsq(X[train],y[train],rcond=None)[0]
                    r=y[test]-X[test]@ct
                    nuisance=np.column_stack([np.ones(test.sum()),d.loc[test,'PS1-g_AV'],m['PS1-i'][test]-16])
                    nc=np.linalg.lstsq(nuisance,r,rcond=None)[0]
                    rows.append(dict(arm=arm,mode=mode,band=b,degree=degree,n=int(mask.sum()),n_train=int(train.sum()),n_test=int(test.sum()),
                        n_train_cells=d.loc[train,'cell'].nunique(),n_test_cells=d.loc[test,'cell'].nunique(),
                        intercept_pivot06=coef[0],slope=coef[1],quadratic=coef[2] if degree==2 else 0,
                        test_mean_residual=float(r.mean()),test_rms=float(np.sqrt(np.mean(r*r))),
                        test_residual_slope_per_PS1_g_extinction=float(nc[1]),test_residual_slope_per_i_mag=float(nc[2])))
                    for idx,res in zip(np.flatnonzero(test),r):
                        pred.append(dict(arm=arm,mode=mode,band=b,degree=degree,row_id=int(idx),cell=d.iloc[idx].cell,residual=float(res)))
                for cat in ['stis_ngsl_v2','calspec23']:
                    offset=design(x[mask],1)@catalogs[cat,b]-y[mask]
                    offsetrows.append(dict(arm=arm,mode=mode,band=b,library=cat,n=int(mask.sum()),
                        mean_synthetic_minus_observed=float(offset.mean()),median_synthetic_minus_observed=float(np.median(offset)),
                        scatter=float(offset.std(ddof=1))))
    for frame,name in [(pd.DataFrame(rows),'observed_fits.csv'),(pd.DataFrame(synthrows),'synthetic_fits.csv'),
                       (pd.DataFrame(offsetrows),'synthetic_offset_estimates.csv'),(pd.DataFrame(pred),'heldout_residuals.csv')]:
        frame.to_csv(out/name,index=False)
    offsets=pd.DataFrame(offsetrows);changes=[]
    for (arm,cat,b),p in offsets.groupby(['arm','library','band']):
        p=p.set_index('mode')
        changes.append(dict(arm=arm,library=cat,band=b,
            delta_mean_correction_stored_minus_zero=float(p.loc['stored_extinction','mean_synthetic_minus_observed']-p.loc['zero_extinction','mean_synthetic_minus_observed'])))
    changes=pd.DataFrame(changes);changes.to_csv(out/'extinction_branch_changes.csv',index=False)
    d['selected_both']=common
    for name,cut in cuts.items():d['selected_'+name]=cut
    d.to_csv(out/'star_selection.csv',index=False)
    report={'scope':'Descriptive DES-minus-PS1 sensitivity to stored stellar extinction vs forced-zero branch. Not a reproduction of calibration MCMC or evidence which branch published release executed.',
        'source_rows':len(d),'common_cut_rows':int(common.sum()),'own_cut_counts':{k:int(v.sum()) for k,v in cuts.items()},
        'selection_limit':'Input AV file was already preselected; zero-extinction cannot recover excluded raw catalog stars.',
        'independence_limit':'Previously used calibration stars; internal sky holdout shares survey zero points, PS1, dust map and stellar SED assumptions. No per-star errors or release execution provenance supplied.',
        'fixed_cohort_changes':changes[changes.arm=='fixed_intersection'].to_dict('records'),
        'synthetic_cut_counts':pd.DataFrame(synthrows).groupby('library').n.first().to_dict(),
        'source_pipeline_differences':'Independent unweighted linear/quadratic fits; no source iterative residual clipping, MCMC prior, white dwarf term, filter shifts, or PS1 anchor offsets. All absolute calibration estimates remain conditional.',
        'rounding':'Stored AV file has ~6-significant-figure magnitudes; do not interpret sub-rounding effects.'}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    inputs=[obsfile,*synfiles,BASE/'dovekie.py',BASE/'scripts/helpers.py',BASE/'DOVEKIE_DEFS.yml']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source).hexdigest(),
        'protocol_sha256':hashlib.sha256(plan).hexdigest(),'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs}},indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
