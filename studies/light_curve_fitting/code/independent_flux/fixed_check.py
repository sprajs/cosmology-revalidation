"""First independent fixed-parameter check; retain disagreements in full."""
from engine import *


def main():
    engine=Engine(); obs,fit,pars=read_inputs()
    accepted=mask_observations(obs)
    accepted.to_csv(OUT/'published_mask_and_folds.csv',index=False,float_format='%.17g')
    pred=pd.read_csv(BUNDLE/'snana_predictions_and_mask.csv',dtype={'CID':str})
    pred=pred[pred.DATA_MODEL==0]
    rows=[]; diagnostics=[]
    for cid,par in pars.items():
        x0,x1,c,t0=par['parameters_x0_x1_c_t0']; p=np.array([np.log(x0),x1,c,t0])
        z,ebv=fit.loc[cid,['zHEL','MWEBV']]
        d=pred[pred.CID==cid].sort_values(['MJD','BAND']).copy()
        d=d[(d.MJD-t0)/(1+z)>=-15]
        d=d[(d.MJD-t0)/(1+z)<=45]
        grids={s:engine.prepare(z,ebv,s) for s in [10.,5.,2.,1.]}
        b,t=d.BAND.values,d.MJD.values
        d['phase']=(t-t0)/(1+z)
        d['sncosmo']=engine.sncosmo_flux(p,b,t,z,ebv)
        for s,g in grids.items():
            d[f'direct_scipy_{s:g}A']=engine.flux(p,b,t,z,g)
        d['direct_sncosmo_interp_10A']=engine.flux(p,b,t,z,grids[10.],interpolation='sncosmo')
        d['direct_sncosmo_interp_1A']=engine.flux(p,b,t,z,grids[1.],interpolation='sncosmo')
        for name in ['sncosmo','direct_scipy_10A','direct_sncosmo_interp_10A']:
            d[name+'_fracdiff']=(d[name]-d.FLUXCAL)/d.FLUXCAL.replace(0,np.nan)
        d['integration_1_vs_2A']=(d['direct_scipy_1A']-d['direct_scipy_2A'])/d['direct_scipy_1A']
        d['integration_1_vs_10A']=(d['direct_scipy_1A']-d['direct_scipy_10A'])/d['direct_scipy_1A']
        d['interpolator_difference']=(d['direct_scipy_1A']-d['direct_sncosmo_interp_1A'])/d['direct_scipy_1A']
        d['sncosmo_direct_integration_difference']=(d.sncosmo-d['direct_sncosmo_interp_1A'])/d.sncosmo
        rows.append(d)
        a=accepted[accepted.CID==cid]
        af,acov=engine.sncosmo_flux(p,a.BAND.values,a.MJD.values,z,ebv,covariance=True)
        ff,fcov,snake,cd=engine.model_covariance(p,a.BAND.values,a.MJD.values,z,grids[10.])
        np.savez_compressed(OUT/f'initial_covariance_{cid}.npz',sncosmo=acov,independent_snana_formula=fcov,
                            flux_sncosmo=af,flux_direct=ff,rows=a.source_phot_row_one_based.values,
                            snake_relative_variance=snake,colour_dispersion=cd)
        diagnostics.append({'CID':cid,'n_accepted':len(a),'n_heldout':int(sum(a.epoch_fold==0)),
                            'model_sigma_ratio_snana_formula_to_sncosmo_median':float(np.median(np.sqrt(np.diag(fcov)/np.diag(acov))))})
    rows=pd.concat(rows,ignore_index=True)
    rows.to_csv(OUT/'fixed_predictions_coarse.csv.gz',index=False,float_format='%.17g')
    numeric=[c for c in rows if 'fracdiff' in c or 'difference' in c or c.startswith('integration_')]
    summary={c:{'median':float(rows[c].median()),'abs_p95':float(rows[c].abs().quantile(.95)),
                'abs_max':float(rows[c].abs().max())} for c in numeric}
    summary['n_model_points']=len(rows)
    summary['zero_reference_flux_points']=rows[rows.FLUXCAL==0][['CID','BAND','phase','FLUXCAL','sncosmo']].to_dict('records')
    summary['per_object']=diagnostics
    summary['per_band']={str(b):{c:float(d[c].median()) for c in numeric} for b,d in rows.groupby('BAND')}
    (OUT/'fixed_summary_coarse.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    provenance('fixed_check',list(OUT.glob('initial_covariance_*.npz'))+[OUT/'fixed_summary_coarse.json',OUT/'fixed_predictions_coarse.csv.gz',OUT/'published_mask_and_folds.csv'],[OUT/'PLAN.json'])
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
