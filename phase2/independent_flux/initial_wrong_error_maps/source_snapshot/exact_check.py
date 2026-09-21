from engine import *


def stats(x):
    x=np.asarray(x); good=np.isfinite(x)
    return {'median':float(np.median(x[good])), 'abs_p95':float(np.quantile(abs(x[good]),.95)),
            'abs_max':float(np.max(abs(x[good]))), 'n_nonfinite':int(sum(~good))}


def main():
    engine=Engine(); rows=[]; diagnostics=[]
    obs,_,_=read_inputs(); accepted=mask_observations(obs)
    for fn in sorted(BUNDLE.glob('objective_*.npz')):
        cid=fn.stem.split('_')[-1]; q=np.load(fn)
        p=q['parameters_x0_x1_c_t0'].copy(); p[0]=np.log(p[0])
        z,ebv=float(q['zHEL'][0]),float(q['MWEBV'][0]); b,t=q['band'],q['MJD']
        gp=engine.prepare(z,ebv,10.); gf=engine.prepare(z,ebv,1.)
        sf,sc=engine.sncosmo_flux(p,b,t,z,ebv,True)
        f,mc,snake,cd=engine.model_covariance(p,b,t,z,gp)
        modelmagerr=np.sqrt(np.diag(mc))/abs(f)*2.5/np.log(10)
        sncosmagerr=np.sqrt(np.diag(sc))/abs(sf)*2.5/np.log(10)
        fine=engine.flux(p,b,t,z,gf,interpolation='sncosmo')
        scipy=engine.flux(p,b,t,z,gf)
        covariance=q['frozen_flux_covariance']
        residual=q['data_flux']-q['model_flux']; residual_ind=q['data_flux']-f
        ownrows=[]
        for band,time in zip(b,t):
            a=obs[(obs.CID==cid)&(obs.BAND==band)]
            assert abs(a.MJD-time).min()<1e-7
            ownrows.append(a.loc[abs(a.MJD-time).idxmin(),'source_phot_row_one_based'])
        assert set(ownrows)==set(accepted[accepted.CID==cid].source_phot_row_one_based)
        d=pd.DataFrame({'CID':cid,'source_phot_row_one_based':ownrows,'BAND':b,'MJD':t,'phase':q['rest_phase'],
                        'snana_flux':q['model_flux'],'sncosmo_flux':sf,'direct_sncosmo_interp_10A':f,
                        'direct_sncosmo_interp_1A':fine,'direct_scipy_1A':scipy,'data_fluxerr':q['data_fluxerr'],
                        'snana_model_magerr':q['model_magerr'],'independent_formula_model_magerr':modelmagerr,
                        'sncosmo_model_magerr':sncosmagerr})
        for k in ['sncosmo_flux','direct_sncosmo_interp_10A','direct_sncosmo_interp_1A','direct_scipy_1A']:
            d[k+'_fracdiff']=(d[k]-d.snana_flux)/d.snana_flux
            d[k+'_sigma_difference']=(d[k]-d.snana_flux)/d.data_fluxerr
        d['model_magerr_formula_fracdiff']=modelmagerr/q['model_magerr']-1
        d['model_magerr_sncosmo_fracdiff']=sncosmagerr/q['model_magerr']-1
        rows.append(d)
        diagnostics.append({'CID':cid,'n':len(t),'model_covariance_diag_ratio_snana_to_sncosmo':float(np.median((np.diag(covariance)-q['data_fluxerr']**2)/np.diag(sc))),
                            'snana_data_chi2':float(residual@np.linalg.solve(covariance,residual)),
                            'independent_fixedparam_data_chi2':float(residual_ind@np.linalg.solve(covariance,residual_ind)),
                            'magerr_formula_ratio_median':float(np.median(modelmagerr/q['model_magerr']))})
    result=pd.concat(rows,ignore_index=True)
    numeric=[c for c in result if 'fracdiff' in c or 'sigma_difference' in c]
    summary={c:stats(result[c]) for c in numeric}
    summary['n']=len(result);summary['per_object']=diagnostics
    summary['mask_agreement']='Exact exported fit rows equal published accepted rows for all 12 objects'
    result.to_csv(OUT/'fixed_predictions_exact.csv.gz',index=False,float_format='%.17g')
    (OUT/'fixed_summary_exact.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    provenance('exact_check',[OUT/'fixed_predictions_exact.csv.gz',OUT/'fixed_summary_exact.json'],[OUT/'PLAN.json'])
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
