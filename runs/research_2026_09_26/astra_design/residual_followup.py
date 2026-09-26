"""Inspect coefficient responses and independently verify predictive score algebra."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular,eigvalsh
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent
DATA=ROOT/'runs/salt_dust_audit/flux_response';x=np.load(DATA/'matrices.npz')
meta=pd.read_csv(DATA/'selected_objects.csv',dtype={'CID':str}).set_index('CID')
gauge=np.array([[1,0,0],[-1,-1,-1],[0,1,0],[0,0,1]],float)
report={'samples':{},'exact_covariance_overlap':{}};blocks=[];all_responses=[];saves={}
for label,p in [('all64',P),('highIa43',P/'residual43')]:
    a=np.load(p/'residual64-arrays.npz');saved=json.loads((p/'residual64.json').read_text())
    objects=pd.read_csv(p/'residual64-objects.csv',dtype={'CID':str})
    m=objects[objects.weighting=='measurement_plus_model'].reset_index(drop=True)
    per=[];checks={};coefs={}
    groups={'object':[(row.CID,np.array([i])) for i,row in enumerate(m.itertuples())],
            'field':[(field,np.flatnonzero(m.field==field)) for field in sorted(m.field.unique())]}
    for arm in ['measurement_plus_model','measurement_only']:
        F=a[arm+'_F'];u=a[arm+'_u'];Fa=F.sum(0);ua=u.sum(0)
        for family,ii in [('observer',np.arange(3)),('rest_phase',np.arange(3,6))]:
            for sigma in [.01,.02,.05]:
                prior=np.eye(3)/sigma**2;Vall=np.linalg.inv(prior+Fa[np.ix_(ii,ii)])
                for cv,grp in groups.items():
                    total=0.
                    for group,ids in grp:
                        fg=F[ids].sum(0)[np.ix_(ii,ii)];ug=u[ids].sum(0)[ii]
                        V=np.linalg.inv(prior+Fa[np.ix_(ii,ii)]-fg)
                        mean=V@(ua[ii]-ug);res=ug-fg@mean
                        # Direct low-rank conditional predictive quadratic.
                        value=mean@ug-.5*mean@fg@mean+.5*res@Vall@res-.5*np.linalg.slogdet(np.eye(3)+V@fg)[1]
                        total+=value
                        if family=='observer' and sigma==.02 and cv=='field':
                            offsets=gauge@mean
                            blocks.append({'sample':label,'weighting':arm,'field':group,'score':value,
                               **dict(zip('griz',offsets))})
                    name=f'{cv}_{family}_{sigma:g}'
                    checks[arm+'_'+name]=float(total-saved['arms'][arm]['candidate_scores'][name])
        assert max(abs(v) for v in checks.values())<1e-9
        sigma=.02;V=np.linalg.inv(np.eye(3)/sigma**2+Fa[:3,:3]);mean=V@ua[:3]
        offsets=gauge@mean
        coefs[arm]={'basis_mean':mean.tolist(),'basis_covariance':V.tolist(),'gray_gauge':'sum(griz)=0',
          'band_order':'griz','band_offset_mean_mag':offsets.tolist(),'band_offset_sd_mag':np.sqrt(np.diag(gauge@V@gauge.T)).tolist(),
          'qualification':'Diagnostic regularized shared residual coefficient, not independently measured calibration.'}
        maps=[]
        for row in m.itertuples():
            R=x[row.CID+'__'+arm+'_response']
            mapping=-R[:,5:9]@gauge
            delta=mapping@mean
            fixed_standard=np.array([1,.16087,-3.11780,0])@delta
            all_responses.append({'sample':label,'weighting':arm,'CID':row.CID,'zHEL':row.zHEL,
              'field':row.field,'pIa':row.pIa,'delta_mB':delta[0],'delta_x1':delta[1],
              'delta_c':delta[2],'delta_t0':delta[3],'delta_fixed_reference_preBBC_standardized_mag':fixed_standard})
            maps.append(mapping)
        saves[label+'_'+arm+'_parameter_response_map']=np.array(maps)
        saves[label+'_'+arm+'_coefficient_covariance']=V
        rr=pd.DataFrame([r for r in all_responses if r['sample']==label and r['weighting']==arm])
        coefs[arm]['response_summary']={key:{'median':float(rr[key].median()),'min':float(rr[key].min()),'max':float(rr[key].max())} for key in ['delta_mB','delta_x1','delta_c','delta_fixed_reference_preBBC_standardized_mag']}
        # Unweighted deterministic sample quartiles are descriptive, not a
        # selected-population evolution estimate or cosmological correction.
        bins=pd.qcut(rr.zHEL,4,labels=False)
        coefs[arm]['descriptive_redshift_quartiles']=[{'n':int((bins==i).sum()),'mean_z':float(rr.loc[bins==i,'zHEL'].mean()),'mean_fixed_reference_preBBC_delta':float(rr.loc[bins==i,'delta_fixed_reference_preBBC_standardized_mag'].mean())} for i in range(4)]
    report['samples'][label]={'coefficient_response':coefs,'max_independent_score_error':max(abs(v) for v in checks.values()),'score_differences':checks}

for cid in ['1280240','1339609']:
    q=np.load(ROOT/f'phase2/official/portable_pilot/objective_{cid}.npz');p=cid+'__'
    dt=abs(q['MJD'][:,None]-x[p+'mjd'][None,:])
    df=abs(q['data_flux'][:,None]-x[p+'flux_observed'][None,:])
    cost=dt/.003+df/np.maximum(q['data_fluxerr'][:,None]*.01,.001)
    cost+=(q['band'][:,None]!=x[p+'band'][None,:])*1e9
    aa,bb=linear_sum_assignment(cost);assert np.array_equal(aa,np.arange(len(aa)))
    assert dt[aa,bb].max()<.003 and df[aa,bb].max()<.005
    actual=q['frozen_flux_covariance'];L=np.linalg.cholesky(actual)
    modeled=x[p+'measurement_plus_model_covariance'][np.ix_(bb,bb)]
    eig=eigvalsh(modeled,actual)
    w=solve_triangular(L,modeled-actual,lower=True)
    w=solve_triangular(L,w.T,lower=True).T
    report['exact_covariance_overlap'][cid]={'epochs':len(aa),'one_to_one_epoch_match':True,
       'max_mjd_difference':float(dt[aa,bb].max()),'max_flux_difference':float(df[aa,bb].max()),
       'model_to_exact_generalized_eigenvalue_range':[float(eig.min()),float(eig.max())],
       'whitened_covariance_difference_frobenius':float(np.linalg.norm(w)),
       'qualification':'Only two overlap objects; native covariance differs from exact SNANA and does not validate64 population noise.'}
pd.DataFrame(all_responses).to_csv(P/'calibration-shaped-responses.csv',index=False)
pd.DataFrame(blocks).to_csv(P/'calibration-shaped-field-coefficients.csv',index=False)
np.savez_compressed(P/'calibration-shaped-shared-response.npz',**saves)
report['response_qualification']='Signed response to REMOVING the fitted observer template from observed flux, with fixed masks, covariance, SALT training and no BBC/selection. Raw parameter deltas are primary; standardized combinations hold reference alpha=.16087,beta=3.11780 and are not complete distance derivatives. Shared covariance is map*V*map^T, not independent diagonal noise.'
coefrows=[]
for sample,v in report['samples'].items():
    for weighting,c in v['coefficient_response'].items():
        for band,mean,sd in zip('griz',c['band_offset_mean_mag'],c['band_offset_sd_mag']):
            coefrows.append({'sample':sample,'weighting':weighting,'band':band,
              'fitted_residual_equivalent_dimming_mag':mean,'conditional_regularized_sd_mag':sd,
              'data_correction_dimming_if_interpreted_as_calibration':-mean,
              'gray_gauge':'sum_griz_zero','ridge_scale_mag':.02})
pd.DataFrame(coefrows).to_csv(P/'calibration-shaped-coefficients.csv',index=False)
report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(P/'residual-followup.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'exact_covariance_overlap':report['exact_covariance_overlap'],
 'summaries':{k:{'max_score_error':v['max_independent_score_error'],'coefficients':{arm:vv['band_offset_mean_mag'] for arm,vv in v['coefficient_response'].items()}} for k,v in report['samples'].items()}},indent=2))
