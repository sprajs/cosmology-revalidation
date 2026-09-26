"""Previously frozen fixed-template transfer, exact native mean/frozen C."""
from pathlib import Path
import sys,json,hashlib,itertools,importlib.util
from types import SimpleNamespace
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import linear_sum_assignment
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).resolve().parent
OUT=BASE/'validation1020';DEST=OUT/'analysis'
spec=importlib.util.spec_from_file_location('flux_response',ROOT/'scripts/salt_dust_audit/flux_response.py');fr=importlib.util.module_from_spec(spec);spec.loader.exec_module(fr)

def project(C,J,f,y,b,rest):
    L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True)
    U,s,V=np.linalg.svd(jw,full_matrices=True);rank=int(np.sum(s>s[0]*1e-10));assert rank==4
    Q=U[:,4:];r=Q.T@solve_triangular(L,y-f,lower=True)
    G=np.column_stack([-fr.K*f*(b==k) for k in 'griz'])
    obs=np.column_stack([G[:,0]-G[:,1],G[:,2]-G[:,1],G[:,3]-G[:,1]])
    T=Q.T@solve_triangular(L,np.column_stack([obs,rest]),lower=True)
    assert max(abs(Q.T@solve_triangular(L,-fr.K*f,lower=True)))<1e-6
    return r,T,float(s[0]/s[-1])

def joint_predictive(u,F,c,V):
    W=np.linalg.inv(np.linalg.inv(V)+F);delta=np.asarray(u)-F@c
    return np.asarray(u)@c-.5*c@F@c+.5*np.einsum('...i,ij,...j->...',delta,W,delta)-.5*np.linalg.slogdet(np.eye(3)+V@F)[1]

def main():
    if DEST.exists():raise RuntimeError('Preserve existing analysis')
    DEST.mkdir();(DEST/'objects').mkdir()
    m=pd.read_csv(OUT/'cohort.csv',dtype={'CID':str});assert len(m)==1020
    dc=np.load(OUT/'frozen-discovery-coefficients.npz');c=dc['basis_mean'];V=dc['basis_covariance']
    dr=np.load(OUT/'frozen-rest-discovery-coefficients.npz');cr=dr['basis_mean'];Vr=dr['basis_covariance']
    points=pd.read_csv(fr.RELEASE/'0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz',sep=r'\s+',comment='#',usecols=['SNID','DATA_MODEL','BAND','MJD','FLUXCAL','FLUXCALERR'],dtype={'SNID':str})
    points=points[(points.DATA_MODEL==1)&points.SNID.isin(m.CID)];pointgroups={str(k):v for k,v in points.groupby('SNID',sort=False)}
    model,bands,paths,zp=fr.build_model();magoff={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    summaries=[];records=[];sufficient={arm:dict(F=[],u=[]) for arm in ['published_mask','first_exact_duplicate']};native_check=[]
    for count,record in enumerate(m.to_dict('records')):
        cid=record['CID'];x=np.load(OUT/'objectives'/f'objective_{cid}.npz');published=pointgroups[cid]
        n=len(x['MJD']);assert n==len(published)==record['expected_epochs']
        cost=abs(x['MJD'][:,None]-published.MJD.to_numpy()[None,:])+(x['band'][:,None]!=published.BAND.to_numpy()[None,:])*1e6+abs(x['data_flux'][:,None]-published.FLUXCAL.to_numpy()[None,:])*1e-4
        aa,bb=linear_sum_assignment(cost);assert np.array_equal(aa,np.arange(n))
        dt=abs(x['MJD']-published.MJD.to_numpy()[bb]);df=abs(x['data_flux']-published.FLUXCAL.to_numpy()[bb])
        assert np.array_equal(x['band'],published.BAND.to_numpy()[bb]) and dt.max()<=.005 and df.max()<.001,(cid,dt.max(),df.max())
        assert np.ptp(x['zHEL'])==np.ptp(x['MWEBV'])==0
        record.update(zip(['x0','x1','c','PKMJD'],x['parameters_x0_x1_c_t0']));record['zHEL']=float(x['zHEL'][0]);row=SimpleNamespace(**record)
        obs=pd.DataFrame({'MJD':x['MJD'],'BAND':x['band'],'FLUXCAL':x['data_flux'],'FLUXCALERR':x['data_fluxerr']});order=obs.sort_values(['MJD','BAND'],kind='stable').index.to_numpy()
        _,a,_,diag=fr.analyze(row,obs,{'MWEBV':x['MWEBV'][0]},model,bands,magoff,False)
        f=x['model_flux'][order];y=x['data_flux'][order];b=x['band'][order];C=x['frozen_flux_covariance'][np.ix_(order,order)]
        J=a['jacobian_flux'].copy();J[:,0]=-fr.K*f
        phase=(a['mjd']-row.PKMJD)/(1+row.zHEL);G=a['nuisance_flux']
        rest=np.column_stack([G[:,2],G[:,1],J[:,2]*np.tanh(phase/20)])
        key=obs.iloc[order].reset_index(drop=True);dup=key.duplicated(subset=['MJD','BAND','FLUXCAL','FLUXCALERR'],keep='first').to_numpy();keep=np.flatnonzero(~dup)
        ngr=int(key.groupby(['MJD','BAND','FLUXCAL','FLUXCALERR'],dropna=False).size().gt(1).sum())
        save=dict(MJD=a['mjd'],band=b,official_flux_model=f,native_flux_model=a['flux_model'],observed_flux=y,quoted_error=x['data_fluxerr'][order],exact_covariance=C,jacobian_flux=J,rest_templates_raw=rest,first_duplicate_keep=keep)
        for arm,ids in [('published_mask',np.arange(n)),('first_exact_duplicate',keep)]:
            r,T,condition=project(C[np.ix_(ids,ids)],J[ids],f[ids],y[ids],b[ids],rest[ids]);u=T.T@r;F=T.T@T
            aa=float(c@u[:3]);I=float(c@F[:3,:3]@c);gain=aa-.5*I
            ar=float(cr@u[3:]);Ir=float(cr@F[3:,3:]@cr);cross=float(c@F[:3,3:]@cr)
            records.append(dict(arm=arm,CID=cid,field=row.field,zHEL=row.zHEL,epoch_count=len(ids),projected_dimension=len(r),projected_chi2=float(r@r),matched_filter=aa,information=I,fixed_prediction_gain=gain,rest_fixed_prediction_gain=ar-.5*Ir,rest_matched_filter=ar,rest_information=Ir,prediction_cross_information=cross,exact_duplicate_extra_epochs=int(dup.sum()),exact_duplicate_groups=ngr,condition=condition))
            sufficient[arm]['u'].append(u);sufficient[arm]['F'].append(F);save[arm+'_r']=r;save[arm+'_T']=T
        np.savez_compressed(DEST/'objects'/f'{cid}.npz',**save)
        native_check.append(dict(CID=cid,matched_epochs=n,max_MJD_match_difference=float(dt.max()),max_flux_match_difference=float(df.max()),native_derivative_halfstep_relative_error=diag['derivative_relative_error'],max_native_minus_official_in_measurement_sigma=float(max(abs((a['flux_model']-f)/a['flux_error'])))))
        if count%50==0:print('processed',count+1,'of',len(m),flush=True)
    rr=pd.DataFrame(records);rr.to_csv(DEST/'object-scores.csv',index=False);pd.DataFrame(native_check).to_csv(DEST/'implementation-checks.csv',index=False)
    fields=sorted(m.field.unique());assert len(fields)==10
    signs=np.array(list(itertools.product([-1.,1.],repeat=len(fields))));N=100000;seed=26092683;rng=np.random.default_rng(seed);results={};arrays={};fieldrows=[]
    for arm,a in sufficient.items():
        u6=np.array(a['u']);F6=np.array(a['F']);u=u6[:,:3];F=F6[:,:3,:3];ua=u.sum(0);Fa=F.sum(0);df=rr[rr.arm==arm]
        ur=u6[:,3:];Fr=F6[:,3:,3:];ura=ur.sum(0);Fra=Fr.sum(0)
        uf=np.array([u[m.field==f].sum(0) for f in fields]);fieldvals=[]
        for field,uu in zip(fields,uf):
            dd=df[df.field==field];r=dict(arm=arm,field=field,objects=len(dd),epochs=int(dd.epoch_count.sum()),matched_filter=float(uu@c),information=float(dd.information.sum()),fixed_prediction_gain=float(dd.fixed_prediction_gain.sum()),projected_chi2=float(dd.projected_chi2.sum()),projected_dimension=int(dd.projected_dimension.sum()),rest_matched_filter=float(dd.rest_matched_filter.sum()),rest_information=float(dd.rest_information.sum()),rest_fixed_prediction_gain=float(dd.rest_fixed_prediction_gain.sum()),prediction_cross_information=float(dd.prediction_cross_information.sum()));fieldrows.append(r);fieldvals.append(r['matched_filter'])
        aval=float(c@ua);I=float(c@Fa@c);gain=aval-.5*I;sign_a=signs@np.array(fieldvals)
        posterior_gain=float(joint_predictive(ua,Fa,c,V));sign_pp=joint_predictive(signs@uf,Fa,c,V)
        eig,vec=np.linalg.eigh(Fa);root=vec*np.sqrt(np.maximum(eig,0));simu=rng.normal(size=(N,3))@root.T;simpp=joint_predictive(simu,Fa,c,V)
        results[arm]={'fixed_template_gain':gain,'matched_filter':aval,'fixed_prediction_information':I,'conditional_Gaussian_signed_Z':aval/np.sqrt(I),'conditional_Gaussian_one_sided_tail':float(norm.sf(aval/np.sqrt(I))),
          'descriptive_transfer_amplitude':aval/I,'conditional_Gaussian_amplitude_se':1/np.sqrt(I),
          'field_sign_exact_tail':float(np.mean(sign_a>=aval-1e-10)),'field_sign_exceed_count':int(np.sum(sign_a>=aval-1e-10)),
          'projected_chi2':float(df.projected_chi2.sum()),'projected_dimension':int(df.projected_dimension.sum()),'epochs':int(df.epoch_count.sum()),
          'joint_discovery_posterior_predictive_gain':posterior_gain,'secondary_joint_predictive_field_sign_tail':float(np.mean(sign_pp>=posterior_gain-1e-10)),
          'secondary_joint_predictive_Gaussian_plus_one_tail':float((1+np.sum(simpp>=posterior_gain))/(N+1)),
          'discovery_posterior_predictive_qualification':'One shared3D discovery posterior integrated once over joint validation likelihood. Frozen conditional Gaussian noise.',
          'objects_with_duplicate_groups':int((df.exact_duplicate_groups>0).sum()),'exact_duplicate_groups':int(df.exact_duplicate_groups.sum()),'extra_exact_duplicate_epochs':int(df.exact_duplicate_extra_epochs.sum()),
          'duplicate_object_fixed_gain':float(df.loc[df.exact_duplicate_groups>0,'fixed_prediction_gain'].sum())}
        # Paired fixed forecasts use their actual covariance, not independence.
        urf=np.array([ur[m.field==f].sum(0) for f in fields]);ar=float(cr@ura);Ir=float(cr@Fra@cr)
        cross=float(c@F6.sum(0)[:3,3:]@cr);difference_I=I+Ir-2*cross;difference_a=aval-ar
        rest_gain=ar-.5*Ir;rest_pp=float(joint_predictive(ura,Fra,cr,Vr))
        sign_ar=signs@(urf@cr);paired_sign=sign_a-sign_ar
        results[arm]['frozen_rest_secondary']={'fixed_gain':rest_gain,'matched_filter':ar,'information':Ir,'conditional_Gaussian_signed_Z':ar/np.sqrt(Ir),'field_sign_exact_tail':float(np.mean(sign_ar>=ar-1e-10)),'joint_discovery_posterior_predictive_gain':rest_pp,
          'paired_observer_minus_rest_gain':gain-rest_gain,'prediction_correlation':cross/np.sqrt(I*Ir),'difference_prediction_information':difference_I,'paired_matched_filter':difference_a,'paired_Gaussian_signed_Z_about_zero_residual':difference_a/np.sqrt(difference_I),'paired_field_sign_exact_tail_about_zero_residual':float(np.mean(paired_sign>=difference_a-1e-10)),
          'qualification':'Paired null calibrations are around zero residual, not a composite test that all rest-frame explanations fail.'}
        arrays[arm+'_u']=u6;arrays[arm+'_F']=F6;arrays[arm+'_field_sign_matched_filter']=sign_a;arrays[arm+'_field_sign_joint_predictive']=sign_pp
    pd.DataFrame(fieldrows).to_csv(DEST/'field-scores.csv',index=False);np.savez_compressed(DEST/'sufficient-arrays.npz',**arrays)
    result={'scope':'Fixed empirical observer-template transfer to1020 highIa DES objects outside64 discovery objects; no physical calibration attribution or cosmology correction.',
      'objects':len(m),'fields':m.field.value_counts().to_dict(),'discovery_coefficient_mean':c.tolist(),'discovery_coefficient_covariance':V.tolist(),'frozen_rest_coefficient_mean':cr.tolist(),'frozen_rest_coefficient_covariance':Vr.tolist(),
      'arms':results,'field_sign_assignments':len(signs),'secondary_Gaussian_draws':N,'seed':seed,
      'limitations':['Selection/clipping and SALT training were conditioned on, not regenerated. Objects are disjoint but observations and trained model are not wholly independent experiments.',
        'Cross-object SALT/calibration covariance is not integrated by per-object SNANA C; common modes can invalidate Gaussian or field-sign significance.',
        'A coherent residual can be a realization of already released shared systematic uncertainty; this test does not reject the fully marginalized released systematic model.',
        'Native independent shape/colour/time derivatives approximate exact SNANA tangent; exact SNANA mean/amplitude used.',
        'High classifier score does not prove all objects are Ia or that classifier selection is independent of tested flux.',
        'Duplicate sensitivity is diagnostic and does not establish which records should be removed.']}
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();result['source_sha256']=sha(__file__)
    result['inputs_sha256']={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'input-manifest.json',OUT/'duplicate-amendment-manifest.json',OUT/'rest-amendment-manifest.json',OUT/'frozen-rest-discovery-coefficients.npz',OUT/'objectives'/'manifest.json',OUT/'frozen-discovery-coefficients.npz',BASE/'validation1020-protocol.md',ROOT/'scripts/salt_dust_audit/flux_response.py']}
    (DEST/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
