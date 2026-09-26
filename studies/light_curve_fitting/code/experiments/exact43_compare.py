"""Frozen implementation-discriminator protocol; preserves earlier residual runs."""
from pathlib import Path
import sys, json, hashlib, itertools, importlib.util
from types import SimpleNamespace
import numpy as np, pandas as pd
from scipy.linalg import solve_triangular, eigh
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).resolve().parent
OUT=BASE/'exact43'/'comparison'
sys.path.insert(0,str(BASE));import residual64 as r64
spec=importlib.util.spec_from_file_location('flux_response',ROOT/'scripts/salt_dust_audit/flux_response.py')
fr=importlib.util.module_from_spec(spec);spec.loader.exec_module(fr)

def main():
    if OUT.exists():raise RuntimeError('Preserve completed comparison')
    OUT.mkdir()
    model,bands,paths,zp=fr.build_model()
    magoff={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    meta=pd.read_csv(r64.DATA/'selected_objects.csv',dtype={'CID':str})
    meta=meta[meta.PROB_SNNV19>.999].copy().reset_index(drop=True)
    rows=pd.read_csv(ROOT/'phase2/hierarchy/data/conditioned-multistart-best/rows.csv',dtype={'CID':str})
    meta=meta.merge(rows[['CID','field']],on='CID',validate='one_to_one')
    old=np.load(r64.DATA/'matrices.npz');names=['official_mean_exactC','native_mean_exactC','native_mean_nativeC']
    arms={k:dict(u=[],F=[],chi2=[],templates=[]) for k in names};archives={};checks=[];epoch=[];coef=[];response=[]
    for count,record in enumerate(meta.to_dict('records')):
        cid=record['CID'];p=cid+'__';x=np.load(BASE/'exact43'/'objectives'/f'objective_{cid}.npz')
        n=len(x['MJD']);assert n==len(old[p+'mjd'])
        cost=abs(x['MJD'][:,None]-old[p+'mjd'][None,:])+(x['band'][:,None]!=old[p+'band'][None,:])*1e6+abs(x['data_flux'][:,None]-old[p+'flux_observed'][None,:])*1e-4
        rr,cc=linear_sum_assignment(cost);assert np.array_equal(rr,np.arange(n))
        assert np.array_equal(x['band'],old[p+'band'][cc]);assert max(abs(x['MJD']-old[p+'mjd'][cc]))<=.005
        assert max(abs(x['data_flux']-old[p+'flux_observed'][cc]))<.001
        assert np.ptp(x['MWEBV'])==0 and np.ptp(x['zHEL'])==0
        record.update(zip(['x0','x1','c','PKMJD'],x['parameters_x0_x1_c_t0']))
        record['zHEL']=float(x['zHEL'][0]);row=SimpleNamespace(**record)
        points=pd.DataFrame({'MJD':x['MJD'],'BAND':x['band'],'FLUXCAL':x['data_flux'],'FLUXCALERR':x['data_fluxerr']})
        order=points.sort_values(['MJD','BAND'],kind='stable').index.to_numpy()
        _,a,_,diag=fr.analyze(row,points,{'MWEBV':x['MWEBV'][0]},model,bands,magoff,False)
        official=x['model_flux'][order];C=x['frozen_flux_covariance'][np.ix_(order,order)]
        y=x['data_flux'][order];difference=a['flux_model']-official
        assert np.array_equal(a['band'],x['band'][order]);assert np.allclose(a['flux_observed'],y)
        # All epochs kept, including repeated timestamp/band multiplicities.
        phase=(a['mjd']-row.PKMJD)/(1+row.zHEL)
        generalized=eigh(C,a['measurement_plus_model_covariance'],eigvals_only=True)
        check={'CID':cid,'epochs':n,'zHEL':row.zHEL,'field':row.field,
          'epoch_match_max_time_difference_day':float(max(abs(x['MJD']-old[p+'mjd'][cc]))),
          'epoch_match_max_flux_difference':float(max(abs(x['data_flux']-old[p+'flux_observed'][cc]))),
          'native_derivative_halfstep_relative_error':diag['derivative_relative_error'],
          'exact_to_native_cov_eigen_min':float(generalized.min()),'exact_to_native_cov_eigen_max':float(generalized.max()),
          'raw_mean_diff_measurement_chi2':float(np.sum((difference/a['flux_error'])**2)),
          'official_full_chi2':float(x['chi2']), 'prior_chi2':float(x['prior_chi2'])}
        for arm in names:
            f=official if arm.startswith('official') else a['flux_model']
            cov=a['measurement_plus_model_covariance'] if arm.endswith('nativeC') else C
            J=a['jacobian_flux'].copy();J[:,0]=-fr.K*f
            G=a['nuisance_flux'].copy();G[:,5:]=np.column_stack([-fr.K*f*(a['band']==b) for b in 'griz'])
            obs=np.column_stack([G[:,5]-G[:,6],G[:,7]-G[:,6],G[:,8]-G[:,6]])
            rest=np.column_stack([G[:,2],G[:,1],J[:,2]*np.tanh(phase/20)])
            L=np.linalg.cholesky(cov);jw=solve_triangular(L,J,lower=True)
            U,s,Vt=np.linalg.svd(jw,full_matrices=True);assert np.sum(s>s[0]*1e-10)==4
            Q=U[:,4:];r=Q.T@solve_triangular(L,y-f,lower=True)
            T=Q.T@solve_triangular(L,np.column_stack([obs,rest]),lower=True)
            d=Q.T@solve_triangular(L,difference,lower=True)
            assert max(abs(Q.T@solve_triangular(L,-fr.K*f,lower=True)))<1e-6
            arms[arm]['u'].append(T.T@r);arms[arm]['F'].append(T.T@T)
            arms[arm]['chi2'].append(float(r@r));arms[arm]['templates'].append(T)
            archives[p+arm+'_T']=T;archives[p+arm+'_r']=r;archives[p+arm+'_mean_difference_projected']=d
            inverseJ=(Vt.T/s)@U[:,:4].T
            archives[p+arm+'_parameter_response_observer']=inverseJ@solve_triangular(L,obs,lower=True)
            check[arm+'_projected_chi2']=float(r@r)
            check[arm+'_projected_mean_difference_chi2']=float(d@d)
        archives.update({p+k:v for k,v in a.items()})
        archives[p+'official_flux_model']=official;archives[p+'exact_covariance']=C
        archives[p+'parameters_x0_x1_c_t0']=x['parameters_x0_x1_c_t0']
        epoch += [dict(CID=cid,band=b,phase=float(t),native_flux=float(nat),official_flux=float(off),observed_flux=float(yy),measurement_error=float(e)) for b,t,nat,off,yy,e in zip(a['band'],phase,a['flux_model'],official,y,a['flux_error'])]
        checks.append(check)
        if count%8==0:print('matched',count+1,'of',len(meta),flush=True)
    for a in arms.values():
        for k in ['u','F','chi2']:a[k]=np.array(a[k])
    np.savez_compressed(OUT/'matched-matrices.npz',**archives)
    pd.DataFrame(checks).to_csv(OUT/'object-checks.csv',index=False)
    pd.DataFrame(epoch).to_csv(OUT/'epoch-predictions.csv',index=False)
    fields=sorted(meta.field.unique());n=len(meta)
    groups={'object':[(row.CID,np.array([i])) for i,row in enumerate(meta.itertuples())],
            'field':[(f,np.flatnonzero(meta.field==f)) for f in fields]}
    seed=26092681;N=10000;rng=np.random.default_rng(seed)
    object_sign=rng.choice([-1.,1.],size=(N,n));field_sign=np.array(list(itertools.product([-1.,1.],repeat=len(fields))))
    field_sign=field_sign[:,np.array([fields.index(f) for f in meta.field])]
    report={};all_per=[];observed=[];signscores=[];fieldscores=[];save={}
    for name,a in arms.items():
        F,u=a['F'],a['u'];spec=r64.prepare_scores(F,groups)
        vals,per=r64.scores(u,spec,True);vals=vals[0]
        all_per += [dict(arm=name,**row) for row in per]
        eigen,V=np.linalg.eigh(F);assert eigen.min()>-1e-8
        root=V*np.sqrt(np.maximum(eigen,0))[:,None,:]
        gaussian=np.einsum('nik,bnk->bni',root,rng.normal(size=(N,n,6)))
        gs,_=r64.scores(gaussian,spec);os,_=r64.scores(object_sign[:,:,None]*u[None],spec)
        fs,_=r64.scores(field_sign[:,:,None]*u[None],spec)
        maximum=float(max(vals));sigma=.02;fall=F.sum(axis=0)[:3,:3];v=np.linalg.inv(np.eye(3)/sigma**2+fall);c=v@u.sum(axis=0)[:3]
        # Sum(griz)=0; contrasts create [c_g, -sum(c), c_i, c_z].
        gauge=np.array([[1,0,0],[-1,-1,-1],[0,1,0],[0,0,1.]])
        for band,value,sd in zip('griz',gauge@c,np.sqrt(np.diag(gauge@v@gauge.T))):coef.append(dict(arm=name,band=band,ridge_sigma=sigma,magnitude_extra_dimming=float(value),conditional_sd=float(sd)))
        # Implementation displacement projected in same observer basis, with same ridge.
        diff_u=np.array([archives[row.CID+'__'+name+'_T'].T@archives[row.CID+'__'+name+'_mean_difference_projected'] for row in meta.itertuples()])
        diffcoef=v@diff_u.sum(axis=0)[:3]
        for row in meta.itertuples():
            dpar=-archives[row.CID+'__'+name+'_parameter_response_observer']@c
            response.append(dict(arm=name,CID=row.CID,zHEL=row.zHEL,field=row.field,delta_mB=float(dpar[0]),delta_x1=float(dpar[1]),delta_c=float(dpar[2]),delta_t0=float(dpar[3]),fixed_reference_preBBC=float(np.array([1,.16087,-3.1178,0])@dpar)))
        stack=np.vstack(a['templates']);angles=np.linalg.svd(np.linalg.qr(stack[:,:3])[0].T@np.linalg.qr(stack[:,3:])[0],compute_uv=False)
        report[name]={'candidate_scores':{s[0]:float(vv) for s,vv in zip(spec,vals)},'max_score':maximum,
          'gaussian_arm_max_tail':r64.tail(gs.max(axis=1),maximum),
          'object_sign_arm_max_tail':r64.tail(os.max(axis=1),maximum),
          'field_sign_arm_max_fraction':float(np.mean(fs.max(axis=1)>=maximum-1e-10)),
          'projected_chi2':float(a['chi2'].sum()),'projected_dimension':sum(len(t) for t in a['templates']),
          'projected_native_minus_official_chi2':float(sum(check[name+'_projected_mean_difference_chi2'] for check in checks)),
          'ridge02_extra_dimming_griz':(gauge@c).tolist(),
          'ridge02_native_minus_official_griz':(gauge@diffcoef).tolist(),
          'subspace_canonical_correlations':angles.tolist()}
        observed.append(vals);signscores.append(os);fieldscores.append(fs)
        save[name+'_u']=u;save[name+'_F']=F;save[name+'_gaussian_scores']=gs
        print(name,maximum,report[name]['field_sign_arm_max_fraction'],flush=True)
    maxobs=float(np.concatenate(observed).max());objectmax=np.concatenate(signscores,axis=1).max(axis=1);fieldmax=np.concatenate(fieldscores,axis=1).max(axis=1)
    ep=pd.DataFrame(epoch);dif=(ep.native_flux-ep.official_flux).to_numpy();frac=dif/ep.official_flux.to_numpy()
    result={'scope':'Secondary implementation discriminator at identical native-objective coordinates and published epochs, not historical runtime reproduction or correction inference.',
      'objects':n,'epochs':len(epoch),'arms':report,'max_score_all_36':maxobs,
      'global_object_sign_tail':r64.tail(objectmax,maxobs),'global_field_sign_exact_fraction':float(np.mean(fieldmax>=maxobs-1e-10)),
      'global_gaussian_Bonferroni_three_arm_tail':min(1.,3*min(r['gaussian_arm_max_tail']['plus_one_tail'] for r in report.values())),
      'native_minus_official_fractional_median':float(np.median(frac)),
      'native_minus_official_absolute_fractional_p95':float(np.quantile(abs(frac),.95)),
      'native_minus_official_absolute_fractional_max':float(max(abs(frac))),
      'native_minus_official_measurement_sigma_p95':float(np.quantile(abs(dif/ep.measurement_error),.95)),
      'native_minus_official_measurement_sigma_max':float(max(abs(dif/ep.measurement_error))),
      'seed':seed,'simulations':N,'field_sign_assignments':len(fieldmax),
      'limitations':['Modern pinned SNANA objective, not unknown historical runtime.',
        'Native shape/colour/time tangent and dust derivatives retained; official amplitude derivative exact.',
        'Frozen covariance/mask and local linear nuisance projection do not regenerate selection, retraining, intrinsic scatter, or BBC.',
        'Shared-field sign exchangeability remains a conditional assumption; related arms are diagnostics, not independent replication.',
        'No measurement-only arm rerun here; adverse earlier measurement-only result is retained.']}
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    result['source_sha256']=sha(__file__);result['inputs_sha256']={str(p.relative_to(ROOT)):sha(p) for p in [BASE/'exact43-protocol.md',BASE/'exact43'/'objectives'/'manifest.json',r64.DATA/'matrices.npz',ROOT/'scripts/salt_dust_audit/flux_response.py',BASE/'residual64.py']}
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(all_per).to_csv(OUT/'block-scores.csv',index=False);pd.DataFrame(coef).to_csv(OUT/'coefficients.csv',index=False);pd.DataFrame(response).to_csv(OUT/'fitted-distance-response.csv',index=False)
    np.savez_compressed(OUT/'score-arrays.npz',**save)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
