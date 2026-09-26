"""Frozen 64-object linearized residual transfer and block-sign diagnostics."""
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,itertools
from scipy.linalg import solve_triangular
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
DATA=ROOT/'runs/salt_dust_audit/flux_response'

def tail(vals,t):
    n=len(vals);k=int(np.sum(vals>=t-1e-10))
    return {'count':k,'draws':n,'plus_one_tail':(k+1)/(n+1),
        'binomial_95_interval':[0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
        1. if k==n else float(beta.ppf(.975,k+1,n-k))]}

def build():
    x=np.load(DATA/'matrices.npz');m=pd.read_csv(DATA/'selected_objects.csv',dtype={'CID':str})
    rows=pd.read_csv(ROOT/'phase2/hierarchy/data/conditioned-multistart-best/rows.csv',dtype={'CID':str})
    m=m.merge(rows[['CID','field']],on='CID',validate='one_to_one',how='left')
    assert len(m)==64 and m.field.notna().all()
    arms={};info=[]
    for weighting in ['measurement_plus_model','measurement_only']:
        us=[];fs=[];chis=[];templates=[]
        for row in m.itertuples():
            p=row.CID+'__';C=x[p+weighting+'_covariance'];L=np.linalg.cholesky(C)
            J=x[p+'jacobian_flux'];G=x[p+'nuisance_flux']
            jw=solve_triangular(L,J,lower=True)
            U,s,V=np.linalg.svd(jw,full_matrices=True);rank=int(np.sum(s>s[0]*1e-10));assert rank==4
            Q=U[:,rank:]
            observed=x[p+'flux_observed'];model=x[p+'flux_model']
            r=Q.T@solve_triangular(L,observed-model,lower=True)
            phase=(x[p+'mjd']-row.PKMJD)/(1+row.zHEL)
            obs=np.column_stack([G[:,5]-G[:,6],G[:,7]-G[:,6],G[:,8]-G[:,6]])
            rest=np.column_stack([G[:,2],G[:,1],J[:,2]*np.tanh(phase/20)])
            T=Q.T@solve_triangular(L,np.column_stack([obs,rest]),lower=True)
            gray=Q.T@solve_triangular(L,-.4*np.log(10)*model,lower=True)
            assert np.max(abs(gray))<1e-6
            u=T.T@r;F=T.T@T
            us.append(u);fs.append(F);chis.append(float(r@r));templates.append(T)
            info.append({'weighting':weighting,'CID':row.CID,'field':row.field,'zHEL':row.zHEL,
              'pIa':row.PROB_SNNV19,'epoch_count':len(r)+rank,'residual_dimension':len(r),
              'projected_chi2':float(r@r),'projected_chi2_per_dim':float(r@r/len(r)),
              'max_gray_response':float(abs(gray).max()),'jacobian_condition':float(s[0]/s[-1])})
        arms[weighting]={'u':np.array(us),'F':np.array(fs),'chi2':np.array(chis),'templates':templates}
    return m,arms,info

def prepare_scores(F,groups):
    total=F.sum(axis=0);spec=[]
    for family,ii in [('observer',np.arange(3)),('rest_phase',np.arange(3,6))]:
        fall=total[np.ix_(ii,ii)]
        for sigma in [.01,.02,.05]:
            prior=np.eye(3)/sigma**2;vall=np.linalg.inv(prior+fall)
            kall=-.5*np.linalg.slogdet(np.eye(3)+sigma**2*fall)[1]
            for cv,grouplist in groups.items():
                folds=[]
                for label,ids in grouplist:
                    ftrain=fall-F[ids].sum(axis=0)[np.ix_(ii,ii)]
                    v=np.linalg.inv(prior+ftrain)
                    k=-.5*np.linalg.slogdet(np.eye(3)+sigma**2*ftrain)[1]
                    folds.append((str(label),ids,v,k))
                spec.append((f'{cv}_{family}_{sigma:g}',ii,vall,kall,folds))
    return spec

def scores(u,spec,return_folds=False):
    u=np.asarray(u);u=u[None] if u.ndim==2 else u
    total=u.sum(axis=1);out=[];per=[]
    for name,ii,vall,kall,folds in spec:
        tt=total[:,ii];la=.5*np.einsum('bi,ij,bj->b',tt,vall,tt)+kall
        gain=np.zeros(len(u))
        for label,ids,v,k in folds:
            train=tt-u[:,ids][:,:,ii].sum(axis=1)
            lt=.5*np.einsum('bi,ij,bj->b',train,v,train)+k
            gain+=la-lt
            if return_folds:per.append({'candidate':name,'block':label,'delta_log_score':float((la-lt)[0]),'max_abs_coefficient':float(abs(v@train[0]).max())})
        out.append(gain)
    return np.array(out).T,per

def main():
    if (OUT/'residual64.json').exists():raise RuntimeError('Preserve completed run')
    m,arms,info=build();n=len(m);fields=sorted(m.field.unique());nf=len(fields)
    groups={'object':[(row.CID,np.array([i])) for i,row in enumerate(m.itertuples())],
            'field':[(f,np.flatnonzero(m.field==f)) for f in fields]}
    N=10000;seed=26092678;rng=np.random.default_rng(seed)
    object_sign=rng.choice([-1.,1.],size=(N,n))
    field_sign=np.array(list(itertools.product([-1.,1.],repeat=nf)))
    field_sign=field_sign[:,np.array([fields.index(f) for f in m.field])]
    report={};all_per=[];observed=[];signscores=[];fieldscores=[];save={};powers={}
    for weighting,a in arms.items():
        F,u=a['F'],a['u'];spec=prepare_scores(F,groups)
        vals,per=scores(u,spec,True);vals=vals[0]
        all_per += [dict(weighting=weighting,**row) for row in per]
        eigen,V=np.linalg.eigh(F);assert eigen.min()>-1e-8
        root=V*np.sqrt(np.maximum(eigen,0))[:,None,:]
        gaussian=np.einsum('nik,bnk->bni',root,rng.normal(size=(N,n,6)))
        gs,_=scores(gaussian,spec);os,_=scores(object_sign[:,:,None]*u[None],spec)
        fs,_=scores(field_sign[:,:,None]*u[None],spec)
        candidate={s[0]:float(v) for s,v in zip(spec,vals)}
        maximum=float(max(vals));threshold=float(np.quantile(gs.max(axis=1),.975,method='higher'))
        stack=np.vstack(a['templates']);so=np.linalg.svd(stack[:,:3],compute_uv=False);sr=np.linalg.svd(stack[:,3:],compute_uv=False)
        qo=np.linalg.qr(stack[:,:3])[0];qr=np.linalg.qr(stack[:,3:])[0]
        angles=np.linalg.svd(qo.T@qr,compute_uv=False)
        report[weighting]={'candidate_scores':candidate,'max_score':maximum,
          'gaussian_arm_max_tail':tail(gs.max(axis=1),maximum),
          'object_sign_arm_max_tail':tail(os.max(axis=1),maximum),
          'field_sign_arm_max_fraction':float(np.mean(fs.max(axis=1)>=maximum-1e-10)),
          'projected_chi2':float(a['chi2'].sum()),'projected_dimension':int(sum(len(t) for t in a['templates'])),
          'observer_singular_values':so.tolist(),'rest_phase_singular_values':sr.tolist(),
          'subspace_canonical_correlations':angles.tolist()}
        observed.append(vals);signscores.append(os);fieldscores.append(fs)
        save[weighting+'_u']=u;save[weighting+'_F']=F;save[weighting+'_gaussian_scores']=gs
        # Prespecified explanatory signal levels; max threshold includes both
        # CV forms, families and scales and a Bonferroni allowance for two arms.
        for idx,label in [(1,'observer_i_minus_r'),(3,'host_RV3p1'),(5,'phase_colour')]:
            amp=.05
            shift=amp*F[:,:,idx]
            gs_alt,_=scores(gaussian+shift[None],spec)
            powers[weighting+'_'+label]={'injected_amplitude':amp,
              'known_direction_snr':float(amp*np.sqrt(F[:,idx,idx].sum())),
              'conservative_global_5percent_power':float(np.mean(gs_alt.max(axis=1)>=threshold))}
        print(weighting,maximum,report[weighting]['gaussian_arm_max_tail'],flush=True)
    maxobs=float(np.concatenate(observed).max());objectmax=np.concatenate(signscores,axis=1).max(axis=1)
    fieldmax=np.concatenate(fieldscores,axis=1).max(axis=1)
    global_gauss=min(1.,2*min(r['gaussian_arm_max_tail']['plus_one_tail'] for r in report.values()))
    result={'scope':'Retrospective 64-object linearized nuisance-projected flux-residual transport; no physical correction/cosmology inference.',
      'objects':n,'epochs':int(sum(len(t)+4 for t in arms['measurement_only']['templates'])),
      'fields':m.field.value_counts().to_dict(),'low_classifier_probability_count':int((m.PROB_SNNV19<=.999).sum()),
      'arms':report,'global_gaussian_Bonferroni_two_arm_tail':global_gauss,
      'max_score_all_24':maxobs,'global_object_sign_tail':tail(objectmax,maxobs),
      'global_field_sign_exact_fraction':float(np.mean(fieldmax>=maxobs-1e-10)),
      'field_sign_assignments':len(fieldmax),'power':powers,'seed':seed,'simulations':N,
      'limitations':['Linearized published-coordinate nuisance refit, not nonlinear likelihood refit.',
       'Model covariance is archived sncosmo variant, not exact official SNANA.',
       'Sign symmetry and independent blocks may fail after selection/clipping or shared model error.',
       'All 64 include imperfect photometric Ia probabilities; no classifier-selected transport claimed.',
       'Ridge scales and small empirical templates do not define complete physical alternatives.']}
    inputs=[DATA/'matrices.npz',DATA/'selected_objects.csv',ROOT/'phase2/hierarchy/data/conditioned-multistart-best/rows.csv',OUT/'residual64-protocol.md']
    sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
    result['inputs_sha256']={str(p.relative_to(ROOT)):sha(p) for p in inputs};result['source_sha256']=sha(__file__)
    (OUT/'residual64.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(all_per).to_csv(OUT/'residual64-block-scores.csv',index=False)
    pd.DataFrame(info).to_csv(OUT/'residual64-objects.csv',index=False)
    np.savez_compressed(OUT/'residual64-arrays.npz',**save)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
