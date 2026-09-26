"""Frozen pilot same-support covariance audit; all original scores/C remain unchanged."""
from pathlib import Path
import sys,json,hashlib,re
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parent/'common_noise_audit'; OUT.mkdir(exist_ok=True)
P=OUT.parent; SIM=ROOT/'runs/research_2026_09_26/simulation_residual_control'; CL=ROOT/'runs/research_2026_09_26/common_classifier_residual'; VAL=P/'validation1020'
sys.path.insert(0,str(ROOT/'scripts/salt_dust_audit'));import flux_response as fr
sys.path.insert(0,str(ROOT/'scripts/phase2/official'));from audit_fits import read_fit
coeff=np.load(VAL/'frozen-discovery-coefficients.npz',allow_pickle=False)['basis_mean']
model,bands,mpaths,zp=fr.build_model();offset={str(q['Filter Name'])[-1]:float(q['Primary Mag']) for q in zp}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest(); hashes={};rows=[]
# Runtime NML echoes prove the covariance separation is enabled, not just defaults.
logs=[VAL/'fit.log']+[SIM/'full256'/a/l/'fit.log' for a in ['P21','G10'] for l in ['approx_minus99','exact_99']]
for p in logs:
 t=p.read_text();assert re.search(r'FUDGE_MAG_COVERR="\s*"',t);assert re.search(r'OPT_COVAR_MWXTERR=1\s*,',t)
 hashes[str(p.relative_to(ROOT))]=sha(p)
real=pd.read_csv(VAL/'analysis/object-scores.csv').query("arm=='published_mask'").copy();real['law']='exact_99';real['population']='real'
sim=pd.read_csv(SIM/'score/object-scores.csv');sim['population']=sim.arm
for row in pd.concat([real,sim],ignore_index=True).itertuples():
 cid=int(row.CID);pop=row.population;law=row.law
 p=(VAL if pop=='real' else SIM/'full256'/pop/law)/'objectives'/f'objective_{cid}.npz'
 a=np.load(p,allow_pickle=False);hashes[str(p.relative_to(ROOT))]=sha(p)
 order=pd.DataFrame({'m':a['MJD'],'b':a['band']}).sort_values(['m','b'],kind='stable').index.to_numpy();n=len(order)
 f=a['model_flux'][order];y=a['data_flux'][order];b=a['band'][order];C=a['frozen_flux_covariance'][np.ix_(order,order)];cp=a['covariance_components'][order];t=a['MJD'][order];phase=a['rest_phase'][order]
 if pop=='real':
  cache=VAL/'analysis/objects'/f'{cid}.npz';q=np.load(cache,allow_pickle=False);hashes[str(cache.relative_to(ROOT))]=sha(cache)
  for key,x in [('MJD',t),('band',b),('official_flux_model',f),('observed_flux',y),('exact_covariance',C)]:assert np.array_equal(q[key],x),(cid,key)
  J=q['jacobian_flux'].copy();J[:,0]=-fr.K*f
 else:
  x0,x1,color,t0=a['parameters_x0_x1_c_t0'];z=float(a['zHEL'][0]);ebv=float(a['MWEBV'][0]);bs=np.array([bands[v] for v in b],dtype=object);scale=np.array([10**(-.4*(.27+offset[v])) for v in b])
  def flux(dx,dc,dt):
   model.set(z=z,x0=x0,x1=x1+dx,c=color+dc,t0=t0+dt,mwebv=ebv,mwrv=3.1,hostebv=0.,hostrv=3.1)
   return scale*model.bandflux(bs,t,zp=27.5,zpsys='ab')
  J=np.column_stack([-fr.K*f,(flux(.001,0,0)-flux(-.001,0,0))/.002,(flux(0,.0001,0)-flux(0,-.0001,0))/.0002,(flux(0,0,.01)-flux(0,0,-.01))/.02])
 L=np.linalg.cholesky(C);wj=solve_triangular(L,J,lower=True);U,s,_=np.linalg.svd(wj,full_matrices=False);assert (s>s[0]*1e-10).sum()==4
 def projected(x):
  w=solve_triangular(L,x,lower=True);return w-U@(U.T@w)
 r=projected(y-f);gray=np.column_stack([-fr.K*f*(b==v) for v in 'griz']);T=projected(np.column_stack([gray[:,0]-gray[:,1],gray[:,2]-gray[:,1],gray[:,3]-gray[:,1]]));v=T@coeff
 chi=float(r@r);M=float(r@v);I=float(v@v)
 assert max(abs(chi-row.projected_chi2),abs(M-row.matched_filter),abs(I-row.information))<1e-7,(pop,law,cid,chi-row.projected_chi2)
 # R4EP sums/products are single precision on the RHS in Fortran; retain
 # this exact diagonal convention, not an unlabelled float64 replacement.
 d=(cp[:,3].astype('f4')**2+cp[:,4].astype('f4')**2).astype('f8');mw=cp[:,0];rem=C-np.diag(d)-np.outer(mw,mw)
 pdroot=projected(np.diag(np.sqrt(d)));pmw=projected(mw);traceD=float(np.sum(pdroot**2));traceMW=float(pmw@pmw);traceR=float(n-4-traceD-traceMW)
 wr=solve_triangular(L,rem,lower=True);wr=solve_triangular(L,wr.T,lower=True).T;mineig=float(np.linalg.eigvalsh((wr+wr.T)/2)[0])
 rec=dict(population=pop,law=law,CID=cid,epochs=n,nu=n-4,chi2=chi,chi2_per_nu=chi/(n-4),M=M,I=I,G=M-I/2,trace_data=traceD,trace_local_MW=traceMW,trace_SALT_remainder=traceR,whitened_remainder_min_eigenvalue=mineig,max_fudge_flux_error=float(abs(cp[:,4]).max()),phase_min=float(phase.min()),phase_max=float(phase.max()))
 for band in 'griz':rec['epochs_'+band]=int((b==band).sum())
 for name,mask in [('early',phase<0),('late',phase>=0),('before_minus10',phase<-10),('after30',phase>30)]:rec['epochs_'+name]=int(mask.sum())
 rows.append(rec)
 if len(rows)%200==0:print('objects',len(rows),flush=True)
obj=pd.DataFrame(rows);obj.to_csv(OUT/'object-noise.csv',index=False,float_format='%.17g')
# Independently construct support and transport; compare all supplied scores.
prob=pd.read_csv(CL/'inference/simulation-probabilities.csv');rp=pd.read_csv(CL/'inference/real-probabilities.csv');good=set(rp.loc[rp.gt999,'CID']);assert len(good)==1006
sim=sim.merge(prob[['arm','CID','gt999']],on=['arm','CID'],validate='many_to_one');sim['snrbin']=np.where(sim.SNRMAX1_archived<15,'lt15','ge15')
real=real[real.CID.isin(good)].copy();rf=read_fit(VAL/'fit.FITRES.TEXT');real=real.merge(rf[['CID','SNRMAX1']].assign(CID=lambda x:x.CID.astype(int)),on='CID',validate='one_to_one');real['snrbin']=np.where(real.SNRMAX1<15,'lt15','ge15');real['cell']=real.field+'_'+pd.cut(real.zHEL,[.05,.2,.35,.5,.65,.8,1.2],include_lowest=True,labels=False).astype(int).astype(str)
tab=pd.read_csv(CL/'score/transport-summary.csv');support=pd.read_csv(CL/'score/supported-object-ledger.csv');boots=pd.read_csv(CL/'score/bootstrap.csv.gz');libs=sorted(set(sim.LIBID));draw=np.random.default_rng(26092691).choice(libs,size=(1000,len(libs)),replace=True)
results=[];max_score_error=0.;max_boot_error=0.
def select(arm,law,stage):
 q=sim[(sim.arm==arm)&(sim.law==law)&sim.archived_basic_quality].copy()
 if stage.startswith('post_classifier'):q=q[q.gt999]
 if stage=='post_classifier_p21_nonnegative' and arm=='P21':q=q[q.support_class=='nonnegative_declared_grid']
 return q
for rr in tab.itertuples():
 stage=rr.stage;ss=stage if rr.support_source=='own_stage' else 'post_classifier_quality';arms={a:select(a,rr.law,ss) for a in ['P21','G10']};target=real.copy();q=select(rr.arm,rr.law,stage)
 for v in [*arms.values(),q,target]:v['key']=v.cell if rr.transport=='field_z' else v.cell+'_'+v.snrbin
 counts={a:v.key.value_counts() for a,v in arms.items()};cells=sorted(set(counts['P21'][counts['P21']>=2].index)&set(counts['G10'][counts['G10']>=2].index)&set(target.key));assert cells==rr.cell_ids.split('|')
 q=q[q.key.isin(cells)].copy();target=target[target.key.isin(cells)].copy();assert len(q)==rr.supported_sim and len(target)==rr.supported_real
 w=(q.key.map(target.key.value_counts())/q.key.map(q.key.value_counts())).to_numpy();assert abs(sum(w)-len(target))<1e-9
 sl=support[(support.arm==rr.arm)&(support.law==rr.law)&(support.transport==rr.transport)&(support.stage==stage)&(support.support_source==rr.support_source)].set_index('CID');assert set(sl.index)==set(q.CID);assert np.max(abs(sl.loc[q.CID].base_weight-w))<1e-10
 s=np.array([w@q.matched_filter,w@q.information,w@q.fixed_gain]);max_score_error=max(max_score_error,np.max(abs(s-[rr.weighted_a,rr.weighted_I,rr.weighted_G])),abs(target.matched_filter.sum()-rr.real_a),abs(target.information.sum()-rr.real_I));assert max_score_error<1e-9
 b=boots[(boots.arm==rr.arm)&(boots.law==rr.law)&(boots.transport==rr.transport)&(boots.stage==stage)&(boots.support_source==rr.support_source)].sort_values('replicate')
 for rep in [0,1,2,500,999]:
  ww=w*q.LIBID.map(pd.Series(draw[rep]).value_counts()).fillna(0).to_numpy();aa=ww@q.matched_filter/(ww@q.information);gg=ww@q.fixed_gain/ww.sum();max_boot_error=max(max_boot_error,abs(aa-b.iloc[rep].amplitude),abs(gg-b.iloc[rep].mean_gain));assert max_boot_error<1e-9
 for pop,base,weight in [(rr.arm,q,w),('real',target,np.ones(len(target)))]:
  x=base[['CID']].merge(obj[(obj.population==pop)&(obj.law==(rr.law if pop!='real' else 'exact_99'))],on='CID',validate='one_to_one');assert np.array_equal(x.CID,base.CID)
  rec=dict(population=pop,comparison_arm=rr.arm,law=rr.law,transport=rr.transport,stage=stage,support_source=rr.support_source,objects=len(x),supported_real=len(target),cells=len(cells),weight_ESS=float(sum(weight)**2/(weight@weight)),weight_max_fraction=float(weight.max()/sum(weight)))
  for col in ['epochs','nu','chi2','M','I','G','trace_data','trace_local_MW','trace_SALT_remainder']+[k for k in x.columns if k.startswith('epochs_')]:rec['sum_'+col]=float(weight@x[col])
  rec['pooled_chi2_per_nu']=rec['sum_chi2']/rec['sum_nu'];rec['amplitude']=rec['sum_M']/rec['sum_I'];rec['phase_min']=float(x.phase_min.min());rec['phase_max']=float(x.phase_max.max())
  # Object quantiles use frozen transport weights (midpoint weighted CDF).
  order=np.argsort(x.chi2_per_nu.to_numpy());cdf=(np.cumsum(weight[order])-.5*weight[order])/sum(weight)
  for pcent in [.25,.5,.75]:rec['object_Q_quantile_'+str(pcent)]=float(np.interp(pcent,cdf,x.chi2_per_nu.to_numpy()[order]))
  results.append(rec)
pd.DataFrame(results).to_csv(OUT/'supported-noise-summary.csv',index=False,float_format='%.17g')
for p in [P/'common-noise-audit-protocol.md',Path(__file__),CL/'score/transport-summary.csv',CL/'score/supported-object-ledger.csv',CL/'score/bootstrap.csv.gz',SIM/'score/object-scores.csv',VAL/'analysis/object-scores.csv',CL/'inference/simulation-probabilities.csv',CL/'inference/real-probabilities.csv',ROOT/'phase2/official/build/SNANA-audit-v3/src/snlc_fit.F90',*mpaths]:hashes[str(p.relative_to(ROOT))]=sha(p)
report={'objects':len(obj),'transport_comparisons':len(tab),'transport_score_max_error':float(max_score_error),'bootstrap_max_error':float(max_boot_error),'min_whitened_SALT_remainder_eigenvalue':float(obj.whitened_remainder_min_eigenvalue.min()),'remainder_negative_below_minus_1e6':int((obj.whitened_remainder_min_eigenvalue< -1e-6).sum()),'max_fudge_flux_error':float(obj.max_fudge_flux_error.max()),'covariance_scope':'Nominal frozen iteration conditional allocation; no covariance changed. Model remainder isolated after source/runtime gates for diagonal data/fudge and full local MW outer product. Any whitened negative eigenvalues retained.','input_sha256':hashes}
(OUT/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'},indent=2))
