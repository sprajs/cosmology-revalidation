"""Independent holdout sufficient-score/support/transport arithmetic."""
from pathlib import Path
import hashlib,json,sys
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;R=P.parents[2];H=R/'runs/research_2026_09_26/simulation_holdout_control';B=R/'runs/research_2026_09_26/simulation_residual_control';CL=R/'runs/research_2026_09_26/common_classifier_residual';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();hashes=0
for d in ['object_scores','transport']:
 m=json.loads((H/d/'manifest.json').read_text())
 for k in ['inputs_sha256','outputs_sha256']:
  for p,h in m[k].items():assert sha(R/p)==h,p;hashes+=1
hold=pd.read_csv(H/'object_scores/holdout-object-scores.csv');s=np.load(H/'object_scores/holdout-sufficient-arrays.npz',allow_pickle=False);c=np.load(P/'validation1020/frozen-discovery-coefficients.npz')['basis_mean'];assert np.array_equal(s['CID'],hold.CID) and np.array_equal(s['arm'],hold.arm)
M=s['u']@c;I=np.einsum('j,njk,k->n',c,s['F'],c);serr=max(float(np.max(abs(M-hold.matched_filter))),float(np.max(abs(I-hold.information))),float(np.max(abs(M-I/2-hold.fixed_gain))));assert serr<1e-10
pilot=pd.read_csv(B/'score/object-scores.csv').query("law=='approx_minus99'").merge(pd.read_csv(CL/'inference/simulation-probabilities.csv')[['arm','CID','gt999']],on=['arm','CID'],validate='one_to_one');pilot['is_holdout']=False;hold['is_holdout']=True;combined=pd.concat([hold,pilot],ignore_index=True);assert not combined.duplicated(['arm','CID']).any();combined['snrbin']=np.where(combined.SNRMAX1_archived<15,'lt15','ge15')
for arm in ['P21','G10']:
 pid=set(pd.read_csv(P/f'simulation_design/{arm}-cohort.csv').CID);hid=set(pd.read_csv(P/f'simulation_holdout_design/{arm}-cohort.csv').CID);assert not pid&hid;assert set(hold.loc[hold.arm==arm,'CID'])<=hid
rp=pd.read_csv(CL/'inference/real-probabilities.csv');real=pd.read_csv(P/'validation1020/analysis/object-scores.csv');real=real[(real.arm=='published_mask')&real.CID.isin(rp.loc[rp.gt999,'CID'])].copy();assert len(real)==1006;real['cell']=real.field+'_'+pd.cut(real.zHEL,[.05,.2,.35,.5,.65,.8,1.2],include_lowest=True,labels=False).astype(int).astype(str);sys.path.insert(0,str(R/'scripts/phase2/official'));from audit_fits import read_fit
rf=read_fit(P/'validation1020/fit.FITRES.TEXT');real=real.merge(rf[['CID','SNRMAX1']].assign(CID=lambda x:x.CID.astype(int)),on='CID');real['snrbin']=np.where(real.SNRMAX1<15,'lt15','ge15')
tab=pd.read_csv(H/'transport/transport-summary.csv');ledger=pd.read_csv(H/'transport/supported-object-ledger.csv');boot=pd.read_csv(H/'transport/bootstrap.csv.gz');libs=sorted(set(combined.LIBID));draw=np.random.default_rng(26092691).choice(libs,size=(1000,len(libs)),replace=True);maxerr=0.;berr=0.
def select(src,arm,stage):
 q=src[src.arm==arm].copy()
 if stage!='all_success':q=q[q.archived_basic_quality]
 if stage.startswith('quality_classifier'):q=q[q.gt999]
 if stage.endswith('_p21_nonnegative') and arm=='P21':q=q[q.support_class=='nonnegative_declared_grid']
 return q
for r in tab.itertuples():
 src=combined[combined.is_holdout].copy() if r.scope=='holdout' else combined.copy();ss=r.stage if r.support_source=='own_stage' else 'quality_classifier';arms={a:select(src,a,ss) for a in ['P21','G10']};target=real.copy();q=select(src,r.arm,r.stage)
 for x in [*arms.values(),target,q]:x['key']=x.cell if r.transport=='field_z' else x.cell+'_'+x.snrbin
 cnt={a:x.key.value_counts() for a,x in arms.items()};cells=sorted(set(cnt['P21'][cnt['P21']>=2].index)&set(cnt['G10'][cnt['G10']>=2].index)&set(target.key));assert cells==r.cell_ids.split('|');q=q[q.key.isin(cells)];target=target[target.key.isin(cells)];assert len(q)==r.supported_sim and len(target)==r.supported_real
 means=q.groupby('key')[['matched_filter','information','fixed_gain','projected_chi2','projected_dimension']].mean();n=target.key.value_counts().reindex(means.index);agg=means.mul(n,axis=0).sum();vals=[agg.matched_filter,agg.information,agg.fixed_gain,agg.projected_chi2,agg.projected_dimension,agg.projected_chi2/agg.projected_dimension,target.matched_filter.sum(),target.information.sum(),target.projected_chi2.sum()/target.projected_dimension.sum()];refs=[r.weighted_a,r.weighted_I,r.weighted_G,r.weighted_chi2,r.weighted_dimension,r.pooled_Q,r.real_a,r.real_I,r.real_pooled_Q];maxerr=max(maxerr,float(max(abs(np.array(vals)-refs))));assert maxerr<1e-8
 w=q.key.map(target.key.value_counts())/q.key.map(q.key.value_counts());mask=lambda d:np.logical_and.reduce([d[k]==getattr(r,k) for k in ['scope','arm','transport','stage','support_source']]);l=ledger[mask(ledger)].set_index('CID');assert set(l.index)==set(q.CID);assert max(abs(l.loc[q.CID].base_weight.to_numpy()-w.to_numpy()))<1e-10;b=boot[mask(boot)].sort_values('replicate')
 for rep in [0,500,999]:
  ww=w.to_numpy()*q.LIBID.map(pd.Series(draw[rep]).value_counts()).fillna(0).to_numpy();vals=[ww@q.matched_filter/(ww@q.information),ww@q.fixed_gain/ww.sum(),ww@q.projected_chi2/(ww@q.projected_dimension)];refs=[b.iloc[rep].amplitude,b.iloc[rep].mean_gain,b.iloc[rep].pooled_Q];berr=max(berr,float(max(abs(np.array(vals)-refs))));assert berr<1e-10
out={'status':'PASS independent arithmetic, cohort-disjointness, support, weights and selected paired-block replicas','hashes_verified':hashes,'holdout_objects':len(hold),'combined_objects':len(combined),'sufficient_score_error':serr,'transport_rows':len(tab),'transport_weighted_sum_error':maxerr,'bootstrap_replicas_checked':len(tab)*3,'bootstrap_error':berr,'scope':'Independent same-catalogue holdout audit; not independent random seed or calibrated physical-null coverage. Mean/tangent and noise-recipe limitations remain.'};(P/'holdout-score-independent.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
