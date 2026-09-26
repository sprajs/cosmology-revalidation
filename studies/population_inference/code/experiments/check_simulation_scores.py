"""Independent sufficient-score, transport and paired block-bootstrap arithmetic."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;D=ROOT/'runs/research_2026_09_26/simulation_residual_control/score';obj=pd.read_csv(D/'object-scores.csv');tab=pd.read_csv(D/'transport-summary.csv');boot=pd.read_csv(D/'bootstrap.csv.gz');a=np.load(D/'sufficient-arrays.npz',allow_pickle=False);c=np.load(P/'validation1020/frozen-discovery-coefficients.npz')['basis_mean'];assert np.array_equal(a['CID'],obj.CID)
M=a['u']@c;I=np.einsum('j,njk,k->n',c,a['F'],c);G=M-I/2;err=max(np.max(abs(M-obj.matched_filter)),np.max(abs(I-obj.information)),np.max(abs(G-obj.fixed_gain)));assert err<1e-11
real=pd.read_csv(P/'validation1020/analysis/object-scores.csv').query("arm=='published_mask'");sys.path.insert(0,str(ROOT/'scripts/phase2/official'));from audit_fits import read_fit
rfit=read_fit(P/'validation1020/fit.FITRES.TEXT');real=real.merge(rfit[['CID','SNRMAX1']].assign(CID=lambda x:x.CID.astype(int)),on='CID',validate='one_to_one');real['cell']=real.field+'_'+pd.cut(real.zHEL,[.05,.2,.35,.5,.65,.8,1.2],include_lowest=True,labels=False).astype(int).astype(str);real['snrbin']=np.where(real.SNRMAX1<15,'lt15','ge15');obj['snrbin']=np.where(obj.SNRMAX1_archived<15,'lt15','ge15')
libs=sorted(set(obj.LIBID));rng=np.random.default_rng(26092691);draw=rng.choice(libs,size=(1000,len(libs)),replace=True);maxtransport=0.;maxboot=0.;rows=[]
def select(arm,law,stage):
 q=obj[(obj.arm==arm)&(obj.law==law)].copy()
 if stage=='archived_quality':q=q[q.archived_basic_quality]
 if stage=='new_quality':q=q[q.new_basic_quality]
 if stage=='p21_grid_nonnegative' and arm=='P21':q=q[q.support_class=='nonnegative_declared_grid']
 return q
for r in tab.itertuples():
 arms={k:select(k,r.law,r.stage) for k in ['P21','G10']};target=real.copy()
 for q in [*arms.values(),target]:q['key']=q.cell if r.transport=='field_z' else q.cell+'_'+q.snrbin
 counts={k:v.key.value_counts() for k,v in arms.items()};support=sorted(set(counts['P21'][counts['P21']>=2].index)&set(counts['G10'][counts['G10']>=2].index)&set(target.key));assert support==r.supported_cell_ids.split('|');t=target[target.key.isin(support)];q=arms[r.arm];q=q[q.key.isin(support)];assert len(q)==r.supported_objects and len(t)==r.supported_real_objects
 # Explicit per-cell averaging times target counts, not implementation's row map.
 cell=q.groupby('key')[['matched_filter','information','fixed_gain']].mean();target_n=t.key.value_counts().reindex(cell.index);sums=cell.mul(target_n,axis=0).sum();amp=sums.matched_filter/sums.information;gain=sums.fixed_gain/len(t);ra=t.matched_filter.sum()/t.information.sum();maxtransport=max(maxtransport,abs(amp-r.transported_amplitude),abs(gain-r.transported_mean_gain),abs(ra-r.real_supported_amplitude));assert maxtransport<1e-10
 b=boot[(boot.arm==r.arm)&(boot.law==r.law)&(boot.stage==r.stage)&(boot.transport==r.transport)].sort_values('replicate');assert len(b)==1000 and b.replicate.tolist()==list(range(1000));weights=q.key.map(t.key.value_counts())/q.key.map(q.key.value_counts())
 for rep in [0,1,2,500,999]:
  multi=pd.Series(draw[rep]).value_counts();w=weights.to_numpy()*q.LIBID.map(multi).fillna(0).to_numpy();aa=float(w@q.matched_filter/(w@q.information));gg=float(w@q.fixed_gain/w.sum());rr=b.iloc[rep];maxboot=max(maxboot,abs(aa-rr.amplitude),abs(gg-rr.mean_gain));assert maxboot<1e-10
 rows.append(dict(arm=r.arm,law=r.law,stage=r.stage,transport=r.transport,weight_ESS=float(weights.sum()**2/(weights@weights)),max_normalized_weight=float(weights.max()/weights.sum()),sim_amplitude=float(amp),real_amplitude=float(ra)))
pd.DataFrame(rows).to_csv(P/'simulation-score-independent-transport.csv',index=False);report={'status':'PASS: numeric sufficient-score arithmetic,32 independently constructed transports and160 independently reconstructed block-bootstrap replicates','objects':len(obj),'numeric_score_max_error':float(err),'transport_max_error':float(maxtransport),'bootstrap_max_error':float(maxboot),'design_limits':'LIBID resampling of one selected catalogue is not calibrated physical-null or independent-seed coverage; unknown real/CC selection and calibrated shared modes remain outside this test.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'object-scores.csv',D/'sufficient-arrays.npz',D/'transport-summary.csv',D/'bootstrap.csv.gz']}};(P/'simulation-score-independent.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k},indent=2))
