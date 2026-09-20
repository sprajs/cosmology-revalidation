"""Read-only audit of released numerical inputs, script syntax and source hashes."""
import ast
import hashlib
import pathlib
import numpy as np
import pandas as pd
from scipy.linalg import cholesky,eigh
from audit import ROOT,SRC,OUT,write

p=pd.read_csv(SRC/'Analysis C1/Pantheon+SH0ES.dat',sep=r'\s+')
z=pd.read_csv(SRC/'Analysis C1/Z_mbcorr.csv')
c=np.load(SRC/'Analysis C1/statsys_mbcorr.npy',allow_pickle=False)
public_path=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov'
public=np.loadtxt(public_path,skiprows=1).reshape(len(p),len(p))
results={'C1':dict(n=len(p),n_unique_CID=int(p.CID.nunique()),column_differences={k:float(abs(z[k]-p[k]).max()) for k in z},cov_shape=list(c.shape),max_cov_asymmetry=float(abs(c-c.T).max()),cov_min_eigenvalue=float(eigh(c,eigvals_only=True,subset_by_index=[0,0])[0]),max_difference_from_public_total_cov=float(abs(c-public).max()),n_main_cut=int(sum((p.zHEL>.00937)&(p.zHEL<.8)))), 'scripts':{}}
for f in SRC.glob('Analysis*/*.py'):
    try:ast.parse(f.read_text());ans='syntax valid'
    except Exception as e:ans=str(e)
    results['scripts'][str(f.relative_to(ROOT))]=ans
write('input-audit.json',results)
f=OUT/'sources/cov_final.npy';c=np.load(f,allow_pickle=False)
ix=np.load(SRC/'index_sorted_lane.npy',allow_pickle=False);d=p.iloc[ix];z=pd.read_csv(SRC/'Analysis C2/Zpan.csv')
cc=dict(matrix_shape=list(c.shape),index_length=len(ix),index_unique=len(set(ix)),min_sorted_z_difference=float(np.diff(d.zHEL).min()),covariance_max_asymmetry=float(abs(c-c.T).max()),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),md5=hashlib.md5(f.read_bytes()).hexdigest(),input_column_max_differences={k:float(abs(z[k]-p[k]).max()) for k in z if k in p},dropped_rows=p.loc[~p.index.isin(ix),['CID','IDSURVEY','zHEL']].to_dict('records'),n_below_point8=int(sum(d.zHEL<.8)),n_main_cut=int(sum((d.zHEL>.00937)&(d.zHEL<.8))))
for k,col in [(0,'mBERR'),(1,'x1ERR'),(2,'cERR')]:
    a=np.diag(c)[k::3];b=d[col].to_numpy()**2
    cc[col]=dict(correlation=float(np.corrcoef(a,b)[0,1]),extra_variance_quantiles=np.quantile(a-b,[0,.5,.95,1]).tolist())
L=cholesky(c,lower=True,check_finite=False);cc.update(cholesky='success',min_cholesky_diagonal=float(np.diag(L).min()))
assert cc['md5']=='09996e2b37009aa7d0c7de13b79c90c7'
write('c2-input-audit.json',cc)
files=list((OUT/'sources').iterdir())
write('source-hashes.json',{str(f.relative_to(ROOT)):dict(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files if f.is_file()})
print('C1 and C2 input audit complete; no pickle or author code executed.')
