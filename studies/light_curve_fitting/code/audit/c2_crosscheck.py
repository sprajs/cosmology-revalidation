"""C2 likelihood checked in standardized-observable coordinates, not source basis."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np
import pandas as pd
from scipy.linalg import cholesky,solve_triangular
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'runs/audit'
SRC=ROOT/'sources/repos/Shin107__Anisotropy-in-Pantheon-Plus'
direction=ROOT/'runs/directional'
fitpath=direction/'c2-fits.json';result_text=fitpath.read_text();fits=json.loads(result_text)
raw=pd.read_csv(SRC/'Analysis C2/Pantheon+SH0ES.dat',sep=r'\s+')
index=np.load(SRC/'index_sorted_lane.npy',allow_pickle=False)
tab=raw.iloc[index].reset_index(drop=True)
mask=(tab.zHEL>.00937)&(tab.zHEL<.8)
blocks=np.nonzero(mask.to_numpy())[0];ind=np.concatenate([np.arange(3*i,3*i+3) for i in blocks])
source_cov=np.load(direction/'sources/cov_final.npy',allow_pickle=False)
cov=source_cov[np.ix_(ind,ind)].copy();del source_cov
d=tab[mask];n=len(d);z=d.zHEL.to_numpy()
xyz=lambda ra,dec:np.column_stack([np.cos(np.deg2rad(dec))*np.cos(np.deg2rad(ra)),np.cos(np.deg2rad(dec))*np.sin(np.deg2rad(ra)),np.sin(np.deg2rad(dec))])
vectors=xyz(d.RA.to_numpy(),d.DEC.to_numpy());axis=xyz(np.array([168.]),np.array([-7.]))[0];cosine=vectors@axis
age=pd.read_csv(SRC/'median_deltaage.csv',header=None).sort_values(0)
spline=CubicSpline(age[0].to_numpy(),age[1].to_numpy())
results={'result_input_sha256':hashlib.sha256(result_text.encode()).hexdigest(),'n':n,
         'row_indices_original_pantheon':index[blocks].tolist(),'checks':[]}
for f in fits:
    assert f['n']==n and f['frame']=='zHEL' and f['reverse_bias']
    for label in ['dipole','isotropic']:
        p=f[label];a=p['alpha'];b=p['beta']
        # Unit determinant T maps m,x,c to m+a*x-b*c,x,c. This diagonalizes the latent covariance.
        transformed=cov.copy()
        transformed[0::3,:]=cov[0::3,:]+a*cov[1::3,:]-b*cov[2::3,:]
        transformed[:,0::3]=transformed[:,0::3]+a*transformed[:,1::3]-b*transformed[:,2::3]
        diag=np.diag_indices(3*n)
        transformed[diag]+=np.tile([p['sigmaM']**2,p['sigmaX']**2,p['sigmaC']**2],n)
        symmetry=float(np.max(abs(transformed-transformed.T)))
        L=cholesky(transformed,lower=True,check_finite=False);del transformed
        y=d[['mB','x1','c']].to_numpy().copy();y[:,0]+=d.biasCor_m_b.to_numpy()+a*y[:,1]-b*y[:,2]
        if f['age_correction']:y[:,0]-=.03*spline(z)
        q=p['q']+(p['qd']*cosine*np.exp(-z/p['S']) if label=='dipole' else 0)
        dl=(299792.458/70)*(z+(1-q)*z*z/2-(1-q-3*q*q+p['j'])*z**3/6)
        y[:,0]-=25+5*np.log10(dl)
        r=solve_triangular(L,y.ravel(),lower=True,check_finite=False)
        design=np.tile(np.eye(3),(n,1));wd=solve_triangular(L,design,lower=True,check_finite=False)
        means=np.linalg.lstsq(wd,r,rcond=None)[0]
        w=r-wd@means;chi=float(w@w);nll=chi+2*np.log(np.diag(L)).sum()+3*n*np.log(2*np.pi)
        published=np.array([p['M'],p['X'],p['C']]);wp=r-wd@published
        results['checks'].append({'age_correction':f['age_correction'],'model':label,
                 'nll2':float(nll),'nll2_reported':p['nll2'],'nll2_difference':float(nll-p['nll2']),
                 'chi2':chi,'chi2_difference':float(chi-p['chi2']),
                 'profiled_M_X_C':means.tolist(),'reported_M_X_C':published.tolist(),
                 'max_mean_difference':float(abs(means-published).max()),'reported_means_minus_profile_chi2':float(wp@wp-chi),
                 'transformed_covariance_max_asymmetry':symmetry,'q0':p['q'],'qd':p['qd'],'S':p['S']})
        assert abs(nll-p['nll2'])<1e-7
        assert abs(means-published).max()<1e-8
        print(results['checks'][-1],flush=True)
# Coordinate-independent hemisphere cross-check via the conventional J2000 rotation.
# Equatorial -> Galactic orthogonal rotation (IAU/Hipparcos convention); transpose is the inverse.
R=np.array([[-.0548755604162154,-.8734370902348850,-.4838350155487132],
            [.4941094278755837,-.4448296299600112,.7469822444972189],
            [-.8676661490190047,-.1980763734312015,.4559837761750669]])
galaxis=xyz(np.array([264.]),np.array([48.]))[0];eqaxis=R.T@galaxis
ra=float(np.rad2deg(np.arctan2(eqaxis[1],eqaxis[0]))%360);dec=float(np.rad2deg(np.arcsin(eqaxis[2])))
hem=raw[(raw.zHEL>.00937)&(raw.zHEL<=.8)];v=xyz(hem.RA.to_numpy(),hem.DEC.to_numpy())
co=v@eqaxis;galco=(v@R.T)@galaxis
def counts(axis):
    dot=v@axis;return [int(np.count_nonzero(dot>=0)),int(np.count_nonzero(dot<0))]
results['coordinates']={'rows':len(v),'matrix_orthogonality_error':float(abs(R@R.T-np.eye(3)).max()),
    'J2000_equatorial_dipole_deg':[ra,dec],
    'galactic_vs_equatorial_dot_max_difference':float(abs(co-galco).max()),
    'proper_counts':counts(eqaxis),'rounded_Ray_counts':counts(xyz(np.array([167.8]),np.array([-7.1]))[0]),
    'rounded_Sah_code_counts':counts(xyz(np.array([168.]),np.array([-7.]))[0]),
    'wrong_Galactic_as_equatorial_counts':counts(galaxis),
    'frame_qualification':'Conventional J2000 rotation; tiny reference-frame differences from Astropy ICRS do not change these counts.'}
ray_sign=v@xyz(np.array([167.8]),np.array([-7.1]))[0]>=0
sah_sign=v@xyz(np.array([168.]),np.array([-7.]))[0]>=0
results['coordinates']['rows_switching_between_rounded_axes']=hem.loc[ray_sign!=sah_sign,['CID','IDSURVEY','RA','DEC','zHEL']].to_dict(orient='records')
assert results['coordinates']['galactic_vs_equatorial_dot_max_difference']<1e-12
output=OUT/'c2-crosscheck.json';output.write_text(json.dumps(results,indent=2)+'\n')
inputs=[fitpath,SRC/'Analysis C2/Pantheon+SH0ES.dat',SRC/'index_sorted_lane.npy',direction/'sources/cov_final.npy',SRC/'median_deltaage.csv',Path(__file__),ROOT/'docs/experiments/independent-audit-plan.md',direction/'sources/2607.20570.txt',direction/'sources/2607.20570.pdf',direction/'sources/2607.20570.html',direction/'sources/2608.02484.txt']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'time_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Independent equivalent-coordinate pointwise C2 and sky-coordinate audit',
          'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'outputs_sha256':{str(output.relative_to(ROOT)):sha(output)},
          'code_imports':'No directional investigator implementation imported','fit_snapshot_age_cases':[f['age_correction'] for f in fits],
          'scope':'pointwise likelihood and means only; no independent optimum search or significance calibration'}
(OUT/'c2-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
