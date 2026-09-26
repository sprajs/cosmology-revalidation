"""Numerical stability follow-up to measured native-versus-Python tangent differences. Same IDs, half steps; no new scientific score or chosen subset."""
from pathlib import Path
import os,sys,re,json,hashlib,subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;O=P/'sim_native_half';E=ROOT/'runs/research_2026_09_26/simulation_residual_control/engineering8'
EXE=ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe';PY=ROOT/'phase2/env-official/bin/python';EXPORT=ROOT/'scripts/research_2026_09_26/export_flux_objectives.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Fixed extremes of engineering8 measured redshift; no score-based selection.
selection={}
for arm in ['P21','G10']:
 m=pd.read_csv(P/f'simulation_design/{arm}-cohort.csv').set_index('CID');ids=list(map(int,(P/f'simulation_design/{arm}-engineering8-cids.txt').read_text().split()));q=m.loc[ids].sort_values('zHEL');selection[arm]=[int(q.index[0]),int(q.index[-1])]
steps={1:.0005,2:.00005,3:.005}
if not O.exists():
 O.mkdir();protocol={'selection':selection,'rule':'Lowest and highest measured zHEL in frozen engineering8 per arm','steps_x1_c_t0':steps,'laws':[-99,99],'native_coordinates':'Nominal free-fit final x0,x1,c,t0 plus/minus only one declared increment; LFIXPAR_ALL=T for derivative measurement','fixed_frozen_inputs':'Same version,KCOR,model,ordinary cuts; require exact accepted band/time/data/error identity with free baseline','covariance':'Nominal free-fit C retained for tangent and projection comparison; varied C not used','comparison':'Python modern F99 tangent as actual Sol scoring approximation, law-matched Python diagnostic, exact native finite derivative; native amplitude exact','scope':'8 object-law combinations spanning measured engineering z range; not a global derivative-error bound or selection null','source_sha256':sha(__file__),'binary_sha256':sha(EXE)};(O/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
else:raise RuntimeError('preserve existing derivatives')
jobs=[]
for arm,ids in selection.items():
 for law in ['approx_minus99','exact_99']:
  baseline=E/arm/law
  for j,h in steps.items():
   for sign in [-1,1]:
    d=O/arm/law/f'p{j}_{"plus" if sign>0 else "minus"}';d.mkdir(parents=True);seed=d/'seed.FITRES';lines=['VARNAMES: CID PKMJD x0 x1 c']
    for cid in ids:
     pars=np.load(baseline/f'objectives/objective_{cid}.npz')['parameters_x0_x1_c_t0'].copy();pars[j]+=sign*h;x0,x1,c,t0=pars;lines.append(f'SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}')
    seed.write_text('\n'.join(lines)+'\n');nml=(baseline/'fit.nml').read_text();nml=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f"    SNCID_LIST_FILE = '{seed}'",nml);nml=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f"    TEXTFILE_PREFIX = '{d/'fit'}'",nml);nml=nml.replace('&FITINP','&FITINP\n    LFIXPAR_ALL = T',1);(d/'fit.nml').write_text(nml);(d/'cids.txt').write_text('\n'.join(map(str,ids))+'\n');jobs.append(d)
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(d):
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,d
 r=subprocess.run([str(PY),str(EXPORT),'--log',str(d/'fit.log'),'--output',str(d/'objectives'),'--expected-cids',str(d/'cids.txt')],cwd=ROOT,env=env,capture_output=True,text=True);(d/'export.log').write_text(r.stdout+r.stderr);assert r.returncode==0,(d,r.stderr)
 return str(d.relative_to(O))
with ThreadPoolExecutor(max_workers=2) as pool:
 for d in pool.map(run,jobs):print('finished',d,flush=True)
# Comparison is separate so a failed mask gate leaves source artifacts intact.
sys.path.insert(0,str(ROOT));from scripts.salt_dust_audit import flux_response as fr
from scripts.salt_dust_audit.snana_extinction import SnanaExtinction
ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
class OldLaw(fr.VariableF99):
 def propagate(self,wave,flux,phase=None):
  ebv,rv=self._parameters;return flux*10**(-.4*ext(wave,rv,ebv,option=-99))
model,bands,paths,zp=fr.build_model();old=fr.sncosmo.Model(source=model.source,effects=[OldLaw(),fr.VariableF99()],effect_names=['mw','host'],effect_frames=['obs','rest']);offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
basis=np.load(P/'validation1020/frozen-discovery-coefficients.npz');coeff=basis['gauge_griz']@basis['basis_mean'];rows=[];arrays={}
for arm,ids in selection.items():
 for law in ['approx_minus99','exact_99']:
  for cid in ids:
   a=np.load(E/arm/law/f'objectives/objective_{cid}.npz');J=np.zeros((len(a['MJD']),4));J[:,0]=-fr.K*a['model_flux'];maxpar=0.
   for j,h in steps.items():
    pp=[]
    for sign in [-1,1]:
     b=np.load(O/arm/law/f'p{j}_{"plus" if sign>0 else "minus"}'/f'objectives/objective_{cid}.npz')
     for key in ['MJD','band','data_flux','data_fluxerr']:assert np.array_equal(a[key],b[key]),(arm,law,cid,j,sign,key)
     expect=a['parameters_x0_x1_c_t0'].copy();expect[j]+=sign*h;err=float(np.max(abs(expect-b['parameters_x0_x1_c_t0'])));maxpar=max(maxpar,err);assert err<1e-12
     pp.append(b['model_flux'])
    J[:,j]=(pp[1]-pp[0])/(2*h)
   C=a['frozen_flux_covariance'];L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True);U,s,V=np.linalg.svd(jw,full_matrices=False);assert (s>s[0]*1e-10).sum()==4
   raw=solve_triangular(L,-fr.K*a['model_flux']*np.array([coeff['griz'.index(b)] for b in a['band']]),lower=True);v=raw-U@(U.T@raw);rw=solve_triangular(L,a['data_flux']-a['model_flux'],lower=True)
   native_a=float(v@rw);native_I=float(v@v);x0,x1,c,t0=a['parameters_x0_x1_c_t0'];bp=np.array([bands[b] for b in a['band']],dtype=object);conv=np.array([10**(-.4*(.27+offsets[b])) for b in a['band']]);z=float(a['zHEL'][0]);ebv=float(a['MWEBV'][0])
   for label,m in [('actual_Python_exact99',model),('law_matched_Python',old if law=='approx_minus99' else model)]:
    def f(th):
     m.set(z=z,t0=t0+th[3],x0=x0*np.exp(-fr.K*th[0]),x1=x1+th[1],c=c+th[2],mwebv=ebv,mwrv=3.1,hostebv=0,hostrv=3.1);return conv*m.bandflux(bp,a['MJD'],zp=27.5,zpsys='ab')
    jj=J.copy()
    for j,h in steps.items():
     e=np.eye(4)[j]*h;jj[:,j]=(f(e)-f(-e))/(2*h)
    w=solve_triangular(L,jj,lower=True);uu,ss,_=np.linalg.svd(w,full_matrices=False);vp=raw-uu@(uu.T@raw)
    row=dict(arm=arm,law=law,CID=cid,zHEL=z,epochs=len(v),comparison=label,parameter_max_error=maxpar,relative_whitened_derivative_norm=float(np.linalg.norm(w-jw)/np.linalg.norm(jw)),projector_spectral_difference=float(np.linalg.norm(U@U.T-uu@uu.T,2)),fixed_vector_difference_norm=float(np.linalg.norm(vp-v)),fixed_vector_relative_difference=float(np.linalg.norm(vp-v)/np.linalg.norm(v)),native_information=native_I,native_matched=native_a,Python_minus_native_information=float(vp@vp-native_I),Python_minus_native_matched=float(vp@rw-native_a));row['Python_minus_native_gain']=row['Python_minus_native_matched']-.5*row['Python_minus_native_information'];rows.append(row)
   arrays[f'{arm}_{law}_{cid}_J_native']=J
pd.DataFrame(rows).to_csv(O/'comparison.csv',index=False);np.savez_compressed(O/'native-jacobians.npz',**arrays)
result={'status':'PASS: exact same accepted rows and coordinates for all native finite differences','checks':rows,'source_sha256':sha(__file__),'protocol_sha256':sha(O/'protocol.json'),'scope':'Finite selected-coordinate tangent control, not full sample error bound'};(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(pd.DataFrame(rows).to_string(index=False),flush=True)
