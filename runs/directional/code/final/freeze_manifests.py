"""Freeze code lineage and complete per-experiment provenance after computations.

C2 began before an unused optional optimizer-bound argument was added to audit.py.
Recover that exact imported helper and verify its previously recorded SHA, rather
than falsely declaring the subsequently edited helper to be what was executed.
"""
import pathlib,hashlib,json,shutil,datetime,platform,subprocess
import numpy,pandas,scipy,astropy

ROOT=pathlib.Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/directional';SRC=ROOT/'sources/repos/Shin107__Anisotropy-in-Pantheon-Plus'
CODE=ROOT/'scripts/directional'
def record(f):
    f=pathlib.Path(f)
    return dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
snap=OUT/'code/final';snap.mkdir(parents=True,exist_ok=True)
for f in CODE.glob('*.py'):shutil.copy2(f,snap/f.name)
shutil.copy2(ROOT/'docs/experiments/directional-plan.md',snap/'directional-plan.md')
old_dir=OUT/'code/c2-executed';old_dir.mkdir(parents=True,exist_ok=True)
old=(CODE/'audit.py').read_text().replace('def fit(self,y,dip=False,scatter=False,start=None,multistart=True,wide=False):','def fit(self,y,dip=False,scatter=False,start=None,multistart=True):').replace('        if wide:bounds[:2]=[(-15,15),(-2000,2000)]\n','')
expected=json.loads((OUT/'initial-results/manifest-data.json').read_text())['inputs']['scripts/directional/audit.py']
assert hashlib.sha256(old.encode()).hexdigest()==expected,'Historical helper reconstruction must match recorded SHA'
(old_dir/'audit.py').write_text(old);shutil.copy2(CODE/'c2.py',old_dir/'c2.py')
versions=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,pandas=pandas.__version__,astropy=astropy.__version__)
base=dict(frozen_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_head_at_freeze=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),environment=versions,blas_threads=1,plan=record(snap/'directional-plan.md'),preregistration_note='Final append-only plan includes dated amendments; original execution record establishes pre-outcome ordering.')
c1inputs=[SRC/'Analysis C1/Z_mbcorr.csv',SRC/'Analysis C1/Pantheon+SH0ES.dat',SRC/'Analysis C1/statsys_mbcorr.npy']
experiments={
 'c1':dict(code=['audit.py'],inputs=c1inputs+[SRC/'median_deltaage.csv'],outputs=['c1-fits.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py data',seed=None,config=dict(cut='.00937<zHEL<.8',frames=['zHEL','zCMB','zHD','zLG'],reverse_bias=[False,True],age_correction=[False,True],likelihood='Gaussian normalized, free additional magnitude scatter, profiled offset',dipole='fixed axis, exponential; S lower bound is first selected frame redshift')),
 'synthetic':dict(code=['audit.py','summarize_boundaries.py'],inputs=c1inputs,outputs=['taylor-errors.json','synthetic-summary.json','synthetic-realizations.csv','synthetic-boundaries.csv'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py synthetic; then .venv/bin/python scripts/directional/summarize_boundaries.py',seed=260609650,config=dict(models=['flat LCDM Om=.3','EdS','flat constant q=0'],zmin=.01,zmax=[.1,.2,.4,.8],frame='zHD',prefactor='zHEL',noise_realizations_per_model_cut=100,covariance='C1 total; no extra scatter',bounds={'q':[-3,3],'J':[-20,20]})),
 'lowz-wide':dict(code=['audit.py','lowz_wide.py'],inputs=c1inputs,outputs=['lowz-wide-summary.json','lowz-wide-realizations.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/lowz_wide.py',seed=260609650,config=dict(cuts=[.1,.2],bounds={'q':[-15,15],'J':[-2000,2000]},purpose='predeclared response to jerk boundary occupancy; same paired seeds')),
 'optimizer-flags':dict(code=['audit.py','check_optimizer_flags.py'],inputs=c1inputs+[OUT/'lowz-wide-realizations.json',OUT/'synthetic-realizations.csv'],outputs=['synthetic-optimizer-flag-check.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/check_optimizer_flags.py',seed=260609650,config=dict(purpose='Independent Nelder-Mead reruns only of six non-success optimizer statuses')),
 'c2':dict(code=['c2.py','audit.py'],inputs=[SRC/'Analysis C2/Pantheon+SH0ES.dat',SRC/'index_sorted_lane.npy',OUT/'sources/cov_final.npy',SRC/'median_deltaage.csv'],outputs=['c2-fits.json','c2-gradient-check.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/c2.py',seed=None,config=dict(cut='.00937<zHEL<.8',frame='zHEL',n=1547,reverse_bias=True,age_correction=[False,True],means='profiled exactly',nuisance='alpha,beta,sigmaM,sigmaX,sigmaC optimized',dipole='RA168 DEC-7, S>=.00938',starts=2,code_lineage='c2-executed/audit.py is verified reconstruction of the actually imported pre-wide helper')),
 'ray-and-coordinates':dict(code=['audit.py','ray_check.py'],inputs=[SRC/'Analysis C1/Pantheon+SH0ES.dat',SRC/'median_deltaage.csv'],outputs=['ray-v2-simple-fit.json','response-coordinate-check.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/ray_check.py',seed=None,config=dict(cut='.00937<zHEL<=.8',J=1,diagonal_sigma_mag=.15,frames=['zHEL','zCMB','zHD'],subsets=['all','dipole','antidipole'],slopes='alpha/beta and offset jointly linear-profiled',age_correction=[False,True])),
 'input-audit':dict(code=['audit.py','input_audit.py'],inputs=c1inputs+[SRC/'Analysis C2/Zpan.csv',SRC/'index_sorted_lane.npy',OUT/'sources/cov_final.npy',ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov'],outputs=['input-audit.json','c2-input-audit.json'],command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/input_audit.py',seed=None,config=dict(pickle='never deserialized',author_python='syntax parsed only; not executed')),
}
for name,e in experiments.items():
    code_dir=old_dir if name=='c2' else snap
    m=dict(**base,experiment=name,command=e['command'],seed=e['seed'],configuration=e['config'],inputs=[record(p) for p in e['inputs']],code=[record(code_dir/n) for n in e['code']],outputs=[record(OUT/n) for n in e['outputs']])
    (OUT/('experiment-'+name+'.json')).write_text(json.dumps(m,indent=2)+'\n')
# Earlier lightweight C2 manifest was emitted at program exit, after the unrelated
# helper edit. Keep explicit actual code lineage consistent in the active manifest.
m=json.loads((OUT/'manifest-c2.json').read_text())
m['inputs'].pop('scripts/directional/audit.py',None)
m['inputs']['runs/directional/code/c2-executed/audit.py']=expected
m['actual_code_snapshot']=[record(old_dir/'audit.py'),record(old_dir/'c2.py')]
m['note']='Actual C2-imported helper differs only by an unused optional C1 optimizer widening added during this run; use verified snapshot for exact historical code.'
(OUT/'manifest-c2.json').write_text(json.dumps(m,indent=2)+'\n')
print('Frozen',len(experiments),'experiment manifests and verified actual C2 code lineage.')
