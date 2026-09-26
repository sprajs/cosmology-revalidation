"""Prespecified two-stage changed-input pilot; amplitudes derived by fixed rule."""
from pathlib import Path
import json,hashlib,shutil,os,subprocess,time,re
import numpy as np
from parse_native import parse
ROOT=Path.cwd();O=Path(__file__).resolve().parent;B=O/'SNANA-v11_04k-output';P=O/'changed-pilot-probes-v3';S=ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-preparation/pilot';CP=S.parent/'protocol.json';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name,h in json.loads(CP.read_text())['inputs_sha256'].items():assert sha(ROOT/name)==h,name
assert not P.exists();P.mkdir();(P/'data').symlink_to(S/'data',target_is_directory=True)
base=(S/'fits/known_Jdw_to_j/fit.nml').read_text();initial=[dict(name='nominal',shift=None,type='fitted'),dict(name='start_minus',shift=-.2,type='fitted'),dict(name='start_plus',shift=.2,type='fitted')]
def prep(a,nml):
 d=P/'fits'/a['name'];d.mkdir(parents=True);(d/'fit.nml').write_text(nml);shutil.copyfile(S/'fits/known_Jdw_to_j/vpec.list',d/'vpec.list');return d
for a in initial:prep(a,base)
inputs=[Path(__file__),O/'parse_native.py',B/'bin/snlc_fit.exe',B/'src/snlc_fit.car',CP,O/'entry-source-proof.json',O/'probes-v3/result.json']+list((P/'fits').rglob('fit.nml'))+list((P/'fits').rglob('vpec.list'))+list((S/'data/known_Jdw_to_j/CSPDR3_RAISIN').glob('*'))
proto={'state':'Frozen before any changed-pilot output. No full42 changed-response execution.','scope':'Only prepared6 Jdw-to-j tokens for2004ef;2005hc exactunchangedcontrol. Native same historical empirical model, fixed shape/AV/headerpeak.','initial_arms':initial,'fixed_arm_rule':'After default succeeds with exact2005hc science/LCPLOT control, for each of the two CIDs create LFIXPAR_ALL=T,INISTP_DLMAG=0,INIVAL_DLMAG=float32(default final native D+offset) with offsets0,−.15,+.15. No other outcome adaptation; actual exported seeddefinesratio.','gates':{'Qclosure_absolute':1e-8,'starts':'Native first-entry shifts±.2 within1e−12; finalD samewithin.001; samephysicalmask','last_two_D':.001,'final_frozen_C_D':.001,'native_mean_scaling_perrow_relative':2e-8,'support':'Every callback phase inside actual−20..70 andshape.7..1.3; nonzero throughput endpoints inside declarednative2900..20000 screen; positive means/C; no clamp. This is numeric support, not empirical adequacy.','control':'2005hc exactscience andLCPLOT bytes versus originalnominal','physical_rows':'Every callback matches baseline physicalrow MJD/flux/error multiset. Bands unchanged or exactly the6knownJ→j. Parent independentledgercheck required.'},'hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs if p.is_file()},'per_invocation_timeout_seconds':30,'total_active_seconds_cap':120,'scope_warning':'Raw processing response, no BBC/cosmology. Large empirical model residuals remain despite numerical closure.'}
(P/'protocol.json').write_text(json.dumps(proto,indent=2)+'\n');print('frozen',sha(P/'protocol.json'),flush=True)
env=os.environ.copy();env.update(SNANA_DIR='/tmp/snana-csp-audit-8c5f0d',SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('CSP_POSTINIT_DSHIFT',None)
used=0;results=[]
def run(a):
 global used
 d=P/'fits'/a['name'];e=env.copy()
 if a['shift'] is not None:e['CSP_POSTINIT_DSHIFT']=str(a['shift'])
 t=time.monotonic()
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(B/'bin/snlc_fit.exe'),'fit.nml'],cwd=d,env=e,stdout=f,stderr=subprocess.STDOUT,timeout=min(30,120-used))
 elapsed=time.monotonic()-t;used+=elapsed;x=dict(name=a['name'],returncode=r.returncode,seconds=elapsed,log_sha256=sha(d/'fit.log'));results.append(x);(d/'execution.json').write_text(json.dumps(x,indent=2)+'\n');(P/'execution-results.json').write_text(json.dumps({'arms':results,'total_seconds':used},indent=2)+'\n');print(x,flush=True)
 if r.returncode:raise RuntimeError('Native failure preserved; later arms not run')
run(initial[0]);d=P/'fits/nominal';orig=ROOT/'runs/research_2026_09_26/csp_native_filter_response/fits/nominal'
control=lambda p:[l for l in p.read_text().splitlines() if l.split() and (l.startswith('VARNAMES:') or len(l.split())>1 and l.split()[1]=='2005hc')]
for name in ['fit.FITRES.TEXT','fit.LCPLOT.TEXT']:assert control(d/name)==control(orig/name),name
last={b['CID']:b for b in parse(d/'fit.log')};assert set(last)=={'2004ef','2005hc'};fixed=[]
for cid,b in last.items():
 for label,offset in [('center',0),('minus',-.15),('plus',.15)]:
  dseed=float(np.float32(b['objective']['D']+offset));nml=re.sub(r"(?m)^(\s*SNCCID_LIST\s*=).*$",lambda m:m[1]+f" '{cid}'",base);nml=nml.replace('&FITINP',f'&FITINP\n LFIXPAR_ALL = T\n INIVAL_DLMAG = {dseed:.17g}\n INISTP_DLMAG = 0.0');a=dict(name=f'{cid}_{label}',shift=None,type='fixed',CID=cid,requested_offset=offset,D_seed_float32=dseed,reference_D=b['objective']['D']);prep(a,nml);fixed.append(a)
seedledger={'rule_source_protocol_sha256':sha(P/'protocol.json'),'arms':initial+fixed,'derived_NML_sha256':{str(q.relative_to(P)):sha(q) for q in (P/'fits').rglob('fit.nml')}};(P/'derived-fixed-probes.json').write_text(json.dumps(seedledger,indent=2)+'\n')
for a in initial[1:]+fixed:run(a)
