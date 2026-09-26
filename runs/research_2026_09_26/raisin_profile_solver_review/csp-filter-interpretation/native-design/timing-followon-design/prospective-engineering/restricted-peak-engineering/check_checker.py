"""Synthetic log-schema tests only; no native calls or simulated observations."""
from pathlib import Path
import importlib.util,json
Q=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('restricted',Q/'run_restricted.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
w=Q/'checker-unit-fixture/noiseless_synthetic';w.mkdir(parents=True,exist_ok=False)
(w/'fit.nml').write_text('NFIT_ITERATION = 12\n')
(w/'execution.json').write_text('{"hard_domain_mode":1}\n')
base=next(x for x in (Q/'normal-1.log').read_text().splitlines() if x.startswith('PROSP_PEAK_BOUNDS ')).split()
b=[]
for i in range(1,13):
 x=base.copy();x[1]='1';x[2]=str(i);b.append(' '.join(x))
summ='PROSP_DOMAIN_SUMMARY 1 1 12 235 1 -20 70 0 57707.80078125 57707.80078125 7.2662189311959082'
support='PROSP_SUPPORT 1 235 0 0 -20 70 0 0.453000009059906 1 1'
log='\n'.join(b+[summ,support])+'\n';(w/'native.log').write_text(log)
by={'1':[{'ITER':i,'objective':{'peak_absolute':57707.80078125}} for i in range(1,13)]}
r=m.domain_review(w,by,['1']);assert r['enabled']
checks=['complete synthetic schema and exact endpoints pass']
for name,change in [('missing_bound',lambda s:s.replace(b[0]+'\n','')),('wrong_call_count',lambda s:s.replace('12 235 1','12 234 1')),('guard_abort',lambda s:s+'PROSP_DOMAIN_ABORT 1 synthetic_test\n')]:
 (w/'native.log').write_text(change(log))
 try:m.domain_review(w,by,['1'])
 except RuntimeError:checks.append(name+' correctly rejected')
 else:raise AssertionError(name)
(w/'native.log').write_text(log)
by['1'][-1]['objective']['peak_absolute']=58000
try:m.domain_review(w,by,['1'])
except RuntimeError:checks.append('out-of-domain final correctly rejected')
else:raise AssertionError('invalid final accepted')
# Preserve this directory as an explicitly synthetic unit fixture, not a native run.
(w/'README.txt').write_text('SYNTHETIC CHECKER FIXTURE. No native process, photons, or fitted outcome. Log strings exercise schema/arithmetic checks only.\n')
r={'pass':True,'native_calls':0,'synthetic_only':True,'checks':checks}
(Q/'checker-unit-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
