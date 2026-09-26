"""Read-only archived-state and explicitly synthetic-log checker tests. No fits."""
from pathlib import Path
import importlib.util,json,hashlib
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('runner',H/'run_start_only.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
original=m.P/'fits-restricted/joint12';by,detail=m.previous.base_block_review(original,117,m.CIDS,True);dom=m.previous.domain_review(original,by,m.CIDS)
# Generalized identity checker must accept complete exact9iteration saved states.
nine=m.P/'fits-restricted/joint9';identity9=m.identity.compare(nine,nine,m.CIDS,117,True,9)
w=H/'executor-unit-fixture';w.mkdir(exist_ok=False)
(w/'README.txt').write_text('SYNTHETIC START-LOG UNIT FIXTURE; no new native process or optimization. Prior/objective tokens come from already saved nominal outputs, and start/readback tokens are constructed solely to test schema/assertions. Never use as a scientific fit.\n')
(w/'execution.json').write_text('{"minuit_peak_shift":2}\n');(w/'fit.nml').write_text('NFIT_ITERATION = 12\n')
log=[]
for tag in ['CSP_ENTRY:','CSP_OBJECTIVE:','PROSP_PEAK_BOUNDS']:log += [' '.join(x) for x in m.records(original/'native.log',tag)]
for c in m.CIDS:
 lo,hi=next(x['absolute_bounds'] for x in dom['details'] if x['CID']==c)
 for it in range(1,13):
  cb=next(x for x in by[c] if x['ITER']==it);source=cb['entry']['peak_entry_absolute'];shift=2 if it==1 else 0;proposed=source+shift
  aa=[c,str(it)]+[format(x,'.17g') for x in [2,shift,source,proposed,source,lo,hi,proposed-lo,hi-proposed]]
  bb=[c,str(it),'1','PKMJD']+[format(x,'.17g') for x in [source,proposed,proposed,source,lo,hi,lo,hi]]
  log += ['PROSP_MNPARM_PEAK '+' '.join(aa),'PROSP_MNPOUT_PEAK '+' '.join(bb)]
text='\n'.join(log)+'\n';(w/'native.log').write_text(text)
r=m.shift_review(w,by,m.CIDS,dom);assert r['active'] and len(r['details'])==8
checks=['exact9iteration saved-state identity checker pass','complete synthetic96readback records pass']
for case in ['missing_actual_shift','changed_common_entry','changed_later_readback']:
 lines=text.splitlines()
 for k,l in enumerate(lines):
  a=l.split()
  if case=='missing_actual_shift' and a[:3]==['PROSP_MNPOUT_PEAK','1','1']:a[7]=a[5];lines[k]=' '.join(a);break
  if case=='changed_common_entry' and a[:3]==['CSP_ENTRY:','1','1']:a[9]=format(float(a[9])+1,'.17g');lines[k]=' '.join(a);break
  if case=='changed_later_readback' and a[:3]==['PROSP_MNPOUT_PEAK','1','2']:a[7]=format(float(a[7])+1e-6,'.17g');lines[k]=' '.join(a);break
 (w/'native.log').write_text('\n'.join(lines)+'\n')
 try:m.shift_review(w,by,m.CIDS,dom)
 except RuntimeError:checks.append(case+' correctly rejected')
 else:raise AssertionError(case)
(w/'native.log').write_text(text)
# NML byte equality for both baseline iteration counts, without creating fit dirs.
for it in [9,12]:
 expected=(m.P/'inputs-v2/fit-template-v3.nml').read_text().format(data_path='../../readme-adapted/ledger',iterations=it,peak_step=2,peak_initializer='INIVAL_PEAKMJD = 57707.80078125',bands='grizJH')
 assert expected==(m.P/'fits-restricted'/('joint'+str(it))/'fit.nml').read_text()
checks.append('nominal NML exact for9and12 before localshift')
out={'pass':True,'native_calls':0,'checks':checks,'identity9_callbacks':identity9['callbacks'],'synthetic_start_records':96,'runner_sha256':hashlib.sha256((H/'run_start_only.py').read_bytes()).hexdigest()}
(H/'executor-unit-result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
