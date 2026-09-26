from pathlib import Path
import hashlib,json,re,difflib
R=Path.cwd(); out=Path(__file__).resolve().parent
q=R/'phase2/pte/restricted-peak-engineering/start-only-hook'
b=R/'phase2/pte/restricted-peak-engineering/build/src'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
expected={'start-only.patch':'6d1b6ce87d632224c411e6d8cb80b167cae93485e3e72db0977d362285d5cd05','source/prosp_minuit_start.h':'38027fe0b92bd803c453c197d91064a4a2596c145ac248c2988ca9d8954506de','source/snana.car':'98783d22cae174ccbc1015c10d68641025abaac95d6dc78262998fc91a221f73'}
assert all(sha(q/k)==v for k,v in expected.items())
assert (q/'snana.base.car').read_bytes()==(b/'snana.car').read_bytes()
assert (q/'genmag_snoopy.base.c').read_bytes()==(b/'genmag_snoopy.c').read_bytes()
old=(q/'snana.base.car').read_text();new=(q/'source/snana.car').read_text()
a=old.index('      SUBROUTINE MNFIT_DRIVER (');z=old.index('      END   ! end of MNFIT_DRIVER',a) if '      END   ! end of MNFIT_DRIVER' in old[a:] else old.index('\n      END\n',a)
an=new.index('      SUBROUTINE MNFIT_DRIVER (');zn=new.index('\n      END\n',an)
assert old[:a]==new[:an] and old[z:]==new[zn:]
ch=new[an:zn]; h=(q/'source/prosp_minuit_start.h').read_text()
assert 'PROSP_MN_INITIAL = INIVAL(IPAR)' in ch
assert 'CALL MNPARM(IPAR, PARNAME(IPAR)\n     &       ,PROSP_MN_INITIAL, INISTP(IPAR)' in ch
assert 'CALL MNPOUT(IPAR,PROSP_MN_NAME,PROSP_MN_VALUE' in ch
# No shared INIVAL assignment added to the driver; local output arguments differ.
assert not re.search(r'^\s*INIVAL\([^\n]*\)\s*=',ch,re.M)
assert 'if(applied!=0.0) *local_start=proposal;' in h
assert 'if(prosp_mn_shift==0.0) return;' in h
assert '(*iteration!=1.0 && *stored!=*source_initial)' in h
records=[]
for case,shift in [('absent',0),('zero',0),('minus',-2),('plus',2)]:
 lines=(q/(case+'.log')).read_text().splitlines()
 cap=[l.split() for l in lines if l.startswith('UNIT_CAPTURE')]
 assert len(cap)==48
 for c in cap:
  it,ip=map(int,c[1:3]); source,local,stored=map(float,c[3:6]);delta=shift if it==1 and ip==3 else 0
  assert local==source+delta and stored==source+delta
 records.append({'case':case,'captured_native_MNPOUT_values':len(cap),'max_delta_error':0})
# Test results belong to peer; independently inspect source hashes, captures and abort logs.
t=json.loads((q/'test-result.json').read_text())
assert t['status']=='PASS'
for k,v in t['source_sha256'].items(): assert sha(q/k)==v
for c in t['checks']: assert c['pass'] and c['returncode']==c['expected']
paths=[q/k for k in expected]+[q/n for n in ['snana.base.car','genmag_snoopy.base.c','test-result.json','unit_driver.f90','test.py','absent.log','zero.log','minus.log','plus.log']]+[b/'snlc_fit.car',b/'minuit.F']
res={'status':'PASS_STATIC_AND_SAVED_TEST_REVIEW','scope':'No native build, fits, photons or rerun of collaborator executables','patch_sha256':expected['start-only.patch'],'checks':records,'saved_peer_cases':len(t['checks']),'science_qualification':'First-iteration objective function and shared prior/covariance initialization are unchanged. ITER1 model-dependent W may differ at different coordinates; later native covariance and recentered prior can differ with trajectory. This is iterative-estimator start stability, not a fixed-objective global-minimum certificate.','inputs':{str(p.relative_to(R)):sha(p) for p in paths}}
(out/'result.json').write_text(json.dumps(res,indent=2)+'\n')
(out/'reviewed-start-only.patch').write_bytes((q/'start-only.patch').read_bytes())
(out/'reviewed-prosp_minuit_start.h').write_bytes((q/'source/prosp_minuit_start.h').read_bytes())
print(json.dumps({'status':res['status'],'patch':res['patch_sha256'],'cases':res['saved_peer_cases']}))
