"""Hash/static-only review; does not compile or invoke SNANA."""
from pathlib import Path
import hashlib,json,shutil,difflib
R=Path(__file__).resolve().parents[6];O=Path(__file__).resolve().parent
P=R/'phase2/pte/restricted-peak-engineering'; B=R/'phase2/pte/fit-support/build/src'
expected={'build/src/snlc_fit.car':'54a013388053d10e4998a356422dda77743d0f90dbd2158e44155960fc5d2236','build/src/genmag_snoopy.c':'e67cc924a2813c3c6813d5acb162639bdf2a90c8e94088c63ce824b75906c38f','build/src/prosp_peak_domain.h':'d0bddc2535fdcb71cd723a947e7ece64fbae2f74b4febbd7aeff9d37bdc81fd7'}
files=[]
for p in ['restricted-peak.patch','static-result.json','snlc_fit.base.car','genmag_snoopy.base.c',*expected]:
 x=P/p;h=hashlib.sha256(x.read_bytes()).hexdigest();assert p not in expected or h==expected[p]
 files.append({'path':str(x.relative_to(R)),'sha256':h,'bytes':x.stat().st_size})
for base,original in [('snlc_fit.base.car','snlc_fit.car'),('genmag_snoopy.base.c','genmag_snoopy.c')]:assert (P/base).read_bytes()==(B/original).read_bytes()
f=(P/'build/src/snlc_fit.car').read_text();c=(P/'build/src/genmag_snoopy.c').read_text();h=(P/'build/src/prosp_peak_domain.h').read_text()
assert f.count('CALL PROSP_BOUND_FITPAR(ITER)')==2
assert f.count('CALL PROSP_PEAK_FCN_GUARD')==1 and f.count('CALL PROSP_PHYSICAL_MEAN_GUARD')==1
assert f.index('CALL PROSP_PEAK_FCN_GUARD') < f.index('USRFUN ( ITER, IFILT_OBS, ZSN, Tobs')
assert f.index('CALL PROSP_PHYSICAL_MEAN_GUARD') < f.index('CALL PROSP_MEAN_SUPPORT')
assert c.endswith('#include "prosp_peak_domain.h"\n')
assert h.count('if(!prosp_domain_on()) return;')==3
assert 's->lo!=lo||s->hi!=hi||*lower!=lo||*upper!=hi' in h
assert 'prosp_domain_abort("physical_mean_outside_table_domain",cid);' in h
shutil.copyfile(P/'restricted-peak.patch',O/'reviewed-restricted-peak.patch')
shutil.copyfile(P/'build/src/prosp_peak_domain.h',O/'reviewed-prosp_peak_domain.h')
r={'verdict':'PASS for bounded private build; native identity/recovery remain unexecuted gates','scope':'Source-only independent review, no builds/native calls/pixel or flux outcomes','verified_base_equals_existing_support_source':True,'expected_source_hashes_match':True,'manual_reviews':['all117 stored R8 times, no flux or truth values in bounds','first hook after overrides, later hook before early return','fixed cadence/z/MJDOFF/tables/bounds asserted across iterations','FCN guard after prior rejection but before physical mean, includes SIGMA_ONLY bypass','USRFUN guard before template/extinction/KCOR and old output-only trace, applies to auxiliary rest calls in their actual phase','absent/0 environment disabled path has no native-array writes','enabled invalid calls abort, no phase clipping or row pruning','native relative/absolute peak conversions and R8/R4 endpoints consistent for current MJDOFF0','summary edge_peak is all-trial contact, final boundary flag must derive from final state'],'files':files}
(O/'patch-review.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
