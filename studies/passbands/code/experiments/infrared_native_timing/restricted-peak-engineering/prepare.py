from pathlib import Path
import difflib,hashlib,json
P=Path(__file__).resolve().parent
b=P/'build/src'
s=(P/'snlc_fit.base.car').read_text()
old='      IF ( .NOT. FIRST_ITER ) RETURN\n\n      SIMCHI2_CHEAT'
new='''      IF ( .NOT. FIRST_ITER ) THEN
         CALL PROSP_BOUND_FITPAR(ITER)
         RETURN
      ENDIF

      SIMCHI2_CHEAT'''
assert s.count(old)==1;s=s.replace(old,new)
old='''c set FITVAL on 1st interation since this array
'''
new='''c Private optional fixed metadata-domain estimator: after all overrides.
      CALL PROSP_BOUND_FITPAR(ITER)

c set FITVAL on 1st interation since this array
'''
assert s.count(old)==1;s=s.replace(old,new)
old='''C ==========================================
+DECK,FITINI_PHOTOZ.'''
new='''C ==========================================
+DECK,PROSP_BOUND_FITPAR.
      SUBROUTINE PROSP_BOUND_FITPAR(ITER)
      IMPLICIT NONE
+CDE,SNDATCOM.
+CDE,SNANAFIT.
+CDE,SNFITCOM.
+CDE,SNLCINP.
+CDE,FILTCOM.
      INTEGER ITER
c C helper returns without any native writes unless explicit mode is 1.
c All stored R8 epochs are used, including bands omitted by a NIR fit.
      CALL PROSP_PEAK_BOUNDS(SNLC_CCID//CHAR(0),ITER,
     & ISNLC_NEPOCH_STORE,SNLC8_MJD,INIVAL(IPAR_zPHOT),
     & INISTP(IPAR_zPHOT),DBLE(MJDOFF),INIVAL(IPAR_PEAKMJD),
     & DBLE(Trange_KCOR(1)),DBLE(Trange_KCOR(2)),
     & DBLE(Zrange_KCOR(2)),FITMODEL_INDEX,MODEL_SNOOPY,
     & INIBND(1,IPAR_PEAKMJD),INIBND(2,IPAR_PEAKMJD))
      RETURN
      END

C ==========================================
+DECK,FITINI_PHOTOZ.'''
assert s.count(old)==1;s=s.replace(old,new)
old='''      if ( LFLAG_PRIOR_ONLY ) RETURN

c -------------------------------------------'''
new='''      if ( LFLAG_PRIOR_ONLY ) RETURN

c Enabled-only fixed-domain assertion before any physical FCN mean.
      CALL PROSP_PEAK_FCN_GUARD(SNLC_CCID//CHAR(0),ITER,
     & PEAKMJD,DBLE(MJDOFF),ZSN)

c -------------------------------------------'''
assert s.count(old)==1;s=s.replace(old,new)
old='''c Output-only prospective support trace. No RNG or model call.
'''
new='''c Enabled-only hard refusal before any physical mean operation.
      CALL PROSP_PHYSICAL_MEAN_GUARD(SNLC_CCID//CHAR(0),
     & ITER,IFILT_OBS,IFITDATA_USRFUN,Trest,Tobs,ZSN,
     & DIST,SHAPE(1),AVHOST,RVHOST,
     & DBLE(Trange_KCOR(1)),DBLE(Trange_KCOR(2)),
     & DBLE(Zrange_KCOR(2)),FITMODEL_INDEX,MODEL_SNOOPY)

c Output-only prospective support trace. No RNG or model call.
'''
assert s.count(old)==1;s=s.replace(old,new)
(b/'snlc_fit.car').write_text(s)
c=(P/'genmag_snoopy.base.c').read_text()+'\n/* Separate default-off supported peak estimator. */\n#include "prosp_peak_domain.h"\n'
(b/'genmag_snoopy.c').write_text(c)
patch=''
for name,base,newfile in [('snlc_fit.car',P/'snlc_fit.base.car',b/'snlc_fit.car'),('genmag_snoopy.c',P/'genmag_snoopy.base.c',b/'genmag_snoopy.c')]:
 patch+=''.join(difflib.unified_diff(base.read_text().splitlines(True),newfile.read_text().splitlines(True),fromfile='base/src/'+name,tofile='restricted/src/'+name))
patch+=''.join(difflib.unified_diff([], (b/'prosp_peak_domain.h').read_text().splitlines(True),fromfile='/dev/null',tofile='restricted/src/prosp_peak_domain.h'))
(P/'restricted-peak.patch').write_text(patch)
