from pathlib import Path
import difflib,hashlib,json
Q=Path(__file__).resolve().parent
s=(Q/'snana.base.car').read_text()
start=s.index('      SUBROUTINE MNFIT_DRIVER (');end=s.index('      END',start)
# Restrict replacements to the routine, never similarly named global arrays.
chunk=s[start:end]
old='''     &  ,ERRSYM      ! local SYMMMETRIC fiterr
'''
new=old+'''      DOUBLE PRECISION PROSP_MN_INITIAL
'''
assert chunk.count(old)==1;chunk=chunk.replace(old,new)
old='''+SELF,IF=MINUIT.
         CALL MNPARM(IPAR, PARNAME(IPAR)
     &       ,INIVAL(IPAR), INISTP(IPAR)
'''
new='''+SELF,IF=MINUIT.
c Separate local optimizer start; shared INIVAL/prior center unchanged.
         PROSP_MN_INITIAL = INIVAL(IPAR)
         IF (PARNAME(IPAR) .EQ. 'PKMJD') THEN
           CALL PROSP_MINUIT_PEAK_START(CID,INIVAL(1),
     &       PARNAME(1)//CHAR(0),INIVAL(IPAR),INISTP(IPAR),
     &       INIBND(1,IPAR),INIBND(2,IPAR),PROSP_MN_INITIAL)
         ENDIF
         CALL MNPARM(IPAR, PARNAME(IPAR)
     &       ,PROSP_MN_INITIAL, INISTP(IPAR)
'''
assert chunk.count(old)==1;chunk=chunk.replace(old,new)
s=s[:start]+chunk+s[end:]
(Q/'source/snana.car').write_text(s)
c=(Q/'genmag_snoopy.base.c').read_text()+'\n/* Private, default-off MNPARM-local peak start check. */\n#include "prosp_minuit_start.h"\n'
(Q/'source/genmag_snoopy.c').write_text(c)
patch=''
for n in ['snana.car','genmag_snoopy.c']:
 patch+=''.join(difflib.unified_diff((Q/(n.split('.')[0]+'.base.'+n.split('.')[1])).read_text().splitlines(True),(Q/'source'/n).read_text().splitlines(True),fromfile='restricted/src/'+n,tofile='start-only/src/'+n))
patch+=''.join(difflib.unified_diff([], (Q/'source/prosp_minuit_start.h').read_text().splitlines(True),fromfile='/dev/null',tofile='start-only/src/prosp_minuit_start.h'))
(Q/'start-only.patch').write_text(patch)
print(hashlib.sha256(patch.encode()).hexdigest())
