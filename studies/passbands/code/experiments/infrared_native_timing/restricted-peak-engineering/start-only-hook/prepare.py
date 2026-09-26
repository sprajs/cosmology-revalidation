from pathlib import Path
import difflib,hashlib,json
Q=Path(__file__).resolve().parent
s=(Q/'snana.base.car').read_text()
start=s.index('      SUBROUTINE MNFIT_DRIVER (');end=s.index('      END',start)
# Restrict replacements to the routine, never similarly named global arrays.
chunk=s[start:end]
old='''     &  ,ERRSYM      ! local SYMMMETRIC fiterr
'''
new=old+'''      DOUBLE PRECISION PROSP_MN_INITIAL,PROSP_MN_VALUE
     & ,PROSP_MN_ERROR,PROSP_MN_LOW,PROSP_MN_HIGH
      INTEGER PROSP_MN_CHECK,PROSP_MN_INDEX
      CHARACTER PROSP_MN_NAME*20
'''
assert chunk.count(old)==1;chunk=chunk.replace(old,new)
old='''+SELF,IF=MINUIT.
         CALL MNPARM(IPAR, PARNAME(IPAR)
     &       ,INIVAL(IPAR), INISTP(IPAR)
'''
new='''+SELF,IF=MINUIT.
c Separate local optimizer start; shared INIVAL/prior center unchanged.
         PROSP_MN_INITIAL = INIVAL(IPAR)
         PROSP_MN_CHECK = 0
         IF (PARNAME(IPAR) .EQ. 'PKMJD') THEN
           CALL PROSP_MINUIT_PEAK_START(CID,INIVAL(1),
     &       PARNAME(1)//CHAR(0),INIVAL(IPAR),INISTP(IPAR),
     &       INIBND(1,IPAR),INIBND(2,IPAR),PROSP_MN_INITIAL,
     &       PROSP_MN_CHECK)
         ENDIF
         CALL MNPARM(IPAR, PARNAME(IPAR)
     &       ,PROSP_MN_INITIAL, INISTP(IPAR)
'''
assert chunk.count(old)==1;chunk=chunk.replace(old,new)
s=s[:start]+chunk+s[end:]
old='''     &       ,IERR   ) 
+SELF.
         IF ( IERR .NE. 0 ) THEN'''
new='''     &       ,IERR   ) 
         IF (IERR.EQ.0 .AND. PROSP_MN_CHECK.EQ.1) THEN
           CALL MNPOUT(IPAR,PROSP_MN_NAME,PROSP_MN_VALUE,
     &       PROSP_MN_ERROR,PROSP_MN_LOW,PROSP_MN_HIGH,
     &       PROSP_MN_INDEX)
           CALL PROSP_MINUIT_PEAK_STORED(CID,INIVAL(1),
     &       INIVAL(IPAR),PROSP_MN_INITIAL,PROSP_MN_VALUE,
     &       INIBND(1,IPAR),INIBND(2,IPAR),PROSP_MN_LOW,
     &       PROSP_MN_HIGH,PROSP_MN_INDEX,
     &       TRIM(PROSP_MN_NAME)//CHAR(0))
         ENDIF
+SELF.
         IF ( IERR .NE. 0 ) THEN'''
assert s.count(old)==1;s=s.replace(old,new)
(Q/'source/snana.car').write_text(s)
c=(Q/'genmag_snoopy.base.c').read_text()+'\n/* Private, default-off MNPARM-local peak start check. */\n#include "prosp_minuit_start.h"\n'
(Q/'source/genmag_snoopy.c').write_text(c)
patch=''
for n in ['snana.car','genmag_snoopy.c']:
 patch+=''.join(difflib.unified_diff((Q/(n.split('.')[0]+'.base.'+n.split('.')[1])).read_text().splitlines(True),(Q/'source'/n).read_text().splitlines(True),fromfile='restricted/src/'+n,tofile='start-only/src/'+n))
patch+=''.join(difflib.unified_diff([], (Q/'source/prosp_minuit_start.h').read_text().splitlines(True),fromfile='/dev/null',tofile='start-only/src/prosp_minuit_start.h'))
(Q/'start-only.patch').write_text(patch)
print(hashlib.sha256(patch.encode()).hexdigest())
