from pathlib import Path
import hashlib,json,difflib
D=Path(__file__).resolve().parent;O=D.parent;S=O.parent/'restricted-peak-engineering/start-only-hook/build/src'
def once(s,a,b):assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)
a=(S/'snana.car').read_text();b=a
b=once(b,'      CHARACTER PROSP_MN_NAME*20\n','''      CHARACTER PROSP_MN_NAME*20
c Output-only MINUIT diagnostics; no shared-state or model changes.
      LOGICAL PROSP_TRACE_ACTIVE,PROSP_STATUS_TRACE
      EXTERNAL PROSP_TRACE_ACTIVE
      INTEGER PROSP_MIN_RETURN
''')
b=once(b,'      IERR = 0\n      NFIXPAR = 0\n      FITCHI2 = 0.0\n','''      IERR = 0
      NFIXPAR = 0
      FITCHI2 = 0.0
      PROSP_STATUS_TRACE = .FALSE.
      IF (CID.EQ.22 .OR. CID.EQ.64) THEN
        PROSP_STATUS_TRACE = PROSP_TRACE_ACTIVE()
      ENDIF
      PROSP_MIN_RETURN = 0
''')
b=once(b,"      PRINT *,' ------------------------------------------------ '\n\nc -------------------------------------------\nc get the fit result", """      IF (PROSP_STATUS_TRACE) PROSP_MIN_RETURN = IERR
      PRINT *,' ------------------------------------------------ '

c -------------------------------------------
c get the fit result""")
b=once(b,'      MNSTAT_COV = ISTAT  ! 3=OK\n','''      MNSTAT_COV = ISTAT  ! 3=OK
c Read the result already obtained; do not invoke any extra MINUIT call.
      IF (PROSP_STATUS_TRACE) THEN
        write(6,'(A,1X,I8,1X,I3,1X,A,1X,L1,1X,4(I6,1X),
     &    3(ES25.17E3,1X))') 'PROSP_MINUIT_STATUS',
     &    CID,INT(INIVAL(1)),TRIM(PARNAME(1)),USE_MINOS,
     &    PROSP_MIN_RETURN,MNSTAT_COV,NPARI,NPARX,
     &    FITCHI2,FEDM,ERRDEF
      ENDIF
''')
(D/'source/snana.car').write_text(b)
c=(S/'snlc_fit.car').read_text();d=c
start=c.index('      SUBROUTINE FCNSNLC(');end=c.index('      END  ! FCNSNLC',start)+len('      END  ! FCNSNLC');part=c[start:end]
part=once(part,'      DOUBLE PRECISION CSP_DIAGW(MXFIT_DATA)\n','''      DOUBLE PRECISION CSP_DIAGW(MXFIT_DATA)
c Private output counters only; no optimizer/model/RNG state writes.
      LOGICAL PROSP_TRACE_ACTIVE,PROSP_TRACE_EMIT
      EXTERNAL PROSP_TRACE_ACTIVE
      INTEGER PROSP_TRACE_ID,PROSP_TRACE_ROWS
      SAVE PROSP_TRACE_ID
      DATA PROSP_TRACE_ID/0/
''')
part=once(part,'     & PEAKMJD,DBLE(MJDOFF),ZSN)\n','''     & PEAKMJD,DBLE(MJDOFF),ZSN)

c Capture existing final-iteration MINUIT physical evaluations only.
      PROSP_TRACE_EMIT = .FALSE.
      PROSP_TRACE_ROWS = 0
      IF ((SNLC_CID.EQ.22 .OR. SNLC_CID.EQ.64) .AND.
     &    ITER.EQ.NFIT_ITERATION .AND.
     &    IFLAG.GE.1 .AND. IFLAG.LE.4) THEN
        IF (PROSP_TRACE_ACTIVE()) THEN
          PROSP_TRACE_ID = PROSP_TRACE_ID + 1
          IF (PROSP_TRACE_ID.LE.4096) THEN
            PROSP_TRACE_EMIT = .TRUE.
            write(6,'(A,1X,A,1X,3(I8,1X),L1,1X,
     &        2(I6,1X),L1,1X,9(ES25.17E3,1X))')
     &        'PROSP_TRIAL_BEGIN',TRIM(SNLC_CCID),ITER,
     &        PROSP_TRACE_ID,IFLAG,LREPEAT_ITER,
     &        NFITDATA,NFITDATA_LOC,USE_FITCOV,
     &        DISTPAR,PEAKMJD+MJDOFF,
     &        INIVAL(IPAR_PEAKMJD)+MJDOFF,ZSN,
     &        SHAPEPAR(1),COLORPAR,RVHOST,MWEBV,DBLE(MJDOFF)
          ELSE IF (PROSP_TRACE_ID.EQ.4097) THEN
            write(6,'(A,1X,A,1X,I8)') 'PROSP_TRIAL_OVERFLOW',
     &        TRIM(SNLC_CCID),PROSP_TRACE_ID
          ENDIF
        ENDIF
      ENDIF
''')
part=once(part,'c convert MAG_ERR into flux-error\n','''c Observe the already returned native mean; never call USRFUN again.
        IF (PROSP_TRACE_EMIT) THEN
          PROSP_TRACE_ROWS = PROSP_TRACE_ROWS + 1
          write(6,'(A,1X,A,1X,4(I8,1X),A1,1X,
     &      8(ES25.17E3,1X))') 'PROSP_TRIAL_ROW',
     &      TRIM(SNLC_CCID),ITER,PROSP_TRACE_ID,ifitdata,
     &      epoch,cfilt,MJD,Trest,flux_model,MAG_ERR,
     &      flux_data,flux_data_errtot,AVwarp,XTMW
        ENDIF

c convert MAG_ERR into flux-error
''')
part=once(part,'c CSP output-only full state; COVMAT2 is inverse when USE_FITCOV.\n','''c Complete existing objective evaluation; no extra likelihood calls.
      IF (PROSP_TRACE_EMIT) THEN
        write(6,'(A,1X,A,1X,3(I8,1X),3(ES25.17E3,1X))')
     &    'PROSP_TRIAL_END',TRIM(SNLC_CCID),ITER,
     &    PROSP_TRACE_ID,PROSP_TRACE_ROWS,
     &    CHI2,CHI2INI,chi2sum_sigma
      ENDIF

c CSP output-only full state; COVMAT2 is inverse when USE_FITCOV.
''')
d=c[:start]+part+c[end:]
d+='''

C ======================================
+DECK,PROSP_TRACE_ACTIVE.
      LOGICAL FUNCTION PROSP_TRACE_ACTIVE()
c Read-only output selector; default off, no caching/shared model state.
      IMPLICIT NONE
      CHARACTER MODE*16
      INTEGER NCHAR,STATUS
      CALL GET_ENVIRONMENT_VARIABLE('PROSP_MINUIT_TRACE',MODE,
     & LENGTH=NCHAR,STATUS=STATUS)
      PROSP_TRACE_ACTIVE = .FALSE.
      IF (STATUS.EQ.1) RETURN
      IF (STATUS.EQ.0 .AND. NCHAR.EQ.1) THEN
        IF (MODE(1:1).EQ.'0') RETURN
        IF (MODE(1:1).EQ.'1') THEN
          PROSP_TRACE_ACTIVE = .TRUE.
          RETURN
        ENDIF
      ENDIF
      write(6,'(A)') 'PROSP_TRACE_ABORT invalid_mode'
      STOP 97
      END
'''
(D/'source/snlc_fit.car').write_text(d)
patch=''
for name,new in [('snana.car',b),('snlc_fit.car',d)]:
 old=(S/name).read_text();patch+=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='base/'+name,tofile='source/'+name))
(D/'output-only.patch').write_text(patch)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(D/'source-hashes.json').write_text(json.dumps({'base':{str(S/n):sha(S/n) for n in ['snana.car','snlc_fit.car']},'new':{n:sha(D/'source'/n) for n in ['snana.car','snlc_fit.car']},'patch':sha(D/'output-only.patch')},indent=2)+'\n')
