"""Isolated output-only v11_04k instrumentation; no native execution."""
from pathlib import Path
import difflib,hashlib,json,shutil

ROOT=Path(__file__).resolve().parents[5]
OUT=Path(__file__).resolve().parent/'instrumentation'
BUILD=OUT/'SNANA-v11_04k-output'
SOURCE=ROOT/'phase2/official/build/SNANA-v11_04k'
PIN=ROOT/'sources/repos/RickKessler__SNANA@v11_04k/src/snlc_fit.car'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not OUT.exists()
OUT.mkdir()
assert (SOURCE/'src/snlc_fit.car').read_bytes()==PIN.read_bytes()
shutil.copytree(SOURCE,BUILD)
p=BUILD/'src/snlc_fit.car';original=p.read_text();s=original

def insert_before(text,anchor,addition):
    assert text.count(anchor)==1,(anchor,text.count(anchor))
    return text.replace(anchor,addition+'\n'+anchor)

# FITPAR_PREP has completed all initialization before returning to MNFIT.
anchor='      END   ! end of FITPAR_PREP'
idx=s.index(anchor);before=s[:idx];after=s[idx:]
point=before.rfind('      RETURN\n')
assert point>0
entry="""c CSP output-only audit: true parameters passed to MNFIT_DRIVER.
      write(6,'(A,1X,A,1X,I3,1X,I6,1X,2L2,1X,
     &  12(ES25.17E3,1X))') 'CSP_ENTRY:',TRIM(SNLC_CCID),
     &  ITER,NFITDATA,USE_FITCOV,LREPEAT_ITER,
     &  INIVAL(IPAR_DLMAG),INIVAL(IPAR_SHAPE),INIVAL(IPAR_AV),
     &  INIVAL(IPAR_PEAKMJD)+MJDOFF,
     &  INISTP(IPAR_DLMAG),INISTP(IPAR_SHAPE),INISTP(IPAR_AV),
     &  INISTP(IPAR_PEAKMJD),
     &  FITVAL(IPAR_DLMAG,ITER-1),FITVAL(IPAR_SHAPE,ITER-1),
     &  FITVAL(IPAR_AV,ITER-1),FITVAL(IPAR_PEAKMJD,ITER-1)+MJDOFF
      CALL FLUSH(6)

"""
s=before[:point]+entry+before[point:]+after
start=s.index('      SUBROUTINE FCNSNLC(');end=s.index('      END  ! FCNSNLC',start)
fc=s[start:end]
fc=insert_before(fc,'      character cfilt*1',
'''c CSP output-only audit: local workspace, not fitted state.
      DOUBLE PRECISION CSP_DIAGW(MXFIT_DATA)''')
anchor='      CHI2 = 1.0E7'
fc=insert_before(fc,anchor,
'''      DO irow=1,MXFIT_DATA
        CSP_DIAGW(irow)=0.0D0
      ENDDO''')
anchor='        R4EP_ALL(epoch,JEP_FLUX_ERRTOT)   = flux_errtot '
addition='''
c Print only final callback of each iteration. No model re-evaluation.
        IF (LFLAG_LAST_MN .and. LFITDATA) THEN
          CSP_DIAGW(ifitdata)=inv_sqsig
          write(6,'(A,1X,A,1X,I3,1X,I6,1X,I6,1X,A,1X,
     &      13(ES25.17E3,1X))') 'CSP_ROW:',TRIM(SNLC_CCID),
     &      ITER,ifitdata,epoch,cfilt,MJD,Trest,flux_model,
     &      mag_err,flux_data,flux_data_errtot,ZSN,MWEBV,
     &      LAMAVG/ZZ,flux_fudge_err,flux_model_err,
     &      DBLE(RESTLAMBDA_USEFIT(1)),DBLE(RESTLAMBDA_USEFIT(2))
        ENDIF
'''
assert fc.count(anchor)==1
fc=fc.replace(anchor,anchor+'\n'+addition)
point=fc.rfind('      RETURN\n')
assert point>0
ending='''c CSP output-only full state; COVMAT2 is inverse when USE_FITCOV.
      IF (LFLAG_LAST_MN) THEN
        write(6,'(A,1X,A,1X,I3,1X,I6,1X,L2,1X,
     &    11(ES25.17E3,1X))') 'CSP_OBJECTIVE:',
     &    TRIM(SNLC_CCID),ITER,NFITDATA,USE_FITCOV,
     &    CHI2,CHI2INI,chi2sum_sigma,
     &    XVAL(IPAR_DLMAG),XVAL(IPAR_SHAPE),XVAL(IPAR_AV),
     &    XVAL(IPAR_PEAKMJD)+MJDOFF,
     &    INIVAL(IPAR_PEAKMJD)+MJDOFF,
     &    DBLE(SNLC_SEARCH_PEAKMJD),DBLE(MJDOFF),
     &    DBLE(OPT_CHI2_SIGMA)
        DO irow=1,NFITDATA
          IF (USE_FITCOV) THEN
            write(6,'(A,1X,A,1X,I3,1X,I6,1X,
     &        *(ES25.17E3,1X))') 'CSP_WROW:',
     &        TRIM(SNLC_CCID),ITER,irow,
     &        (COVMAT2(irow,icol),icol=1,NFITDATA)
          ELSE
            write(6,'(A,1X,A,1X,I3,1X,I6,1X,ES25.17E3)')
     &        'CSP_WDIAG:',TRIM(SNLC_CCID),ITER,irow,
     &        CSP_DIAGW(irow)
          ENDIF
        ENDDO
        CALL FLUSH(6)
      ENDIF

'''
fc=fc[:point]+ending+fc[point:]
s=s[:start]+fc+s[end:]
p.write_text(s)
patch=''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='original/src/snlc_fit.car',tofile='instrumented/src/snlc_fit.car'))
(OUT/'output-only.patch').write_text(patch)
added=[x[1:] for x in patch.splitlines() if x.startswith('+') and not x.startswith('+++')]
assert all(len(x)<=72 or x.startswith('c') for x in added),[(len(x),x) for x in added if len(x)>72]
# Reuse every unmodified build object/library, rebuilding only the affected
# Fortran target through the existing pinned recipe. Shared build untouched.
script='''#!/usr/bin/env bash
set -euo pipefail
task_root="/home/szymon/Documents/ChatGPT/supernova"
task_build="'''+str(BUILD)+'''"
export SNANA_DIR="$task_build"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SNANA_DIR/src"
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make -j2 ../bin/snlc_fit.exe FFC="$task_compiler" \\
  LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" \\
  EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
'''
(OUT/'build.sh').write_text(script)
protocol={
 'state':'Frozen before compiling or native instrumented outputs. Output-only port, no response fits.',
 'original_source_sha256':sha(PIN),'original_binary_sha256':sha(SOURCE/'bin/snlc_fit.exe'),
 'instrumented_source_sha256':sha(p),'patch_sha256':sha(OUT/'output-only.patch'),
 'prepare_source_sha256':sha(Path(__file__)),'build_script_sha256':sha(OUT/'build.sh'),
 'build':'Private52MB copy of alreadycompiled v11_04k; same compiler/libraries/flags. 180s buildcap. Shared original untouched.',
 'changes':'Only terminalwrites and newlocal diagonalW array. No USRFUN calls added, objective/C/acceptedmask modifications, startup overrides or newparameters.',
 'schemas':{
   'CSP_ENTRY':'CID ITER NFITDATA USE_FITCOV LREPEAT_ITER D shape AV absolutePeak; fourINISTP; fourpreviousFITVAL',
   'CSP_ROW':'CID ITER acceptedIndex sourceEpoch band; MJD restphase modelF modelMagErr dataF dataErr z MWEBV restMeanLambda fudgeFluxErr modelFluxErr restLambdaMin restLambdaMax',
   'CSP_OBJECTIVE':'CID ITER NFITDATA USE_FITCOV; totalQ priorQ sigmaQ D shape AV absolutePeak priorCenter searchPeak MJDOFF OPT_CHI2_SIGMA',
   'CSP_WROW':'CID ITER row followed by fullinverseC row (native USE_FITCOV)',
   'CSP_WDIAG':'CID ITER row inversevariance (native diagonal objective)'},
 'verification':'Root runs instrumented nominal identicalpilot inputs; require originalFITRES science rows/LCPLOT equality before arithmetic. Export every callbackblock in order, never overwrite earlier sameCID/ITER states; require exactrows for each block. Qclosure residualTWres+priorQ+sigmaQ=totalQ. Positiveamplitude kernel onlyif sigmaQ absent/constant.',
 'means':'After equivalence, standaloneall-fixed NML Dref and Dref±.15, with nativefloat32 seedaccounting; samephysicalrows, modelmean scaling2e−8 normalized tolerance. C differences diagnostic, originalstateC frozen forprofile. No modernoracle substitution.',
 'starts':'CSP_ENTRY records truepostinitialization vector. Distinct INIVAL_DLMAG NML requests may be erased; this output-only amendment doesnot force starts. If erased, report untested; any postinitoverride needs separateamendment.',
 'iterations':'Retain eachfinalcallback D/C/mask and bothiteration2/3. Need finaltwoDgap≤.001mag; Cmetricreported, not claimed state-independent objective. Repeatediterations separatelyidentified by blocksequence.',
 'support':'CSP_ROW restMeanLambda is nativefilter-mean support screen, not broadbandendpoint proof. KCORFilterTrans restwavelength extent and empiricalshape/phase grid checked offline; firstiterationrowsincluded andclipping/nonpositive mean flagged.'}
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps({'protocol_sha256':sha(OUT/'protocol.json'),'build_script':str(OUT/'build.sh'),'build':str(BUILD)}))
