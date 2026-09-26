"""Create an isolated opt-in native mean oracle; do not modify upstream build."""
from pathlib import Path
import shutil,subprocess,json,hashlib
O=Path(__file__).resolve().parent;ROOT=O.parents[3]
SRC=ROOT/'phase2/official/build/SNANA-audit-v3';OUT=O/'fixed-c-profile';BUILD=OUT/'native-build'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    if BUILD.exists():assert sha(BUILD/'src/snlc_fit.F90')==sha(SRC/'src/snlc_fit.F90')
    else:shutil.copytree(SRC,BUILD,symlinks=True)
    f=BUILD/'src/snlc_fit.F90';original=f.read_text()
    decl='''    ! RAISIN optional mean-only oracle; original objective/mask/C are untouched.
    logical, save :: raisin_stream_done = .false.
    integer raisin_status, raisin_id, raisin_row, raisin_epoch, raisin_filter
    character(len=16) raisin_environment
    double precision raisin_pars(4), raisin_shape(2), raisin_means(MXFIT_DATA)
    double precision raisin_tobs, raisin_wavecor, raisin_avwarp
    double precision raisin_magkcor(2), raisin_xtav, raisin_xtmw, raisin_magerr
'''
    marker='    character cfilt*1\n'
    start=original.index('    SUBROUTINE FCNSNLC(');end=original.index('  END SUBROUTINE FCNSNLC',start)
    section=original[start:end];assert section.count(marker)==1
    modified=original[:start]+section.replace(marker,marker+decl)+original[end:]
    body='''    IF (LFLAG_LAST_MN .and. ITER==NFIT_ITERATION .and. .not.raisin_stream_done) THEN
      raisin_stream_done=.true.
      call get_environment_variable('RAISIN_PROBE_STREAM',raisin_environment,status=raisin_status)
      IF (raisin_status==0 .and. trim(raisin_environment)=='1') THEN
        write(6,'(A,1X,I6)') 'RAISIN_READY:',NFITDATA
        call flush(6)
        DO
          read(5,*,iostat=raisin_status) raisin_id,raisin_pars
          IF (raisin_status/=0) EXIT
          IF (raisin_id<0) EXIT
          raisin_shape(1)=raisin_pars(2)
          raisin_shape(2)=XVAL(IPAR_SHAPE2)
          DO raisin_row=1,NFITDATA
            raisin_epoch=EPLIST_FIT(raisin_row)
            raisin_filter=I4EP_ALL(raisin_epoch,IEP_IFILT_OBS)
            ifitdata_usrfun=raisin_row
            raisin_tobs=R8EP_MJD(raisin_row)-raisin_pars(4)
            raisin_wavecor=dble(R4EP_ALL(raisin_epoch,JEP_WAVECOR))
            raisin_means(raisin_row)=USRFUN(ITER,raisin_filter,ZSN,raisin_tobs, &
              raisin_wavecor,raisin_shape,raisin_pars(1),raisin_pars(3),RVHOST,MWEBV, &
              .false.,raisin_avwarp,raisin_magkcor,raisin_xtav,raisin_xtmw,raisin_magerr)
          ENDDO
          write(6,'(A,1X,I10,1X,I6,1X,*(ES25.17E3,1X))') 'RAISIN_MEAN:', &
            raisin_id,NFITDATA,raisin_pars,(raisin_means(raisin_row),raisin_row=1,NFITDATA)
          call flush(6)
        ENDDO
        write(6,'(A)') 'RAISIN_STREAM_END'
        call flush(6)
      ENDIF
    ENDIF
'''
    marker='    RETURN\n  END SUBROUTINE FCNSNLC'
    assert modified.count(marker)==1;modified=modified.replace(marker,body+marker)
    f.write_text(modified)
    (OUT/'snlc_fit.original.F90').write_text(original)
    patch=subprocess.run(['diff','-u',str(OUT/'snlc_fit.original.F90'),str(f)],text=True,capture_output=True)
    (OUT/'mean-stream.patch').write_text(patch.stdout)
    with(OUT/'build.log').open('w') as log:
      r=subprocess.run(['make','-C',str(BUILD/'src'),'snlc_fit'],stdout=log,stderr=subprocess.STDOUT)
    assert r.returncode==0,'Read build.log; originals unchanged'
    result=dict(status='Built isolated optional mean-only oracle; numerical gates pending',
      original_source_sha256=sha(SRC/'src/snlc_fit.F90'),modified_source_sha256=sha(f),
      source_builder_sha256=sha(__file__),binary_sha256=sha(BUILD/'bin/snlc_fit.exe'),patch_sha256=sha(OUT/'mean-stream.patch'),
      build_log_sha256=sha(OUT/'build.log'),original_binary_sha256=sha(SRC/'bin/snlc_fit.exe'))
    (OUT/'build-manifest.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
