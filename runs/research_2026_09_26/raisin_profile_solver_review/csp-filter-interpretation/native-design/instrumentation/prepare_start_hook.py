from pathlib import Path
import json,hashlib,shutil,difflib
O=Path(__file__).resolve().parent;S=O/'SNANA-v11_04k-output/src/snlc_fit.car';B=O/'SNANA-v11_04k-output/bin/snlc_fit.exe';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
shutil.copyfile(S,O/'output-only-v2-snlc_fit.car');shutil.copyfile(B,O/'output-only-v2-snlc_fit.exe');s=S.read_text()
a='      REAL*4 Z4, LMIN4';assert s.count(a)==2;s=s.replace(a,a+'''\n      REAL*8 CSP_DSHIFT, CSP_DENTRY
      INTEGER CSP_ENV_STATUS, CSP_READ_STATUS
      CHARACTER CSP_ENV*64''',1)
a='c CSP output-only audit: true parameters passed to MNFIT_DRIVER.';assert s.count(a)==1
s=s.replace(a,'''c Explicit diagnostic start shift, only after native initialization.
c Absent environment leaves all fitted state unchanged.
      IF (FIRST_ITER) THEN
        CSP_ENV=' '
        CALL GET_ENVIRONMENT_VARIABLE('CSP_POSTINIT_DSHIFT',
     &    CSP_ENV,STATUS=CSP_ENV_STATUS)
        CSP_DSHIFT=0.0D0
        IF (CSP_ENV_STATUS .EQ. -1) STOP 72
        IF (LEN_TRIM(CSP_ENV) .GT. 0) THEN
          READ(CSP_ENV,*,IOSTAT=CSP_READ_STATUS) CSP_DSHIFT
          IF (CSP_READ_STATUS .NE. 0) STOP 73
          IF (ABS(CSP_DSHIFT) .GT. 1.0D0) STOP 74
          IF (CSP_DSHIFT .NE. CSP_DSHIFT) STOP 75
          CSP_DENTRY=INIVAL(IPAR_DLMAG)
          INIVAL(IPAR_DLMAG)=CSP_DENTRY+CSP_DSHIFT
          write(6,'(A,1X,A,1X,I3,1X,3(ES25.17E3,1X))')
     &     'CSP_START_SHIFT:',TRIM(SNLC_CCID),ITER,
     &     CSP_DENTRY,CSP_DSHIFT,INIVAL(IPAR_DLMAG)
        ENDIF
      ENDIF

'''+a)
S.write_text(s);old=(O/'output-only-v2-snlc_fit.car').read_text();patch=''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='output-v2/src/snlc_fit.car',tofile='start-hook-v3/src/snlc_fit.car'));(O/'start-hook-v3.patch').write_text(patch)
assert all(len(x[1:])<=72 or x.startswith('+c') for x in patch.splitlines() if x.startswith('+') and not x.startswith('+++'))
p={'reason':'Parent authorized true post-initialization ±0.2mag starts; ordinary requested seeds can be overwritten. First iteration only offset after native prep, default absent leaves state untouched. No objective/C/model code changes.','before_source_sha256':sha(O/'output-only-v2-snlc_fit.car'),'before_binary_sha256':sha(O/'output-only-v2-snlc_fit.exe'),'after_source_sha256':sha(S),'patch_sha256':sha(O/'start-hook-v3.patch'),'script_sha256':sha(Path(__file__)),'gate':'Default nominal FITRES science rows/LCPLOT must match original again; active shift must export exactly baseline entry ±.2. Fixed final masks, finalD agreement≤.001, iterations/amplitude gate same.','start_design_amendment':'Explicitly ±.2 around actual post-initialization D (2004ef35.7428562716,2005hc36.5285993367), not archived finalD. Original plan archived-centering is amended before changed-filter responses; root authorized. Actual center versus archived finalD reported; tests initialization sensitivity, not different model.','other_probes':'All-fixed native Dref,Dref±.15 using NML float32 seeds; actual exports define ratio. No changed-filter outputs.','build_cap_seconds':120,'resource_authorization':'Parent allows up to360s total build cap; native600s remains.'}
(O/'start-hook-v3-protocol.json').write_text(json.dumps(p,indent=2)+'\n');print(sha(O/'start-hook-v3-protocol.json'))
