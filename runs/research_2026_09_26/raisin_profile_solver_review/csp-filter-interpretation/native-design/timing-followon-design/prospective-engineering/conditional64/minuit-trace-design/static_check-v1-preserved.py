from pathlib import Path
import re,json,hashlib,difflib
D=Path(__file__).resolve().parent;S=D.parent.parent/'restricted-peak-engineering/start-only-hook/build/src'
result={}
for name in ['snana.car','snlc_fit.car']:
 old=(S/name).read_text();new=(D/'source'/name).read_text()
 def calls(s):return re.findall(r'^      \s*CALL\s+(\w+)',s,re.M|re.I)
 a,b=calls(old),calls(new);added=b.copy()
 for x in a:added.remove(x)
 assert [x for x in b if x!='GET_ENVIRONMENT_VARIABLE']==a
 result[name]={'base_sha256':hashlib.sha256(old.encode()).hexdigest(),'new_sha256':hashlib.sha256(new.encode()).hexdigest(),'added_CALL_statements':added,'all_original_CALL_sequence_exact':True}
 addedlines=[x[1:] for x in difflib.unified_diff(old.splitlines(),new.splitlines()) if x.startswith('+') and not x.startswith('+++')]
 assert not any(re.search(r'\b(USRFUN|FCNSNLC|MNSTAT|MNMIGR|MNPARM|MNEXCM|FITINI_COV|RAN\w*)\s*\(',x,re.I) for x in addedlines if not x.lstrip().startswith(('c','C','+')))
 result[name]['added_executable_max_line_length']=max(len(x) for x in addedlines if x and x[0] not in 'cC+')
 assert result[name]['added_executable_max_line_length']<=72
part=(D/'source/snlc_fit.car').read_text();a=part.index('      SUBROUTINE FCNSNLC(');b=part.index('      END  ! FCNSNLC',a);fc=part[a:b]
assert fc.index('CALL PROSP_PEAK_FCN_GUARD')<fc.index("'PROSP_TRIAL_BEGIN'")<fc.index('USRFUN ( ITER')<fc.index("'PROSP_TRIAL_ROW'")<fc.index("'PROSP_TRIAL_END'")
assert 'IFLAG.GE.1 .AND. IFLAG.LE.4' in fc and 'ITER.EQ.NFIT_ITERATION' in fc and 'PROSP_TRACE_ID.LE.4096' in fc
result['source_only_static_pass']=True;result['native_builds_or_runs']=0
(D/'static-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
