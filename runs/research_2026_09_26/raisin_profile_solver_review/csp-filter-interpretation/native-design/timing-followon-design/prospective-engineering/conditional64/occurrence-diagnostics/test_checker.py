from pathlib import Path
import importlib.util,json,hashlib
import checker as c
D=Path(__file__).resolve().parent;O=D.parent
fixture_results=[]
def blocks(seq,flags=None):return [{'ITER':i,'entry':{'LREPEAT_ITER':flags[j] if flags is not None else False}} for j,i in enumerate(seq)]
for name,seq,flags,should_pass in [
 ('normal',list(range(1,13)),None,True),
 ('one_native_final_retry',list(range(1,13))+[12],[False]*12+[True],True),
 ('missing_iteration',[1,2,4]+list(range(5,13)),None,False),
 ('arbitrary_repeat',[1,2,2]+list(range(3,13)),None,False),
 ('unflagged_final_retry',list(range(1,13))+[12],None,False),
 ('double_final_retry',list(range(1,13))+[12,12],[False]*12+[True,True],False),
 ('spurious_repeat_flag',list(range(1,13)),[False]*11+[True],False)]:
 try:c.occurrence_pattern(blocks(seq,flags),12);ok=True
 except RuntimeError:ok=False
 assert ok==should_pass,name;fixture_results.append({'case':name,'accepted':ok,'expected':should_pass})
by,result=c.review_saved(O/'fits/joint12',[str(x) for x in range(1,65)])
c.save(D/'saved-joint12-check.json',result)
c.save(D/'test-result.json',{'fixtures':fixture_results,'saved_original_joint12_checks_pass':True,'callbacks':sum(len(x) for x in by.values()),'MINIMIZE4_callback_count':len(result['native']['abnormal_minimize_returns']),'no_native_calls':True})
print(json.dumps({'fixtures_pass':len(fixture_results),'saved_callbacks':sum(len(x) for x in by.values()),'warnings':len(result['native']['abnormal_minimize_returns'])}))
