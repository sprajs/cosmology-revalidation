"""Frozen second-tier ±180-day CAOM metadata expansion after zero exact headers."""
import hashlib,json,time
from pathlib import Path
from query_caom import filters,request,count_from

BASE=Path(__file__).resolve().parent
P=json.loads((BASE/'protocol.json').read_text())
H=json.loads((BASE/'header-results.json').read_text())
assert all(z['eligible_count']<4 for z in H['visits'].values())
started=time.monotonic(); transferred=H['total_metadata_bytes']; requests=[]; visits={}
for visit,mid in P['science_visit_midpoints_mjd'].items():
    f=filters(P,mid,180,True)
    label=f'{visit}_180d_dark_target_expanded'
    payload={'service':'Mast.Caom.Filtered','format':'json','params':{'columns':'COUNT_BIG(*)','filters':f,'obstype':'all'}}
    answer,rec=request(payload,label+'_count',P['limits'],started,transferred)
    requests.append(rec);transferred+=rec['bytes']; n=count_from(answer)
    rows=[]; page=1
    while len(rows)<n:
        payload={'service':'Mast.Caom.Filtered','format':'json','params':{'columns':'*','filters':f,'obstype':'all'},'pagesize':500,'page':page}
        answer,rec=request(payload,label+f'_page{page}',P['limits'],started,transferred)
        requests.append(rec);transferred+=rec['bytes'];got=answer.get('data',[])
        if not got:break
        rows.extend(got);page+=1
    visits[visit]={'count':n,'rows':rows}
out={'status':'metadata_only','visits':visits,'requests':requests,'elapsed_seconds':time.monotonic()-started,'cumulative_metadata_bytes':transferred,'header_result_sha256':hashlib.sha256((BASE/'header-results.json').read_bytes()).hexdigest()}
(BASE/'expanded-caom-results.json').write_text(json.dumps(out,indent=2)+'\n')
