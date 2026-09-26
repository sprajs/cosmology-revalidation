"""Prepare, but do not execute, frozen CSP native filter interventions."""
from pathlib import Path
from decimal import Decimal as D
from collections import Counter
import csv,json,hashlib,re,shutil
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_native_filter_response'
OUT=BASE/'cohort-preparation'
DESIGN=ROOT/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design'
REL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=False)
    p=json.loads((DESIGN/'protocol.json').read_text())
    for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h,n
    assert json.loads((BASE/'initial-baseline-result.json').read_text())['initial_baseline_pass']
    cohort=list(csv.DictReader((DESIGN/'cohort42.csv').open()))
    changes=list(csv.DictReader((DESIGN/'known-relabel-rows.csv').open()))
    expected=Counter()
    for r in changes:expected[(r['CID'],D(r['absolute_MJD']),r['old_band'])]+=int(r['multiplicity'])
    assert sum(expected.values())==188
    files=[Path(__file__),DESIGN/'protocol.json',DESIGN/'cohort42.csv',DESIGN/'known-relabel-rows.csv',BASE/'initial-baseline-result.json',BASE/'execution-protocol.json']
    jobs=[];ledger=[]
    for scope,ids in [('pilot',p['pilot_CIDs']),('full',[r['CID'] for r in cohort])]:
        for arm in ['nominal','nominal_copy','known_Jdw_to_j']:
            data=OUT/scope/'data'/arm/'CSPDR3_RAISIN';data.mkdir(parents=True)
            names=[];changed=0
            for c in cohort:
                if c['CID'] not in ids:continue
                source=ROOT/c['photometry_path'];original=source.read_text();cid=c['CID'];seen=Counter();lines=[]
                for i,l in enumerate(original.splitlines(keepends=True)):
                    tokens=l.split()
                    if tokens and tokens[0]=='OBS:':
                        key=(cid,D(tokens[1]),tokens[2])
                        if key in expected:
                            seen[key]+=1
                            if arm=='known_Jdw_to_j':
                                new,n=re.subn(r'^(OBS:\s+\S+\s+)J(?=\s)',r'\1j',l);assert n==1
                                assert new.split()[:2]+new.split()[3:]==l.split()[:2]+l.split()[3:]
                                l=new;changed+=1;ledger.append(dict(scope=scope,CID=cid,source_line=i+1,MJD=tokens[1],old_band='J',new_band='j'))
                    lines.append(l)
                assert seen==Counter({k:v for k,v in expected.items() if k[0]==cid})
                dest=data/source.name;dest.write_text(''.join(lines));names.append(source.name);files.append(dest)
                if arm!='known_Jdw_to_j' or not seen:assert dest.read_bytes()==source.read_bytes()
            q=data/'CSPDR3_RAISIN.LIST';q.write_text(''.join(n+'\n' for n in names));files.append(q)
            q=data/'CSPDR3_RAISIN.README';shutil.copyfile(REL/'photometry/RAISIN/CSPDR3_RAISIN'/q.name,q);files.append(q)
            fit=OUT/scope/'fits'/arm;fit.mkdir(parents=True)
            text=(BASE/'fits/nominal/fit.nml').read_text();text,n=re.subn(r"(?m)^(\s*PRIVATE_DATA_PATH\s*=).*$",lambda m:m[1]+f" '../../data/{arm}'",text);assert n==1
            nml=fit/'fit.nml';nml.write_text(text);files.append(nml)
            q=fit/'vpec.list';shutil.copyfile(BASE/'fits/nominal/vpec.list',q);files.append(q)
            jobs.append(dict(scope=scope,arm=arm,objects=len(names),changed_rows=changed,nml=str(nml.relative_to(ROOT))))
    with (OUT/'changed-row-ledger.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=ledger[0]);w.writeheader();w.writerows(ledger)
    files.append(OUT/'changed-row-ledger.csv')
    result=dict(status='Inputs frozen before any changed-filter outcome. Prepared inputs are not executed or scientifically certified.',parent_design_sha256=sha(DESIGN/'protocol.json'),parent_pilot_pass=True,jobs=jobs,change='Only exact physical-unique Jdw rows change J→j. Other bytes preserved. Four ambiguous2006kf rows retainedJ. Optionalstress not prepared yet.',gates='Requires source-instrumentation/optimizer/support gates and fullcohort baseline before final42summary.600sactive native budget including1.331s parentinitialpilot.',inputs_sha256={str(f.relative_to(ROOT)):sha(f) for f in files})
    (OUT/'protocol.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes());print(json.dumps(jobs,indent=2));print('protocol',sha(OUT/'protocol.json'))
if __name__=='__main__':main()
