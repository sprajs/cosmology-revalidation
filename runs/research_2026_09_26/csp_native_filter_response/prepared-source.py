"""Execute source-gated native CSP baseline first; no relabel until reviewed."""
from pathlib import Path
import csv,gzip,hashlib,json,os,re,shutil,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/research_2026_09_26/csp_native_filter_response'
DESIGN=ROOT/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design'
REL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
PHOT=REL/'photometry/RAISIN/CSPDR3_RAISIN'
BUILD=ROOT/'phase2/official/build/SNANA-v11_04k'
BIN=BUILD/'bin/snlc_fit.exe'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def parse(text,tag):
    cols=None;rows=[]
    for line in text.splitlines():
        a=line.split()
        if not a:continue
        if a[0]=='VARNAMES:':cols=a[1:]
        elif a[0]==tag:
            assert cols and len(a[1:])==len(cols)
            rows.append(dict(zip(cols,a[1:])))
    return rows
def prepare():
    OUT.mkdir(exist_ok=False)
    p=json.loads((DESIGN/'protocol.json').read_text())
    for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h,n
    data=OUT/'data/CSPDR3_RAISIN';data.mkdir(parents=True)
    names={}
    for n in (PHOT/'CSPDR3_RAISIN.LIST').read_text().splitlines():
        text=(PHOT/n).read_text();cid=re.findall(r'(?m)^SNID:\s*(\S+)',text)[0];names[cid]=n
    files=[Path(__file__),DESIGN/'protocol.json',DESIGN/'CSP-nominal.FITRES.gz',BIN]
    for cid in p['pilot_CIDs']:
        q=data/names[cid];shutil.copyfile(PHOT/names[cid],q);files.append(q)
    q=data/'CSPDR3_RAISIN.LIST';q.write_text(''.join(names[c]+'\n' for c in p['pilot_CIDs']));files.append(q)
    q=data/'CSPDR3_RAISIN.README';shutil.copyfile(PHOT/q.name,q);files.append(q)
    nml=REL/'lcfitting/REFAC_CSP_RAISIN_nir_sys.nml';files.append(nml)
    raw='  &SNLCINP'+nml.read_text().split('&SNLCINP',1)[1]
    substitutions={'PRIVATE_DATA_PATH':"'../../data'",'KCOR_FILE':f"'{REL/'kcor/kcor_CSPDR3_BD17.fits'}'",'SNTABLE_LIST':"'FITRES(text:key) LCPLOT(text:key)'",'TEXTFILE_PREFIX':"'fit'",'HEADER_OVERRIDE_FILE':"'vpec.list'",'FITMODEL_NAME':f"'{REL/'model/snoopy.B18'}'"}
    for key,value in substitutions.items():
        raw,n=re.subn(rf'(?m)^(\s*{key}\s*=).*$',lambda m:m[1]+' '+value,raw);assert n==1,key
    for arm in ['nominal','nominal_copy']:
        d=OUT/'fits'/arm;d.mkdir(parents=True)
        f=d/'fit.nml';f.write_text(raw);files.append(f)
        f=d/'vpec.list';shutil.copyfile(REL/'vpec/vpec_baseline_raisin.list',f);files.append(f)
    files.extend(f for f in (REL/'model/snoopy.B18').rglob('*') if f.is_file())
    files.append(REL/'kcor/kcor_CSPDR3_BD17.fits')
    dump(OUT/'execution-protocol.json',dict(status='Initial pilot nominal+identical-copy only, frozen before any CSP native fit. No changed-filter response authorized by this runner until gates and separate instrumentation amendment.',science_design_sha256=sha(DESIGN/'protocol.json'),pilot_CIDs=p['pilot_CIDs'],active_seconds_cap=600,inputs_sha256={str(f.relative_to(ROOT)):sha(f) for f in files},baseline_gates=p['baseline_gates'],optimizer_gate=p['optimizer_gate']))
    (OUT/'prepared-source.py').write_bytes(Path(__file__).read_bytes())
    print('frozen',sha(OUT/'execution-protocol.json'))
def run():
    p=json.loads((OUT/'execution-protocol.json').read_text())
    for n,h in p['inputs_sha256'].items():assert sha(ROOT/n)==h,n
    env=os.environ.copy();env.update(SNANA_DIR=str(BUILD),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    used=0
    for arm in ['nominal','nominal_copy']:
        d=OUT/'fits'/arm;assert not(d/'fit.log').exists()
        start=time.monotonic()
        with (d/'fit.log').open('x') as handle:
            process=subprocess.run([str(BIN),'fit.nml'],cwd=d,env=env,stdout=handle,stderr=subprocess.STDOUT,timeout=max(1,600-used))
        elapsed=time.monotonic()-start;used+=elapsed
        dump(d/'execution.json',dict(returncode=process.returncode,seconds=elapsed,cumulative_active_seconds=used,protocol_sha256=sha(OUT/'execution-protocol.json'),log_sha256=sha(d/'fit.log')))
        print(arm,process.returncode,elapsed,flush=True)
        if process.returncode:raise RuntimeError('Preserved native execution failure; no response run')
def review():
    p=json.loads((OUT/'execution-protocol.json').read_text())
    archived={r['CID']:r for r in parse(gzip.open(DESIGN/'CSP-nominal.FITRES.gz','rt').read(),'SN:')}
    results={};texts={}
    for arm in ['nominal','nominal_copy']:
        f=OUT/'fits'/arm/'fit.FITRES.TEXT';assert f.exists(),f
        text=f.read_text();texts[arm]='\n'.join(l for l in text.splitlines() if l.startswith('SN:'))
        results[arm]={r['CID']:r for r in parse(text,'SN:')}
    assert set(results['nominal'])==set(p['pilot_CIDs'])
    copies_equal=texts['nominal']==texts['nominal_copy']
    lc1=OUT/'fits/nominal/fit.LCPLOT.TEXT';lc2=OUT/'fits/nominal_copy/fit.LCPLOT.TEXT'
    lc_equal=lc1.read_bytes()==lc2.read_bytes()
    rows=[]
    for cid in p['pilot_CIDs']:
        r=results['nominal'][cid];a=archived[cid];d=float(r['DLMAG'])-float(a['DLMAG']);q=float(r['FITCHI2'])-float(a['FITCHI2']);nd=r['NDOF']==a['NDOF']
        fixed=all(float(r[k])==0 for k in ['STRETCHERR','AVERR','PKMJDERR']) and float(r['STRETCH'])==1 and float(r['AV'])==0
        gates=dict(D=abs(d)<=.001,Q=abs(q)<=max(.01,.001*abs(float(a['FITCHI2']))),NDOF=nd,ERRFLAG=int(r['ERRFLAG_FIT'])==0,fixed=fixed)
        rows.append(dict(CID=cid,delta_DLMAG=d,delta_FITCHI2=q,native=r,archived=a,gates=gates,pass_initial_archive_gate=all(gates.values())))
    result=dict(pilot=rows,identical_science_rows=copies_equal,identical_LCPLOT=lc_equal,initial_baseline_pass=copies_equal and lc_equal and all(r['pass_initial_archive_gate'] for r in rows),still_required='Actual fullmask/support/iteration-state/amplitude checks and real post-initialization starts. No changed-filter response executed.')
    dump(OUT/'initial-baseline-result.json',result)
    dump(OUT/'initial-manifest.json',{str(f.relative_to(OUT)):sha(f) for f in OUT.rglob('*') if f.is_file() and f.name!='initial-manifest.json'})
    print(json.dumps({k:v for k,v in result.items() if k!='pilot'},indent=2));print([(r['CID'],r['delta_DLMAG'],r['delta_FITCHI2'],r['gates']) for r in rows])
if __name__=='__main__':{'prepare':prepare,'run':run,'review':review}[sys.argv[1]]()
