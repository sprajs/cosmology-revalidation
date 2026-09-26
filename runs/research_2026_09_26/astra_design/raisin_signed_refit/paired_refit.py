"""Frozen source-derived SNooPy epoch-restoration experiment. No source edits."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,re,shutil,subprocess,time
import numpy as np

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
REL=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
AUD=ROOT/'runs/research_2026_09_26/raisin_flux_sign_audit'
BUILD=ROOT/'phase2/official/build/SNANA-current';EXE=BUILD/'bin/snlc_fit.exe'
DOC=ROOT/'docs/research-2026-09-26/raisin-signed-refit-design.md'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')
def csvsave(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def phot(p):
    rows=[];cols=None;head={}
    for li,line in enumerate(p.read_text().splitlines(),1):
        t=line.split()
        if not t:continue
        if t[0]=='VARLIST:':cols=t[1:]
        elif t[0]=='OBS:':
            d=dict(zip(cols,t[1:]));d.update(_line=li,_index=len(rows),_text=line)
            rows.append(d)
        elif t[0].endswith(':'):head[t[0][:-1]]=t[1:]
    return head,rows
def assemble(original,obs):
    lines=[]
    for line in original.splitlines():
        if line.startswith('OBS:') or line.strip()=='END:':continue
        if line.startswith('NOBS:'):line='NOBS: '+str(len(obs))
        lines.append(line)
    return '\n'.join(lines+obs+['END:'])+'\n'
def commonobs(r):return f"OBS: {r['MJD']} {r.get('FLT',r.get('BAND'))} NULL {r['FLUXCAL']} {r['FLUXCALERR']}"

def prepare():
    assert not (OUT/'input-manifest.json').exists()
    hand=json.loads((AUD/'lineage-handoff.json').read_text())
    for p,h in hand['inputs_sha256'].items():
        path=ROOT/'runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json' if p=='code-tree.json' else AUD/p
        assert sha(path)==h,(p,'lineage hash')
    acq=json.loads((AUD/'author-acquisition.json').read_text())
    for r in acq['records']:assert sha(AUD/'author-source'/r['path'])==r['sha256']
    gates=list(csv.DictReader((AUD/'lineage-object-gates.csv').open()))
    cids=sorted(r['CID'] for r in gates if r['year']=='2016');assert len(cids)==10
    aliases=json.loads((AUD/'protocol.json').read_text())['des_aliases']
    inputpaths=[DOC,Path(__file__),AUD/'lineage-handoff.json',AUD/'lineage-object-gates.csv',AUD/'author-negative-candidates.csv',AUD/'author-acquisition.json',
       REL/'lcfitting/REFAC_DES_RAISIN_optical_sys.nml',REL/'kcor/kcor_DES_NIR.fits',REL/'kcor/kcor_DES_NIR.input',REL/'model/snoopy.B18/SNooPy_B18.fits',REL/'model/snoopy.B18/snoopy.info',REL/'vpec/vpec_baseline_raisin.list',
       EXE,BUILD/'src/snlc_fit.F90',BUILD/'src/snana.F90',BUILD/'src/genmag_snoopy.c',BUILD/'src/sntools_modelgrid_read.c',OUT/'author-optical-FITOPT000.FITRES.gz',OUT/'author-fitres-acquisition.json']
    cohort=[];ledger=[]
    for arm in ['R','A','B','H','Rcopy']:
        d=OUT/'data'/f'RSR_{arm}';d.mkdir(parents=True,exist_ok=True)
        (d/f'RSR_{arm}.LIST').write_text(''.join(c+'.snana.dat\n' for c in cids))
        (d/f'RSR_{arm}.README').write_text('DOCUMENTATION:\n  PURPOSE: Frozen conditional RAISIN signed-epoch refit\nDOCUMENTATION_END:\n')
    for cid in cids:
        released=REL/'photometry/RAISIN/DES_RAISIN'/f'{cid}.snana.dat'
        raw=AUD/'author-source/data/Photometry/DES'/f"des_real_{int(aliases[cid]):08d}.dat"
        inputpaths.extend([released,raw]);h,r=phot(released);_,s=phot(raw)
        s=[x for x in s if x.get('FLT',x.get('BAND')) in 'griz']
        for x in s:assert np.isfinite([float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])]).all() and float(x['FLUXCALERR'])>0
        neg=[x for x in s if float(x['FLUXCAL'])<=0];pos=[x for x in s if float(x['FLUXCAL'])>0]
        assert all(not any(y['FLT']==x['FLT'] and abs(float(y['MJD'])-float(x['MJD']))<=.00055 for y in r) for x in neg),'hybrid ambiguous collision'
        nir=[x['_text'] for x in r if x['FLT'] not in 'griz'];pobs=[commonobs(x) for x in pos];nobs=[commonobs(x) for x in neg]
        texts={'R':released.read_text(),'Rcopy':released.read_text(),'A':assemble(released.read_text(),nir+pobs),'B':assemble(released.read_text(),nir+pobs+nobs)}
        # H preserves every existing released line/order except NOBS and END location.
        old=released.read_text();new=[]
        for line in old.splitlines():
            if line.startswith('NOBS:'):line='NOBS: '+str(len(r)+len(neg))
            if line.strip()=='END:':new.extend(nobs)
            new.append(line)
        texts['H']='\n'.join(new)+'\n'
        for arm,txt in texts.items():
            p=OUT/'data'/f'RSR_{arm}'/released.name;p.write_text(txt)
            hh,rr=phot(p)
            assert all(hh.get(k)==h.get(k) for k in ['SNID','MWEBV','REDSHIFT_HELIO','REDSHIFT_FINAL','PEAKMJD','RA','DEC'])
            assert int(hh['NOBS'][0])==len(rr)
            if arm in ['R','Rcopy']:assert p.read_bytes()==released.read_bytes()
            if arm=='H':assert [x['_text'] for x in rr[:len(r)]]==[x['_text'] for x in r]
        # Remove appended nonpositives from B: same native A observation records exactly.
        assert phot(OUT/'data/RSR_B'/released.name)[1][:-len(neg)]==phot(OUT/'data/RSR_A'/released.name)[1] if neg else True
        cohort.append(dict(CID=cid,source_CID=aliases[cid],released_path=str(released.relative_to(ROOT)),raw_path=str(raw.relative_to(ROOT)),released_optical=sum(x['FLT'] in 'griz' for x in r),author_positive=len(pos),author_nonpositive=len(neg),peak_header=float(h['PEAKMJD'][0]),zHEL=float(h['REDSHIFT_HELIO'][0])))
        for x in s:ledger.append(dict(CID=cid,source_path=str(raw.relative_to(ROOT)),source_line=x['_line'],source_index=x['_index'],band=x['FLT'],MJD=x['MJD'],FLUXCAL=x['FLUXCAL'],FLUXCALERR=x['FLUXCALERR'],PHOTFLAG=x.get('PHOTFLAG',''),included_A=float(x['FLUXCAL'])>0,included_B=True,appended_H=float(x['FLUXCAL'])<=0))
    csvsave(OUT/'cohort.csv',cohort);csvsave(OUT/'source-row-ledger.csv',ledger)
    source=(REL/'lcfitting/REFAC_DES_RAISIN_optical_sys.nml').read_text().split('  &SNLCINP',1)[1]
    source='  &SNLCINP'+source
    replacements={'PRIVATE_DATA_PATH':str(OUT/'data'),'KCOR_FILE':str(REL/'kcor/kcor_DES_NIR.fits'),'HEADER_OVERRIDE_FILE':str(OUT/'vpec_baseline_documented.list'),'FITMODEL_NAME':str(REL/'model/snoopy.B18')}
    for key,v in replacements.items():source=re.sub(r'(?m)^\s*'+key+r'\s*=.*$',f" {key} = '{v}'",source)
    source=source.replace("SNTABLE_LIST = 'FITRES LCPLOT(text:key)'","SNTABLE_LIST = 'FITRES(text:key) LCPLOT(text:key)'")
    source=source.replace('  &SNLCINP','  &SNLCINP\n OPT_SETPKMJD = 0\n OPT_MWCOLORLAW = 99\n MXLC_PLOT = 1000')
    for stage in ['engineering','engineering_copy','full']:
        for peak in ['free','fixed']:
            for arm in (['R'] if stage=='engineering' else ['Rcopy'] if stage=='engineering_copy' else ['R','A','B','H']):
                d=OUT/'fits'/stage/peak/arm;d.mkdir(parents=True,exist_ok=True)
                n=source
                n=re.sub(r"VERSION_PHOTOMETRY\s*=.*",f"VERSION_PHOTOMETRY = 'RSR_{arm}'",n)
                n=re.sub(r'TEXTFILE_PREFIX\s*=.*',f"TEXTFILE_PREFIX = '{d}/fit'",n)
                if stage!='full':n=re.sub(r"SNCCID_LIST\s*=.*",f"SNCCID_LIST = '{cids[0]}'",n)
                if peak=='fixed':n=n.replace('  &FITINP','  &FITINP\n INISTP_PEAKMJD = 0.0')
                (d/'fit.nml').write_text(n)
    inputpaths.extend(p for p in (OUT/'data').rglob('*') if p.is_file())
    inputpaths.extend((OUT/'fits').rglob('fit.nml'));inputpaths.extend([OUT/'cohort.csv',OUT/'source-row-ledger.csv'])
    save(OUT/'input-manifest.json',dict(status='Frozen before native fits',cohort=cids,arms=['R','A','B','H'],timing=['free','fixed'],inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputpaths},lineage_error_bridge='Not exact; raw author error retained in A/B, hybrid in H',native_command=[str(EXE),'absolute fit.nml'],resource='one process one thread; engineering runtime gate first'))
    print(json.dumps({'prepared':cids,'negative_candidates':sum(r['author_nonpositive'] for r in cohort),'manifest_sha256':sha(OUT/'input-manifest.json')}))

def run(stage,peak,arm):
    m=json.loads((OUT/'input-manifest.json').read_text())
    for p,h in m['inputs_sha256'].items():assert sha(ROOT/p)==h,p
    d=OUT/'fits'/stage/peak/arm;assert not (d/'fit.log').exists()
    env=os.environ.copy();env.update(SNANA_DIR=str(BUILD),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    start=time.monotonic()
    with (d/'fit.log').open('x') as f:
        proc=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=900)
    out=dict(returncode=proc.returncode,seconds=time.monotonic()-start,manifest_sha256=sha(OUT/'input-manifest.json'),log_sha256=sha(d/'fit.log'))
    save(d/'execution.json',out);print(json.dumps(out))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('--stage',default='engineering');p.add_argument('--peak',default='free');p.add_argument('--arm',default='R');a=p.parse_args()
    if a.action=='prepare':prepare()
    else:run(a.stage,a.peak,a.arm)
