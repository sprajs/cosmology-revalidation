"""Post-primary, prespecified numerical stability controls; preserve all initial fits."""
from pathlib import Path
import argparse,json,os,re,subprocess,time
import numpy as np
from astropy.io import fits
import paired_refit as p
O=p.OUT

def prepare():
    assert not (O/'stability-protocol.json').exists()
    jobs=[]
    for start in [.85,1.15]:
        for peak in ['free','fixed']:
            for arm in ['R','A','B','H']:
                d=O/'fits/multistart'/str(start)/peak/arm;d.mkdir(parents=True)
                s=(O/'fits/full'/peak/arm/'fit.nml').read_text()
                s=re.sub(r'(?m)^\s*INIVAL_SHAPE\s*=.*$',f' INIVAL_SHAPE = {start}',s)
                s=s.replace(str(O/'fits/full'/peak/arm/'fit'),str(d/'fit'))
                (d/'fit.nml').write_text(s)
                jobs.append(dict(kind='multistart',start=start,peak=peak,arm=arm,nml=str((d/'fit.nml').relative_to(p.ROOT)),sha256=p.sha(d/'fit.nml')))
    grid=fits.getdata(p.REL/'model/snoopy.B18/SNooPy_B18.fits','LUMI-GRID').field(0).astype(float)
    # IDs named from primary source implementation concern, not from distance-shift sign.
    from check_results import table
    fr={r['CID']:r for r in table(O/'fits/full/free/B/fit.FITRES.TEXT')}
    for cid in ['DES16C1cim','DES16S1agd','DES16E1dcx']:
        value=float(fr[cid]['STRETCH']);k=int(np.argmin(abs(grid-value)));step=(grid[-1]-grid[0])/(len(grid)-1)
        for offset in [-1.,-.25,-.02,0.,.02,.25,1.]:
            shape=float(grid[k]+offset*step)
            if not(grid[0]<=shape<=grid[-1]):continue
            d=O/'fits/profile'/cid/str(offset);d.mkdir(parents=True)
            s=(O/'fits/engineering/free/R/fit.nml').read_text().replace("'RSR_R'","'RSR_B'")
            s=s.replace("SNCCID_LIST = 'DES16C1cim'",f"SNCCID_LIST = '{cid}'")
            s=re.sub(r'(?m)^\s*INIVAL_SHAPE\s*=.*$',f' INIVAL_SHAPE = {shape:.17g}',s)
            s=s.replace('  &FITINP','  &FITINP\n INISTP_SHAPE = 0.0')
            s=s.replace(str(O/'fits/engineering/free/R/fit'),str(d/'fit'))
            (d/'fit.nml').write_text(s)
            jobs.append(dict(kind='profile',CID=cid,grid_index=k,offset_cells=offset,fixed_shape=shape,nml=str((d/'fit.nml').relative_to(p.ROOT)),sha256=p.sha(d/'fit.nml')))
    p.save(O/'stability-protocol.json',dict(status='Frozen after primary results reveal initialization dependence; before new stability fits',primary_manifest_sha256=p.sha(O/'input-manifest.json'),source_sha256=p.sha(Path(__file__)),jobs=jobs,decision_rule='No outcome-sign based minimum selection. Compare objectives only with identical native accepted row masks; otherwise preserve as differing-mask solutions. Flag nonquadratic/grid-knot profiles. Local profile is not proof of global convergence.',resource='one process, oneBLAS thread; two complete alternate starts plus21 one-object profile fits'))
    print(len(jobs),'jobs prepared')

def run(kind):
    pr=json.loads((O/'stability-protocol.json').read_text());assert p.sha(Path(__file__))==pr['source_sha256']
    env=os.environ.copy();env.update(SNANA_DIR=str(p.BUILD),SNDATA_ROOT=str(p.ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(p.ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for j in pr['jobs']:
        if j['kind']!=kind:continue
        n=p.ROOT/j['nml'];assert p.sha(n)==j['sha256'];d=n.parent;assert not(d/'fit.log').exists()
        t=time.monotonic()
        with (d/'fit.log').open('x') as f:r=subprocess.run([str(p.EXE),str(n)],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=900)
        result=dict(returncode=r.returncode,seconds=time.monotonic()-t,protocol_sha256=p.sha(O/'stability-protocol.json'),log_sha256=p.sha(d/'fit.log'));p.save(d/'execution.json',result)
        print(json.dumps(dict(job=j['nml'],**result)),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','multistart','profile']);a=ap.parse_args()
    if a.action=='prepare':prepare()
    else:run(a.action)
