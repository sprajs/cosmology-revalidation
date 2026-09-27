#!/usr/bin/env python3
"""Rebuild documented native OPT_SETPKMJD=16 clump timing on public SMP flux.

The caller supplies an already built SNANA tree; executable/output are isolated.
No original scientific product is edited. Compiler/object identities are recorded.
"""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from common import ROOT,WORK,RESULTS,sha

def main():
    p=argparse.ArgumentParser();p.add_argument('--snana-tree',type=Path,default=ROOT/'.work/survey-physics/SNANA');p.add_argument('--sndataroot',type=Path,default=ROOT/'.work/survey-physics/SNDATA_ROOT');p.add_argument('--library-path',default=str(Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/official/build/sysroot/usr/lib'));a=p.parse_args()
    folder=WORK/'native-dataprep';(folder/'obj').mkdir(parents=True,exist_ok=True);(folder/'bin').mkdir(exist_ok=True)
    exe=folder/'bin/snana.exe';src=a.snana_tree.resolve()
    if not exe.exists():
        for obj in (src/'obj').glob('*.o'):shutil.copy2(obj,folder/'obj'/obj.name)
        args=['make','-C',str(src/'src'),'snana','SRC='+str(src/'src'),'OBJ='+str(folder/'obj'),'BIN='+str(folder/'bin')]
        with (folder/'build.log').open('w') as log:subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True)
    nml=folder/'dataprep.nml';nml.write_text(f"""&SNLCINP
 OPT_SETPKMJD = 16
 SNTABLE_LIST = 'SNANA(text:key)'
 TEXTFILE_PREFIX = 'dataprep'
 OPT_YAML = 1
 PRIVATE_DATA_PATH = '{WORK/'release/0_DATA'}'
 VERSION_PHOTOMETRY = 'DES-SN5YR_DES'
 PHOTFLAG_MSKREJ = 1016
&END
""")
    env=os.environ.copy();env.update(SNANA_DIR=str(src),SNDATA_ROOT=str(a.sndataroot.resolve()),LD_LIBRARY_PATH=a.library_path,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    with (folder/'dataprep.log').open('w') as log:r=subprocess.run([str(exe),str(nml)],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    log=(folder/'dataprep.log').read_text();ok=r.returncode==0 and 'ENDING PROGRAM GRACEFULLY' in log
    record={'status':'passed' if ok else 'failed','code_sha256':sha(__file__),'command_arguments':vars(a)|{'snana_tree':str(a.snana_tree),'sndataroot':str(a.sndataroot)},'returncode':r.returncode,'dependencies':{str(q.relative_to(ROOT)) if q.is_relative_to(ROOT) else str(q):sha(q) for q in [exe,nml,src/'src/snana.F90',src/'src/Makefile',folder/'build.log',folder/'dataprep.log']},'native_object_sha256':{str(q.relative_to(ROOT)):sha(q) for q in (folder/'obj').glob('*.o')},'output_sha256':{str(q.relative_to(ROOT)):sha(q) for q in folder.glob('dataprep.*') if q.suffix in ['.TEXT','.YAML']},'scope':'Native recordedcurrentSNANA clumpestimator frompublicSMP; exactauthorhistoricaldataprepbinary/outputnotassumed.'}
    (RESULTS/'dataprep.json').write_text(json.dumps(record,indent=2)+'\n');print('dataprep',ok);assert ok
if __name__=='__main__':main()
