from pathlib import Path
import shutil,hashlib,json,datetime
R=Path.cwd();O=Path(__file__).resolve().parent;A=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d';S=A/'build';D=O/'build-v11_04d';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not D.exists();shutil.copytree(S,D,symlinks=True)
script='''#!/usr/bin/env bash
set -euo pipefail
task_root=/home/szymon/Documents/ChatGPT/supernova
task_build="'''+str(D)+'''"
export SNANA_DIR="$task_build"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset CERN_DIR ROOT_DIR SNANA_PYTHON_DIR
cd "$task_build/src"
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make -j1 ../bin/snlc_sim.exe FFC="$task_compiler" LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
'''
(O/'build.sh').write_text(script)
sourcefiles={str(p.relative_to(R)):sha(p) for p in sorted(S.rglob('*')) if p.is_file() and not p.is_symlink()}
for p,h in sourcefiles.items():assert sha(D/Path(p).relative_to(S.relative_to(R)))==h
x={'status':'Pre-build freeze; root authorizes generator implementation readiness only, no photons/fits.','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':'10ec91297e4482d593cb5d3d055b10d4aa915071','private_copy':str(D.relative_to(R)),'build_cap_seconds':180,'command':['bash',str((O/'build.sh').relative_to(R))],'parallelism':1,'build_script_sha256':sha(O/'build.sh'),'copied_input_files_sha256':sourcefiles,'protected_original_fitter_sha256':sha(S/'bin/snlc_fit.exe'),'upstream_build_result_sha256':sha(A/'build-result.json'),'source_changes':'None beyond source flags already present and reviewed in copied fitter build; generator uses same source and existing compiled library/objects by hash.'}
(O/'build-protocol.json').write_text(json.dumps(x,indent=2)+'\n');print('buildprotocol',sha(O/'build-protocol.json'),'files',len(sourcefiles))
