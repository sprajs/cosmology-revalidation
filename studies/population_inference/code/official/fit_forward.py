#!/usr/bin/env python3
"""Fit one completed forward simulation with the immutable nominal DES config.

Example: python fit_forward.py --version PH2_pilot02_P21 --label forward_p21
Use --prepare-only to write the frozen namelist without launching a fit.
"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official'
ap=argparse.ArgumentParser();ap.add_argument('--version',required=True);ap.add_argument('--label',required=True);ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
assert all(c.isalnum() or c in '_-' for c in args.label+args.version)
data=ROOT/'phase2/literature/simulations/outputs';version=data/args.version
assert version.is_dir(),version
template=OUT/'inputs/snana_mock0001.nml'
s=template.read_text().replace(str(ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS/SNIa_SIMULATIONS'),str(data)).replace('PIP_D5YR_SIM_V2_DATADESSIM_4D_P21-0001',args.version).replace('snana_mock0001','snana_'+args.label)
config=OUT/f'inputs/snana_{args.label}.nml'
if config.exists():assert config.read_text()==s,'Existing different configuration must not be overwritten'
else:config.write_text(s)
log=OUT/f'diagnostics/snana-{args.label.replace("_","-")}.log'
if args.prepare_only:print(config);raise SystemExit
assert not log.exists(),'Preserve existing run; choose a new label'
env=os.environ.copy();env.update(SNANA_DIR=str(OUT/'build/SNANA-current'),SNDATA_ROOT=str(OUT/'inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(OUT/'build/sysroot/usr/lib'))
with log.open('w') as f:
    proc=subprocess.run([str(OUT/'build/SNANA-current/bin/snlc_fit.exe'),str(config)],stdout=f,stderr=subprocess.STDOUT,env=env)
record={'version':args.version,'label':args.label,'returncode':proc.returncode,'config':str(config.relative_to(ROOT)),'config_sha256':hashlib.sha256(config.read_bytes()).hexdigest(),'executable_sha256':hashlib.sha256((OUT/'build/SNANA-current/bin/snlc_fit.exe').read_bytes()).hexdigest(),'log':str(log.relative_to(ROOT)),'graceful':'ENDING PROGRAM GRACEFULLY.' in log.read_text()}
(OUT/f'diagnostics/forward-fit-{args.label}.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2));raise SystemExit(proc.returncode or (0 if record['graceful'] else 1))
