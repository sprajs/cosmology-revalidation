#!/usr/bin/env python3
"""Freeze author-accepted-mask arm; distinct from reconstructing clipping."""
from pathlib import Path
import re,json,hashlib
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official'
data=ROOT/'phase2/assumptions/conditional_mask'
parent=OUT/'inputs/conditional';parent.mkdir(exist_ok=True)
version='DES-SN5YR_DES_CONDITIONAL';link=parent/version
if not link.exists():link.symlink_to(data,target_is_directory=True)
files=[]
for i in range(4):
    label=f'des_conditioned{i}'
    s=(OUT/f'inputs/snana_des_mask32{i}.nml').read_text()
    s=s.replace('PHOTFLAG_MSKREJ = 32','PHOTFLAG_MSKREJ = 536870912')
    s=s.replace(str(ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA'),str(parent))
    s=s.replace("VERSION_PHOTOMETRY = 'DES-SN5YR_DES'",f"VERSION_PHOTOMETRY = '{version}'")
    s=s.replace(f'snana_des_mask32{i}',f'snana_{label}')
    s=re.sub(r'(?m)^(\s*DELCHI2_REJECT\s*=).*$',r'\1 1.0E9 ! no further clipping in author-mask arm',s)
    s=re.sub(r'(?m)^(\s*FITWIN_TREST\s*=).*$',r'\1 -999.0, 999.0 ! freeze epoch membership',s)
    s=re.sub(r'(?m)^(\s*CUTWIN_TREST\s*=).*$',r'\1 -999.0, 999.0 ! no new phase selection',s)
    # Extra plot output gives an independent accepted-epoch count/mask check.
    s=s.replace("SNTABLE_LIST      = 'FITRES(text:host)'","SNTABLE_LIST      = 'FITRES(text:host) LCPLOT(text:col)'")
    s=s.replace('&SNLCINP','&SNLCINP\n MXLC_PLOT = 1000',1)
    p=OUT/f'inputs/snana_{label}.nml'
    if p.exists():assert p.read_text()==s
    else:p.write_text(s)
    files.append(p)
record={'purpose':'Conditional fit at published accepted epoch mask, default starts. Tests objective/minimizer conditional on author selection, not independent clipping reproduction.','no_new_phase_mask':True,'no_new_chi2_clipping':True,'reject_flag':536870912,'source_manifest':str((data/'manifest.json').relative_to(ROOT)),'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(OUT/'inputs/conditioned-mask-contract.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
