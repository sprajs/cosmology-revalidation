#!/usr/bin/env python3
"""Prepare noise-truth misspecification independently of the fitting error model."""
from pathlib import Path
import json,hashlib,re,argparse
R=Path(__file__).resolve().parents[3];O=R/'phase2/literature/simulations';D=O/'inputs'
def keep(p,t):
 if p.exists() and p.read_text()!=t:raise RuntimeError(f'immutable {p}')
 if not p.exists():p.write_text(t)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
src=R/'phase2/official/inputs/SNDATA_ROOT/simlib/DES/DES-SN5YR_DES_FLUXERRMODEL_SIM.DAT'
text=src.read_text();scaled=text.replace('MAPNAME: FLUXERR_SCALE','MAPNAME: FLUXERR_SCALE\nSCALE_FLUXERR_TRUE: 1.2\nSCALE_FLUXERR_DATA: 1.0')
mapfile=D/'DES_FLUXERRMODEL_TRUE120_DATA100.DAT';keep(mapfile,scaled)
parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot01');args=parser.parse_args()
base=(D/f'PH2_{args.tag}_P21.input').read_text();runs=[]
for suffix,extra,noise in [('P21_rho000','FLUXERRMODEL_REDCOV: NONE\n',False),('P21_rho090','FLUXERRMODEL_REDCOV: g:0.9,r:0.9,i:0.9,z:0.9\n',False),('P21_noisetrue120','',True)]:
 version=f'PH2_{args.tag}_'+suffix
 inp=D/(version+'.input');b=base.replace(f'PH2_{args.tag}_P21',version)+extra
 if noise:b=re.sub(r'^FLUXERRMODEL_FILE:.*$',f'FLUXERRMODEL_FILE: {mapfile}',b,flags=re.M)
 keep(inp,b);runs.append({'version':version,'input':str(inp.relative_to(R)),'sha256':sha(inp)})
manifest={'purpose':'Misspecified noise injections, fitted with unchanged measured error convention and same SALT pipeline.','runs':runs,'source_map':str(src.relative_to(R)),'source_sha256':sha(src),'scaled_map':str(mapfile.relative_to(R)),'scaled_map_sha256':sha(mapfile),'definitions':{'rho000':'Disable added-fudge interepoch correlation; template common-mode remains.','rho090':'Set same-band added-fudge interepoch correlation0.9, total flux correlation is smaller/different.','noisetrue120':'Multiply true noise scale by1.2 while leaving reported DATA-error scaling1.0. This creates purposeful underreported uncertainty; not a proposed calibrated nominal model.'},'source_convention':'snlc_sim.c24303-24313 addedvarianceF=(scale²-1)*searchvariance;24506-24512 covariance rho*sqrt(FiFj);24654-24674 independentsearch+commonbandtemplate+fudge.','inference':'Perturbation values are stress brackets, not empirically estimated posterior uncertainties. Compare them to independent real-noise diagnostics before using as scientific uncertainty prior.','coverage':'Single pilot realizes misspecification and selection effect; frequentist coverage requires repeated independently seeded realizations and confidence intervals on coverage.'}
keep(D/f'preregistration-noise-{args.tag}.json',json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
