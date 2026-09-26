#!/usr/bin/env python3
import hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
files=[
'runs/research_2026_09_26/astra_design/csp_optical_filter_identification/Krisciunas2017.pdf',
'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316/README',
'phase2/official/inputs/SNDATA_ROOT/filters/CSP/CSP_TAMU_20180316/J_DUP_TAM_scan_atm.dat',
'phase2/official/inputs/SNDATA_ROOT/filters/CSP/CSP_TAMU_20180316/Jrc2_SWO_TAM_scan_atm.dat',
'phase2/official/inputs/SNDATA_ROOT/standards/bd_17d4708_stisnic_007.dat',
'phase2/official/inputs/SNDATA_ROOT/snsed/Hsiao07.dat',
'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.input',
'phase2/official/inputs/SNDATA_ROOT/kcor/Pantheon+/calib_Pantheon+_CSPDR3.input',
'sources/repos/RickKessler__SNANA@v11_04k/src/kcor.c',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/snpy-acquisition.json',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/snpy/snpy/filters/filters/LCO/Dupont/filters.dat',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/snpy/snpy/filters/filters/LCO/Dupont/J_DUP_TAM_scan_atm.dat',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/snpy/snpy/filters/filters/LCO/Swope/filters.dat',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/snpy/snpy/filters/filters/LCO/Swope/J_SWO_TAM_scan_atm.dat',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/bd17-photometry-protocol.json',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/bd17_photometry.py',
'runs/research_2026_09_26/csp_wirc_operator_feasibility/bd17-photometry-result.json',
]
m={'purpose':'CSP physical WIRC-J source and numeric-feasibility audit only', 'files':[]}
for f in files:
 b=(root/f).read_bytes()
 m['files'].append({'path':f,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
(out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
