#!/usr/bin/env python3
"""Check the inherited age model's foreground errors and printed dust contract.

Uses a stub map/law only to test preserved function data flow; synthetic changes
are not measurements of real map errors. The dust-law comparison is evaluated
directly from the pinned FSPS formula/defaults.
"""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"runs/assumption_audit/pantheon"
SRC=OUT/"sources"
UPSTREAM=ROOT/"sources/repos/benjaminrose__mc-age/util.py"


def main():
    tree=ast.parse(UPSTREAM.read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="correct_dust")
    class StubMap:
        def __init__(self,path):pass
        def ebv(self,ra,dec):return .04
    ns=dict(np=np,pd=pd,SFDMap=StubMap,fitzpatrick99=lambda wavelength,av:av*(5500/wavelength))
    exec(compile(ast.Module(body=[function],type_ignores=[]),str(UPSTREAM),"exec"),ns)
    table=pd.DataFrame([dict(ra=1.,dec=2.,u=20.,g=19.,r=18.,i=17.,z=16.,
                             err_u=.1,err_g=.08,err_r=.05,err_i=.05,err_z=.07)])
    out=ns["correct_dust"](table)
    changed=[col for col in table if not np.array_equal(table[col],out[col])]
    assert changed==["u","g","r","i","z"]
    assert np.array_equal(table[["u","g","r","i","z"]].iloc[0],[20,19,18,17,16])
    wavelengths=np.array([3543.,4770.,5500.,6231.,7625.,9134.])
    d2=.3;d1=2*d2;ratio=wavelengths/5500
    printed_tau=d1*ratio**(-.7)
    executed_tau=d1*ratio**(-1.)+d2*ratio**(-.7)
    factor=2.5/np.log(10.)
    pd.DataFrame(dict(wavelength_A=wavelengths,printed_young_A_mag=factor*printed_tau,
                     implemented_young_A_mag=factor*executed_tau,
                     difference_mag=factor*(executed_tau-printed_tau))).to_csv(OUT/"host-dust-contract.csv",index=False)
    result=dict(foreground_test="Preserved upstream function with synthetic deterministic map/law stubs",
        changed_columns=changed,photometric_error_columns_changed=False,
        covariance_returned=False,original_table_preserved=True,
        paper_young_optical_depth="dust1*(lambda/5500)^(-0.7)",
        source_young_optical_depth="dust1*(lambda/5500)^(-1)+dust2*(lambda/5500)^(-0.7)",
        dust1=d1,dust2=d2,
        Vband_printed_A_mag=float(factor*d1),Vband_implemented_A_mag=float(factor*(d1+d2)),
        interpretation="Documented-vs-implemented host attenuation discrepancy. No estimate of a SN extinction error or resulting age bias.")
    (OUT/"host-dust-contract.json").write_text(json.dumps(result,indent=2)+"\n")
    inputs=[UPSTREAM,Path(__file__),SRC/"fsps-ae31b2f-add_dust.f90",SRC/"fsps-ae31b2f-attn_curve.f90",SRC/"fsps-ae31b2f-sps_vars.f90"]
    (OUT/"host-dust-contract-manifest.json").write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
