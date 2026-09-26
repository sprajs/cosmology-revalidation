#!/usr/bin/env python3
"""Constructed broadband check of low-RV F99 extrapolation, not a data fit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from flux_response import build_model, ROOT


def main():
    out=ROOT/'runs/salt_dust_audit/f99_support'
    out.mkdir(parents=True,exist_ok=True)
    model,bands,paths,_=build_model()
    rows=[]; unsupported=[]
    for z in [.06,.1,.2,.5,.8,1.1]:
        for phase in [-10.,0.,20.]:
            model.set(z=z,t0=0,x0=1e-4,x1=0,c=0,mwebv=0,hostebv=0,hostrv=3.1)
            t=phase*(1+z)
            active={name:b for name,b in bands.items() if bool(model.bandoverlap(b))}
            unsupported.extend({'z':z,'phase':phase,'band':name} for name in bands if name not in active)
            base={name:float(model.bandflux(b,t)) for name,b in active.items()}
            for rv in [.3,.4,.5,.66,.8,1.,1.5,2.,3.1]:
                model.set(hostebv=.1,hostrv=rv)
                for name,b in active.items():
                    reddened=float(model.bandflux(b,t))
                    assert base[name]>0 and reddened>0
                    effective_rest=b.wave_eff/(1+z)
                    rows.append(dict(z=z,phase=phase,RV=rv,EBV=.1,band=name,
                                     effective_rest_wavelength_A=effective_rest,
                                     within_nominal_mean_wavelength_cut=bool(3500<=effective_rest<=8000),
                                     extinction_mag=-2.5*np.log10(reddened/base[name])))
    frame=pd.DataFrame(rows)
    frame.to_csv(out/'broadband_grid.csv',index=False)
    selected=frame[frame.within_nominal_mean_wavelength_cut]
    example=frame[(frame.z==.1)&(frame.phase==0)&(frame.RV==.4)&(frame.band=='i')].iloc[0]
    assert example.extinction_mag < -.02
    # Standard RV examples remain attenuating in this bounded grid.
    assert (frame.loc[frame.RV>=2,'extinction_mag']>0).all()
    result={'status':'Constructed SALT3 broadband consequence of low-RV F99 extrapolation; not an observed flux excess or population bias estimate.',
            'example':example.to_dict(),
            'grid_points_inside_diagnostic_mean_wavelength_cut':len(selected),
            'negative_extinction_grid_points_inside_cut':int((selected.extinction_mag<0).sum()),
            'unsupported_band_phase_redshift_combinations_excluded':unsupported,
            'warning':'Grid counts are not event fractions. sncosmo wave_eff implements its own effective wavelength; the cutoff flag is diagnostic, not an exact reconstruction of SNANA selection. RV<2 is extrapolation; F99 is not a validated physical absorption law there.',
            'provenance':'Uses released SALT3.DES5YR and DES filters with independent extinction.fitzpatrick99. Source-only checks separately compare this function with SNANA.',
            'next_test':'Regenerate and refit identical-seed populations with physically nonnegative extinction-law support, recompute selection/BBC, and propagate resulting distance differences; simple truncation changes the population and is not an automatic correction.'}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    paths += [Path(__file__).resolve(),ROOT/'scripts/salt_dust_audit/flux_response.py']
    sha=lambda p:hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
    (out/'manifest.json').write_text(json.dumps({'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in paths],
        'outputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in [out/'summary.json',out/'broadband_grid.csv']]},indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
