"""Explicitly select equation-reconstruction variants; never label author arrays."""
from core import ROOT,manifest
import pandas as pd

p=ROOT/'runs/mapping/csfh-dtd-curves.csv'
d=pd.read_csv(p)
out=ROOT/'runs/cosmology/templates';out.mkdir(parents=True,exist_ok=True)
outputs=[];config=[]
for cos,dtd,stat,label in [
    ('lcdm_H70','C14_smooth','median','c14-lcdm70-median'),
    ('son_cpl_H63p6','C14_smooth','median','c14-cpl63-median'),
    ('son_cpl_H63p6','C14_smooth','mean','c14-cpl63-mean'),
    ('son_cpl_H63p6','W26_cut40','mean','cut40-cpl63-mean')]:
    v=d[(d.cosmology==cos)&(d.csfh=='B13')&(d.dtd==dtd)]
    assert len(v)>100 and v.z.max()>=2.3
    col=f'correction_subtracted_mag_{stat}_delay_gyr'
    dest=out/f'{label}.csv'
    v[['z',col]].rename(columns={col:'delta_mu'}).to_csv(dest,index=False)
    outputs.append(dest);config.append({'file':str(dest.relative_to(ROOT)),'cosmology':cos,'csfh':'B13','dtd':dtd,'statistic':stat})
manifest(out,'Fixed-cosmology sensitivity templates from independent equations; not exact Son correction',[p],
         {'variants':config,'slope_mag_per_Gyr':.03,'cosmology_recomputed_during_fit':False},outputs)
