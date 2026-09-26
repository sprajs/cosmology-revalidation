#!/usr/bin/env python
# D. Jones - 4/26/21

import glob
import snana
import numpy as np
import astropy.table as at
from txtobj import txtobj
import os

_outdir="lcs/DEHVILS_DR1_ATLAS/"
_outdir_pkmjd="lcs/DEHVILS_DR1_PKMJD/"

def add_pkmjds():

    for indir,outdir,fitresfile in zip([
            'data/Photometry/CSPDR3_RAISIN',
            'data/Photometry/PS1_RAISIN',
            'data/Photometry/DES_RAISIN'],
            ['data/Photometry/CSPDR3_RAISIN_PKMJD',
             'data/Photometry/PS1_RAISIN_PKMJD',
             'data/Photometry/DES_RAISIN_PKMJD'],
            ['output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
             'output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
             'output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']):
    
        if 'CSP' in indir: files = glob.glob(f"{indir}/*PKMJD.DAT")
        else: files = glob.glob(f"{indir}/*snana.dat")
        fr = txtobj(fitresfile,fitresheader=True)

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        for f in files:
            sn = snana.SuperNova(f)
            #sn.PEAKMJD = sn.SEARCH_PEAKMJD
            iPkMJD = fr.CID == sn.SNID
            if not len(fr.CID[iPkMJD]): continue

            with open(f"{outdir}/{f.split('/')[-1]}","w") as fout, open(f) as fin:
                for line in fin:
                    if line.startswith('SEARCH_PEAKMJD:') or line.startswith('PEAKMJD'):
                        print(f'PEAKMJD: {fr.PKMJD[iPkMJD][0]:.1f}',file=fout)
                    else:
                        print(line.replace('\n',''),file=fout)

    
if __name__ == "__main__":
    add_pkmjds()
