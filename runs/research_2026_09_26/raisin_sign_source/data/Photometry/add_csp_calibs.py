#!/usr/bin/env python
# D. Jones - 12/17/21
"""put Cepheid calibs with CSP data into the correct format"""

import pylab as plt
import numpy as np
import os
import glob
import astropy.table as at

massdict = {'2012fr':10.3991,'2012ht':8.39097,'2015F':10.6429}
masserrdict = {'2012fr':0.01,'2012ht':0.01,'2015F':0.01}
radict = {'2012fr':53.3999583,'2012ht':163.3447917,'2015F':114.0656667}
decdict = {'2012fr':-36.1271389,'2012ht':16.7763611,'2015F':-69.5063889}
ebvdict = {'2012fr':0.018,'2012ht':0.025,'2015F':0.175}
zdict = {'2012fr':0.005476,'2012ht':0.003555,'2015F':0.004846}
zcmbdict = {'2012fr':5.1514123081166652E-003,'2012ht':4.6858315838131936E-003,'2015F':5.2405429582595975E-003}
pkmjddict = {'2012fr':56242.0469,'2012ht':56295.6992,'2015F':57107.1055}
nobsdict = {'2012fr':516,'2012ht':362,'2015F':351}

def main():
    for snid in ['2012fr','2012ht','2015F']:
        
        snanaheader = f"""SURVEY:   CSP
SNID:     {snid}
FILTERS:  ugriBnmoJjYyH
RA:         {radict[snid]}  # deg
DEC:       {decdict[snid]}  # deg
REDSHIFT_HELIO: {zdict[snid]} +- 0.0001  # Helio
REDSHIFT_FINAL: {zcmbdict[snid]} +- 0.0001  # CMB
MWEBV:     {ebvdict[snid]:.3f} +-  {ebvdict[snid]*0.05:.3f}        # SFD98 * 0.86
PEAKMJD: {pkmjddict[snid]}
HOSTGAL_LOGMASS: {massdict[snid]} +- {masserrdict[snid]}
FAKE: 0

# -----------------------------------
NOBS: {nobsdict[snid]:.0f}
NVAR: 7
VARLIST: MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG   MAGERR
"""
        with open(f'CSPDR3_RAISIN/CSPDR3_{snid}.PKMJD.DAT','w') as fout:
            print(snanaheader,file=fout)
        
            # for 2012fr, we have to parse optical and NIR files
            if snid == '2012fr':
                opt = at.Table.read('2012fr_contreras_opt.txt',format='ascii')
                nir = at.Table.read('2012fr_contreras_NIR.txt',format='ascii')

                for d in opt:
                    for filt in 'ugriBV':
                        if d[filt] < 0: continue
                        flux = 10**(-0.4*(d[filt]-27.5))
                        fluxerr = flux*0.4*np.log(10)*d[filt+'err']/1000.
                        print(f"OBS: {d['JD']-2400000.5:.3f} {filt.replace('V','n')} VOID {flux:.4f} {fluxerr:.4f} {d[filt]:.4f} {d[filt+'err']:.4f}",file=fout)
                for d in nir:
                    for filt in 'YJH':
                        if d[filt] < 0: continue
                        flux = 10**(-0.4*(d[filt]-27.5))
                        fluxerr = flux*0.4*np.log(10)*d[filt+'err']/1000.
                        print(f"OBS: {d['JD']-2400000.5:.3f} {filt.replace('V','n')} VOID {flux:.4f} {fluxerr:.4f} {d[filt]:.4f} {d[filt+'err']:.4f}",file=fout)
                        
            # for 2012ht, 2015F, these are in a more normal format
            if snid == '2012ht' or snid == '2015F':
                data = at.Table.read('csp_calibs.txt',format='ascii',delimiter='|')
                for d in data:
                    if d['SN'] == 'SN'+snid:
                        if d['omag'] < 0: continue
                        flux = 10**(-0.4*(d['omag']-27.5))
                        fluxerr = flux*0.4*np.log(10)*d['_omag']
                        print(f"OBS: {d['MJD']:.3f} {d['f'].replace('V','n')} VOID {flux:.4f} {fluxerr:.4f} {d['omag']:.4f} {d['_omag']:.4f}",file=fout)
            print('END:',file=fout)
            
if __name__ == "__main__":
    main()
