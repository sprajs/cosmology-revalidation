#!/usr/bin/env python
# D. Jones - 12/10/20
# update raisin photometry data files with the latest
# and greatest

import numpy as np
import glob
import os
from txtobj import txtobj
import snana

lameffdict = {'J':12516.293,
              'H':15392.331}

def main():

    rp = txtobj('../phot/raisin_photometry.txt')
    files = glob.glob('data/Photometry/PS1_RAISIN/*.snana.dat')
    files = np.append(files,glob.glob('data/Photometry/DES_RAISIN/*.snana.dat'))
    yf = txtobj('data/YtmplFlux.txt')
    
    for f in files:
        newfile = f.replace('.snana.dat','_latest.snana.dat')
        sysfile = f.replace('.snana.dat','_sys.snana.dat')
        if 'PS1_RAISIN' in sysfile:
            sysfile = sysfile.replace('/PS1_RAISIN/','/PS1_RAISIN_TMPLSYS/')
        elif 'DES_RAISIN' in sysfile:
            sysfile = sysfile.replace('/DES_RAISIN/','/DES_RAISIN_TMPLSYS/')
        else: raise RuntimeError('not sure which survey!')
        print(f)
        sn = snana.SuperNova(f)
        with open(f) as fin, open(newfile,'w') as fout, open(sysfile,'w') as fout2:
            for line in fin:
                if not line.startswith('OBS:'):
                    print(line.replace('\n',''),file=fout)
                    print(line.replace('\n',''),file=fout2)
                else:
                    filt = line.split()[2]
                    if filt == 'H' or filt == 'J':
                        # then figure out the new phot
                        iPhot = np.where(np.abs(rp.MJD-float(line.split()[1])) < 0.002)[0]
                        if len(iPhot) != 1:
                            raise RuntimeError('halp meh')
                        flux = 10**(-0.4*(rp.diffmag_aper[iPhot][0]-27.5))
                        flux_err = flux*0.4*np.log(10)*rp.diffmagerr_aper[iPhot][0]
                        mjd = rp.MJD[iPhot][0]

                        # correction needed?
                        lameff_divide = 13954.312
                        if filt == 'H' and lameffdict['H']/(1+sn.z) < lameff_divide:
                            apply_correction = True
                        elif filt == 'J' and lameffdict['J']/(1+sn.z) < lameff_divide:
                            apply_correction = True
                        else: apply_correction = False
                        
                        # apply correction?
                        if apply_correction:
                            iSN = sn.SNID == yf.snid
                            print(sn.SNID)
                            flux_factor = 1+10**(-0.4*yf.mags[iSN][0])
                            flux *= flux_factor
                            flux_err *= flux_factor
                            print(sn.SNID,flux_factor)
                        
                            # sys. correction
                            flux_factor_sys = 1+10**(-0.4*yf.mags_min[iSN][0])
                            flux_sys = flux/flux_factor*flux_factor_sys
                            flux_err_sys = flux_err/flux_factor*flux_factor_sys
                        else:
                            flux_sys = flux
                            flux_err_sys = flux_err
                            
                        newline = f"OBS:  {mjd:.3f}  {filt} NULL   {flux:.3f} {flux_err:.3f}"
                        print(newline,file=fout)
                        newline_sys = f"OBS:  {mjd:.3f}  {filt} NULL   {flux_sys:.3f} {flux_err_sys:.3f}"
                        print(newline_sys,file=fout2)
                    else:
                        print(line.replace('\n',''),file=fout)
                        print(line.replace('\n',''),file=fout2)
                        
        #import pdb; pdb.set_trace()
        os.system(f"mv {newfile} {f}")
        
if __name__ == "__main__":
    main()
