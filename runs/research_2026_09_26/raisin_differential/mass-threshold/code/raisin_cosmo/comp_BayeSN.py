#!/usr/bin/env python
# D. Jones - 4/23/21
# see what the typical difference between BayeSN and SNooPy
# distances is for the purpose of systematics

import numpy as np
from txtobj import txtobj
import cosmo
import pylab as plt; plt.ion()

def main():
    frs_lowz = txtobj('output/fit_nir_sys_oldpkmjds/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_raisin = txtobj('output/fit_nir_sys_oldpkmjds/PS1_RAISIN/PS1_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_des = txtobj('output/fit_nir_sys_oldpkmjds/DES_RAISIN/DES_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    for k in frs_raisin.__dict__.keys():
        frs_raisin.__dict__[k] = np.append(frs_raisin.__dict__[k],frs_des.__dict__[k])
    
    frb_lowz = txtobj('output/BayeSN/CSP_RAISIN_BayeSN_NIR_theta_-1_dists_jones_snoopy_tmax.txt',fitresheader=False)
    frb_raisin = txtobj('output/BayeSN/RAISIN_BayeSN_theta_-1_NIR_dists_jones_tmax.txt',fitresheader=False)

    # low-z

    muls,mulb,mulserr,mulberr,z = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    for j,i in enumerate(frs_lowz.CID):
        iB = np.where(frb_lowz.sn == i)[0]
        if not len(iB): continue
        muls = np.append(muls,frs_lowz.DLMAG[j])#-cosmo.mu(frs_lowz.zHD))
        mulb = np.append(mulb,frb_lowz.__dict__['mu_mu+eta'][iB])
        mulserr = np.append(mulserr,frs_lowz.DLMAGERR[j])
        mulberr = np.append(mulberr,frb_lowz.__dict__['muerr_mu+eta'][iB])
        z = np.append(z,frs_lowz.zHD[j])
        
    # raisin
    for j,i in enumerate(frs_raisin.CID):
        iB = np.where(frb_raisin.sn == i)[0]
        if not len(iB): continue
        muls = np.append(muls,frs_raisin.DLMAG[j])#-cosmo.mu(frs_raisin.zHD))
        mulb = np.append(mulb,frb_raisin.__dict__['mu_mu+eta'][iB])
        mulserr = np.append(mulserr,frs_raisin.DLMAGERR[j])
        mulberr = np.append(mulberr,frb_raisin.__dict__['muerr_mu+eta'][iB])
        z = np.append(z,frs_raisin.zHD[j])

    iZ1 = z < 0.1
    iZ2 = z > 0.1
    muls -= np.average(muls-mulb,weights=1/(mulserr**2.+mulberr**2.))
    print(np.average(muls[iZ1]-mulb[iZ1],weights=1/(mulserr[iZ1]**2.+mulberr[iZ1]**2.)))
    print(np.average(muls[iZ2]-mulb[iZ2],weights=1/(mulserr[iZ2]**2.+mulberr[iZ2]**2.)))
    
    import pdb; pdb.set_trace()

def fig():

    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    for ax in [ax1,ax2]:
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    
    frs_lowz = txtobj('output/fit_nir_sys_oldpkmjds/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_raisin = txtobj('output/fit_nir_sys_oldpkmjds/PS1_RAISIN/PS1_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    frs_des = txtobj('output/fit_nir_sys_oldpkmjds/DES_RAISIN/DES_RAISIN/FITOPT000.FITRES.gz',fitresheader=True)
    for k in frs_raisin.__dict__.keys():
        frs_raisin.__dict__[k] = np.append(frs_raisin.__dict__[k],frs_des.__dict__[k])
    
    frb_lowz = txtobj('output/BayeSN/CSP_RAISIN_BayeSN_NIR_theta_-1_dists_jones_snoopy_tmax.txt',fitresheader=False)
    frb_raisin = txtobj('output/BayeSN/RAISIN_BayeSN_theta_-1_NIR_dists_jones_tmax.txt',fitresheader=False)

    muls,mulb,mulserr,mulberr,z = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    for j,i in enumerate(frs_lowz.CID):
        iB = np.where(frb_lowz.sn == i)[0]
        if not len(iB): continue
        muls = np.append(muls,frs_lowz.DLMAG[j])#-cosmo.mu(frs_lowz.zHD))
        mulb = np.append(mulb,frb_lowz.__dict__['mu_mu+eta'][iB])
        mulserr = np.append(mulserr,frs_lowz.DLMAGERR[j])
        mulberr = np.append(mulberr,frb_lowz.__dict__['muerr_mu+eta'][iB])
        z = np.append(z,frs_lowz.zHD[j])

    # raisin
    for j,i in enumerate(frs_raisin.CID):
        iB = np.where(frb_raisin.sn == i)[0]
        if not len(iB): continue
        muls = np.append(muls,frs_raisin.DLMAG[j])#-cosmo.mu(frs_raisin.zHD))
        mulb = np.append(mulb,frb_raisin.__dict__['mu_mu+eta'][iB])
        mulserr = np.append(mulserr,frs_raisin.DLMAGERR[j])
        mulberr = np.append(mulberr,frb_raisin.__dict__['muerr_mu+eta'][iB])
        z = np.append(z,frs_raisin.zHD[j])

    
    ax1.set_ylabel(r'BayeSN: HR $-$ ${\rm \overline{HR}}$ (mag)')
    ax1.set_xlabel(r'SNooPy: HR $-$ ${\rm \overline{HR}}$ (mag)')

    iLowz = z < 0.1
    iHighz = z > 0.1
    ax1.plot([-0.8,0.8],[-0.8,0.8],color='k',lw=2)
    ax1.errorbar(muls[iLowz]-cosmo.mu(z[iLowz])-np.median(muls-cosmo.mu(z)),
                 mulb[iLowz]-cosmo.mu(z[iLowz])-np.median(mulb-cosmo.mu(z)),
                 xerr=mulserr[iLowz],yerr=mulberr[iLowz],fmt='.',color='b',
                 label='low-$z$')
    ax1.errorbar(muls[iHighz]-cosmo.mu(z[iHighz])-np.median(muls-cosmo.mu(z)),
                 mulb[iHighz]-cosmo.mu(z[iHighz])-np.median(mulb-cosmo.mu(z)),
                 xerr=mulserr[iHighz],yerr=mulberr[iHighz],fmt='.',color='r',
                 label='RAISIN')
    ax1.set_xlim([-0.8,0.8])
    ax1.set_ylim([-0.8,0.8])
    ax1.legend(loc='upper left')
    
    ax2.set_xlabel('$z_{CMB}$',fontsize=15)
    ax2.set_ylabel(r'SNooPy $-$ BayeSN: $\Delta\mu_{S-B} - \overline{\Delta\mu_{S-B}}$ (mag)')

    ax2.axhline(0,color='k',lw=2)
    ax2.errorbar(z,muls-mulb-np.median(muls-mulb),yerr=np.sqrt(mulserr**2.+mulberr**2.),fmt='.',color='k')
    ax2.set_xscale('log')
    ax2.set_xlim([0.01,1.0])
    
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    #main()
    fig()
