#!/usr/bin/env python
# D. Jones - 7/14/20

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
plt.ion()
import numpy as np
from txtobj import txtobj
from scipy.optimize import minimize
import scipy.stats
import os
import cosmo
import astropy.table as at
import getmu
plt.rcParams['figure.figsize'] = (11,6)
#plt.subplots_adjust(hspace=0)

#_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_CIDS.LIST',unpack=True,dtype=str),
#                           #np.loadtxt('output/goodcids/CfA_CIDS.LIST',unpack=True,dtype=str),
#                           np.loadtxt('output/goodcids/PS1_CIDS.LIST',unpack=True,dtype=str),
#                           np.loadtxt('output/goodcids/DES_CIDS.LIST',unpack=True,dtype=str)))
_goodcids = np.concatenate((np.loadtxt('output/goodcids/CSP_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/PS1_GOODCIDS_LATEST.LIST',unpack=True,dtype=str),
                            np.loadtxt('output/goodcids/DES_GOODCIDS_LATEST.LIST',unpack=True,dtype=str)))

nirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                     '$RAISIN_ROOT/cosmo/output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                     '$RAISIN_ROOT/cosmo/output/fit_nir/DES_RAISIN.FITRES.TEXT']
nirdatafitreslist = [os.path.expandvars(filepath) for filepath in nirdatafitreslist]
nirdatafitreslist = ['output/fitres_cosmo/CSP.FITRES',
                     'output/fitres_cosmo/PS1.FITRES',
                     'output/fitres_cosmo/DES.FITRES']

nirstretchdatafitreslist = ['output/fit_nir/CSP_RAISIN_NIR_SHAPE.FITRES.TEXT',
                            'output/fit_nir/PS1_RAISIN_NIR_SHAPE.FITRES.TEXT',
                            'output/fit_nir/DES_RAISIN_NIR_SHAPE.FITRES.TEXT']


opticalnirdatafitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_optnir.FITRES.TEXT',
                            '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_optnir.FITRES.TEXT',
                            '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_optnir.FITRES.TEXT']
opticalnirdatafitreslist = [os.path.expandvars(filepath) for filepath in opticalnirdatafitreslist]

salt2fitreslist = ['$RAISIN_ROOT/cosmo/output/fit_optical/CSP_RAISIN_SALT3.FITRES.TEXT',
                   '$RAISIN_ROOT/cosmo/output/fit_optical/PS1_RAISIN_SALT3.FITRES.TEXT',
                   '$RAISIN_ROOT/cosmo/output/fit_optical/DES_RAISIN_SALT3.FITRES.TEXT']
salt2fitreslist = [os.path.expandvars(filepath) for filepath in salt2fitreslist]

def checknewmasses():
    hmn = txtobj('/Users/David/Downloads/lephare_finalout_tmp.delme')
    hmo = txtobj('hosts/host_properties_raisin.txt')
    for j,i in enumerate(hmo.snid):
        if hmo.logmass[j] > 0:
            print(i)#,hmn.logmass[hmn.snid == i])

def add_ps1hosts():
    import glob
    import snana
    from coordutils import GetSexigesimalString

    eventlist = []
    ps1files = glob.glob(os.path.expandvars('$SNDATA_ROOT/lcmerge/Pantheon_PS1MD_TEXT/*txt'))
    for f in ps1files:
        sn = snana.SuperNova(f)
        ra,dec = GetSexigesimalString(sn.RA,sn.DECL)
        print(f'YSEmasterlist.py -n PSc%06i {ra} {dec} --type SNIa -z {sn.REDSHIFT_FINAL.split()[0]}'%sn.SNID)
        eventlist += [sn.SNID]
        

    print(','.join(eventlist))


def add_hosts():
    import glob
    import snana
    from coordutils import GetSexigesimalString

    eventlist = []
    ps1files = glob.glob('data/Photometry/PS1_RAISIN/*dat')
    for f in ps1files:
        sn = snana.SuperNova(f)
        ra,dec = GetSexigesimalString(sn.RA,sn.DEC)
        print(f'YSEmasterlist.py -n {sn.SNID} {ra} {dec} --type SNIa -z {sn.REDSHIFT_FINAL.split()[0]}')
        eventlist += [sn.SNID]
        
    desfiles = glob.glob('data/Photometry/DES_RAISIN/*dat')
    for f in desfiles:
        sn = snana.SuperNova(f)
        ra,dec = GetSexigesimalString(sn.RA,sn.DEC)
        print(f'YSEmasterlist.py -n {sn.SNID} {ra} {dec} --type SNIa -z {sn.REDSHIFT_FINAL.split()[0]}')
        eventlist += [sn.SNID]

    print(','.join(eventlist))
        
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)
        
def lnlikefunc(x,p_iae=None,mu_i=None,sigma_i=None,sigma=None,z=None,survey=None):

    if sigma or sigma == 0.0:
        # fix the dispersion
        x[2] = sigma; x[3] = sigma

    p_iae[np.where(p_iae == 0)] == 1e-4
    return -np.sum(np.logaddexp(-(mu_i-x[0])**2./(2.0*(sigma_i**2.+x[2]**2.)) +\
                np.log((1-0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[2]**2.+sigma_i**2.))),
                -(mu_i-x[1])**2./(2.0*(sigma_i**2.+x[3]**2.)) +\
                np.log((0.01*p_iae)/(np.sqrt(2*np.pi)*np.sqrt(x[3]**2.+sigma_i**2.)))))

def neglnlikefunc(x,p_lm=None,mu_i=None,sigma_i=None,sigma=None,survey=None,z=None):

    iCSP = survey == 5
    iPS1_midz = (survey == 15) & (z < 0.4306)
    iPS1_highz = (survey == 15) & (z >= 0.4306)
    iDES_midz = (survey == 10) & (z < 0.4306)
    iDES_highz = (survey == 10) & (z >= 0.4306)
    
    mu_lowz,mu_midz,mu_highz = x[0],x[1],x[2]
    #mu_midz = mu_highz
    #mu_lowz = mu_highz
    sigint_csp,sigint_ps1,sigint_des = x[3],x[4],x[5]
    mass_step = x[6]
    
    # sigint split by sample, but distance split by redshift
    # each one with a low-mass and high-mass component
    loglike_csp = -np.sum(np.logaddexp(-(mu_i[iCSP]+mass_step-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((1-0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))),
                                       -(mu_i[iCSP]-mu_lowz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.)))))
            
    loglike_ps1_midz = -np.sum(np.logaddexp(-(mu_i[iPS1_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((1-0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))),
                                            -(mu_i[iPS1_midz]-mu_midz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.)))))
    
    loglike_ps1_highz = -np.sum(np.logaddexp(-(mu_i[iPS1_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((1-0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))),
                                             -(mu_i[iPS1_highz]-mu_highz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.)))))
            
    loglike_des_midz = -np.sum(np.logaddexp(-(mu_i[iDES_midz]+mass_step-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((1-0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))),
                                            -(mu_i[iDES_midz]-mu_midz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.)))))
            
    loglike_des_highz = -np.sum(np.logaddexp(-(mu_i[iDES_highz]+mass_step-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((1-0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))),
                                             -(mu_i[iDES_highz]-mu_highz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.)))))

    return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz


def neglnlikefunc_cosmo(x,p_lm=None,mu_i=None,sigma_i=None,sigma=None,survey=None,z=None):

    iCSP = survey == 5
    iPS1_midz = (survey == 15) & (z < 0.4306)
    iPS1_highz = (survey == 15) & (z >= 0.4306)
    iDES_midz = (survey == 10) & (z < 0.4306)
    iDES_highz = (survey == 10) & (z >= 0.4306)
    
    mu_allz = x[0]
    sigint_csp,sigint_ps1,sigint_des = x[1],x[2],x[3]
    mass_step = x[4]
    
    # sigint split by sample, but distance split by redshift
    # each one with a low-mass and high-mass component
    loglike_csp = -np.sum(np.logaddexp(-(mu_i[iCSP]+mass_step-mu_allz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((1-0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.))),
                                       -(mu_i[iCSP]-mu_allz)**2./(2.0*(sigma_i[iCSP]**2.+sigint_csp**2.)) + \
                                       np.log((0.01*p_lm[iCSP])/(np.sqrt(2*np.pi)*np.sqrt(sigint_csp**2.+sigma_i[iCSP]**2.)))))
            
    loglike_ps1_midz = -np.sum(np.logaddexp(-(mu_i[iPS1_midz]+mass_step-mu_allz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((1-0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.))),
                                            -(mu_i[iPS1_midz]-mu_allz)**2./(2.0*(sigma_i[iPS1_midz]**2.+sigint_ps1**2.)) + \
                                            np.log((0.01*p_lm[iPS1_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_midz]**2.)))))
    
    loglike_ps1_highz = -np.sum(np.logaddexp(-(mu_i[iPS1_highz]+mass_step-mu_allz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((1-0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.))),
                                             -(mu_i[iPS1_highz]-mu_allz)**2./(2.0*(sigma_i[iPS1_highz]**2.+sigint_ps1**2.)) + \
                                             np.log((0.01*p_lm[iPS1_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_ps1**2.+sigma_i[iPS1_highz]**2.)))))
            
    loglike_des_midz = -np.sum(np.logaddexp(-(mu_i[iDES_midz]+mass_step-mu_allz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((1-0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.))),
                                            -(mu_i[iDES_midz]-mu_allz)**2./(2.0*(sigma_i[iDES_midz]**2.+sigint_des**2.)) + \
                                            np.log((0.01*p_lm[iDES_midz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_midz]**2.)))))
            
    loglike_des_highz = -np.sum(np.logaddexp(-(mu_i[iDES_highz]+mass_step-mu_allz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((1-0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.))),
                                             -(mu_i[iDES_highz]-mu_allz)**2./(2.0*(sigma_i[iDES_highz]**2.+sigint_des**2.)) + \
                                             np.log((0.01*p_lm[iDES_highz])/(np.sqrt(2*np.pi)*np.sqrt(sigint_des**2.+sigma_i[iDES_highz]**2.)))))

    return loglike_csp + loglike_ps1_midz + loglike_ps1_highz + loglike_des_midz + loglike_des_highz



def apply_all_cuts(fr,fropt,restrict_to_good_list=False):

    # AV
    #iGoodAV = np.zeros(len(fr.CID),dtype=bool)
    #for j,i in enumerate(fr.CID):
    #   if i in fropt.CID and fropt.AV[fropt.CID == i] < 0.3*fropt.RV[fropt.CID == i]:
    #       iGoodAV[j] = True

    # reasonable stretch
    #iGoodSt = np.zeros(len(fr.CID),dtype=bool)
    #for j,i in enumerate(fr.CID):
    #   if i in fropt.CID and fropt.STRETCH[fropt.CID == i] > 0.8 and fropt.STRETCH[fropt.CID == i] < 1.3:
    #       iGoodSt[j] = True

    #for k in fr.__dict__.keys():
    #   fr.__dict__[k] = fr.__dict__[k][iGoodAV & iGoodSt]

    #restrict_to_good_list = False
    if restrict_to_good_list:
        # then case-by-case removal of bad stuff
        # because some RAISIN things have bad photometry
        iGood = np.array([],dtype=int)
        for j,i in enumerate(fr.CID):
            if i in _goodcids: iGood = np.append(iGood,j)
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

    return fr


def main(boundary=10.0):

    fig = plt.figure()#constrained_layout=True)
    plt.subplots_adjust(top=0.95)
    gs = GridSpec(3, 3, figure=fig)
    gs.update(wspace=0.0, hspace=0.4)

    #plt.tick_params(left)
    
    axmain = fig.add_subplot(gs[1:, :])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])


    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,ax,title in zip(
            nirdatafitreslist,opticalnirdatafitreslist,[ax1,ax2,ax3],['CSP','PS1','DES']):

        fr = txtobj(frfile,fitresheader=True)
        fropt = txtobj(froptfile,fitresheader=True)
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2. + 0.02**2.)
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
        
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        ax.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
        ax.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

        # light shading for the dispersion
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

        ax.axvline(boundary,ls='--',lw=2,color='0.2')
        ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
                    yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
        sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
                        s=30,zorder=9,cmap='RdBu_r')
        ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        
        ax.set_ylim([-0.7,0.7])
        ax.set_xlim([7,13])
        ax.set_title(title)
        ax.set_xlabel('log(M/M$_{\odot}$)')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.xaxis.set_ticks([8,9,10,11,12])
        #import pdb; pdb.set_trace()
    ax1.set_ylabel('Hubble Resid (mag)')
    ax2.yaxis.set_ticklabels([])
    ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    ax3.yaxis.tick_right()
    ax3.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)


    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))

    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])

    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[1]
    resid_full[iHighz] -= md.x[2]
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.3f}$\pm${steperr:.3f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7,13])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    plt.savefig('figs/raisin_massstep.png',dpi=200)
    
        
    import pdb; pdb.set_trace()

def main_empiricalstretch(boundary=10.44):

    fig = plt.figure()
    plt.subplots_adjust(top=0.95)
    gs = GridSpec(3, 3, figure=fig)
    gs.update(wspace=0.0, hspace=0.4)

    
    axmain = fig.add_subplot(gs[1:, :])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    st = txtobj('st_av_corr_mags.txt')
    hmcorr = -0.487*st.STRETCH + 0.503
    hmcorrmod = (-0.487+0.142)*st.STRETCH + 0.503
    
    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,ax,title in zip(
            nirdatafitreslist,opticalnirdatafitreslist,[ax1,ax2,ax3],['CSP','PS1','DES']):

        fr = txtobj(frfile,fitresheader=True)
        fropt = txtobj(froptfile,fitresheader=True)
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2. + 0.02**2.)
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.

        # resid correction for stretch
        print(np.std(fr.resid))
        for i in range(len(fr.CID)):
            fr.resid[i] -= hmcorr[st.CID == fr.CID[i]][0]
        print(np.std(fr.resid))
#            import pdb; pdb.set_trace()
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        ax.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
        ax.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

        # light shading for the dispersion
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

        ax.axvline(boundary,ls='--',lw=2,color='0.2')
        ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
                    yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
        sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
                        s=30,zorder=9,cmap='RdBu_r')
        ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        
        ax.set_ylim([-0.7,0.7])
        ax.set_xlim([7,13])
        ax.set_title(title)
        ax.set_xlabel('log(M/M$_{\odot}$)')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.xaxis.set_ticks([8,9,10,11,12])
        #import pdb; pdb.set_trace()
    ax1.set_ylabel('Hubble Resid (mag)')
    ax2.yaxis.set_ticklabels([])
    ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    ax3.yaxis.tick_right()
    ax3.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)


    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))

    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])

    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[1]
    resid_full[iHighz] -= md.x[2]
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.3f}$\pm${steperr:.3f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7,13])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    plt.savefig('figs/raisin_massstep.png',dpi=200)
    
        
    import pdb; pdb.set_trace()

    
def main_stretch(boundary=10.0):

    fig = plt.figure()#constrained_layout=True)
    plt.subplots_adjust(top=0.95)
    gs = GridSpec(3, 3, figure=fig)
    gs.update(wspace=0.0, hspace=0.4)

    #plt.tick_params(left)
    
    axmain = fig.add_subplot(gs[1:, :])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])


    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,ax,title in zip(
            nirstretchdatafitreslist,opticalnirdatafitreslist,[ax1,ax2,ax3],['CSP','PS1 (RAISIN1)','DES (RAISIN2)']):

        fr = txtobj(frfile,fitresheader=True)
        fropt = txtobj(froptfile,fitresheader=True)
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

        # need to apply PV terms
        peczerr=0.00083
        zerr = peczerr*5.0/np.log(10)*(1.0+fr.zHD)/(fr.zHD*(1.0+fr.zHD/2.0))
        fr.DLMAGERR = np.sqrt(fr.DLMAGERR**2. + zerr**2. + 0.055**2.*fr.zHD**2.)
            
        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2. + 0.02**2.)
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
        print(np.std(fr.resid))
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        ax.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
        ax.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

        # light shading for the dispersion
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

        ax.axvline(boundary,ls='--',lw=2,color='0.2')
        ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
                    yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
        sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
                        s=30,zorder=9,cmap='RdBu_r')
        ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        
        ax.set_ylim([-0.8,0.8])
        ax.set_xlim([7,13])
        ax.set_title(title)
        ax.set_xlabel('log(M/M$_{\odot}$)')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.xaxis.set_ticks([8,9,10,11,12])
    ax1.set_ylabel('Hubble Resid (mag)')
    ax2.yaxis.set_ticklabels([])
    ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    ax3.yaxis.tick_right()
    ax3.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)


    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])
    import pdb; pdb.set_trace()
    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[1]
    resid_full[iHighz] -= md.x[2]
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.3f}$\pm${steperr:.3f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.8,0.8])
    axmain.set_xlim([7,13])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    plt.savefig('figs/raisin_massstep.png',dpi=200)
    
        
    import pdb; pdb.set_trace()

def main_cosmo(boundary=10.44,w=-1.01,biascor=True):

    fig = plt.figure()#constrained_layout=True)
    plt.subplots_adjust(top=0.95)
    gs = GridSpec(3, 3, figure=fig)
    gs.update(wspace=0.0, hspace=0.4)

    #plt.tick_params(left)
    
    axmain = fig.add_subplot(gs[1:, :])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])


    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,ax,title in zip(
            nirstretchdatafitreslist,opticalnirdatafitreslist,[ax1,ax2,ax3],['CSP','PS1 (RAISIN1)','DES (RAISIN2)']):

        fr = txtobj(frfile,fitresheader=True)
        if biascor:
            z,bc,bce = np.loadtxt(f'{title.split()[0]}_biascor_stretch.txt',unpack=True)
            biascor_interp = np.interp(fr.zHD,z[bc == bc],bc[bc == bc])
#            import pdb; pdb.set_trace()
            fr.DLMAG -= biascor_interp
            
        fropt = txtobj(froptfile,fitresheader=True)
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD,w0=w)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2. + 0.02**2.)
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        ax.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
        ax.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

        # light shading for the dispersion
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

        ax.axvline(boundary,ls='--',lw=2,color='0.2')
        ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
                    yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
        sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
                        s=30,zorder=9,cmap='RdBu_r')
        ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        
        ax.set_ylim([-0.8,0.8])
        ax.set_xlim([7,13])
        ax.set_title(title)
        ax.set_xlabel('log(M/M$_{\odot}$)')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.xaxis.set_ticks([8,9,10,11,12])
    ax1.set_ylabel('Hubble Resid (mag)')
    ax2.yaxis.set_ticklabels([])
    ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    ax3.yaxis.tick_right()
    ax3.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)


    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

    md = minimize(neglnlikefunc_cosmo,(0,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #import pdb; pdb.set_trace()
    #md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
    #              args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #import pdb; pdb.set_trace()
    step,steperr = md.x[4],np.sqrt(md.hess_inv[4,4])
    #md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
    #              args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])

    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[0]
    resid_full[iHighz] -= md.x[0]
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.3f}$\pm${steperr:.3f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7,13])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    #plt.savefig('figs/raisin_massstep.png',dpi=200)
    
        
    import pdb; pdb.set_trace()

def main_highz(boundary=10.44,w=-1.00,biascor=True):

    fig = plt.figure()#constrained_layout=True)
    plt.subplots_adjust(top=0.95)
    gs = GridSpec(3, 3, figure=fig)
    gs.update(wspace=0.0, hspace=0.4)

    #plt.tick_params(left)
    
    axmain = fig.add_subplot(gs[1:, :])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])


    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,ax,title in zip(
            nirstretchdatafitreslist[1:],opticalnirdatafitreslist[1:],[ax1,ax2,ax3][1:],['CSP','PS1 (RAISIN1)','DES (RAISIN2)'][1:]):

        fr = txtobj(frfile,fitresheader=True)
        if biascor:
            z,bc,bce = np.loadtxt(f'{title.split()[0]}_biascor_stretch.txt',unpack=True)
            biascor_interp = np.interp(fr.zHD,z[bc == bc],bc[bc == bc])
            fr.DLMAG -= biascor_interp
            
        fropt = txtobj(froptfile,fitresheader=True)
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zHD,w0=w)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]

        fr.HOST_LOGMASS_ERR = np.sqrt(fr.HOST_LOGMASS_ERR**2. + 0.02**2.)
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        ax.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
        ax.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

        # light shading for the dispersion
        ax.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
        ax.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

        ax.axvline(boundary,ls='--',lw=2,color='0.2')
        ax.errorbar(fr.HOST_LOGMASS,fr.resid,xerr=fr.HOST_LOGMASS_ERR,
                    yerr=fr.DLMAGERR,color='0.6',fmt='',ls='None')
        sc = ax.scatter(fr.HOST_LOGMASS,fr.resid,c=100-fr.p_hm,
                        s=30,zorder=9,cmap='RdBu_r')
        ax.text(0.05,0.95,f"$\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        
        ax.set_ylim([-0.8,0.8])
        ax.set_xlim([7,13])
        ax.set_title(title)
        ax.set_xlabel('log(M/M$_{\odot}$)')
        ax.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
        ax.xaxis.set_ticks([8,9,10,11,12])
    ax1.set_ylabel('Hubble Resid (mag)')
    ax2.yaxis.set_ticklabels([])
    ax3.tick_params(top="on",bottom="on",left="off",right="on",direction="inout",length=8, width=1.5)
    ax3.yaxis.tick_right()
    ax3.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)


    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

    md = minimize(neglnlikefunc_cosmo,(0,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #import pdb; pdb.set_trace()
    #md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
    #              args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #import pdb; pdb.set_trace()
    step,steperr = md.x[4],np.sqrt(md.hess_inv[4,4])
    #md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
    #              args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])

    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[0]
    resid_full[iHighz] -= md.x[0]
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.6,zorder=100)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='salmon',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.3f}$\pm${steperr:.3f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7,13])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    #plt.savefig('figs/raisin_massstep.png',dpi=200)
    
        
    import pdb; pdb.set_trace()

    
    
def main_opt():
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    main_snoopy_opt(axmain=ax1)
    main_salt2(axmain=ax2)
    import pdb; pdb.set_trace()
    
def main_snoopy_opt(boundary=10,axmain=None):

    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    av_full,stretch_full = np.array([]),np.array([])
    
    for frfile in opticalnirdatafitreslist:

        fr = txtobj(frfile,fitresheader=True)
        fr = apply_all_cuts(fr,None,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zCMB)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]
        
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.HOST_LOGMASS_ERR[i] = np.sqrt(fr.HOST_LOGMASS_ERR[i]**2. + 0.01**2.)
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
            #if fr.p_hm[i] != fr.p_hm[i]: import pdb; pdb.set_trace()
        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)
        av_full = np.append(av_full,fr.AV)
        stretch_full = np.append(stretch_full,fr.STRETCH)

    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))

        
    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    #import pdb; pdb.set_trace()
    
    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)
    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7.9,12.1])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_title('SNooPy Optical',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])

    
def main_salt2(boundary=10,axmain=None):

    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

    for frfile in salt2fitreslist:

        fr = txtobj(frfile,fitresheader=True)
        fr = apply_all_cuts(fr,None,restrict_to_good_list=True)
        fr = getmu.getmu(fr,sigint=0.0)
        #fr.mu = fr.mB+0.14*fr.x1-3.1*fr.c+19.36
        fr.resid = fr.mu - cosmo.mu(fr.zCMB)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]
        
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                np.sqrt(float(fr.HOST_LOGMASS_ERR[i])**2.+0.02**2.))*100.
        #import pdb; pdb.set_trace()        
        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.muerr)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)

    
    md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                  args=(mp_full,resid_full,residerr_full,None))

    resid_iaa,resid_iae = md.x[0],md.x[1]
    scat_iaa,scat_iae = md.x[2],md.x[3]
    residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
    covar = np.sqrt(np.abs(md.hess_inv[1,0]))
    #step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)
    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))
    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])
    
    axmain.plot(np.arange(boundary-10,boundary,0.001),
                [resid_iae]*len(np.arange(boundary-10,boundary,0.001)),
                lw=2,color='0.2')
    axmain.plot(np.arange(boundary,boundary+10,0.001),
                [resid_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                boundary,10,lw=2,color='0.2')
        
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+residerr_iae]*len(np.arange(boundary-10,boundary,0.001)),color='blue')
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+residerr_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',alpha=0.4)

    # light shading for the dispersion
    axmain.fill_between(np.arange(boundary-10,boundary,0.001),
                        [resid_iae-scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        [resid_iae+scat_iae]*len(np.arange(boundary-10,boundary,0.001)),
                        color='lightblue',zorder=1,alpha=0.4)
    axmain.fill_between(np.arange(boundary,boundary+10,0.001),
                        [resid_iaa-scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        [resid_iaa+scat_iaa]*len(np.arange(boundary,boundary+10,0.001)),
                        color='red',zorder=1,alpha=0.6)

    axmain.axvline(boundary,ls='--',lw=2,color='0.2')
    axmain.errorbar(mass_full,resid_full,xerr=masserr_full,
                    yerr=residerr_full,color='0.6',fmt='',ls='None')
    sc = axmain.scatter(mass_full,resid_full,c=100-mp_full,
                    s=30,zorder=9,cmap='RdBu_r')
    axmain.text(0.02,0.95,f"Total $\Delta_M$ = {step:.2f}$\pm${steperr:.2f}",
                va='top',ha='left',transform=axmain.transAxes,
                bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.7},
                zorder=100,fontsize=15)

    axmain.set_ylim([-0.5,0.5])
    axmain.set_xlim([7.9,12.1])
    axmain.set_ylabel('Hubble Resid (mag)',fontsize=15)
    axmain.set_title('SALT2 Optical',fontsize=15)
    axmain.set_xlabel('log(M/M$_{\odot}$)',fontsize=15)
    axmain.tick_params(top="on",bottom="on",left="on",right="on",direction="inout",length=8, width=1.5)
    axmain.xaxis.set_ticks([8,9,10,11,12])
    import pdb; pdb.set_trace()
    
def add_masses():

    hosts = at.Table.read('hosts/host_properties_raisin.txt',format='ascii')
    cspfiles = np.loadtxt('data/Photometry/CSPDR3_RAISIN/CSPDR3_RAISIN.LIST',unpack=True,dtype=str)
    cspfiles = [f"data/Photometry/CSPDR3_RAISIN/{c}" for c in cspfiles]
    ps1files = np.loadtxt('data/Photometry/PS1_RAISIN/PS1_RAISIN.LIST',unpack=True,dtype=str)
    ps1files = [f"data/Photometry/PS1_RAISIN/{p}" for p in ps1files]
    desfiles = np.loadtxt('data/Photometry/DES_RAISIN/DES_RAISIN.LIST',unpack=True,dtype=str)
    desfiles = [f"data/Photometry/DES_RAISIN/{d}" for d in desfiles]

    for files in [cspfiles,ps1files,desfiles]:
        for f in files:
            hostmass = None
            with open(f) as fin, open(f.replace('.DAT','_LOGMASS.DAT').replace('.dat','_LOGMASS.DAT'),'w') as fout:
                for line in fin:
                    line = line.replace('\n','')
                    if line.startswith('SNID:'):
                        snid = line.split()[1]
                        if snid in hosts['snid']:
                            hostmass = hosts['logmass'][hosts['snid'] == snid][0]
                            if hostmass > 10: hostmasserr =  hosts['logmass'][hosts['snid'] == snid][0] -  hosts['logmass_low'][hosts['snid'] == snid][0]
                            else: hostmasserr =  hosts['logmass_high'][hosts['snid'] == snid][0] -  hosts['logmass'][hosts['snid'] == snid][0]
                        else: hostmass = -99.0
                        print(line,file=fout)
                    elif line.startswith('HOSTGAL_LOGMASS:') and hostmass > 0:
                        print(f'HOSTGAL_LOGMASS: {hostmass:.3f} +- {hostmasserr:.3f}',file=fout)
                    else:
                        print(line,file=fout)
            os.system(f"mv {f.replace('.DAT','_LOGMASS.DAT').replace('.dat','_LOGMASS.DAT')} {f}")

def shapecolor(boundary=10):
    """
    figure out how much the shape- and color-corrections 
    would be expected to change the measured mass step(s)
    """

    mp_full,mass_full,masserr_full,resid_full,residerr_full,survey_full,z_full = \
        np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
    x1,c,sbv,av = np.array([]),np.array([]),np.array([]),np.array([])
    
    for frfile,froptfile,frsalt2file in zip(
            nirdatafitreslist,opticalnirdatafitreslist,salt2fitreslist):

        fr = txtobj(frfile,fitresheader=True)
        fropt = txtobj(froptfile,fitresheader=True)
        frsalt2 = txtobj(frsalt2file,fitresheader=True)
        
        fr = apply_all_cuts(fr,fropt,restrict_to_good_list=True)
        fr.resid = fr.DLMAG - cosmo.mu(fr.zCMB)
        iGood = np.where(fr.HOST_LOGMASS_ERR < 5)[0]
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]
        iGoodOpt = np.array([],dtype=int)
        for i in fr.CID:
            if i in fropt.CID:
                iGoodOpt = np.append(iGoodOpt,np.where(fropt.CID == i)[0][0])            
        for k in fropt.__dict__.keys():
            fropt.__dict__[k] = fropt.__dict__[k][iGoodOpt]

        iGoodSALT2 = np.array([],dtype=int)
        for i in fr.CID:
            if i in frsalt2.CID:
                iGoodSALT2 = np.append(iGoodSALT2,np.where(frsalt2.CID == i)[0][0])            
        for k in frsalt2.__dict__.keys():
            frsalt2.__dict__[k] = frsalt2.__dict__[k][iGoodSALT2]

        iGoodOpt = np.array([],dtype=int)
        for i in frsalt2.CID:
            if i in fropt.CID:
                iGoodOpt = np.append(iGoodOpt,np.where(fropt.CID == i)[0][0])            
        for k in fropt.__dict__.keys():
            fropt.__dict__[k] = fropt.__dict__[k][iGoodOpt]
        iGood = np.array([],dtype=int)
        for i in frsalt2.CID:
            if i in fr.CID:
                iGood = np.append(iGood,np.where(fr.CID == i)[0][0])            
        for k in fr.__dict__.keys():
            fr.__dict__[k] = fr.__dict__[k][iGood]
            
            
        fr.p_hm = np.zeros(len(fr.CID))
        for i in range(len(fr.CID)):
            fr.p_hm[i] = scipy.stats.norm.cdf(
                boundary,float(fr.HOST_LOGMASS[i]),
                float(fr.HOST_LOGMASS_ERR[i]))*100.
        
        md = minimize(lnlikefunc,(0.0,0.0,0.1,0.1),
                      args=(fr.p_hm,fr.resid,fr.DLMAGERR,None))

        resid_iaa,resid_iae = md.x[0],md.x[1]
        scat_iaa,scat_iae = md.x[2],md.x[3]
        residerr_iaa,residerr_iae = np.sqrt(md.hess_inv[0,0]),np.sqrt(md.hess_inv[1,1])
        covar = np.sqrt(np.abs(md.hess_inv[1,0]))
        step,steperr = resid_iae-resid_iaa,np.sqrt(residerr_iae**2.+residerr_iaa**2.-2*covar**2.)

        mp_full = np.append(mp_full,fr.p_hm)
        mass_full = np.append(mass_full,fr.HOST_LOGMASS)
        masserr_full = np.append(masserr_full,fr.HOST_LOGMASS_ERR)
        resid_full = np.append(resid_full,fr.resid)
        residerr_full = np.append(residerr_full,fr.DLMAGERR)
        survey_full = np.append(survey_full,fr.IDSURVEY)
        z_full = np.append(z_full,fr.zHD)

        sbv = np.append(sbv,fropt.STRETCH)
        av = np.append(av,fropt.AV)
        x1 = np.append(x1,frsalt2.x1)
        c = np.append(c,frsalt2.c)

    md = minimize(neglnlikefunc,(0,0.01,0.02,0.09,0.1,0.11,0.1),
                  args=(mp_full,resid_full,residerr_full,None,survey_full,z_full))

    step,steperr = md.x[6],np.sqrt(md.hess_inv[6,6])

    iCSP = survey_full == 5
    iMidz = ((survey_full == 15) & (z_full < 0.4306)) | ((survey_full == 10) & (z_full < 0.4306))
    iHighz = ((survey_full == 15) & (z_full >= 0.4306)) | ((survey_full == 10) & (z_full >= 0.4306))
    resid_full[iCSP] -= md.x[0]
    resid_full[iMidz] -= md.x[1]
    resid_full[iHighz] -= md.x[2]
    # high-mass: s_BV = 1.00923, A_V = 0.0546
    # low-mass: s_BV = 1.113075, A_V = -0.0227
    # stretch: -0.048 in Y, 0.043 in i, 0.093 in J, 0.06 in H
    # color: 0.0546+0.0227 = 0.077/3.1 = 0.025 = E(B-V)
    # iYJH: 0.0463, 0.0273, 0.0201, 0.0129
    # without shape/color, mass step is zero
    # low-mass has higher stretch, correcting for this makes low-mass fainter by ~0.05 mag (Dm~ 0.05 +/- 0.05)
    # high-mass has redder color, correcting for this makes high-mass brighter by 0.03 (Dm ~ 0.08 +/- 0.05)
    import pdb; pdb.set_trace()


    
if __name__ == "__main__":
    #add_hosts()
    #main_salt2()
    #add_masses()
    #main()

    #main_empiricalstretch()

    main_stretch()
    #main_opt()
    #shapecolor()
    #checknewmasses()
    #main_cosmo()
    #main_highz()
