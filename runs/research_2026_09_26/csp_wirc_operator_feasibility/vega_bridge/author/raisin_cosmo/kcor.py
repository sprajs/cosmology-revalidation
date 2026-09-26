#!/usr/bin/env python
# D. Jones - 11/21/19
"""see how different SN Ia models affect distance measurements via k-corrections"""

import os
import numpy as np
from txtobj import txtobj
from util.common import *
import pylab as plt
plt.ion()
import sncosmo
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from astropy.io import fits

def savitzky_golay(y, window_size=5, order=3, deriv=0):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except: # ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')


class interp_kcor:

    def __init__(self):
        self.nircompfiles = ['kcor/snsed/NIR_composite_N=3_Nspec=6_phase=m9.92.txt',
                             'kcor/snsed/NIR_composite_N=5_Nspec=16_phase=m3.01.txt',
                             'kcor/snsed/NIR_composite_N=4_Nspec=8_phase=p6.93.txt',
                             'kcor/snsed/NIR_composite_N=2_Nspec=3_phase=p15.22.txt',
                             'kcor/snsed/NIR_composite_N=2_Nspec=4_phase=p38.73.txt']
        self.nirphases = [-9.92,-3.01,6.93,15.22,38.73]
        self.niroutphases = np.concatenate((np.linspace(self.nirphases[0]+0.01,self.nirphases[1],20),
                                            np.linspace(self.nirphases[1]+0.01,self.nirphases[2],20),
                                            np.linspace(self.nirphases[2]+0.01,self.nirphases[3],20),
                                            np.linspace(self.nirphases[3]+0.01,self.nirphases[4],20)))

        self.fecompfiles = ['kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110826_SpeX.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110828_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110831_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110903_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110907_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110910_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110913_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110918_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110922_GNIRS.fits',
                            'kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110927_GNIRS.fits']
        self.fephases = np.array([55799.2713604993,55801.2509995606,55804.2339751990,55807.2384183044,
                                  55811.1063987008,55814.2100968323,55817.2389874040,55822.1314159553,
                                  55826.2153242983,55831.2127987367])-55814.895
        self.feoutphases = np.concatenate((np.linspace(self.fephases[0]+0.01,self.fephases[1],20),
                                           np.linspace(self.fephases[1]+0.01,self.fephases[2],20),
                                           np.linspace(self.fephases[2]+0.01,self.fephases[3],20),
                                           np.linspace(self.fephases[3]+0.01,self.fephases[4],20),
                                           np.linspace(self.fephases[4]+0.01,self.fephases[5],20),
                                           np.linspace(self.fephases[5]+0.01,self.fephases[6],20),
                                           np.linspace(self.fephases[6]+0.01,self.fephases[7],20),
                                           np.linspace(self.fephases[7]+0.01,self.fephases[8],20),
                                           np.linspace(self.fephases[8]+0.01,self.fephases[9],20)))
        
        self.hphase,self.hwave,self.hflux = np.loadtxt('kcor/snsed/Hsiao07.dat',unpack=True)
        self.hflux = self.hflux.reshape([len(np.unique(self.hphase)),len(np.unique(self.hwave))])
        self.hphase = np.unique(self.hphase)
        self.hwave = np.unique(self.hwave)

    def main(self):

        flux1out = np.zeros([50,len(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])]),len(self.hwave)])
        flux1out_tmp = np.zeros([50,len(self.niroutphases),len(self.hwave)])
        
        for h in range(50):
            print(h)
            flux1in = np.zeros([len(self.nirphases),len(self.hwave)])
            for n,p,i in zip(self.nircompfiles,self.nirphases,range(len(self.nirphases))):
                wave,flux1 = \
                    np.loadtxt(n,unpack=True,usecols=[0,h+1])
                flux1 = np.interp(self.hwave,wave,flux1)
                flux1in[i,:] = flux1
            for i,w in enumerate(self.hwave):
                for p1,p2,j in zip(self.nirphases[:-1],self.nirphases[1:],range(len(self.nirphases))):
                    # hsiao from start to end points
                    phase = np.linspace(p1,p2,20)
                    hflux = np.interp(phase,self.hphase,self.hflux[:,i])
                    hflux *= ((flux1in[j+1,i] + flux1in[j,i])/2.)/np.median(hflux)*2

                    slope_old = (hflux[phase == p2]-hflux[phase == p1])/(p2-p1)
                    slope_new = (flux1in[j+1,i] - flux1in[j,i])/(p2-p1)
                    phi = np.arctan(-(slope_old-slope_new)/(1+slope_old*slope_new))
                    phase_new = phase*np.cos(phi) - hflux*np.sin(phi)
                    hflux_new = phase*np.sin(phi) + hflux*np.cos(phi)
                    hflux_new2 = hflux_new + (flux1in[j,i]-hflux_new[phase == p1])

                    
                    ifluxout = np.where((self.niroutphases > p1) & (self.niroutphases <= p2))[0]
                    flux1out_tmp[h,ifluxout,i] = hflux_new2

                        
                    
            # HEAVILY smooth telluric regions
            for j in range(np.shape(flux1out_tmp)[1]):
                flux1out_tmp[h,j,(self.hwave > 13200) & (self.hwave < 14200)] = \
                    savitzky_golay(flux1out_tmp[h,j,(self.hwave > 13200) & (self.hwave < 14200)],window_size=91,order=1)
                flux1out_tmp[h,j,(self.hwave > 18000) & (self.hwave < 19000)] = \
                    savitzky_golay(flux1out_tmp[h,j,(self.hwave > 18000) & (self.hwave < 19000)],window_size=91,order=1)
                
            int1d = interp1d(self.niroutphases,flux1out_tmp,axis=1)
            
            flux1out[h,:,:] = int1d(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])])[h,:,:]

            # kcors complain if optical/UV is way off, so just use Hsiao for this
            phaseidx = np.where((self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1]))[0]
            for j in range(np.shape(flux1out)[1]):
                iJoin = (self.hwave > 8200) & (self.hwave < 8300)
                scale = np.median(flux1out[h,j,iJoin])/np.median(self.hflux[phaseidx[j],iJoin])
                flux1out[h,j,self.hwave < 8250] = self.hflux[phaseidx[j],self.hwave < 8250]*scale

            #import pdb; pdb.set_trace()
        
        with open('kcor/snsed/hsiao_bootstrapped.txt','w') as fout:
            for i,p in enumerate(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])]):
                for j,w in enumerate(self.hwave):
                    print(f"{p:.1f} {w:.1f} {np.nanmean(flux1out[:,i,j]):.5f} {np.nanstd(flux1out[:,i,j]):.5f}",file=fout)
        with open('kcor/snsed/hsiao_bootstrapped_sys.txt','w') as fout:
            for i,p in enumerate(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])]):
                for j,w in enumerate(self.hwave):
                    print(f"{p:.1f} {w:.1f} {np.nanmean(flux1out[:,i,j])+np.nanstd(flux1out[:,i,j]):.5f}",file=fout)
        with open('kcor/snsed/hsiao_errors.txt','w') as fout:
            for i,p in enumerate(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])]):
                for j,w in enumerate(self.hwave):
                    print(f"{p:.1f} {w:.1f} {self.hflux[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])][i,j]} {np.nanstd(flux1out[:,i,j])*self.hflux[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])][i,j]/np.nanmean(flux1out[:,i,j])}",file=fout)
        with open('kcor/snsed/hsiao_final_sys.txt','w') as fout:
            for i,p in enumerate(self.hphase[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])]):
                for j,w in enumerate(self.hwave):
                    print(f"{p:.1f} {w:.1f} {self.hflux[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])][i,j]+  np.nanstd(flux1out[:,i,j])*self.hflux[(self.hphase > self.niroutphases[0]) & (self.hphase < self.niroutphases[-1])][i,j]/np.nanmean(flux1out[:,i,j])}",file=fout)

    def fetmpl(self):

        flux1out = np.zeros([len(self.hphase[(self.hphase > self.feoutphases[0]) & (self.hphase < self.feoutphases[-1])]),len(self.hwave)])
        flux1out_tmp = np.zeros([len(self.feoutphases),len(self.hwave)])
        
        flux1in = np.zeros([len(self.fephases),len(self.hwave)])
        for n,p,i in zip(self.fecompfiles,self.fephases,range(len(self.fephases))):
            a = fits.open(n)
            wave,flux1 = a[0].data[0]*10000,a[0].data[1]
            flux1 = np.interp(self.hwave,wave,flux1)
            flux1in[i,:] = flux1
        for i,w in enumerate(self.hwave):
            for p1,p2,j in zip(self.fephases[:-1],self.fephases[1:],range(len(self.fephases))):
                # hsiao from start to end points
                phase = np.linspace(p1,p2,20)
                hflux = np.interp(phase,self.hphase,self.hflux[:,i])
                hflux *= ((flux1in[j+1,i] + flux1in[j,i])/2.)/np.median(hflux)*2

                slope_old = (hflux[phase == p2]-hflux[phase == p1])/(p2-p1)
                slope_new = (flux1in[j+1,i] - flux1in[j,i])/(p2-p1)
                phi = np.arctan(-(slope_old-slope_new)/(1+slope_old*slope_new))
                phase_new = phase*np.cos(phi) - hflux*np.sin(phi)
                hflux_new = phase*np.sin(phi) + hflux*np.cos(phi)
                hflux_new2 = hflux_new + (flux1in[j,i]-hflux_new[phase == p1])



                ifluxout = np.where((self.feoutphases > p1) & (self.feoutphases <= p2))[0]
                flux1out_tmp[ifluxout,i] = hflux_new2


        int1d = interp1d(self.feoutphases,flux1out_tmp,axis=0)

        flux1out[:,:] = int1d(self.hphase[(self.hphase > self.feoutphases[0]) & (self.hphase < self.feoutphases[-1])])[:,:]

        # kcors complain if optical/UV is way off, so just use Hsiao for this
        phaseidx = np.where((self.hphase > self.feoutphases[0]) & (self.hphase < self.feoutphases[-1]))[0]
        for j in range(np.shape(flux1out)[0]):
            iJoin = (self.hwave > 8200) & (self.hwave < 8300)
            scale = np.median(flux1out[j,iJoin])/np.median(self.hflux[phaseidx[j],iJoin])
            flux1out[j,self.hwave < 8250] = self.hflux[phaseidx[j],self.hwave < 8250]*scale

        #import pdb; pdb.set_trace()
        with open('kcor/snsed/2011fe_sequence.txt','w') as fout:
            for i,p in enumerate(self.hphase[(self.hphase > self.feoutphases[0]) & (self.hphase < self.feoutphases[-1])]):
                for j,w in enumerate(self.hwave):
                    print(f"{p:.1f} {w:.1f} {flux1out[i,j]}",file=fout)

        #array([-0.17221, -0.10029,  0.05265, -0.06227, -0.10707, -0.14443,
        #-0.11017, -0.20619, -0.11256, -0.06934, -0.08742, -0.10396,
        #-0.11662, -0.11504, -0.03108, -0.1229 ,  0.09538, -0.0183 ,
        # 0.01113, -0.04181,  0.00148])
                    
                    
    def plot(self):

        # 2011fe spectra:
        # SN2011fe_20110907_GNIRS.fits: -4 days
        # SN2011fe_20110918_GNIRS.fits: +7 dats
        # SN2011fe_20110927_GNIRS.fits: +16 days
        a = fits.open('kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110907_GNIRS.fits')
        fewave1,feflux1 = a[0].data[0]*10000,a[0].data[1]
        a = fits.open('kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110918_GNIRS.fits')
        fewave2,feflux2 = a[0].data[0]*10000,a[0].data[1]
        a = fits.open('kcor/snsed/SN2011fe_NIR_spec/SN2011fe_20110927_GNIRS.fits')
        fewave3,feflux3 = a[0].data[0]*10000,a[0].data[1]
        
        ax = plt.axes()
        phase,wave,flux,fluxerr = np.loadtxt('kcor/snsed/hsiao_bootstrapped.txt',unpack=True)
        hphase,hwave,hflux,hfluxerr = np.loadtxt('kcor/snsed/hsiao_errors.txt',unpack=True)
        phase = np.unique(phase)
        wave = np.unique(wave)
        flux = flux.reshape([len(phase),len(wave)])
        fluxerr = fluxerr.reshape([len(phase),len(wave)])

        hphase = np.unique(hphase)
        hwave = np.unique(hwave)
        hflux = hflux.reshape([len(hphase),len(hwave)])
        hfluxerr = hfluxerr.reshape([len(hphase),len(hwave)])

        
        for mphase,offset in zip([-5,5,15,25],[0,1,2,3]):

            if offset == 0:
                fewave,feflux = fewave1,feflux1
            elif offset == 1:
                fewave,feflux = fewave2,feflux2
            elif offset == 2:
                fewave,feflux = fewave3,feflux3
            
            sort_flux = np.sort(flux[phase == mphase][0][wave > 8000])
            n_pix = len(wave[wave > 8000])
            minval = sort_flux[round(n_pix*0.05)]
            maxval = sort_flux[round(n_pix*0.95)]
            minval = minval*0.5
            maxval = maxval*1.1
            scale = maxval - minval

            sort_hflux = np.sort(self.hflux[self.hphase == mphase,:][0][self.hwave > 8000])
            n_pix = len(self.hwave[self.hwave > 8000])
            hminval = sort_hflux[round(n_pix*0.05)]
            hmaxval = sort_hflux[round(n_pix*0.95)]
            hminval = hminval*0.5
            hmaxval = hmaxval*1.1
            hscale = hmaxval - hminval

            sort_feflux = np.sort(feflux[fewave > 8000])
            n_pix = len(fewave[fewave > 8000])
            feminval = sort_feflux[round(n_pix*0.05)]
            femaxval = sort_feflux[round(n_pix*0.95)]
            feminval = feminval*0.5
            femaxval = femaxval*1.1
            fescale = femaxval - feminval

            
            if offset == 0:
                ax.plot(wave, (flux[phase == mphase,:][0]-minval)/scale + offset,
                        color='k',label='bootstrapped composite spectra')
                ax.plot(self.hwave,(self.hflux[self.hphase == mphase,:][0]-hminval)/hscale + offset,color='r',label='Hsiao+07')
                ax.plot(fewave,(feflux-feminval)/fescale + offset,color='C0',label='2011fe (Hsiao+13)')
            else:
                ax.plot(wave, (flux[phase == mphase,:][0]-minval)/scale + offset,
                        color='k')
                ax.plot(self.hwave,(self.hflux[self.hphase == mphase,:][0]-hminval)/hscale + offset,color='r')
                if offset != 3: ax.plot(fewave,(feflux-feminval)/fescale + offset,color='C0')
            ax.fill_between(wave, (flux[phase == mphase,:][0]-fluxerr[phase == mphase,:][0]-minval)/scale + offset,
                            (flux[phase == mphase,:][0]+fluxerr[phase == mphase,:][0]-minval)/scale+offset,
                            color='k',alpha=0.2)
            ax.fill_between(hwave, (hflux[hphase == mphase,:][0]-hfluxerr[hphase == mphase,:][0]-hminval)/hscale + offset,
                            (hflux[phase == mphase,:][0]+hfluxerr[phase == mphase,:][0]-hminval)/hscale+offset,
                            color='r',alpha=0.2)
            #import pdb; pdb.set_trace()
            ax.text(0.99,offset/5+0.1,f"phase = ${mphase:.0f}$",ha='right',transform=ax.transAxes)
        ax.set_xlim([8000,20000])
        ax.set_ylim([0,5])
        ax.set_xlabel(r'Wavelength (${\mathrm{\AA}}$)',fontsize=15)
        ax.set_ylabel('Flux',fontsize=15)
        ax.legend()
        import pdb; pdb.set_trace()
        #mphase = 10
        #plt.plot(wave[phase == mphase],flux[phase == mphase]*np.median(hflux[(hphase == mphase) & (hwave > 1000) & (hwave < 12000)])/np.median(flux[(phase == mphase) & (wave > 10000) & (wave < 12000)]))
        #plt.plot(hwave[hphase == mphase],hflux[hphase == mphase])

    def plot_dist(self):

        from scipy.stats import binned_statistic
        from txtobj import txtobj
        
        ax = plt.axes()
        frfiles = ['output/fit_nir/CSP_RAISIN.FITRES.TEXT',
                   'output/fit_nir/PS1_RAISIN.FITRES.TEXT',
                   'output/fit_nir/DES_RAISIN.FITRES.TEXT']
        frfiles_h07 = ['output/fit_nir/CSP_RAISIN_H07.FITRES.TEXT',
                       'output/fit_nir/PS1_RAISIN_H07.FITRES.TEXT',
                       'output/fit_nir/DES_RAISIN_H07.FITRES.TEXT']
        goodcids = ['output/goodcids/CSP_GOODCIDS_LATEST.LIST',
                    'output/goodcids/PS1_GOODCIDS_LATEST.LIST',
                    'output/goodcids/DES_GOODCIDS_LATEST.LIST']
        
        dlmag_base,dlmag_base_err,dlmag_h07,dlmag_h07_err,z = \
            np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
        for f,fh,gc in zip(frfiles,frfiles_h07,goodcids):
            fr = txtobj(f,fitresheader=True)
            frh = txtobj(fh,fitresheader=True)
            goodcidlist = np.loadtxt(gc,unpack=True,dtype=str)
            for j,i in enumerate(fr.CID):
                if i in goodcidlist:
                    dlmag_base = np.append(dlmag_base,fr.DLMAG[j])
                    dlmag_h07 = np.append(dlmag_h07,frh.DLMAG[frh.CID == i][0])
                    dlmag_base_err = np.append(dlmag_base_err,fr.DLMAGERR[j])
                    dlmag_h07_err = np.append(dlmag_h07_err,frh.DLMAGERR[frh.CID == i][0])

                    z = np.append(z,fr.zHD[j])

        delmag = dlmag_base - dlmag_h07
        delmagerr = np.sqrt(dlmag_base_err**2. + dlmag_h07_err**2.)
        
        def weighted_avg(idx):
            if len(idx) < 2: return np.nan
            average = np.average(delmag[idx], weights=1/delmagerr[idx]**2.)
            return average
        
        def weighted_avg_err(idx):
            if len(idx) < 2: return np.nan
            average = np.average(delmag[idx], weights=1/delmagerr[idx]**2.)
            variance = np.average((delmag[idx]-average)**2, weights=1/delmagerr[idx]**2.)  # Fast and numerically precise
            return np.sqrt(variance/len(idx))
        
                    
        zbins = np.linspace(0,0.7,15)
        dmu_binned = binned_statistic(z,range(len(dlmag_base)),bins=zbins,statistic=weighted_avg).statistic
        dmuerr_binned = binned_statistic(z,range(len(dlmag_base)),bins=zbins,statistic=weighted_avg_err).statistic
        ax.errorbar((zbins[1:]+zbins[:-1])/2.,dmu_binned,yerr=dmuerr_binned,fmt='o-',label='Composite $-$ Hsiao+07')
        ax.set_xlabel('$z_{CMB}$',fontsize=15)
        ax.set_ylabel('\Delta \mu',fontsize=15)
        ax.axhline(0,color='k',lw=2)
        import pdb; pdb.set_trace()
        
class kcor:

    def __init__(self):

        self.snsed = ['kcor/snsed/Hsiao07.dat',
                      'kcor/snsed/sn91t_flux.v1.1.dat',
                      #'kcor/snsed/salt2_m0_P18.dat',
                      'kcor/snsed/snflux_1a_Nugent2002.dat']

        self.kcor_input_ps1 = 'kcor/kcor_PS1MD_NIR.input'
        self.kcor_input_des = 'kcor/kcor_DES_NIR.input'
        self.kcor_input_hst = ''

        self.nml_ps1 = 'fit/PS1_RAISIN.nml'
        self.nml_des = 'fit/DES_RAISIN.nml'
        self.nml_hst = ''
        
    def mk_kcor(self):
        for ia_template in self.snsed[1:2]:
            ia_template_name = ia_template.split('/')[-1].split('.')[0]

            for kcor in [self.kcor_input_ps1,self.kcor_input_des]:
                with open(kcor,'r') as fin:
                    lines = fin.readlines()
                new_kcor = kcor.replace('.input','_%s.input'%ia_template_name)
                with open(new_kcor,'w') as fout:
                    for line in lines:
                        if not line.startswith('SN_SED') and not line.startswith('OUTFILE'):
                            print(line.replace('\n',''),file=fout)
                        elif line.startswith('SN_SED'):
                            print('SN_SED: %s'%ia_template,file=fout)
                        elif line.startswith('OUTFILE'):
                            print('OUTFILE: %s'%kcor.replace('.input','_%s.fits'%ia_template_name),file=fout)

                os.system('kcor.exe %s'%new_kcor)
            
    def run_lcfit(self):
        for ia_template in self.snsed:
            ia_template_name = ia_template.split('/')[-1].split('.')[0]
            for nml in [self.nml_ps1,self.nml_des]:
                with open(nml,'r') as fin:
                    lines = fin.readlines()
                new_nml = nml.replace('.nml','_%s.nml'%ia_template_name)
                with open(new_nml,'w') as fout:
                    for line in lines:
                        if 'TEXTFILE_PREFIX' not in line and 'KCOR_FILE' not in line:
                            print(line.replace('\n',''),file=fout)
                        elif 'TEXTFILE_PREFIX' in line:
                            print("  TEXTFILE_PREFIX = 'output/kcor_sys/%s'"%nml.split('/')[-1].replace('.nml','_%s'%ia_template_name),file=fout)
                        elif 'KCOR_FILE' in line:
                            print("  KCOR_FILE = %s"%line.split()[-1].replace('.fits','_%s.fits'%ia_template_name),file=fout)

                print('snlc_fit.exe %s'%new_nml)
                import pdb; pdb.set_trace()
                os.system('snlc_fit.exe %s'%new_nml)

    def comp_distances(self):
        hsiao_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_Hsiao07.FITRES.TEXT'; hsiao_frfile_des = 'output/kcor_sys/DES_RAISIN_Hsiao07.FITRES.TEXT'
        nugent_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_snflux_1a_Nugent2002.FITRES.TEXT'; nugent_frfile_des = 'output/kcor_sys/DES_RAISIN_snflux_1a_Nugent2002.FITRES.TEXT'
        p18_frfile_ps1 = 'output/kcor_sys/PS1_RAISIN_sn91t_flux.FITRES.TEXT'; p18_frfile_des = 'output/kcor_sys/DES_RAISIN_sn91t_flux.FITRES.TEXT'

        frbase = txtobj(hsiao_frfile_ps1,fitresheader=True); frbase2 = txtobj(hsiao_frfile_des,fitresheader=True)
        frn = txtobj(nugent_frfile_ps1,fitresheader=True); frn2 = txtobj(nugent_frfile_des,fitresheader=True)
        frp = txtobj(p18_frfile_ps1,fitresheader=True); frp2 = txtobj(p18_frfile_des,fitresheader=True)
        frbase = concat_simple(frbase,frbase2)
        frn = concat_simple(frn,frn2)
        frp = concat_simple(frp,frp2)
        
        frbase = mkraisincuts(frbase)
        frn = mkraisincuts(frn)
        frp = mkraisincuts(frp)
        
        iSort = np.argsort(frbase.zHD)
        frbase = sort_fitres(frbase,iSort); frn = sort_fitres(frn,iSort); frp = sort_fitres(frp,iSort)

        plt.clf()
        plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='Nugent02 - Hsiao+07')
        #plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='Nugent91T - Hsiao+07')
        plt.xlabel('$z_{CMB}$',fontsize=15)
        plt.ylabel('$\Delta\mu$',fontsize=15)
        plt.axhline(0,color='k',lw=2)
        plt.legend()
        import pdb; pdb.set_trace()

    def comp_distances_new(self):

        fe_frfile_ps1 = 'output/fit_nir/PS1_RAISIN_2011fe.FITRES.TEXT'; fe_frfile_des = 'output/fit_nir/DES_RAISIN_2011fe.FITRES.TEXT'
        fecomp_frfile_ps1 = 'output/fit_nir/PS1_RAISIN_2011fecomp.FITRES.TEXT'; fecomp_frfile_des = 'output/fit_nir/DES_RAISIN_2011fecomp.FITRES.TEXT'
        fe_frfile_csp = 'output/fit_nir/CSP_RAISIN_2011fe.FITRES.TEXT'; fecomp_frfile_csp = 'output/fit_nir/CSP_RAISIN_2011fecomp.FITRES.TEXT'
        
        frbase = txtobj(fecomp_frfile_ps1,fitresheader=True); frbase2 = txtobj(fecomp_frfile_des,fitresheader=True); frbase3 = txtobj(fecomp_frfile_csp,fitresheader=True)
        frbase = concat_simple(frbase,frbase2)
        frbase = concat_simple(frbase,frbase3)
        frn = txtobj(fe_frfile_ps1,fitresheader=True); frn2 = txtobj(fe_frfile_des,fitresheader=True); frn3 = txtobj(fe_frfile_csp,fitresheader=True)
        frn = concat_simple(frn,frn2)
        frn = concat_simple(frn,frn3)
        #frp = concat_simple(frp,frp2)
        
        frbase = mkraisincuts(frbase)
        frn = mkraisincuts(frn)
        #frp = mkraisincuts(frp)
        
        iSort = np.argsort(frbase.zHD)
        frbase = sort_fitres(frbase,iSort); frn = sort_fitres(frn,iSort) #; frp = sort_fitres(frp,iSort)

        plt.clf()
        plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='2011fe - Hsiao+07')
        #plt.plot(frbase.zCMB,frn.DLMAG-frbase.DLMAG,'o-',label='Nugent91T - Hsiao+07')
        plt.xlabel('$z_{CMB}$',fontsize=15)
        plt.ylabel('$\Delta\mu$',fontsize=15)
        plt.axhline(0,color='k',lw=2)
        plt.legend()
        import pdb; pdb.set_trace()

        #(Pdb) np.median(frn.DLMAG[frn.zHD < 0.1]-cosmo.mu(frn.zHD[frn.zHD < 0.1]))
        #-0.05808668625499536
        #(Pdb) np.median(frn.DLMAG[frn.zHD > 0.1]-cosmo.mu(frn.zHD[frn.zHD > 0.1]))
        #-0.04264189178691069
        #(Pdb) np.median(frbase.DLMAG[frbase.zHD < 0.1]-cosmo.mu(frbase.zHD[frbase.zHD < 0.1]))
        #-0.06454616341290631
        #(Pdb) np.median(frbase.DLMAG[frbase.zHD > 0.1]-cosmo.mu(frbase.zHD[frbase.zHD > 0.1]))
        
        
if __name__ == "__main__":

    kc = kcor()
    #kc.mk_kcor()
    #kc.run_lcfit()
    #kc.comp_distances_new()

    kc = interp_kcor()
    #kc.main()
    kc.plot()
    #kc.plot_dist()
    #kc.fetmpl()

