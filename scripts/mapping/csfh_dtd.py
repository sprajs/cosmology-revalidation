"""Equation reconstruction; no survey selection or fitted cosmology posterior.

B13 Eq F1/Table 6; MD14 Eq 15. C14 Eq 3 with alpha=20 from Park A3.
Radiation-free, flat CPL; Gyr ages. All corrections are zero at z=0.
"""
from common import ROOT, OUT, manifest
import json
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from astropy.cosmology import Flatw0waCDM
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HUBBLE_GYR = 977.7922216807892  # (Mpc/km) / seconds-per-Gyr
COSMOS = {
    'lcdm_H70': dict(H0=70., Om0=.3, w0=-1., wa=0.),
    'son_cpl_H70': dict(H0=70., Om0=.353, w0=-.42, wa=-1.75),
    'son_cpl_H63p6': dict(H0=63.6, Om0=.353, w0=-.42, wa=-1.75),
}
DTDS = {'C14_smooth': (.3, -1., 20.), 'cut300': (.3, -1., None), 'W26_cut40': (.04, -1.13, None)}

def csfh(z, form):
    if form == 'B13':
        return .180*np.exp(-np.logaddexp(np.log(10)*(-.997)*(z-1.243), np.log(10)*.241*(z-1.243)))
    if form == 'MD14':
        return .015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
    raise ValueError(form)

def clock(c, n=60001):
    # Integrate dt/d(ln a)=1/H; matter-era boundary is negligible at a=e^-16.
    la = np.linspace(-16., 0., n)
    a = np.exp(la)
    de = (1-c['Om0'])*a**(-3*(1+c['w0']+c['wa']))*np.exp(-3*c['wa']*(1-a))
    einv = 1/np.sqrt(c['Om0']*a**-3+de)
    age = (2*a[0]**1.5/(3*np.sqrt(c['Om0'])) + cumulative_trapezoid(einv, la, initial=0))*HUBBLE_GYR/c['H0']
    return PchipInterpolator(la, age), PchipInterpolator(age, 1/a-1, extrapolate=False)

def calculate(c, sfh, dtd, n=12001, clock_n=60001):
    age_of_lna, z_of_age = clock(c, clock_n)
    tp, s, alpha = DTDS[dtd]
    rows = []
    for z in np.linspace(0, 2.5, 251):
        t = float(age_of_lna(-np.log1p(z)))
        # Separate continuous domain for truncated models avoids grid misplacement at cutoff.
        lo = 1e-7 if alpha is not None else tp
        delay = np.geomspace(lo, t*(1-1e-9), n)
        formed_time = t-delay
        sf = csfh(z_of_age(formed_time), sfh)
        sf = np.nan_to_num(sf) # negligible tail earlier than first cosmic-time grid point
        if alpha is not None:
            x = np.log(delay/tp)
            phi = np.exp(alpha*x-np.logaddexp((alpha-s)*x, 0))
        else:
            phi = delay**s
        pdf = sf*phi
        cdf = cumulative_trapezoid(pdf, delay, initial=0)
        area = cdf[-1]
        mean = np.trapezoid(pdf*delay, delay)/area
        median = np.interp(area/2, cdf, delay)
        # Galaxy integrated formed-mass and surviving-mass means, NOT SN host weighted.
        dg = np.geomspace(1e-7, t*(1-1e-9), n)
        sg = np.nan_to_num(csfh(z_of_age(t-dg), sfh))
        survive = 1-.046*np.log(dg/.000276+1) # C14 A2: 0.276 Myr -> Gyr
        gm = np.trapezoid(sg*dg, dg)/np.trapezoid(sg, dg)
        gms = np.trapezoid(sg*survive*dg, dg)/np.trapezoid(sg*survive, dg)
        rows.append(dict(z=z, cosmic_age_gyr=t, mean_delay_gyr=mean, median_delay_gyr=median,
                         formed_mass_mean_gyr=gm, surviving_mass_mean_gyr=gms))
    df = pd.DataFrame(rows)
    for col in ['mean_delay_gyr', 'median_delay_gyr', 'formed_mass_mean_gyr', 'surviving_mass_mean_gyr']:
        df['delta_'+col] = df.loc[0, col]-df[col]
        df['correction_subtracted_mag_'+col] = .030*df['delta_'+col]
    return df

def main():
    rows, tests = [], {}
    for cname, c in COSMOS.items():
        af, _ = clock(c)
        za = np.array([0., .1, .5, 1., 2., 10.])
        astropy = Flatw0waCDM(**c, Tcmb0=0).age(za).value
        err = np.max(np.abs(af(-np.log1p(za))-astropy))
        tests[cname+'_astropy_age_max_abs_gyr'] = float(err)
        assert err < 2e-6, (cname, err)
        for sfh in ['B13', 'MD14']:
            for dtd in DTDS:
                d = calculate(c, sfh, dtd)
                d['cosmology'], d['csfh'], d['dtd'] = cname, sfh, dtd
                rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    fine = calculate(COSMOS['son_cpl_H63p6'], 'B13', 'C14_smooth', n=24001, clock_n=120001)
    coarse = df.query("cosmology=='son_cpl_H63p6' and csfh=='B13' and dtd=='C14_smooth'").reset_index(drop=True)
    for col in ['mean_delay_gyr', 'median_delay_gyr', 'formed_mass_mean_gyr']:
        tests['double_resolution_max_abs_'+col] = float(np.max(np.abs(fine[col]-coarse[col])))
        assert tests['double_resolution_max_abs_'+col] < .0001
    analytic_age0 = 2/(3*np.sqrt(.7))*np.arcsinh(np.sqrt(.7/.3))*HUBBLE_GYR/70
    tests['lcdm_analytic_age0_gyr'] = float(analytic_age0)
    assert abs(df.query("cosmology=='lcdm_H70'").iloc[0].cosmic_age_gyr-analytic_age0) < 2e-6
    low_h=calculate(dict(COSMOS['son_cpl_H70'],H0=69.3),'B13','C14_smooth')
    high_h=calculate(dict(COSMOS['son_cpl_H70'],H0=70.7),'B13','C14_smooth')
    for col in ['delta_mean_delay_gyr','delta_median_delay_gyr']:
        tests['d_log_'+col+'_d_log_H0_at70_z1']=float((np.log(high_h.loc[100,col])-np.log(low_h.loc[100,col]))/np.log(70.7/69.3))
    full = OUT/'csfh-dtd-curves.csv'
    df.to_csv(full, index=False)
    summary = df[np.isclose(df.z, 0) | np.isclose(df.z, 1) | np.isclose(df.z, 2)]
    summary_path = OUT/'csfh-dtd-summary.csv'
    summary.to_csv(summary_path, index=False)
    check = OUT/'csfh-dtd-checks.json'
    check.write_text(json.dumps(tests, indent=2)+'\n')
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for dtd in DTDS:
        q = df.query("cosmology=='son_cpl_H63p6' and csfh=='B13' and dtd==@dtd")
        ax[0].plot(q.z, q.mean_delay_gyr, label=dtd+' mean')
        ax[0].plot(q.z, q.median_delay_gyr, ls='--', label=dtd+' median')
        ax[1].plot(q.z, q.correction_subtracted_mag_median_delay_gyr, label=dtd+' median')
        ax[1].plot(q.z, q.correction_subtracted_mag_mean_delay_gyr, ls=':', label=dtd+' mean')
    ax[0].set(ylabel='Progenitor delay (Gyr)', xlabel='Redshift', title='Cosmic-volume equation reconstruction')
    ax[1].set(ylabel='Subtracted distance modulus (mag)', xlabel='Redshift', title='Fixed 0.030 mag/Gyr coefficient')
    for a in ax:
        a.legend(fontsize=7); a.grid(alpha=.2)
    figpath=OUT/'csfh-dtd.png'; fig.savefig(figpath,dpi=160); plt.close(fig)
    inputs = [__file__, ROOT/'scripts/mapping/common.py', ROOT/'docs/experiments/mapping-plan.md',
              ROOT/'uv.lock', ROOT/'papers/text/son2025-published.txt', ROOT/'papers/text/childress2014-published.txt',
              ROOT/'papers/text/wiseman2026-published.txt',
              ROOT/'papers/text/2605.12596v1.txt', ROOT/'papers/text/2503.14738v3.txt',
              ROOT/'runs/mapping/sources-2026-09-20/behroozi2013.pdf', ROOT/'runs/mapping/sources-2026-09-20/madau2014.pdf']
    manifest('csfh-dtd', dict(cosmologies=COSMOS, dtds=DTDS, redshift_grid='0:0.01:2.5', delay_points=12001,
        radiation=False, selection=False, age_slope_mag_gyr=-.030, seed=None,
        interpretation='Equation-level reconstruction; not exact S25 author configuration or posterior'), inputs,
        [full,summary_path,check,figpath])
    print(summary.query('z==1')[['cosmology','csfh','dtd','delta_mean_delay_gyr','delta_median_delay_gyr']].to_string(index=False))
    print(json.dumps(tests,indent=2))

if __name__=='__main__': main()
