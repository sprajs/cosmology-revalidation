"""Independent numerical controls for the six RAISIN/CSP clean workflows.

Run from any directory with the clean environment. Input parsing, photon integration,
GLS fits and covariance projections intentionally do not call workflow helpers.
Bootstrap intervals are conditional descriptive resampling, not systematic-error CIs.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import numpy as np
from astropy.io import fits
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import log_ndtr
from scipy.stats import multivariate_normal

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
INPUTS = set()
CHECKS = []


def checked(name, actual, expected, atol=1e-10):
    good = bool(np.allclose(actual, expected, atol=atol, rtol=1e-10))
    CHECKS.append({'name': name, 'passed': good, 'absolute_gap': float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))})
    if not good:
        raise AssertionError(f'{name}: {actual} != {expected}')


def read_text(path):
    INPUTS.add(path)
    return path.read_text()


def table(path):
    INPUTS.add(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        lines = f.read().splitlines()
    names = next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    result = {}
    for line in lines:
        if not line.startswith('SN:'):
            continue
        tokens = line.split()[1:]
        assert len(tokens) == len(names)
        row = dict(zip(names, tokens))
        cid = row.pop('CID')
        assert cid not in result
        result[cid] = {}
        for k, v in row.items():
            try:
                result[cid][k] = float(v)
            except ValueError:
                result[cid][k] = v
    return result


def phot(path):
    header, rows, columns = {}, [], None
    for line in read_text(path).splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'VARLIST:':
            columns = p[1:]
        elif p[0] == 'OBS:':
            assert len(p)-1 == len(columns)
            rows.append(dict(zip(columns, p[1:])))
        elif p[0].endswith(':'):
            header[p[0][:-1]] = ' '.join(p[1:])
    return header, rows


def numeric(path):
    INPUTS.add(path)
    return np.loadtxt(path)


def csvread(path):
    return list(csv.DictReader(read_text(path).splitlines()))


def summary(results, name):
    return json.loads(read_text(results[name] / 'summary.json'))


def interval(x):
    return [float(t) for t in np.quantile(x, [0.025, .5, .975])]


def raisin_check(results, rng):
    base = DATA / 'raisin/distances/w'
    tabs = {b: table(base / f'{b}_dist/RAISIN_combined_FITOPT000.FITRES') for b in ('nir','optical','opticalnir')}
    ids = list(tabs['nir'])
    assert all(list(t) == ids for t in tabs.values())
    assert all([t[c]['zHD'] for c in ids] == [tabs['nir'][c]['zHD'] for c in ids] for t in tabs.values())
    z = np.array([tabs['nir'][c]['zHD'] for c in ids])
    low, high = z < .1, z > .2
    checked('raisin cohort', [len(ids), sum(low), sum(high)], [79,42,37])
    w = high/sum(high)-low/sum(low)
    published = summary(results, 'raisin')
    comparisons = []
    for branch in ('optical','opticalnir'):
        values = {}
        for variant, bc, mc in [('released',0,0),('without_mass',0,-1),('without_bias',1,0),('without_bias_or_mass',1,-1)]:
            d = np.array([(tabs[branch][c]['DLMAG']+bc*tabs[branch][c]['DLMAG_biascor']+mc*tabs[branch][c]['MASS_CORR'])-(tabs['nir'][c]['DLMAG']+bc*tabs['nir'][c]['DLMAG_biascor']+mc*tabs['nir'][c]['MASS_CORR']) for c in ids])
            value = float(np.mean(d[high])-np.mean(d[low]))
            expected = next(r for r in published['contrasts'] if r['branch_minus_nir']==branch and r['variant']==variant)
            checked(f'{branch} {variant} contrast', value, expected['high_minus_low_mag'])
            se = np.sqrt(d[high].var(ddof=1)/sum(high)+d[low].var(ddof=1)/sum(low))
            checked(f'{branch} {variant} SE', se, expected['descriptive_paired_se_mag'])
            # Resample physical objects within each redshift group, preserving paired branches.
            boot = rng.choice(d[high], (10000,sum(high)), replace=True).mean(axis=1)-rng.choice(d[low], (10000,sum(low)), replace=True).mean(axis=1)
            loo = []
            for k in range(len(ids)):
                include = np.arange(len(ids)) != k
                loo.append(float(np.mean(d[high & include])-np.mean(d[low & include])))
            values[variant] = value
            comparisons.append({'branch':branch,'variant':variant,'high_minus_low_mag':value,'paired_descriptive_se_mag':float(se),'within_redshift_object_bootstrap_95_percentile_mag':interval(boot),'leave_one_object_out_range_mag':[min(loo),max(loo)]})
        checked(f'{branch} additive terms closure',values['without_bias_or_mass'],values['without_bias']+values['without_mass']-values['released'])
    covariance = []
    for branch in tabs:
        folder = base / f'{branch}_syst'
        raw = numeric(folder/'RAISIN_all.covmat')
        n = int(raw[0]); assert n==79 and len(raw)==1+n*n
        off = raw[1:].reshape(n,n)
        alltab = numeric(folder/'RAISIN_all_lcparams_cosmosis.txt')
        stat = numeric(folder/'RAISIN_stat_lcparams_cosmosis.txt')
        checked(branch+' covariance row z identity', alltab[:,1], z, 5e-7)
        checked(branch+' covariance all/stat z identity', alltab[:,1], stat[:,1])
        C = off+np.diag(alltab[:,5]**2-stat[:,5]**2)
        checked(branch+' covariance symmetry', C,C.T)
        # Independently remove intercept using centering by row/column means.
        target = C-C.mean(axis=0)[None,:]-C.mean(axis=1)[:,None]+C.mean()
        vectors=[]
        for k in range(1,29):
            other=table(base/f'{branch}_dist/RAISIN_combined_FITOPT{k:03d}.FITRES')
            v=np.array([other[c]['DLMAG']-tabs[branch][c]['DLMAG'] for c in ids]); v-=v.mean()
            vectors.append(np.outer(v,v).ravel())
        X=np.array(vectors).T
        coef,_,rank,_=np.linalg.lstsq(X,target.ravel(),rcond=None)
        residual=np.linalg.norm(target.ravel()-X@coef)/np.linalg.norm(target)
        # Literal six-decimal error rounding contributes to the restored diagonal.
        error_bound=2*(alltab[:,5]+stat[:,5])*.5e-6+2*(.5e-6)**2
        tokens=read_text(folder/'RAISIN_all.covmat').split()[1:]
        entry_bound=np.array([float(Decimal(10)**Decimal(Decimal(x).as_tuple().exponent))/2 for x in tokens]).reshape(n,n)
        np.fill_diagonal(entry_bound,error_bound)
        # Symmetric matrix perturbation operator norm is bounded by maximum absolute row sum.
        rounding_bound=float(entry_bound.sum(axis=1).max())
        min_eigen=float(np.linalg.eigvalsh(C)[0])
        constrained=next(r for r in published['covariance_reconstruction'] if r['branch']==branch)['relative_covariance_shape_residual']
        assert residual <= constrained+1e-9
        covariance.append({'branch':branch,'systematic_min_eigenvalue':min_eigen,'conservative_literal_rounding_operator_bound':rounding_bound,'negative_systematic_eigenvalue_within_rounding_bound':bool(min_eigen>=-rounding_bound),'total_covariance_min_eigenvalue':float(np.linalg.eigvalsh(C+np.diag(stat[:,5]**2))[0]),'exported_covariance_diagonal_max_abs':float(abs(np.diag(off)).max()),'contrast_systematic_variance':float(w@C@w),'variant_outer_product_rank':int(rank),'unconstrained_span_relative_residual':float(residual),'nonnegative_relative_residual':constrained,'meaning':'Even allowing negative coefficients cannot reduce residual below the unconstrained span residual; this diagnoses this centered raw-DLMAG FITOPT 1-28 vector representation, not a defective physical covariance or every possible author source reconstruction.'})
    seen={}; total=Counter(); header_differences=[]
    for path in sorted((DATA/'raisin/photometry/RAISIN').rglob('*')):
        if not path.is_file() or path.suffix.lower()!='.dat': continue
        header,rows=phot(path); cid=header['SNID'].split()[0]
        if cid in seen:
            assert cid not in ids
            continue
        seen[cid]=header
        for r in rows:
            f,e=float(r['FLUXCAL']),float(r['FLUXCALERR'])
            total.update(rows=1,negative=int(f<0),zero=int(f==0),nonfinite=int(not np.isfinite(f)),nonpositive_error=int(e<=0))
        if cid in ids:
            hp=float(header['PEAKMJD'].split()[0]); fp=tabs['nir'][cid]['PKMJD']
            if abs(hp-fp)>.001:
                header_differences.append({'CID':cid,'header_peak':hp,'FITRES_peak':fp,'gap_days':hp-fp,'header_float32':float(np.float32(hp)),'FITRES_float32':float(np.float32(fp))})
    for k,v in total.items(): checked('photometry '+k,v,published['photometry_totals'][k])
    return {'contrasts':comparisons,'covariance':covariance,'photometry_totals':dict(total),'peak_header_gaps_above_0_001_day':header_differences,'bootstrap_scope':'Descriptive object resampling conditional on these releases, fixed selection and fixed redshift strata. Shared calibration/model systematics and population selection are not incorporated.'}


def timing_check(results):
    nir=table(DATA/'timing/nir.FITRES.gz'); joint=table(DATA/'timing/optnir.FITRES.gz')
    common=sorted(set(nir)&set(joint)); old=summary(results,'timing')
    fields=['SIM_PKMJD','SIM_DLMAG','SIM_AV','SIM_STRETCH','SIM_RV']
    assert all(nir[c][f]==joint[c][f] for c in common for f in fields)
    selected=np.array([joint[c]['AV']<.3*joint[c]['RV'] and .75<joint[c]['STRETCH']<1.185 for c in common])
    checked('timing membership',[len(nir),len(joint),sum(selected)],[30000,29995,25606])
    output={'missing_joint_CIDs':sorted(set(nir)-set(joint)),'common_n':len(common),'selected_n':int(sum(selected)),'simulation_description':'The NIR fitted peak is exactly the stored initializer; the 0.0101-day scatter is therefore not recovered timing precision.'}
    for branch,t in [('nir',nir),('joint',joint)]:
        dt=np.array([t[c]['PKMJD']-t[c]['SIM_PKMJD'] for c in common])
        checked(branch+' timing SD', dt.std(), old[branch+'_peak_minus_truth']['sd'])
        checked(branch+' timing mean', dt.mean(), old[branch+'_peak_minus_truth']['mean'])
        dmu=np.array([t[c]['DLMAG']-t[c]['SIM_DLMAG'] for c in common])
        output[branch]={'peak_residual_mean_day':float(dt.mean()),'peak_residual_population_sd_day':float(dt.std()),'distance_residual_mean_mag':float(dmu.mean()),'selected_distance_residual_mean_mag':float(dmu[selected].mean()),'rejected_distance_residual_mean_mag':float(dmu[~selected].mean()),'distance_timing_Pearson_r':float(np.corrcoef(dt,dmu)[0,1]),'exact_initializer_matches':sum(t[c]['PKMJD']==t[c]['PKMJDINI'] for c in common)}
    return output


def dec(x): return Decimal(x)


def lineage_check(results):
    raw=[]; snpy=[]; nirs={'Y','Ydw','J','Jrc2','Jdw','H','Hdw'}; merge={'Jdw':'J','Hdw':'H'}
    path=DATA/'csp/photometry.tgz'; INPUTS.add(path)
    with tarfile.open(path) as a:
        for line in a.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
            p=line.split()
            if p[1] in nirs: raw.append((p[0].lower(),p[1],*map(dec,p[2:])))
        for m in a.getmembers():
            if m.isfile() and m.name.endswith('_snpy.txt'):
                lines=a.extractfile(m).read().decode().splitlines(); cid=lines[0].split()[0].lower();band=None
                for line in lines[1:]:
                    p=line.split()
                    if not p or p[0].startswith('#'):continue
                    if p[0]=='filter':band=p[1]
                    elif band in nirs:snpy.append((cid,band,*map(dec,p)))
    rc=Counter((n,merge.get(b,b),t,m,e) for n,b,t,m,e in raw);sc=Counter(snpy)
    assert rc==sc
    lookup=defaultdict(list); physical=defaultdict(set)
    for n,b,t,m,e in snpy:lookup[n].append((b,t,m,e))
    for n,b,t,m,e in raw:physical[n,merge.get(b,b),t].add(b)
    labels={'J':'J','Jrc2':'j','Y':'Y','Ydw':'y','H':'H'}; rev={v:k for k,v in labels.items()}
    base=DATA/'raisin/photometry/RAISIN/CSPDR3_RAISIN'; mismatch=[];matched=0;matched_nonempty=0;raw_label_count=0;listed=0
    cosmology_ids=set(table(DATA/'raisin/distances/w/nir_dist/RAISIN_combined_FITOPT000.FITRES'))
    for filename in read_text(base/'CSPDR3_RAISIN.LIST').splitlines():
        if not filename.strip() or filename.startswith('#'): continue
        header,obs=phot(base/filename); cid=header['SNID'].split()[0]; key='sn'+cid.lower(); listed+=1
        predictions=[]
        for b,t,m,e in lookup[key]:
            f=10**(-.4*(float(m)-27.5)); er=f*np.expm1(.4*np.log(10)*float(e))
            predictions.append((labels[b],t,m,e,dec(f'{f:.5e}'),dec(f'{er:.5e}')))
        actual=[(r['FLT'],dec(r['MJD'])-53000,dec(r['MAG']),dec(r['MAGERR']),dec(r['FLUXCAL']),dec(r['FLUXCALERR'])) for r in obs if r['FLT'] in rev]
        pc,ac=Counter(predictions),Counter(actual)
        for b,t,*_ in actual:raw_label_count+=physical[key,rev[b],t]=={'Jdw'}
        if pc==ac:
            matched+=1
            matched_nonempty+=bool(sum(ac.values()))
        else:mismatch.append({'CID':cid,'archived_SNpy_rows':sum(pc.values()),'released_rows':sum(ac.values()),'source_object_missing_in_this_DR3_archive':key not in lookup or len(lookup[key])==0,'in_79_object_cosmology_cohort':cid in cosmology_ids,'predicted_only':sum((pc-ac).values()),'released_only':sum((ac-pc).values())})
    old=summary(results,'csp-lineage')
    checked('CSP lineage counts',[len(raw),len(snpy),listed,matched,raw_label_count],[old['raw_NIR_rows'],old['snpy_NIR_rows'],old['listed_objects'],old['converter_exact_objects'],old['released_unique_Jdw_rows']])
    return {'raw_rows':len(raw),'intermediate_rows':len(snpy),'listed_objects':listed,'exact_objects':matched,'exact_nonempty_objects':matched_nonempty,'empty_object_identities':matched-matched_nonempty,'mismatches':mismatch,'raw_Jdw_rows':sum(r[1]=='Jdw' for r in raw),'released_unambiguous_Jdw_rows':raw_label_count,'interpretation':'All three converter nonmatches lack source-object rows in this archived DR3 input; they are source-coverage failures, not demonstrated converter errors. Their later dates (2012, 2015) distinguish them from the 2004-2009 DR3 data.'}


def integrate_polynomial(wave,sed,band):
    lo,hi=band[0,0],band[-1,0]; assert wave[0]<=lo and hi<=wave[-1]
    knots=np.unique(np.concatenate([band[:,0],wave[(wave>lo)&(wave<hi)]])); h=np.diff(knots); L=knots[:-1]
    T=np.interp(knots,band[:,0],band[:,1]); S=np.interp(knots,wave,sed)
    a,b=T[:-1],np.diff(T);c,d=S[:-1],np.diff(S)
    # Analytic integral of (L+h*u)*(a+b*u)*(c+d*u) over 0<=u<=1.
    return float(np.sum(h*(L*a*c+(h*a*c+L*(a*d+b*c))/2+(h*(a*d+b*c)+L*b*d)/3+h*b*d/4)))


def passbands_check(results):
    path=DATA/'raisin/kcor/kcor_CSPDR3_BD17.fits';INPUTS.add(path)
    with fits.open(path) as f:
        h=f['SN SED'].header;w=h['LMIN']+np.arange(h['NBL'])*h['LBIN'];times=h['TMIN']+np.arange(h['NBT'])*h['TBIN'];S=f['SN SED'].data.field(0).astype(float).reshape(len(times),len(w));rw=f['PrimarySED'].data.field(0).astype(float);rf=f['PrimarySED'].data['BD17'].astype(float)
    pairs={'J_WIRC_minus_RC1':('Jrc1_SWO_TAM_scan_atm.dat','J_DUP_TAM_scan_atm.dat'),'H_WIRC_minus_RetroCam':('H_SWO_TAM_scan_atm.dat','H_DUP_TAM_scan_atm.dat'),'Y_WIRC_minus_RetroCam':('Y_SWO_TAM_scan_atm.dat','Y_DUP_TAM_scan_atm.dat'),'J_RC2_minus_RC1':('Jrc1_SWO_TAM_scan_atm.dat','Jrc2_SWO_TAM_scan_atm.dat')}
    bands={n:numeric(DATA/'csp/filters'/n) for n in set(sum(pairs.values(),()))}
    ref={n:integrate_polynomial(rw,rf,b) for n,b in bands.items()};by_pair=defaultdict(list);maxgap=0
    rows=csvread(results['csp-passbands']/'passband_grid.csv')
    for row in rows:
        z=float(row['z']);phase=float(row['phase']);ix=np.flatnonzero(times==phase).item();a,b=pairs[row['pair']]
        ca=integrate_polynomial(w*(1+z),S[ix]/(1+z),bands[a]);cb=integrate_polynomial(w*(1+z),S[ix]/(1+z),bands[b]);dm=-2.5*np.log10(cb/ref[b])+2.5*np.log10(ca/ref[a]);gap=abs(dm-float(row['delta_mag_equal_BD17']));maxgap=max(maxgap,gap); assert gap<1e-11
        by_pair[row['pair']].append(row)
    return {'analytic_piecewise_cubic_max_gap_mag':maxgap,'rows':len(rows),'pairs':{k:{'equal_BD17_mag_range':[min(float(r['delta_mag_equal_BD17']) for r in v),max(float(r['delta_mag_equal_BD17']) for r in v)],'max_abs_phase_contrast_mag':max(abs(float(r['phase_contrast_vs_phase0'])) for r in v),'max_abs_redshift_contrast_mag':max(abs(float(r['redshift_contrast_vs_z0'])) for r in v)} for k,v in by_pair.items()},'scope':'Exact integration identity for tabulated piecewise-linear curves. This validates numerical integration, not spectral realism, physical calibration, or a cosmological distance correction.'}


def baseline_check(results,rng):
    cohort=json.loads(read_text(DATA/'signed/cohort.json'));old=summary(results,'signed-baseline');output={}
    for cut in (180,365):
        groups=[];objectq=[];objectdf=[]
        for case in cohort:
            h,raw=phot(DATA/'signed'/case['raw_file']);assert h['SNID'].split()[0]==case['source_CID']; rows=[r for r in raw if float(r['MJD'])<case['peak_header']-cut]
            sq,sd=0,0
            for band in 'griz':
                chosen=[r for r in rows if r['FLT']==band];y=np.array([float(r['FLUXCAL']) for r in chosen]);e=np.array([float(r['FLUXCALERR']) for r in chosen]);assert np.isfinite(y).all() and np.isfinite(e).all() and (e>0).all();v=1/e
                # QR projection solves the one-column weighted linear model independently.
                qmat,_=np.linalg.qr(v[:,None]);residual=y/e-qmat@(qmat.T@(y/e));Q=float(residual@residual);df=len(y)-1;sq+=Q;sd+=df
                groups.append({'CID':case['CID'],'band':band,'n':len(y),'negative':int(sum(y<0)),'Q':Q,'dof':df})
            objectq.append(sq);objectdf.append(sd)
        Q=sum(g['Q'] for g in groups);df=sum(g['dof'] for g in groups);checked(f'baseline Q cut {cut}',Q,old['cuts'][str(cut)]['Q'],1e-8);checked(f'baseline rows cut {cut}',sum(g['n'] for g in groups),old['cuts'][str(cut)]['rows'])
        q,d=np.array(objectq),np.array(objectdf);idx=rng.integers(0,len(cohort),size=(10000,len(cohort)));boot=q[idx].sum(axis=1)/d[idx].sum(axis=1)
        output[str(cut)]={'rows':sum(g['n'] for g in groups),'negative':sum(g['negative'] for g in groups),'groups':len(groups),'Q':Q,'dof':df,'Q_per_dof':Q/df,'whole_object_bootstrap_Q_per_dof_95_percentile':interval(boot),'leave_one_object_out_Q_per_dof_range':[float(min((Q-q)/(df-d))),float(max((Q-q)/(df-d)))],'object_Q_per_dof':[{'CID':case['CID'],'Q_per_dof':float(a/b)} for case,a,b in zip(cohort,q,d)]}
    output['scope']='Whole-object bootstrap preserves within-object dependence, but ten selected objects are few and may not be independent across fields. The two nested time cuts are sensitivity analyses, not independent experiments. No Gaussian p-value or empirical error renormalization is inferred.'
    return output


def signs_check(results):
    run=json.loads(read_text(results['sign-selection']/'run.json'));cfg=run['config'];rows=csvread(results['sign-selection']/'recovery.csv');N=cfg['epochs'];rng=np.random.default_rng(cfg['seed']);template=np.exp(-np.linspace(-2,2,N)**2/2);C=np.eye(N)+.06**2*np.ones((N,N));W=np.linalg.inv(C);observed=defaultdict(dict)
    for row in rows:observed[int(row['replicate'])][row['method']]=float(row['amplitude'])
    glsgap=0;naivegap=0;quadgap=0; selected=[]
    for rep in range(cfg['replicates']):
        y=cfg['true_amplitude']*template+.06*rng.normal()+rng.normal(size=N);keep=y>0
        exact=float(template@W@y/(template@W@template));glsgap=max(glsgap,abs(exact-observed[rep]['all_signed']))
        wk=np.linalg.inv(C[np.ix_(keep,keep)]);t=template[keep]; naive=float(t@wk@y[keep]/(t@wk@t));naivegap=max(naivegap,abs(naive-observed[rep]['positive_naive']))
        if rep % 11==0:
            # Direct latent integral, conditioned on retained observations via dense covariance.
            K=C[np.ix_(keep,keep)];inv=np.linalg.inv(K);cross=.06*np.ones(sum(keep));postvar=1-float(cross@inv@cross)
            def objective(amplitude):
                mean=amplitude*template;res=y[keep]-mean[keep];postmean=float(cross@inv@res);lp=multivariate_normal.logpdf(y[keep],mean=mean[keep],cov=K)
                f=lambda z: np.exp(-z*z/2-.5*np.log(2*np.pi)+np.sum(log_ndtr(-mean[~keep]-.06*(postmean+np.sqrt(postvar)*z))))
                probability,err=quad(f,-12,12,epsabs=1e-13,epsrel=1e-11)
                assert probability>0 and err<max(1e-11,1e-8*probability)
                return -(lp+np.log(probability))
            fit=minimize_scalar(objective,bounds=(-3,3),method='bounded',options={'xatol':1e-9});assert fit.success;gap=abs(fit.x-observed[rep]['censored_with_known_schedule']);quadgap=max(quadgap,gap);selected.append(rep)
    checked('all signed GLS amplitude control',glsgap,0,1e-6);checked('positive-only GLS amplitude control',naivegap,0,1e-6);checked('adaptive-quadrature censor fit control',quadgap,0,1e-6)
    methods=list(observed[0]);out={}
    for method in methods:
        v=np.array([observed[r][method] for r in range(cfg['replicates'])]);paired=v-np.array([observed[r]['all_signed'] for r in range(cfg['replicates'])]);mean=float(paired.mean());se=float(paired.std(ddof=1)/np.sqrt(len(v)))
        out[method]={'mean':float(v.mean()),'bias':float(v.mean()-cfg['true_amplitude']),'paired_mean_difference_from_all_signed':mean,'paired_difference_mcse':se,'paired_normal_95_interval':[mean-1.96*se,mean+1.96*se]}
    return {'replicates':cfg['replicates'],'all_signed_closed_GLS_max_amplitude_gap':glsgap,'positive_only_closed_GLS_max_amplitude_gap':naivegap,'adaptive_conditional_latent_integral_fit_max_gap':quadgap,'adaptive_checked_replicates':selected,'methods':out,'scope':'Simulation-only. Pairing uses identical random draws and reduces Monte Carlo error in method comparison. Neither the generating Gaussian errors nor the known omitted-epoch schedule has been established for the real survey.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--results',default='baseline');parser.add_argument('--output',type=Path,default=ROOT/'validation/reports/raisin.json');args=parser.parse_args()
    names=['raisin','timing','csp-lineage','csp-passbands','signed-baseline','sign-selection'];results={n:ROOT/'results'/f'{args.results}/{n}' for n in names}
    for n,path in results.items():assert path.is_dir(),(n,path)
    rng=np.random.default_rng(26092631)
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'workflow_result_directory':args.results,'audit_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'status':'passed','raisin':raisin_check(results,rng),'timing':timing_check(results),'csp_lineage':lineage_check(results),'csp_passbands':passbands_check(results),'signed_baseline':baseline_check(results,rng),'sign_selection':signs_check(results),'checks':CHECKS}
    report['inspected_inputs_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUTS)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':'passed','checks':len(CHECKS),'input_files':len(INPUTS),'report':str(args.output)}))

if __name__=='__main__':main()
