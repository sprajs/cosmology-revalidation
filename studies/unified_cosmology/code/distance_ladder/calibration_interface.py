"""Pinned Pantheon+SH0ES Gaussian interface and calibration-decomposition audit.

This is an alternative SN compilation, not a factor to multiply by Dovekie or
an additional SH0ES H0 prior. No cosmological fit is performed by this script.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace
import urllib.request

import numpy as np
import pandas as pd
from scipy import linalg

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/calibration-interface'
LADDER = ROOT / '.work/unified-cosmology/distance-ladder'
OUT = ROOT / 'studies/unified_cosmology/results/distance_ladder/calibration-interface.json'
COMMIT = 'c447f0fea703fcd0fff57de5000947b5ca81286b'
BASE = f'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/{COMMIT}/'
SOURCES = {
    'Pantheon+SH0ES.dat': ('Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat', '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8'),
    'Pantheon+SH0ES_STAT+SYS.cov': ('Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov', 'abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc'),
    'Pantheon+SH0ES_STATONLY.cov': ('Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STATONLY.cov', '9f177129a332735d3637affd20054080d5260815f3ca0809120c05b2c902297f'),
    'distances-README': ('Pantheon+_Data/4_DISTANCES_AND_COVAR/README', 'e2b0d262757f01c1794a938c78d32600a21e289b2a0320e5c660c4c6fc9aa87e'),
    'Pantheon+SH0ES_cosmosis_likelihood.py': ('Pantheon+_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon+SH0ES_cosmosis_likelihood.py', '345fac3781a5cb930b95e91c1c07eb17dcf99b441703bb5e449477519240a59d'),
}
HOST_ALIASES = {'N0105': 'N105A', 'N0976': 'N976A'}
SN_ALIASES = {'2005df_ANU': '2005df', '2008fv_comb': '2008fv'}
R22_SHA = '572e5d2c0f32bf6e80a118eec13ed0107e1e2830f9ebe90b5bc4e9d1b7fa5df1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def acquire(work=WORK, allow_download=False):
    work.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, (upstream, expected) in SOURCES.items():
        target = work / name
        previous = ROOT / '.work/unified-cosmology/inference/pantheon' / Path(upstream).name
        if not target.exists() and previous.exists() and sha(previous) == expected:
            shutil.copyfile(previous, target)
        if not target.exists():
            if not allow_download:
                raise FileNotFoundError(f'{target}: rerun with --acquire')
            with urllib.request.urlopen(BASE + upstream, timeout=120) as response:
                body = response.read()
            if hashlib.sha256(body).hexdigest() != expected:
                raise ValueError(f'Upstream identity mismatch: {upstream}')
            target.write_bytes(body)
        if sha(target) != expected:
            raise ValueError(f'Input identity mismatch: {target}')
        records[rel(target)] = {'url': BASE + upstream, 'sha256': expected, 'bytes': target.stat().st_size}
    return records


def covariance(path):
    with path.open() as stream:
        n = int(stream.readline())
        matrix = np.loadtxt(stream).reshape(n, n)
    assert np.isfinite(matrix).all() and np.max(abs(matrix-matrix.T)) < 3.01e-8
    return matrix


class ReleasedCalibration:
    """Released covariance, calibrators and a single improper flat M measure.

    Supply D_A(zHD) in Mpc at ``z_hd_noncalibrator`` in the given row order.
    Repeated observations remain separate and retain their full covariance.
    The absolute evidence normalization of the flat-M measure is undefined.
    This class never adds Cepheid, H0, peculiar-velocity or diagonal errors.
    """
    def __init__(self, work=WORK):
        self.input_records = acquire(Path(work))
        original = pd.read_csv(Path(work) / 'Pantheon+SH0ES.dat', sep=r'\s+')
        self.mask = (original.zHD.to_numpy() > .01) | original.IS_CALIBRATOR.to_numpy().astype(bool)
        self.original_indices = np.flatnonzero(self.mask)
        self.data = original.loc[self.mask].reset_index(drop=True)
        self.calibrator = self.data.IS_CALIBRATOR.to_numpy().astype(bool)
        self.z_hd_noncalibrator = self.data.zHD.to_numpy()[~self.calibrator]
        self.z_hel_noncalibrator = self.data.zHEL.to_numpy()[~self.calibrator]
        self.C_raw = covariance(Path(work) / 'Pantheon+SH0ES_STAT+SYS.cov')[np.ix_(self.mask, self.mask)]
        self.raw_covariance_antisymmetry_max = float(np.max(abs(self.C_raw-self.C_raw.T)))
        # The printed release differs from its transpose by <=3e-8 mag^2.
        # A proper Gaussian needs a symmetric covariance; raw bytes are retained.
        self.C = (self.C_raw+self.C_raw.T)/2
        self.chol = linalg.cholesky(self.C, lower=True)
        self.one_white = linalg.solve_triangular(self.chol, np.ones(len(self.data)), lower=True)
        self.flat_m_precision = float(self.one_white @ self.one_white)
        self.logdet_covariance = float(2 * np.log(np.diag(self.chol)).sum())
        self.log_normalization = float(self.logdet_covariance + np.log(self.flat_m_precision) + (len(self.data)-1)*np.log(2*np.pi))

    def theory_mu(self, angular_diameter_distance_mpc):
        da = np.asarray(angular_diameter_distance_mpc, dtype=float)
        if da.shape != self.z_hd_noncalibrator.shape or not np.all(np.isfinite(da) & (da > 0)):
            raise ValueError('Require positive finite D_A for every noncalibrator row in original selection order')
        mu = self.data.CEPH_DIST.to_numpy().copy()
        mu[~self.calibrator] = 5*np.log10((1+self.z_hd_noncalibrator)*(1+self.z_hel_noncalibrator)*da)+25
        return mu

    def evaluate(self, angular_diameter_distance_mpc):
        residual = self.data.m_b_corr.to_numpy() - self.theory_mu(angular_diameter_distance_mpc)
        rw = linalg.solve_triangular(self.chol, residual, lower=True)
        mhat = float(self.one_white @ rw / self.flat_m_precision)
        projected = rw-self.one_white*mhat
        chi2 = float(projected @ projected)
        return {'loglike': -.5*(chi2+self.log_normalization), 'chi2': chi2,
                'M_conditional_mean': mhat, 'M_conditional_sigma': self.flat_m_precision**-.5,
                'logdet_covariance': self.logdet_covariance, 'flat_m_precision': self.flat_m_precision,
                'data_rows': len(self.data), 'covariance_replaced': False}


def published_hosts(work):
    pdf = LADDER/'Riess2022v3.pdf'
    if sha(pdf) != R22_SHA:
        raise ValueError('R22 source changed')
    textfile = work/'Riess2022v3-layout.txt'
    subprocess.run(['pdftotext', '-layout', str(pdf), str(textfile)], check=True)
    entries = []
    for line in textfile.read_text().splitlines():
        t = line.split()
        if len(t) != 12 or not t[0].isdigit() or not t[1].startswith(('M','N','U')) or not 1 <= int(t[0]) <= 42:
            continue
        try:
            entry = {'table6_number': int(t[0]), 'paper_host': t[1], 'paper_SN': t[2],
                     'table6_SN_free_mu_rounded': float(t[9]), 'table6_SN_free_sigma_rounded': float(t[10])}
        except ValueError:
            continue
        entry['host'] = HOST_ALIASES.get(t[1], t[1])
        entries.append(entry)
    assert [e['table6_number'] for e in entries] == list(range(1,43))
    assert len({e['paper_SN'].lower() for e in entries}) == 42
    return entries, textfile


def method(path, class_name, method_name, globals_dict):
    """Execute only the unmodified inspected method body, avoiding framework setup."""
    cls = next(x for x in ast.parse(path.read_text()).body if isinstance(x, ast.ClassDef) and x.name == class_name)
    fun = next(x for x in cls.body if isinstance(x, ast.FunctionDef) and x.name == method_name)
    namespace = dict(globals_dict)
    exec(compile(ast.Module(body=[fun], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[method_name]


def numerical_checks(model, work):
    import cobaya
    from cobaya.likelihoods.sn.pantheonplusshoes import PantheonPlusShoes
    from cobaya.likelihoods.sn.pantheonplus import PantheonPlus
    cp = Path(inspect.getfile(PantheonPlusShoes))
    pp = Path(inspect.getfile(PantheonPlus))
    raw = method(cp, 'PantheonPlusShoes', 'alpha_beta_logp', {'np': np})
    profile = method(pp, 'PantheonPlus', '_marginalize_abs_mag', {'np': np})
    cosmosis = work/'Pantheon+SH0ES_cosmosis_likelihood.py'
    theory = method(cosmosis, 'PantheonLikelihood', 'extract_theory_points', {'np': np, 'names': SimpleNamespace(supernova_params='supernova_params')})
    # Execute official methods with the original nonsymmetric printed matrix's
    # full inverse, independently of our symmetric-average Cholesky.
    inv = linalg.solve(model.C_raw, np.eye(len(model.data)), assume_a='gen')
    cal = model.calibrator
    o = SimpleNamespace(mag=model.data.m_b_corr.to_numpy(), invcov=inv.copy(), is_calibrator=cal,
                        ceph_dist=model.data.CEPH_DIST.to_numpy(), use_abs_mag=True)
    p = SimpleNamespace(**vars(o)); p.invcov=inv.copy(); p.use_abs_mag=False; profile(p)
    symmetric_profile=SimpleNamespace(**vars(o))
    symmetric_profile.invcov=linalg.cho_solve((model.chol,True),np.eye(len(model.data)))
    symmetric_profile.use_abs_mag=False;profile(symmetric_profile)
    s = SimpleNamespace(zCMB=model.data.zHD.to_numpy(), zHEL=model.data.zHEL.to_numpy(), is_calibrator=cal,
                        cepheid_distance=model.data.CEPH_DIST.to_numpy(), x_section='distances', x_name='z',
                        y_section='distances', y_name='D_A', kind='linear')
    z = np.unique(model.data.zHD.to_numpy())
    records=[]
    for k in range(12):
        # Fixed synthetic distance curves, not fitted cosmologies or data-driven parameter estimates.
        da = (2800+50*k)*z/((1+z)*(1+(.12+.01*k)*z))
        at_rows=np.interp(model.z_hd_noncalibrator,z,da)
        mu=model.theory_mu(at_rows)
        M=-19.4+.025*k
        block={('distances','z'):z, ('distances','D_A'):da, ('supernova_params','M'):M}
        official=theory(s,block)
        r=model.data.m_b_corr.to_numpy()-mu-M
        rw=linalg.solve_triangular(model.chol,r,lower=True)
        chisq=float(rw@rw)
        cobaya_fixed=-2*raw(o,mu.copy()-25,Mb=M)
        cobaya_profile=-2*raw(p,mu.copy()-25)
        cobaya_symmetric_profile=-2*raw(symmetric_profile,mu.copy()-25)
        check=model.evaluate(at_rows)
        # Independently profile the actual raw quadratic. Its antisymmetric
        # precision part cancels; use the stable residual at the scalar optimum.
        v=model.data.m_b_corr.to_numpy()-mu;one=np.ones(len(v))
        raw_m=float((one@inv@v+v@inv@one)/(2*(one@inv@one)))
        raw_profile=float((v-raw_m)@inv@(v-raw_m))
        records.append({'case':k,'cosmosis_mean_max_abs_mag':float(np.max(abs(official-(mu+M)))),
                        'cobaya_fixed_chi2_abs_error':abs(cobaya_fixed-chisq),
                        'cobaya_flat_M_raw_covariance_chi2_abs_error':abs(cobaya_profile-check['chi2']),
                        'cobaya_flat_M_symmetric_covariance_chi2_abs_error':abs(cobaya_symmetric_profile-check['chi2']),
                        'profiled_raw_quadratic_chi2_abs_error':abs(raw_profile-check['chi2'])})
    write(work/'official-method-checks.json', records)
    assert max(x['cosmosis_mean_max_abs_mag'] for x in records)<1e-12
    assert max(x['cobaya_fixed_chi2_abs_error'] for x in records)<1e-7
    assert max(x['cobaya_flat_M_symmetric_covariance_chi2_abs_error'] for x in records)<1e-7
    assert max(x['profiled_raw_quadratic_chi2_abs_error'] for x in records)<1e-7
    # Independent integration identity including correlated host nuisance determinant.
    rng=np.random.default_rng(273601);errs=[]
    for k in range(16):
        n,h=9,3;U=rng.normal(size=(n,n)); C=U@U.T+np.eye(n)
        V=rng.normal(size=(h,h)); H=V@V.T+np.eye(h)
        J=np.eye(h)[rng.integers(h,size=n)];r=rng.normal(size=n)
        S=C+J@H@J.T;P=linalg.solve(C,np.eye(n),assume_a='pos');Q=linalg.solve(H,np.eye(h),assume_a='pos')
        B=np.column_stack((np.ones(n),J));K=B.T@P@B;K[1:,1:]+=Q
        b=B.T@P@r
        integrated=float(r@P@r-b@linalg.solve(K,b,assume_a='pos')+np.linalg.slogdet(C)[1]+np.linalg.slogdet(H)[1]+np.linalg.slogdet(K)[1])
        T=linalg.solve(S,np.eye(n),assume_a='pos');one=np.ones(n);aa=one@T@one
        direct=float(r@T@r-(one@T@r)**2/aa+np.linalg.slogdet(S)[1]+np.log(aa))
        shifted=r+53;shift=float(shifted@T@shifted-(one@T@shifted)**2/aa+np.linalg.slogdet(S)[1]+np.log(aa))
        errs.append({'host_integral_error':abs(integrated-direct),'offset_error':abs(shift-direct)})
    assert max(x['host_integral_error'] for x in errs)<1e-10
    assert max(x['offset_error'] for x in errs)<1e-9
    return {'official_method_cases':records,'synthetic_host_integration_cases':16,
            'literal_default_Cobaya_raw_covariance_projection_agrees_at_1e_7':False,
            'raw_projection_mismatch_explanation':'Printed C has <=3e-8 antisymmetry; Cobaya rank-one projection assumes symmetric inverse. Symmetric Gaussian and explicit-M raw quadratic agree; literal raw projection mismatch retained.',
            'host_integral_max_abs_error':max(x['host_integral_error'] for x in errs),
            'global_53mag_offset_max_abs_error':max(x['offset_error'] for x in errs),
            'Cobaya_version':cobaya.__version__, 'source_sha256':{rel(q):sha(q) for q in [cp,pp,cp.with_suffix('.yaml')]}}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--acquire',action='store_true');ap.add_argument('--work',type=Path,default=WORK);ap.add_argument('--output',type=Path,default=OUT)
    args=ap.parse_args();work=args.work;records=acquire(work,args.acquire)
    factor_path=ROOT/'studies/unified_cosmology/results/distance_ladder/cepheid-factor.json'
    full_path=LADDER/'gls-result.npz'
    factor=json.loads(factor_path.read_text());hosts=factor['host_order'];H=np.array(factor['covariance_distance_modulus']);mu=np.array(factor['mean_distance_modulus'])
    dependencies_before={rel(q):sha(q) for q in [factor_path,full_path,LADDER/'Riess2022v3.pdf']}
    full=np.load(full_path)['q'][:37]
    table,textpath=published_hosts(work);by_SN={v['paper_SN'].lower():v for v in table}
    data=pd.read_csv(work/'Pantheon+SH0ES.dat',sep=r'\s+');cal=data.IS_CALIBRATOR.to_numpy().astype(bool);cc=np.flatnonzero(cal)
    ledger=[]
    for i in cc:
        row=data.iloc[i];physical=SN_ALIASES.get(row.CID,row.CID).lower();entry=by_SN[physical]
        ledger.append({'original_row_zero_based':int(i),'CID':row.CID,'IDSURVEY':int(row.IDSURVEY),'physical_SN':physical,
                       **entry,'host_index':hosts.index(entry['host']),'CEPH_DIST':float(row.CEPH_DIST),
                       'SN_RA':float(row.RA),'SN_DEC':float(row.DEC),'PKMJD':float(row.PKMJD),
                       'host_coordinates_missing':bool(row.HOST_RA==-999 and row.HOST_DEC==-999)})
    assert len(ledger)==77 and len({v['physical_SN'] for v in ledger})==42 and {v['host'] for v in ledger}==set(hosts)
    alias=data[data.CID.isin(['2005df','2005df_ANU'])];assert len(alias)==2 and alias.RA.nunique()==alias.DEC.nunique()==1 and np.ptp(alias.PKMJD)<.2
    idx=np.array([v['host_index'] for v in ledger]);sn=np.array([v['physical_SN'] for v in ledger])
    shared=[]
    for j,h in enumerate(hosts):
        obs=[v['CEPH_DIST'] for v in ledger if v['host']==h];assert len(set(obs))==1
        shared.append({'host':h,'released_mu':obs[0],'SN_free_mu':float(mu[j]),'full_ladder_mu':float(full[j]),
                       'release_minus_SN_free_mag':float(obs[0]-mu[j]),'release_minus_full_ladder_mag':float(obs[0]-full[j])})
    C_raw=covariance(work/'Pantheon+SH0ES_STAT+SYS.cov');C=(C_raw+C_raw.T)/2
    S=covariance(work/'Pantheon+SH0ES_STATONLY.cov');sc=S[np.ix_(cc,cc)]
    pairs=[]
    for i in range(37):
        for j in range(i+1,37):
            block=sc[np.ix_(idx==i,idx==j)];assert np.ptp(block)==0
            pairs.append({'host1':hosts[i],'host2':hosts[j],'STATONLY_shared_cov_mag2':float(block[0,0]),'SN_free_cov_mag2':float(H[i,j]),'difference_mag2':float(block[0,0]-H[i,j])})
    siblings=[]
    for j,h in enumerate(hosts):
        ii=np.flatnonzero(idx==j);vv=[sc[a,b] for a in ii for b in ii if sn[a]!=sn[b]]
        if vv:
            assert np.ptp(vv)==0
            siblings.append({'host':h,'physical_SN':sorted(set(sn[ii])),'STATONLY_sibling_cov_mag2':float(vv[0]),'SN_free_host_variance_mag2':float(H[j,j]),'ratio':float(vv[0]/H[j,j])})
    selected=(data.zHD.to_numpy()>.01)|cal;hf=np.flatnonzero(selected&~cal)
    cross=C[np.ix_(cc,hf)];corr=cross/np.sqrt(np.diag(C)[cc,None]*np.diag(C)[hf][None,:])
    minus=C.copy();minus[np.ix_(cc,cc)]-=H[np.ix_(idx,idx)]
    model=ReleasedCalibration(work);checks=numerical_checks(model,work)
    for name,values in [('calibrator-host-ledger.csv',ledger),('host-distance-comparison.csv',shared),('distinct-host-covariance.csv',pairs)]:
        pd.DataFrame(values).to_csv(work/name,index=False)
    delta=np.array([r['release_minus_SN_free_mag'] for r in shared]);d_full=np.array([r['release_minus_full_ladder_mag'] for r in shared]);d_cov=np.array([r['difference_mag2'] for r in pairs])
    unusual=np.flatnonzero(np.any(S[np.ix_(cc,np.flatnonzero(~cal))]!=0,axis=0));urows=np.flatnonzero(~cal)[unusual]
    result={'status':'passed_released_Gaussian_interface_checks_no_cosmological_fit','official_commit':COMMIT,
            'released_target_usable':True,'independent_host_factor_replacement_certified':False,
            'selection':{'all_rows':len(data),'selected_rows':int(selected.sum()),'calibrator_rows':int(cal.sum()),'noncalibrator_rows':len(hf),'raw_calibrator_CID_strings':data[cal].CID.nunique(),'physical_calibrator_SN':42,'calibrator_hosts':37,'calibrators_above_z_cut':int(np.sum(cal&(data.zHD.to_numpy()>.01))),'calibrator_host_coordinates_all_missing':all(v['host_coordinates_missing'] for v in ledger)},
            'mapping':{'source':'Riess2022v3 Table6, printed page30; SN identities only, never brightness matching','host_aliases':HOST_ALIASES,'SN_aliases':SN_ALIASES,'SN_2005df_alias_geometry_identical':True,'exact_ladder_SN_measurement_row_mapping':False},
            'distance_comparison':{'release_minus_SN_free_min_mag':float(delta.min()),'max_mag':float(delta.max()),'median_mag':float(np.median(delta)),'rms_mag':float(np.sqrt(np.mean(delta**2))),'rms_after_common_offset_mag':float(np.std(delta)),'release_minus_full_ladder_rms_mag':float(np.sqrt(np.mean(d_full**2))),'release_minus_full_ladder_max_abs_mag':float(abs(d_full).max()),'offset_cause_established':False},
            'covariance':{'dimension':len(C),'raw_antisymmetry_max_mag2':float(abs(C_raw-C_raw.T).max()),'raw_antisymmetric_ordered_entries':int(np.count_nonzero(C_raw-C_raw.T)),'Gaussian_covariance_rule':'(C_raw+C_raw.T)/2','total_min_eigenvalue_mag2':float(linalg.eigh(C,eigvals_only=True,subset_by_index=[0,0])[0]),'selected_min_eigenvalue_mag2':float(linalg.eigh(model.C,eigvals_only=True,subset_by_index=[0,0])[0]),'calibrator_noncalibrator_nonzero':int(np.count_nonzero(cross)),'crossblock_elements':cross.size,'crossblock_max_abs_correlation':float(abs(corr).max()),'STATONLY_distinct_host_pairs':len(pairs),'nonzero_distinct_host_pairs':sum(p['STATONLY_shared_cov_mag2']!=0 for p in pairs),'host_pair_block_spread_max':0.0,'zero_crosscovariance_host':['N1365'],'distinct_host_difference_from_SN_free_rms_mag2':float(np.sqrt(np.mean(d_cov**2))),'distinct_host_difference_max_abs_mag2':float(abs(d_cov).max()),'sibling_host_entries':siblings,'host_diagonals_not_identifiable_from_distinct_SN_pairs':37-len(siblings),'STATONLY_selected_calibrator_noncalibrator_nonzero':int(np.count_nonzero(S[np.ix_(cc,hf)])),'STATONLY_correlated_excluded_noncalibrators':data.iloc[urows][['CID','IS_CALIBRATOR','zHD','CEPH_DIST']].to_dict('records'),'total_minus_our_candidate_host_cov_min_eigenvalue_mag2':float(linalg.eigh(minus,eigvals_only=True,subset_by_index=[0,0])[0]),'positive_residual_certifies_component_identity':False},
            'likelihood':{'observed_column':'m_b_corr','calibrator_mean':'CEPH_DIST + M','other_mean':'5log10[(1+zHD)(1+zHEL)D_A(zHD)/Mpc]+25+M','covariance':'selected symmetric-average STAT+SYS; no extra errors','M':'one improper flat measure dM in magnitudes, analytically integrated','loglike':'-0.5[chi2_projected+logdetC+log(1^T C^-1 1)+(N-1)log(2pi)]','extra_SH0ES_or_Cepheid_factor':False,'logdetC':model.logdet_covariance,'flat_M_precision':model.flat_m_precision,'normalized_evidence_claim':False},
            'validation':checks,'input_sources':records,
            'dependencies_sha256':{rel(q):sha(q) for q in [factor_path,full_path,LADDER/'Riess2022v3.pdf',textpath]},
            'source_sha256':{rel(Path(__file__)):sha(__file__)},
            'output_sha256':{rel(work/n):sha(work/n) for n in ['calibrator-host-ledger.csv','host-distance-comparison.csv','distinct-host-covariance.csv']},
            'limitations':['Published release defines a conditional Gaussian target; these checks do not validate its original covariance construction.','No exact embedded Cepheid covariance or construction program was identified in the pinned official release.','STATONLY contains host-shared terms and is not a verified SN-only replacement.','Differences in means/covariances versus the independently reconstructed ladder are preserved; origin unresolved.','Use as an alternative to Dovekie/Pantheon-only, never multiply overlapping SN compilations or add SH0ES H0 as independent data.']}
    # No mutable upstream or local producer input may change during this audit.
    for path,record in records.items():assert sha(ROOT/path)==record['sha256']
    for path,digest in dependencies_before.items():assert sha(ROOT/path)==digest
    write(args.output,result)
    print(json.dumps({k:result[k] for k in ['status','selection','distance_comparison','covariance','validation']},indent=2))


if __name__=='__main__':
    main()
