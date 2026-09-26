"""Metadata-only BayeSN distance-identification design; no SN outcome fits.

Run with phase2/env-official/bin/python. No package installation or BayeSN import.
"""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def main():
    membership_path = ROOT / 'runs/research_2026_09_26/raisin_differential/frozen-membership.csv'
    rows = list(csv.DictReader(membership_path.open()))
    assert len(rows) == 79 and len({r['CID'] for r in rows}) == 79
    # Frozen without using dust, fitted distance, colour, or residual outcomes.
    pilot = []
    for survey in ['CSP', 'PS1MD', 'DES']:
        group = sorted((r for r in rows if r['survey'] == survey), key=lambda r: (float(r['zHEL']), r['CID']))
        assert len(group) == {'CSP': 42, 'PS1MD': 19, 'DES': 18}[survey]
        for numerator in [1, 2]:
            idx = numerator * (len(group) - 1) // 3
            r = dict(group[idx])
            r.update(rank_1based=idx + 1, survey_N=len(group), quantile_numerator=numerator,
                     photometry_sha256=digest(ROOT / r['photometry_path']))
            pilot.append(r)
    with (OUT / 'pilot-membership.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(pilot[0])); w.writeheader(); w.writerows(pilot)
    protocol = {
        'status': 'design frozen before new BayeSN likelihood or posterior evaluation; no long fits run',
        'membership_source': str(membership_path.relative_to(ROOT)),
        'membership_sha256': digest(membership_path),
        'pilot_rule': 'For each survey CSP/PS1MD/DES, sort by (zHEL,CID); zero-based floor((N-1)/3), floor(2*(N-1)/3).',
        'pilot_ids': [r['CID'] for r in pilot],
        'full_estimand_membership': {'low': '42 CSP', 'high': '19 PS1MD +18 DES; same79 released objects'},
        'model': {'name': 'M20', 'numpyro_code_commit': '08eef9e54188f601506ef9f7f47fc65f9f52a946',
                  'model_file_commit': '1d77e0be8203f09fff0a6b725abee15011a4d7b9',
                  'fixed': ['W0', 'W1', 'L_Sigma_epsilon', 'M0=-19.5mag', 'sigma0=0.088mag'],
                  'per_SN': ['theta~Normal(0,1)', 'epsilon_white[42]~Normal(0,I)', 'AV~Exponential(scale=tauA)',
                             'RV~Normal(muR,sigmaR) normalized on RV>=0.5', 'D=mu+deltaM'],
                  'group': ['muR~Uniform(1,5)', 'sigmaR~HalfNormal(2)', 'tauA~HalfCauchy(1)'],
                  'groups': ['CSP', 'PS1MD+DES'],
                  'common_RV_secondary': 'RV per group~Uniform(1,6), no sigmaR; must be shared outside per-SN plate'},
        'distance_arms': {
            'A': 'mu~Normal(mu_ext,sigma_ext); H0=73.24 km/s/Mpc, Om0=.28, flat LCDM; retain grey deltaM~N(0,.088^2)',
            'B': 'independent mu~Uniform(20,50)mag for every SN, retain same grey deltaM; no redshift-conditioned amplitude prior',
            'B_sensitivity': ['Uniform(15,55)mag', 'Uniform(25,50)mag'],
            'B_measure': 'Integrate D with the exact convolution of Uniform(mu) and Normal(deltaM); never flat a',
            'frame_gate': 'Separate zHEL spectral/time transform from corrected zHD external distance. Freeze exact mu_ext/error vectors after preprocessing closure.'},
        'unchanged_between_arms': ['physical objects', 'rows and bandpasses', 'FLUXCAL/error/covariance', 'MW extinction',
                                  'trained model and intrinsic priors', 'dust hierarchy and hyperpriors', 'phase/tmax treatment',
                                  'selection conditioning', 'sampler settings'],
        'primary_estimands': ['posterior of muR_high-muR_low', 'posterior of tauA_high-tauA_low (mag)',
                               'paired arm change of each posterior mean, SD, equal-tail68/95% intervals; MCSE'],
        'secondary_diagnostics': ['per-SN AV and RV posterior changes', 'posterior D (identified amplitude), not independently identified mu and deltaM',
                                  'RV likelihood information as AV tends to0', 'RV lower-bound sensitivity reported, no outcome-driven retuning'],
        'no_claims': ['no new posterior in this task', 'no percentage of colour information from a posterior variance ratio',
                      'no cosmology fit', 'no claim broad prior identifies arbitrary grey evolution',
                      'no intrinsic population evolution inference without selection normalization and training transport checks'],
        'amplitude_integration_gate': {'fixed_measurement_C': True, 'log_integral_absolute_error': 1e-6,
                                       'gradient_scaled_error': 1e-5, 'checks': ['adaptive reference vs doubled quadrature resolution',
                                        'direct-D sampling subset', 'include low-SNR a->0 tails', 'all normalization factors and Jacobians',
                                        'b<=0, negative flux, and remote prior boundaries']},
        'forward_gate': ['same M20 assets; source-linked filter/zero-point/epoch/peak mapping',
                         'legacy-vs-modern and refined broadband integral discrepancy<0.1 quoted error per retained row',
                         'reference-spectrum calibration closure<1mmag; preserve any known source offsets rather than rescale to fit',
                         'no flux S/N cut, clamp, SALT covariance, or error floor unless exact paper preprocessing requires it',
                         'tmax behavior must be established or explicitly labelled reconstruction with same behavior in both arms'],
        'recovery_gate': {'injected_AV_mag': [0.0, 0.3, 0.8], 'injected_RV': [2.0, 3.1],
                          'theta': 0.0, 'epsilon': 'zero for deterministic closure; then draws from fixed trained prior',
                          'grey_stress_mag': [0.0, 0.15],
                          'tests': ['zero-noise likelihood closure', 'amplitude invariance after profiling under multiplicative flux rescaling',
                                    'RV unidentifiable at AV=0 in likelihood-only diagnostic',
                                    'ten fixed-seed synthetic six-SN datasets screen bias/coverage and prior-driven recovery; not coverage certification']},
        'convergence_gate': {'chains': 4, 'initial_warmup_per_chain': 1000, 'initial_samples_per_chain': 1000,
                             'Rhat_max': 1.01, 'bulk_ESS_min': 400, 'tail_ESS_min': 200, 'divergences': 0,
                             'E_BFMI_min': 0.3, 'max_treedepth_fraction_max': 0.01,
                             'mean_MCSE_over_posterior_SD_max': 0.05,
                             'seed_replication': 'differences compatible with combined MCSE; no tuning priors to induce convergence',
                             'distance_boundary_probability_max': 0.001,
                             'prior_sensitivity': 'Report all changes; boundary/width dependence fails the claim that broad arm removed distance information'},
        'resource_caps': {'CPU': 2, 'GPU': 'not required or benchmarked', 'preflight_wall_minutes': 15,
                          'compile_plus100_gradient_evaluations_minutes': 10,
                          'six_object_two_arm_pilot_wall_minutes': 60,
                          'stop_rule': 'No silent relaxation of gates. Stop with timing/checkpoint if budget exceeded; extrapolate ESS/hour before larger run.'},
        'blockers_before_real_posterior': ['exact 2024 Stan population executable/config/chains not located',
                                         '2024 epoch masks/filter calibration/peak/error-floor behavior not linked',
                                         'exact mu_ext/sigma_ext inputs not exported as paper execution assets',
                                         'NumPyro environment missing direct imports; no install performed'],
    }
    save('protocol.json', protocol)

    # M20 paper Table1 metadata only; no distances or fit outcomes are retained.
    paper = (OUT / 'Mandel-M20-paper.txt').read_text()
    names = re.findall(r'^\s+(SN\d{4}[A-Za-z]+)\s+(CfA|CSP|K04a|K04b|K03|St07|L09|K07|P08)\s+', paper, re.M)
    assert len(names) == 79 and len(set(n for n, _ in names)) == 79, (len(names), names)
    roster = {n[2:].lower(): (n, source) for n, source in names}
    matches = []
    for row in rows:
        match = roster.get(row['CID'].lower())
        if match:
            matches.append({'CID': row['CID'], 'survey': row['survey'], 'paper_name': match[0],
                            'training_data_source': match[1], 'match': 'unique official SN designation, case normalized; no coordinates in Table1'})
    with (OUT / 'M20-training-name-overlap.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(matches[0])); w.writeheader(); w.writerows(matches)
    save('M20-training-overlap.json', {'paper_table1_objects': 79, 'raisin_matches': len(matches),
                                     'survey_counts': {s: sum(m['survey'] == s for m in matches) for s in ['CSP','PS1MD','DES']},
                                     'status': 'paper table input proxy; exact executed M20 training roster/photometry not obtained',
                                     'paper_sha256': digest(OUT / 'Mandel-M20-paper.pdf')})

    pkgs = ['numpy', 'scipy', 'pandas', 'matplotlib', 'h5py', 'sncosmo', 'arviz', 'astropy', 'extinction',
            'jax', 'jaxlib', 'numpyro', 'ruamel.yaml', 'tqdm', 'bayesn', 'pystan', 'cmdstanpy']
    source = 'import json,importlib.metadata as m,sys,os\nr={}\n'
    source += f'for p in {pkgs!r}:\n try:r[p]=m.version(p)\n except m.PackageNotFoundError:r[p]=None\n'
    source += 'print(json.dumps(dict(executable=sys.executable,python=sys.version,cpu_count=os.cpu_count(),packages=r)))'
    runtime = []
    for env in ['phase2/env-official/bin/python', 'phase2/hierarchy/.venv/bin/python', 'phase2/classification/env/bin/python']:
        if (ROOT / env).exists():
            r = subprocess.run([str(ROOT / env), '-c', source], capture_output=True, text=True, check=True)
            runtime.append(json.loads(r.stdout))
    save('runtime-inventory.json', {'environments': runtime, 'packages_installed_this_task': [],
                                    'BayeSN_imported': False, 'CmdStan_home_exists': (Path.home()/'.cmdstan').exists()})

    # Algebra test only: deliberately artificial vectors; never presented as SN fit.
    f = np.array([1., 0.7, 0.2]); y = np.array([1.2, 0.5, -0.1]); C = np.array([[.04,.005,0],[.005,.09,0],[0,0,.16]])
    w = np.linalg.inv(C); q=float(f@w@f); b=float(f@w@y); c=float(y@w@y)
    ahat=max(0.,b/q); profile=float((y-ahat*f)@w@(y-ahat*f))
    lo, hi = .1, 3.; m=b/q
    analytic=np.exp(-.5*(c-b*b/q))*np.sqrt(2*np.pi/q)*(ndtr(np.sqrt(q)*(hi-m))-ndtr(np.sqrt(q)*(lo-m)))/(hi-lo)
    numeric=quad(lambda a: np.exp(-.5*(y-a*f)@w@(y-a*f))/(hi-lo),lo,hi,epsabs=1e-13)[0]
    assert abs(profile-(c-max(b,0.)**2/q)) < 1e-12
    assert abs(analytic-numeric) < 1e-12
    K=.4*np.log(10); lowmu,highmu=20.,50.; sigma0=.088; ref=35.
    density=lambda D:(ndtr((D-lowmu)/sigma0)-ndtr((D-highmu)/sigma0))/(highmu-lowmu)
    norm=quad(density,lowmu-10*sigma0,highmu+10*sigma0,points=[lowmu,highmu],epsabs=1e-12)[0]
    assert abs(norm-1)<1e-12
    like=lambda D: np.exp(-.5*(y-np.exp(-K*(D-ref))*f)@w@(y-np.exp(-K*(D-ref))*f))
    direct=quad(lambda D: like(D)*density(D),lowmu-10*sigma0,highmu+10*sigma0,
                points=[lowmu,ref,highmu],epsabs=1e-12)[0]
    nested=quad(lambda d: quad(lambda mu: like(mu+d)/(highmu-lowmu),lowmu,highmu,
                               points=[ref-d],epsabs=1e-11)[0]*np.exp(-.5*(d/sigma0)**2)/(np.sqrt(2*np.pi)*sigma0),
                -8*sigma0,8*sigma0,epsabs=1e-11)[0]
    assert abs(direct-nested)<1e-10
    save('amplitude-algebra-check.json', {'status':'pure Gaussian-vector algebra; no SN data/model substitution',
                                         'profile_chi2':profile, 'flat_amplitude_integral_error':abs(analytic-numeric),
                                         'convolved_uniform_mu_density_normalization':norm,
                                         'uniform_mu_plus_grey_integral_error':abs(direct-nested),
                                         'improper_flat_mu_tail_likelihood_limit':float(np.exp(-.5*c))})
    print(json.dumps({'pilot': protocol['pilot_ids'], 'training_overlap': len(matches),
                      'protocol_sha256': digest(OUT/'protocol.json'), 'algebra_checks':'pass'}, indent=2))


if __name__ == '__main__':
    main()
