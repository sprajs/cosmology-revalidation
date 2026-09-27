#!/usr/bin/env python3
"""Independent equation, optimization, coverage and CPU throughput checks."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json, time, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from model import Library, Grid, age_bound, feasible_fit, fieller

OUT = Path(__file__).resolve().parents[2]/'results/physical_ages'


def primal_check(matrix, flux, error, ages, threshold, bound, maximize, case):
    """Check a conic endpoint with independently initialized original-mass fits.

    The original column-scaled fractional formulation is retained as a diagnostic.
    The adopted SLSQP coordinates are mass fractions (a simplex) and total formed
    mass divided by the independent BVLS best-fit mass.  These solve the original
    flux problem, not the Charnes--Cooper conic problem.  No conic witness or
    injected true mixture is used to construct any starting point.
    """
    sign = -1 if maximize else 1
    initial, minimum_chi2 = feasible_fit(matrix, flux, error)
    scales = np.linalg.norm(matrix/error[:, None], axis=0)
    ap = matrix/scales
    numerator, denominator = ages/scales, 1/scales

    def old_objective(v):
        return sign*(numerator@v)/max(denominator@v, 1e-100)

    def old_gradient(v):
        total = max(denominator@v, 1e-100)
        return sign*(numerator*total-(numerator@v)*denominator)/total**2

    old = minimize(old_objective, initial*scales, method='SLSQP', jac=old_gradient,
        bounds=[(0, None)]*len(ages), constraints=[{
            'type': 'ineq',
            'fun': lambda v: (threshold-np.sum(((ap@v-flux)/error)**2))/threshold,
            'jac': lambda v: -2*ap.T@((ap@v-flux)/error**2)/threshold}],
        options={'ftol': 1e-12, 'maxiter': 3000})
    old_mass = old.x/scales
    old_age = float(ages@old_mass/old_mass.sum())
    legacy = dict(case=case, extremum='maximum' if maximize else 'minimum',
        success=bool(old.success), status=int(old.status), message=str(old.message),
        iterations=int(old.nit), age_Gyr=old_age,
        difference_gyr=old_age-bound['age'],
        chi2=float(np.sum(((matrix@old_mass-flux)/error)**2)),
        fractional_objective_gradient_norm=float(np.linalg.norm(old_gradient(old.x))))

    n = len(ages)
    mass_scale = float(initial.sum())
    assert mass_scale > 0 and minimum_chi2 < threshold
    aw, yw = matrix*mass_scale/error[:, None], flux/error

    def objective(x):
        return sign*(ages@x[:n])

    def gradient(x):
        return np.r_[sign*ages, 0.]

    def constraint(x):
        return (threshold-np.sum((x[n]*(aw@x[:n])-yw)**2))/threshold

    def constraint_gradient(x):
        residual = x[n]*(aw@x[:n])-yw
        return -2/threshold*np.r_[x[n]*(aw.T@residual), (aw@x[:n])@residual]

    def coordinates(mass):
        return np.r_[mass/mass.sum(), mass.sum()/mass_scale]

    # Six fixed start types; their mixture targets do not depend on the endpoint.
    # A bisection in original masses puts each alternative strictly inside the
    # same ellipsoid as BVLS.  Convexity makes this an independent feasibility
    # construction, not an optimization using a conic answer.
    random = np.random.default_rng(20260927+case)
    targets = [('uniform', np.ones(n)/n), ('youngest', np.eye(n)[0]),
               ('oldest', np.eye(n)[-1]),
               ('dirichlet_sparse', random.dirichlet(np.full(n, .3))),
               ('dirichlet_dense', random.dirichlet(np.ones(n)))]
    starts = [('bvls', initial)]
    interior_limit = minimum_chi2 + .8*(threshold-minimum_chi2)
    for name, fraction in targets:
        spectrum = matrix@fraction/error
        amplitude = max(float(spectrum@(flux/error)/(spectrum@spectrum)), 1e-100)
        target = fraction*amplitude
        low, high = 0., .9
        for _ in range(55):
            middle = (low+high)/2
            mass = (1-middle)*initial+middle*target
            if np.sum(((matrix@mass-flux)/error)**2) <= interior_limit:
                low = middle
            else:
                high = middle
        starts.append((name, (1-low)*initial+low*target))

    x0 = coordinates(initial)
    step = 1e-6
    finite = np.array([(constraint(x0+step*d)-constraint(x0-step*d))/(2*step)
                       for d in np.eye(n+1)])
    derivative_error = float(np.max(abs(finite-constraint_gradient(x0)) /
                                   np.maximum(1., abs(constraint_gradient(x0)))))
    attempts = []
    for name, mass in starts:
        result = minimize(objective, coordinates(mass), method='SLSQP', jac=gradient,
            bounds=[(0, 1)]*n+[(0, None)], constraints=[
                {'type': 'eq', 'fun': lambda x: np.sum(x[:n])-1,
                 'jac': lambda x: np.r_[np.ones(n), 0.]},
                {'type': 'ineq', 'fun': constraint, 'jac': constraint_gradient}],
            options={'ftol': 1e-12, 'maxiter': 3000})
        recovered = result.x[:n]*result.x[n]*mass_scale
        total = float(recovered.sum())
        val = float(ages@recovered/total) if total > 0 else None
        residual = float(np.sum(((matrix@recovered-flux)/error)**2))
        simplex_error = float(abs(result.x[:n].sum()-1))
        feasible = (total > 0 and residual <= threshold+1e-7 and
                    simplex_error <= 1e-9 and float(result.x.min()) >= -1e-10)
        attempts.append(dict(start=name, uses_conic_witness=False,
            success=bool(result.success), status=int(result.status),
            message=str(result.message), iterations=int(result.nit),
            age_Gyr=val, difference_gyr=None if val is None else val-bound['age'],
            chi2=residual, simplex_error=simplex_error,
            minimum_coordinate=float(result.x.min()), feasible=bool(feasible),
            start_chi2=float(np.sum(((matrix@mass-flux)/error)**2))))
    candidates = [x for x in attempts if x['success'] and x['feasible']]
    best = min(candidates, key=lambda x: sign*x['age_Gyr']) if candidates else None
    # This is stricter than the previous 2e-4 Gyr endpoint test.  Exit flags alone
    # never certify a bound; original-coordinate feasibility is checked first.
    passed = (best is not None and abs(best['difference_gyr']) < 1e-5 and
              derivative_error < 1e-6)
    return legacy, dict(case=case, extremum=legacy['extremum'], passed=bool(passed),
        conic_outward_bound_Gyr=bound['age'], conic_primal_Gyr=bound['primal_value'],
        conic_duality_gap=bound['duality_gap'], independent_best=best,
        constraint_gradient_max_relative_error=derivative_error, attempts=attempts)


def main():
    start=time.perf_counter(); lib=Library(); grid=Grid(lib); rng=np.random.default_rng(27192026)
    native=[]; dense=[]
    for iz,z in enumerate(lib.data['native_redshifts']):
        reference=lib.data['native_sdss_mags'][:,:,iz,:].transpose(0,2,1)
        reference=reference-reference[:,2:3,:]
        for mode, destination in [('native',native),('dense',dense)]:
            predicted=-2.5*np.log10(lib.photometry(lib.spectra,z,0,integration=mode))
            predicted-=predicted[:,2:3,:]
            destination.extend((predicted-reference).ravel().tolist())
    assert max(abs(np.array(native)))<1e-4
    assert max(abs(np.array(dense)))<.003
    quadrature_checks=[];old_des_checks=[]
    for family,path in [('sdss',None),('des',Path(__file__).resolve().parents[4]/'.work/host-transport/filters/des-deep-responses.npz')]:
        instrument=Library(response_path=path)
        for z in [0.,.1,.3,.6,1.05]:
            for tau in [0.,1.]:
                precise=instrument.photometry(instrument.spectra,z,.03,tau=tau)
                dense_flux=instrument.photometry(instrument.spectra,z,.03,tau=tau,integration='dense')
                delta=-2.5*np.log10(precise/dense_flux);delta-=delta[:,2:3,:]
                quadrature_checks.append(dict(family=family,z=z,tau=tau,max_colour_difference_mag=float(abs(delta).max())))
                if family=='des':
                    old=instrument.photometry(instrument.spectra,z,.03,tau=tau,integration='native')
                    delta=-2.5*np.log10(old/precise);delta-=delta[:,2:3,:]
                    old_des_checks.append(float(abs(delta).max()))
        flat=Grid(instrument).kernels(.3,0).sum(axis=1)
        assert np.allclose(flat,1.3,rtol=0,atol=1e-13)
    assert max(x['max_colour_difference_mag'] for x in quadrature_checks)<3e-6
    # Compare optimized native kernels with independent interpolation/trapezoids.
    matrix_checks=[]
    for z in [.001,.12,.374]:
        ages, matrices, dn, mass=grid.at_redshift(z,.035)
        ia=np.argmin(abs(ages-1)); norm=None
        for j in [0,7,23,67]:
            pars=grid.parameters[j]; m=np.where(lib.data['metallicity']==pars['metallicity'])[0][0]
            aa,spectra,_=lib.at_redshift(z)
            reference=lib.photometry(spectra[m],z,.035,pars['tau'],pars['power'])
            reference/=reference[2,ia]; candidate=matrices[j]/matrices[j,2,ia]
            matrix_checks.append(float(np.max(abs(candidate/reference-1))))
        assert ages.max()<=13.5 and mass.min()>0 and mass.max()<=1
    assert max(matrix_checks)<2e-12
    # Distinct nonnegative physical mixtures; no model fitting used to make truth.
    ages, models, _, _=grid.at_redshift(.08,.023)
    threshold=float(chi2.ppf(.95,5)); coverage=[]; region=[]; solver_fail=0
    nondetection=[age_bound(models[0],np.full(5,-.1),np.ones(5),ages,threshold,m) for m in [False,True]]
    assert nondetection[0]['age']==min(ages) and nondetection[1]['age']==max(ages)
    negative_debug=[]
    assert age_bound(models[0],np.full(5,-5.),np.ones(5),ages,threshold,diagnostics=negative_debug) is None
    assert negative_debug[-1]['status']=='PrimalInfeasible'
    independent=[]; legacy=[]; widths=[]
    for k in range(500):
        a=models[int(rng.integers(len(models)))]; weight=rng.dirichlet(np.ones(len(ages))*.3)
        mean=a@weight; error=np.maximum(mean*.04,mean.max()*.008)
        y=mean+rng.normal(size=5)*error; truth=float(ages@weight)
        region.append(bool(np.sum(((y-mean)/error)**2)<=threshold))
        lo=age_bound(a,y,error,ages,threshold)
        hi=age_bound(a,y,error,ages,threshold,True)
        if lo is None or hi is None:
            fit, chisq=feasible_fit(a,y,error)
            solver_fail+=int(chisq<=threshold-1e-6)
            coverage.append(False); continue
        coverage.append(lo['age']-1e-5<=truth<=hi['age']+1e-5)
        widths.append(hi['age']-lo['age'])
        if k<12:
            for maximize, bound in [(False,lo),(True,hi)]:
                old, checked=primal_check(a,y,error,ages,threshold,bound,maximize,k)
                legacy.append(old); independent.append(checked)
    assert solver_fail==0
    assert all(c or not r for c,r in zip(coverage,region))
    good=[x['independent_best'] for x in independent if x['independent_best'] is not None]
    independent_pass=len(independent)==24 and all(x['passed'] for x in independent)
    print('Independent primal validation:',json.dumps(dict(passed=independent_pass,
        verified_endpoints=sum(x['passed'] for x in independent),
        legacy_single_start=legacy)),flush=True)
    # Correlated Gaussian red/blue means: Fieller retains covariance explicitly.
    cov=np.array([[.04,.008],[.008,.01]])
    draws=rng.multivariate_normal([1.6,1.],cov,size=20000)
    cover=[]; unbounded=0
    for red,blue in draws:
        interval=fieller(red,blue,cov[0,0],cov[1,1],cov[0,1])
        unbounded+=int(interval is None)
        if interval is not None: cover.append(interval[0]<=1.6<=interval[1])
    # Timings include real conversions; do not infer GPU advantage from hardware.
    times={}
    for name,operation in [('literal_photometry_68_models',lambda:lib.photometry(grid.spectra,.1,.02)),('compiled_BLAS_photometry_68_models',lambda:grid.at_redshift(.1,.02))]:
        operation(); t=time.perf_counter()
        for _ in range(10):operation()
        times[name]=(time.perf_counter()-t)/10
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),passed=independent_pass,
        native_sampling_AB_colour_max_difference_mag=float(max(abs(np.array(native)))),dense_vs_native_FSPS_max_colour_difference_mag=float(max(abs(np.array(dense)))),
        final_gauss_vs_dense_checks=quadrature_checks,discarded_DES_native_sampling_max_colour_error_mag=max(old_des_checks),
        final_quadrature='Four-point Gauss-Legendre on union of passband and redshifted stellar-spectrum knots. Native FSPS discretization agreement is a separate check, not the adopted integration rule.',
        nondetection='Weak signed nondetection has full age support; incompatible strongly negative data are retained as an explicit model failure.',
        vectorized_operator_max_relative_difference=max(matrix_checks),
        conic_coverage={'draws':500,'seed':27192026,'coverage':float(np.mean(coverage)),'input_region_coverage':float(np.mean(region)),
            'solver_failures_on_feasible_inputs':solver_fail,'median_width_Gyr':float(np.median(widths)),
            'meaning':'Conditional on included SSP basis, fixed true nuisance cell and Gaussian diagonal flux errors. Not coverage under incorrect stellar physics.'},
        independent_primal={'endpoints':len(independent),
            'attempts':sum(len(x['attempts']) for x in independent),
            'converged':sum(a['success'] for x in independent for a in x['attempts']),
            'verified_endpoints':sum(x['passed'] for x in independent),
            'maximum_best_difference_Gyr':max(abs(x['difference_gyr']) for x in good) if good else None,
            'criteria':{'age_difference_Gyr':1e-5,'maximum_chi2_excess':1e-7,
                'maximum_simplex_error':1e-9,'minimum_coordinate':-1e-10,
                'require_successful_independent_start_for_every_endpoint':True},
            'method':'SLSQP original flux likelihood; simplex formed-mass fractions and independently scaled amplitude. Six predetermined start types per endpoint, constructed from BVLS without conic witnesses or injected truth. Best feasible successful extremum checked against outward conic bound. Failed and suboptimal attempts retained.',
            'multistart_seed_base':20260927,'records':independent,
            'legacy_single_start':{'method':'Original column-scaled masses with fractional age objective; diagnostic only, never used to validate new endpoints.',
                'records':legacy,'converged':sum(x['success'] for x in legacy),
                'maximum_converged_difference_Gyr':max(abs(x['difference_gyr']) for x in legacy if x['success'])}},
        fieller={'draws':20000,'bounded':len(cover),'unbounded':unbounded,'bounded_coverage':float(np.mean(cover))},
        single_CPU_seconds=times,operator_speedup=times['literal_photometry_68_models']/times['compiled_BLAS_photometry_68_models'],
        elapsed_seconds=time.perf_counter()-start,code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('model.py')]})
    (OUT/'numerical-validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    assert independent_pass, 'Independent primal validation failed; all attempts saved in numerical-validation.json.'

if __name__=='__main__':main()
