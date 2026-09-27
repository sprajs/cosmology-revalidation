"""Bounded gradient diagnostic of frozen synthetic optical BayeSN densities.

No posterior sampling, held-out outcomes, observed photometry, or model edits.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=Path(__file__).resolve().parent/'bayesn_heldout'
DESIGN=Path(__file__).with_name('bayesn-synthetic-gradient-design.json')
NAMES=['D','log_AV','logit_RV','theta']+[f'epsilon_white[{i}]'for i in range(42)]+['logit_tau']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(path,obj):
    assert not path.exists(),f'Preserve previous execution: {path}'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def run(args):
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
        JAX_PLATFORM_NAME='cpu',JAX_ENABLE_X64='true',
        XLA_FLAGS='--xla_cpu_multi_thread_eigen=false --xla_force_host_platform_device_count=1')
    os.sched_setaffinity(0,[max(os.sched_getaffinity(0)) if args.cpu is None else args.cpu])
    opened=set()
    def guard(event,arguments):
        if event!='open':return
        p=str(arguments[0]);low=p.lower()
        if 'synthetic-nir.json' in low or 'nir-sealed' in low or '/photometry/' in low or low.endswith('truth.json'):
            raise RuntimeError('Gradient audit cannot open held-out/truth/source photometry: '+p)
        if isinstance(arguments[0],(str,bytes,os.PathLike)):
            opened.add(str(Path(p).resolve()))
    sys.addaudithook(guard)
    import numpy as np
    import jax
    import jax.numpy as jnp
    from numpyro.infer.util import log_density,potential_energy
    from scipy.special import expit,log_expit,log_ndtr
    from scipy.stats import norm,expon,uniform
    sys.path.insert(0,str(PACKAGE))
    import model
    design=json.loads(DESIGN.read_text());start=time.monotonic()
    case=args.work/'case-00';optical=case/'optical.json';config=args.work/'release-filter-config.yaml'
    worker=json.loads((case/'LCDM/worker-start.json').read_text())
    current={p.name:sha(p)for p in PACKAGE.iterdir()if p.suffix in ('.py','.json')}
    assert current==worker['source_sha256'] and sha(optical)==worker['optical_sha256']
    assert sha(config)==worker['filter_config_sha256']
    data=json.loads(optical.read_text())
    assert data['scope']=='synthetic' and all(r['role']=='optical'for r in data['rows'])
    times=np.array([r['trigger_rest_time']for r in data['rows']]);error=np.array(data['errors']);flux=np.array(data['flux'])
    inputs={str(p.relative_to(ROOT)):sha(p)for p in [DESIGN,optical,config,case/'LCDM/worker-start.json']}
    selected=[]
    for spec in design['retained_states']:
        folder=case/spec['arm'];recordpath=folder/f"chain-{spec['chain']}.json"
        record=json.loads(recordpath.read_text());path=folder/record['archive']
        assert sha(path)==record['sha256']
        inputs.update({str(p.relative_to(ROOT)):sha(p)for p in [recordpath,path]})
        with np.load(path)as raw:
            index=spec['draw'] if isinstance(spec['draw'],int)else int(np.flatnonzero(raw['stat_diverging'])[0])
            p=np.r_[raw['D'][index],raw['AV'][index],raw['RV'][index],raw['theta'][index],raw['epsilon_white'][index],raw['tau'][index]]
            assert p.shape==(47,) and np.isfinite(p).all()
            selected.append(dict(kind='retained_state',arm=spec['arm'],chain=spec['chain'],draw=index,
                divergence_associated_retained_state=bool(raw['stat_diverging'][index]),physical_coordinates=p.tolist()))
    base=np.array(selected[0]['physical_coordinates']);phase=times-base[-1]
    knotrow=int(np.argmin(abs(phase-np.rint(phase))));knot=float(np.rint(phase[knotrow]));boundary=float(times[knotrow]-knot)
    for offset in [-1e-4,0.,1e-4]:
        p=base.copy();p[-1]=boundary+offset;assert -10<p[-1]<20
        selected.append(dict(kind='phase_control',arm='LCDM',control_delta_tau_days=offset,
            knot_optical_row=knotrow,knot_phase=knot,physical_coordinates=p.tolist()))
    # Freeze exact states and hashes before constructing/evaluating the forward operator.
    args.output.mkdir(parents=True,exist_ok=True)
    execution=dict(design_sha256=sha(DESIGN),source_sha256=sha(__file__),input_sha256=inputs,
        reviewed_package_sha256=current,states=selected,CPU_affinity=sorted(os.sched_getaffinity(0)),
        scope='Synthetic optical states only; failed/incomplete parent sampling does not qualify posterior inference.')
    save(args.output/'execution-design.json',execution)
    kernel=model.Kernel(args.archive,config,data['metadata'],data['rows'],2400)
    def sites(p):return dict(D=p[0],AV=p[1],RV=p[2],theta=p[3],epsilon_white=p[4:46],tau=p[46])
    def jconstrain(u):
        return u.at[1].set(jnp.exp(u[1])).at[2].set(1.2+4.8*jax.nn.sigmoid(u[2])).at[46].set(-10+30*jax.nn.sigmoid(u[46]))
    def jacobian(u):
        return u[1]+jnp.log(4.8)+jax.nn.log_sigmoid(u[2])+jax.nn.log_sigmoid(-u[2])+jnp.log(30.)+jax.nn.log_sigmoid(u[46])+jax.nn.log_sigmoid(-u[46])
    def unconstrain(p):
        u=p.copy();u[1]=np.log(p[1]);u[2]=np.log((p[2]-1.2)/(6-p[2]));u[46]=np.log((p[46]+10)/(20-p[46]));return u
    def constrain(u):
        p=u.copy();p[1]=np.exp(u[1]);p[2]=1.2+4.8*expit(u[2]);p[46]=-10+30*expit(u[46]);return p
    def independent_physical(p,arm):
        if not(p[1]>0 and 1.2<p[2]<6 and -10<p[46]<20):return -np.inf
        d=p[0]
        if arm=='LCDM':ld=norm.logpdf(d,data['metadata']['mu_LCDM'],np.hypot(data['metadata']['sigma_external'],.088))
        else:
            a,b=((d-20)/.088,(d-50)/.088)if d<35 else((50-d)/.088,(20-d)/.088)
            la,lb=log_ndtr(a),log_ndtr(b);ld=la+np.log(-np.expm1(lb-la))-np.log(30.)
        pred=np.asarray(kernel.forward(jnp.array(p)))
        return float(ld+expon.logpdf(p[1],scale=.329)+uniform.logpdf(p[2],loc=1.2,scale=4.8)
            +norm.logpdf(p[3])+np.sum(norm.logpdf(p[4:46]))+uniform.logpdf(p[46],loc=-10,scale=30)
            +np.sum(norm.logpdf(flux,loc=pred,scale=error)))
    def independent_unconstrained(u,arm):
        jac=u[1]+np.log(4.8)+log_expit(u[2])+log_expit(-u[2])+np.log(30)+log_expit(u[46])+log_expit(-u[46])
        return independent_physical(constrain(u),arm)+jac
    objectives={};physical={};native={}
    for arm in ['LCDM','broad']:
        likelihood=model.numpyro_model(kernel,data,arm)
        lp=lambda p,likelihood=likelihood:log_density(likelihood,(),{},sites(p))[0]
        physical[arm]=jax.jit(jax.value_and_grad(lp))
        fun=lambda u,lp=lp:lp(jconstrain(u))+jacobian(u)
        objectives[arm]=jax.jit(jax.value_and_grad(fun))
        native[arm]=jax.jit(jax.value_and_grad(lambda u,likelihood=likelihood:-potential_energy(likelihood,(),{},sites(u))))
    steps=np.array([1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,3e-6,1e-6,3e-7,1e-7])
    rows=[]
    for number,spec in enumerate(selected):
        p=np.array(spec['physical_coordinates']);u=unconstrain(p);arm=spec['arm']
        val,gradient=objectives[arm](jnp.array(u));value=float(val);grad=np.asarray(gradient)
        nv,ng=native[arm](jnp.array(u));nativegrad=np.asarray(ng)
        pv,pg=physical[arm](jnp.array(p));pvalue=float(pv);pgrad=np.asarray(pg)
        independent=independent_unconstrained(u,arm)
        assert np.isfinite(value) and np.isfinite(grad).all() and abs(value-independent)<1e-8
        assert abs(float(nv)-value)<1e-8
        coordinate=[];all_step_records=[]
        for j,name in enumerate(NAMES):
            fd=[];cross=[]
            for relative_step in steps:
                h=float(relative_step*max(1,abs(u[j])));left=u.copy();right=u.copy();left[j]-=h;right[j]+=h
                minus=independent_unconstrained(left,arm);plus=independent_unconstrained(right,arm)
                derivative=(plus-minus)/(2*h)
                crosses=bool(np.any(np.floor(times-constrain(left)[-1])!=np.floor(times-constrain(right)[-1])))
                normalized=abs(derivative-grad[j])/(1+abs(grad[j]))
                fd.append(derivative);cross.append(crosses)
                all_step_records.append(dict(coordinate=name,relative_step=float(relative_step),actual_step=h,
                    finite_difference=float(derivative),autodiff=float(grad[j]),normalized_error=float(normalized),
                    crosses_integer_phase=crosses))
            errors=abs(np.array(fd)-grad[j])/(1+abs(grad[j]));ok=(errors<1e-4)&~np.array(cross)
            paired=bool(np.any(ok[:-1]&ok[1:]));smooth=np.flatnonzero(~np.array(cross))
            exactcontrol=spec['kind']=='phase_control'and spec['control_delta_tau_days']==0 and j==46
            coordinate.append(dict(coordinate=name,autodiff=float(grad[j]),
                status='intentional_phase_knot_nondifferentiability'if exactcontrol else('smooth_agreement'if paired else('no_noncrossing_steps'if not len(smooth)else'smooth_disagreement')),
                adjacent_steps_agree=paired,smallest_step_error=float(errors[-1]),
                best_noncrossing_normalized_error=float(errors[smooth].min())if len(smooth)else None,
                noncrossing_step_count=len(smooth)))
        tau=[]
        for h in [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]:
            pm=p.copy();pp=p.copy();pm[-1]-=h;pp[-1]+=h
            vm=independent_physical(pm,arm);vp=independent_physical(pp,arm)
            tau.append(dict(step_days=h,left_derivative=(pvalue-vm)/h,right_derivative=(vp-pvalue)/h,
                central_derivative=(vp-vm)/(2*h),autodiff=float(pgrad[-1]),
                crosses_integer_phase=bool(np.any(np.floor(times-pm[-1])!=np.floor(times-pp[-1])))))
        raw=args.output/f'state-{number:02d}.json'
        record=dict(state=spec,unconstrained_coordinates=u.tolist(),logdensity=value,independent_density_error=abs(value-independent),
            native_potential_value_error=abs(float(nv)-value),native_potential_gradient_max_error=float(np.max(abs(nativegrad-grad))),
            nearest_integer_phase_distance=float(np.min(abs(times-p[-1]-np.rint(times-p[-1])))),
            coordinate_summary=coordinate,finite_differences=all_step_records,physical_tau_sweep=tau)
        save(raw,record)
        rows.append(dict(state_index=number,selection={k:v for k,v in spec.items()if k!='physical_coordinates'},
            raw_path=str(raw.relative_to(ROOT)),raw_sha256=sha(raw),nearest_integer_phase_distance=record['nearest_integer_phase_distance'],
            independent_density_error=record['independent_density_error'],native_potential_value_error=record['native_potential_value_error'],
            native_potential_gradient_max_error=record['native_potential_gradient_max_error'],
            coordinate_summary=coordinate,physical_tau_sweep=tau))
        print(json.dumps(dict(state=number,statuses={status:sum(v['status']==status for v in coordinate)for status in sorted({v['status']for v in coordinate})},
            seconds=time.monotonic()-start)),flush=True)
    assert model.source_hashes()==current
    assert {p:sha(ROOT/p)for p in inputs}==inputs
    archive_inputs={str(Path(p).relative_to(args.archive)):sha(p)for p in sorted(opened)
        if Path(p).is_file()and Path(p).is_relative_to(args.archive)and '.venv'not in Path(p).parts}
    failures=[{'state':r['state_index'],'coordinate':c['coordinate'],'status':c['status']}for r in rows for c in r['coordinate_summary']if c['status']in ['smooth_disagreement','no_noncrossing_steps']]
    result=dict(status='completed_bounded_gradient_diagnostic',source_sha256=sha(__file__),design_sha256=sha(DESIGN),
        execution_design_sha256=sha(args.output/'execution-design.json'),input_sha256=inputs,reviewed_package_sha256=current,
        archive_opened_input_sha256=archive_inputs,states=rows,smooth_diagnostic_failures=failures,
        all_predeclared_smooth_checks_passed=not failures,seconds=time.monotonic()-start,
        versions={'numpy':np.__version__,'jax':jax.__version__},native_spectra_or_cosmology_calls=0,
        constraints=['Synthetic optical only; audit hook rejected held-out/truth/source-photometry opens.',
            'No new chains or changes to the frozen package; every saved finite-difference step retained.',
            'Retained divergence states are not failed leapfrog locations. A phase derivative jump is not proof that it caused a divergence.',
            'This local check cannot establish global derivative correctness, posterior adequacy, or convergence.'])
    save(args.output/'summary.json',result)
    print(json.dumps(dict(status=result['status'],smooth_diagnostic_failures=failures,seconds=result['seconds'])),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--archive',type=Path,default=Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85');p.add_argument('--cpu',type=int)
    args=p.parse_args();args.work=args.work.resolve();args.output=args.output.resolve();args.archive=args.archive.resolve();run(args)
