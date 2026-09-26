"""Independent saved-artifact verification; launches no native fits."""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import json
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit'
GATE=BASE/'native-profile-gate-v2'
CID='DES16C1cim'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def text_table(path):
    names=None;rows=[]
    for line in Path(path).read_text().splitlines():
        q=line.split()
        if not q:continue
        if q[0]=='VARNAMES:':names=q[1:]
        if q[0] in ['SN:','ROW:','OBS:']:
            assert names is not None and len(q[1:])==len(names)
            rows.append(dict(zip(names,q[1:])))
    return rows

def phot(path):
    columns=None;rows=[]
    for lineno,line in enumerate(Path(path).read_text().splitlines(),1):
        q=line.split()
        if q and q[0]=='VARLIST:':columns=q[1:]
        elif q and q[0]=='OBS:':
            row=dict(zip(columns,q[1:]));row['source_line']=lineno;rows.append(row)
    return rows

def parse_log(path):
    lines=Path(path).read_text().splitlines()
    assert 'ENDING PROGRAM GRACEFULLY.' in '\n'.join(lines)
    records={x:[] for x in ['FLUX','COVINV','OBJECTIVE','SETUP']}
    for line in lines:
        q=line.split()
        if not q:continue
        for k in records:
            if q[0]=='PHASE2_'+k+':':
                assert q[1]==CID
                records[k].append(q[2:])
    # The native driver invokes its final callback twice. Require complete,
    # exactly identical copies rather than silently taking the last record.
    repeats=len(records['OBJECTIVE'])
    assert repeats==len(records['SETUP'])==2
    objective=records['OBJECTIVE'][0];n=int(objective[0])
    for label in records:
        count=n if label in ['FLUX','COVINV'] else 1
        assert len(records[label])==repeats*count
        assert records[label][:count]==records[label][count:]
        records[label]=records[label][:count]
    assert [int(q[0]) for q in records['FLUX']]==list(range(1,n+1))
    assert [int(q[0]) for q in records['COVINV']]==list(range(1,n+1))
    values=np.array([[float(x) for x in q[2:]] for q in records['FLUX']])
    W=np.array([[float(x) for x in q[1:]] for q in records['COVINV']])
    return dict(n=n,repeated_identical_callbacks=repeats,band=np.array([q[1] for q in records['FLUX']]),v=values,W=W,
                objective=np.array([float(x) for x in objective[1:]]),
                setup=np.array([float(x) for x in records['SETUP'][0]]))

def main():
    protocol=json.loads((GATE/'protocol.json').read_text())
    manifest=json.loads((GATE/'manifest.json').read_text())
    input_ok={p:sha(ROOT/p)==h for p,h in protocol['input_hashes'].items()}
    output_ok={p:sha(GATE/p)==h for p,h in manifest['all_files_sha256'].items()}
    assert all(input_ok.values()) and all(output_ok.values())
    original=next(r for r in text_table(BASE/'fits/full/fixed/B/fit.FITRES.TEXT') if r['CID']==CID)
    audit=text_table(GATE/'reference/fit.FITRES.TEXT')[0]
    assert original.keys()==audit.keys()
    fitres_diff={k:[original[k],audit[k]] for k in original if original[k]!=audit[k]}
    assert not fitres_diff,fitres_diff
    original_lc=[r for r in text_table(BASE/'fits/full/fixed/B/fit.LCPLOT.TEXT') if r['CID']==CID]
    audit_lc=text_table(GATE/'reference/fit.LCPLOT.TEXT')
    assert original_lc==audit_lc
    lcflags=Counter(r['DATAFLAG'] for r in original_lc)
    assert lcflags['1']==74
    native_manifest=json.loads((BASE/'input-manifest.json').read_text())['inputs_sha256']
    original_native_paths=['phase2/official/build/SNANA-current/bin/snlc_fit.exe','phase2/official/build/SNANA-current/src/snlc_fit.F90']
    original_binary_ok={p:sha(ROOT/p)==native_manifest[p] for p in original_native_paths}
    assert all(original_binary_ok.values())
    logs={n:parse_log(GATE/n/'fit.log') for n in ['reference','fixed','distance_plus','distance_minus','sorted','reversed']}
    exports={n:dict(np.load(GATE/n/'native.npz')) for n in logs}
    checks={}
    for name,log in logs.items():
        v=log['v'];w=log['W'];obj=log['objective'];z=exports[name]
        execution=json.loads((GATE/name/'execution.json').read_text())
        assert execution['returncode']==0
        assert execution['log_sha256']==sha(GATE/name/'fit.log')
        assert execution['nml_sha256']==sha(GATE/name/'fit.nml')
        raw={'band':log['band'],'MJD':v[:,0],'rest_phase':v[:,1],'model_flux':v[:,2],
             'model_magerr':v[:,3],'data_flux':v[:,4],'data_fluxerr':v[:,5],
             'zHEL':v[:,6],'MWEBV':v[:,7],'W':w,'parameters':obj[2:6],
             'setup':log['setup'],'total_chi2':obj[0],'prior_chi2':obj[1]}
        for k,a in raw.items():assert np.array_equal(a,z[k]),(name,k)
        np.linalg.cholesky(w)
        # Independent inversion route: Cholesky solve rather than np.linalg.inv.
        c=cho_solve(cho_factor(w,lower=True),np.eye(log['n']))
        residual=v[:,4]-v[:,2]
        q=float(residual@w@residual)
        reconstructed=q+obj[1]
        checks[name]={'epochs':log['n'],'npz_raw_log_exact':True,'identical_final_callbacks':log['repeated_identical_callbacks'],
                      'C_cholesky_vs_export_max_abs':float(np.max(abs(c-z['C']))),
                      'W_symmetry_max_abs':float(np.max(abs(w-w.T))),
                      'C_min_eigenvalue':float(np.linalg.eigvalsh(c).min()),
                      'quadratic_plus_prior_minus_native_FCN':reconstructed-obj[0],
                      'prior_chi2':float(obj[1])}
        assert log['n']==74 and abs(reconstructed-obj[0])<1e-7

    ref=exports['reference'];fixed=exports['fixed']
    key=lambda a:list(zip(a['band'],a['MJD']))
    keys=key(fixed);assert len(set(keys))==74
    raw_b=phot(BASE/'data/RSR_B'/f'{CID}.snana.dat')
    raw_a=phot(BASE/'data/RSR_A'/f'{CID}.snana.dat')
    bkeys=Counter((r['FLT'],float(r['MJD'])) for r in raw_b)
    for k in keys:assert bkeys[k]==1
    bmap={(r['FLT'],float(r['MJD'])):r for r in raw_b}
    amap={(r['FLT'],float(r['MJD'])):r for r in raw_a}
    accepted=[]
    for j,k in enumerate(keys):
        row=bmap[k]
        assert float(np.float32(row['FLUXCAL']))==fixed['data_flux'][j]
        assert float(np.float32(row['FLUXCALERR']))==fixed['data_fluxerr'][j]
        if k in amap:
            assert row['FLUXCAL']==amap[k]['FLUXCAL'] and row['FLUXCALERR']==amap[k]['FLUXCALERR']
        assert (k in amap)==(float(row['FLUXCAL'])>0)
        accepted.append({'index':j,'band':k[0],'MJD':k[1],'B_source_line':row['source_line'],
                         'A_positive_member':k in amap,'data_flux':float(fixed['data_flux'][j]),
                         'data_fluxerr':float(fixed['data_fluxerr'][j])})

    probes={}
    for name,offset in [('fixed',0.),('distance_plus',.15),('distance_minus',-.15),('sorted',0.),('reversed',0.)]:
        z=exports[name];zkeys=key(z);assert Counter(keys)==Counter(zkeys)
        mapping=np.array([zkeys.index(k) for k in keys])
        for field in ['data_flux','data_fluxerr','MWEBV','zHEL']:
            assert np.array_equal(fixed[field],z[field][mapping]),(name,field)
        requested=ref['parameters'].copy();requested[0]+=offset
        rounded=requested.astype('f4').astype(float)
        assert np.array_equal(z['parameters'],rounded)
        delta=float(z['parameters'][0]-fixed['parameters'][0])
        predicted=fixed['model_flux']*10**(-.4*delta)
        diff=z['model_flux'][mapping]-predicted
        probes[name]={'actual_DLMAG_offset':delta,'float32_seed_exact':True,
                      'pointwise_mean_scaling_max_relative':float(np.max(abs(diff/predicted))),
                      'mean_scaling_max_quoted_sigma':float(np.max(abs(diff/fixed['data_fluxerr']))),
                      'C_fractional_Frobenius_change_from_reference':float(np.linalg.norm(z['C'][np.ix_(mapping,mapping)]-ref['C'])/np.linalg.norm(ref['C']))}
        assert probes[name]['pointwise_mean_scaling_max_relative']<1e-8

    # Independently solve a one-column whitened least-squares problem.
    L=np.linalg.cholesky(ref['C'])
    white_h=solve_triangular(L,fixed['model_flux'],lower=True)
    white_y=solve_triangular(L,ref['data_flux'],lower=True)
    a=float(np.linalg.lstsq(white_h[:,None],white_y,rcond=None)[0][0])
    D=float(fixed['parameters'][0]-2.5*np.log10(a))
    Q=float(np.sum((white_y-a*white_h)**2))
    source=OUT.parent/'fixed_covariance_profile.py'
    spec=importlib.util.spec_from_file_location('reviewed_profile_kernel',source)
    kernel=importlib.util.module_from_spec(spec);spec.loader.exec_module(kernel)
    kr=kernel.profile_amplitude(ref['data_flux'],ref['W'],fixed['model_flux'],fixed['parameters'][0],other_prior=float(ref['prior_chi2']))
    reported=json.loads((GATE/'result.json').read_text())['analytic_positive_amplitude']
    amplitude={'independent_amplitude':a,'independent_DLMAG':D,'independent_data_quadratic':Q,
               'kernel_amplitude_difference':kr['amplitude']-a,'kernel_DLMAG_difference':kr['DLMAG']-D,
               'kernel_data_quadratic_difference':kr['data_quadratic']-Q,
               'reported_amplitude_difference':reported['a']-a,'reported_DLMAG_difference':reported['DLMAG']-D,
               'reported_quadratic_difference':reported['data_chi2']-Q}
    assert max(abs(amplitude[k]) for k in amplitude if k.endswith('difference'))<1e-10

    A=np.array([j for j,k in enumerate(keys) if k in amap]);N=np.array([j for j,k in enumerate(keys) if k not in amap])
    C=ref['C'];W=ref['W'];CA=C[np.ix_(A,A)]
    marginal_precision=cho_solve(cho_factor(CA,lower=True),np.eye(len(A)))
    schur=W[np.ix_(A,A)]-W[np.ix_(A,N)]@np.linalg.solve(W[np.ix_(N,N)],W[np.ix_(N,A)])
    restriction={'A_positive_rows':len(A),'B_signed_rows':74,'removed_nonpositive_rows':len(N),
                 'shared_raw_flux_and_errors_exact':True,'A_C_min_eigenvalue':float(np.linalg.eigvalsh(CA).min()),
                 'inverse_marginalC_vs_precision_Schur_max_abs':float(np.max(abs(marginal_precision-schur))),
                 'incorrect_W_principal_block_fractional_difference':float(np.linalg.norm(W[np.ix_(A,A)]-marginal_precision)/np.linalg.norm(marginal_precision)),
                 'rule':'A covariance is C_B[A,A], then invert. W_B[A,A] is conditional precision and must not be used as marginal precision.'}
    assert restriction['inverse_marginalC_vs_precision_Schur_max_abs']<1e-10
    np.savez_compressed(OUT/'restriction-check.npz',A_indices=A,N_indices=N,C_A=CA,W_A=marginal_precision,W_A_via_Schur=schur)
    (OUT/'accepted-rows.json').write_text(json.dumps(accepted,indent=2)+'\n')
    result={'status':'Independent saved-native proof check passes; no new fits',
            'hash_checks':{'input_count':len(input_ok),'output_count':len(output_ok),'all_pass':True},
            'FITRES_all_columns_identical':len(original),'FITRES_differences':fitres_diff,
            'original_binary_source_match_refit_manifest':original_binary_ok,
            'LCPLOT_all_target_rows_identical':len(original_lc),'LCPLOT_DATAFLAG_counts':dict(lcflags),
            'runs':checks,'probes':probes,'amplitude':amplitude,'A_B_restriction':restriction,
            'scope':'One fixed-peak modern native state only. No global profile or observed sign-correction inferred.',
            'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [GATE/'manifest.json',GATE/'protocol.json',source,Path(__file__),BASE/'data/RSR_A'/f'{CID}.snana.dat']}}
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
