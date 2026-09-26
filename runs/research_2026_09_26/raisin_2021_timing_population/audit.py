#!/usr/bin/env python3
"""Read-only metadata audit of the coherent author 2021 DES RAISIN FITRES pair."""
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
P = json.loads((HERE/'protocol.json').read_text())
assert P['status'].startswith('Frozen')


def checked(path, sha):
    p = ROOT/path
    assert hashlib.sha256(p.read_bytes()).hexdigest() == sha, p
    return p


def read_fitres(path):
    names = None
    rows = {}
    header = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#') and names is None:
                header.append(line.strip())
            elif line.startswith('VARNAMES:'):
                assert names is None
                names = line.split()[1:]
            elif line.startswith('SN:'):
                assert names is not None
                vals = line.split()[1:]
                assert len(vals) == len(names), (path, len(vals), len(names))
                d = dict(zip(names, vals))
                cid = d['CID']
                assert cid not in rows, (path, cid)
                rows[cid] = d
    assert rows and names
    return rows, header, names


def stats(x):
    a = np.asarray(x, dtype=float)
    assert len(a) and np.isfinite(a).all()
    median = float(np.median(a))
    q = [0.001,0.01,0.05,0.25,0.5,0.75,0.95,0.99,0.999]
    return {'n':int(len(a)), 'mean':float(np.mean(a)),
            'sample_sd':float(np.std(a,ddof=1)) if len(a)>1 else None,
            'median':median, 'scaled_mad':float(1.4826*np.median(np.abs(a-median))),
            'min':float(np.min(a)), 'max':float(np.max(a)),
            'quantiles':{str(p):float(np.quantile(a,p)) for p in q}}


def number(row, key):
    return float(row[key])


def sample_metrics(cids, nir, joint):
    nr = [nir[c] for c in cids]
    jr = [joint[c] for c in cids]
    metrics = {
        'nir_peak_minus_initializer_day': [number(r,'PKMJD')-number(r,'PKMJDINI') for r in nr],
        'initializer_minus_sim_peak_day': [number(r,'PKMJDINI')-number(r,'SIM_PKMJD') for r in nr],
        'nir_peak_minus_sim_peak_day': [number(r,'PKMJD')-number(r,'SIM_PKMJD') for r in nr],
        'joint_peak_minus_sim_peak_day': [number(r,'PKMJD')-number(r,'SIM_PKMJD') for r in jr],
        'joint_peak_minus_nir_peak_day': [number(j,'PKMJD')-number(n,'PKMJD') for n,j in zip(nr,jr)],
        'zHEL': [number(r,'zHEL') for r in nr],
        'SIM_ZCMB': [number(r,'SIM_ZCMB') for r in nr],
    }
    out = {k:stats(v) for k,v in metrics.items()}
    x=np.abs(np.asarray(metrics['initializer_minus_sim_peak_day']))
    out['initializer_near_truth_counts']={str(t):int(np.sum(x<=t)) for t in [0.01,0.03,0.05]}
    edges=np.array([0,.2,.3,.4,.5,.6,.7,.8,1.,2.])
    z=np.asarray(metrics['zHEL'])
    out['zHEL_bins']={'edges':edges.tolist(), 'counts':np.histogram(z,bins=edges)[0].tolist(),
                      'underflow':int(np.sum(z<edges[0])), 'overflow':int(np.sum(z>=edges[-1]))}
    return out


def main():
    src=P['inputs']
    nirpath=checked(src['nir'],src['nir_sha256'])
    jointpath=checked(src['optnir'],src['optnir_sha256'])
    memberpath=checked(src['plotted500'],src['plotted500_sha256'])
    nmlpaths=[HERE/'author/fit/sim/REFAC_DES_RAISIN_NIR.nml',HERE/'author/fit/sim/REFAC_DES_RAISIN_optnir.nml']
    for p,sha in zip(nmlpaths,[src['nir_nml_sha256'],src['optnir_nml_sha256']]):
        assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
    nir,nh,nn=read_fitres(nirpath)
    joint,jh,jn=read_fitres(jointpath)
    for col in ['PKMJDINI','PKMJD','SIM_PKMJD','zHEL','SIM_ZCMB','SIM_LIBID','ERRFLAG_FIT']:
        assert col in nn and col in jn
    assert any('v11_04d' in s for s in nh) and any('v11_04d' in s for s in jh)
    with memberpath.open(newline='') as f: member=list(csv.DictReader(f))
    fixed500=[r['CID'] for r in member]
    assert len(fixed500)==len(set(fixed500))==500
    first8=[str(i) for i in range(1,9)]
    common=set(nir)&set(joint)
    cohorts={'all_common':sorted(common,key=int), 'plotted500':fixed500, 'first8':first8}
    for name,cids in cohorts.items():
        missing=[c for c in cids if c not in common]
        assert not missing,(name,missing)
    meta=['PKMJDINI','SIM_PKMJD','SIM_LIBID','SIM_ZCMB','SIM_STRETCH','SIM_AV','SIM_RV','zHEL','zCMB','zHD']
    mismatches={k:[c for c in cohorts['all_common'] if nir[c][k]!=joint[c][k]] for k in meta}
    summaries={name:sample_metrics(ids,nir,joint) for name,ids in cohorts.items()}
    nir_all=sorted(nir,key=int)
    summaries['all_nir_unpaired']={
       'nir_peak_minus_initializer_day':stats([number(nir[c],'PKMJD')-number(nir[c],'PKMJDINI') for c in nir_all]),
       'initializer_minus_sim_peak_day':stats([number(nir[c],'PKMJDINI')-number(nir[c],'SIM_PKMJD') for c in nir_all]),
       'nir_peak_minus_sim_peak_day':stats([number(nir[c],'PKMJD')-number(nir[c],'SIM_PKMJD') for c in nir_all]),
       'zHEL':stats([number(nir[c],'zHEL') for c in nir_all]),
    }
    flags={label:{k:dict(sorted(Counter(r[k] for r in tab.values()).items()))
                  for k in ['ERRFLAG_FIT','CUTFLAG_SNANA','TYPE','SIM_TYPE_INDEX']}
           for label,tab in [('nir',nir),('optnir',joint)]}
    missing_joint=sorted(set(nir)-set(joint),key=int)
    joint_only=sorted(set(joint)-set(nir),key=int)
    assert len(missing_joint)==5 and not joint_only
    missfields=['CID','ERRFLAG_FIT','CUTFLAG_SNANA','TYPE','FIELD','zHEL','SIM_ZCMB','SIM_LIBID','PKMJDINI','SIM_PKMJD','PKMJD','NDOF','SNRMAX1','SNRMAX2','SNRMAX3']
    with (HERE/'missing-joint.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=missfields);w.writeheader()
        for c in missing_joint:w.writerow({k:nir[c][k] for k in missfields})
    bymember={r['CID']:r for r in member}
    firstfields=['CID','zHEL','SIM_LIBID','n_J','n_H','nir_peak_minus_initializer_day','initializer_minus_sim_peak_day','joint_peak_minus_sim_peak_day','joint_minus_nir_peak_day']
    with (HERE/'first8.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=firstfields);w.writeheader()
        for c in first8:
            n,j=nir[c],joint[c]
            w.writerow({'CID':c,'zHEL':n['zHEL'],'SIM_LIBID':n['SIM_LIBID'],
                        'n_J':bymember[c]['n_J'],'n_H':bymember[c]['n_H'],
                        'nir_peak_minus_initializer_day':number(n,'PKMJD')-number(n,'PKMJDINI'),
                        'initializer_minus_sim_peak_day':number(n,'PKMJDINI')-number(n,'SIM_PKMJD'),
                        'joint_peak_minus_sim_peak_day':number(j,'PKMJD')-number(j,'SIM_PKMJD'),
                        'joint_minus_nir_peak_day':number(j,'PKMJD')-number(n,'PKMJD')})
    cadence={
      'n_J':stats([int(r['n_J']) for r in member]),
      'n_H':stats([int(r['n_H']) for r in member]),
      'total_J':sum(int(r['n_J']) for r in member),
      'total_H':sum(int(r['n_H']) for r in member),
      'n_J_plus_H_matches_nir_NDOF_plus_1':sum(int(r['n_J'])+int(r['n_H'])==int(float(nir[r['CID']]['NDOF']))+1 for r in member),
    }
    result={
      'source_revision':'2021-11-11 author FITRES pair at aaa709ead7a7339d56a2a4604d78991329d9368f',
      'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
      'source_headers':{'nir':nh,'optnir':jh},
      'counts':{'nir':len(nir),'optnir':len(joint),'common':len(common),'nir_only':len(missing_joint),'optnir_only':len(joint_only),'fixed_plot_subset':len(fixed500),'first8':len(first8)},
      'metadata_mismatch_counts':{k:len(v) for k,v in mismatches.items()},
      'metadata_mismatch_cids':{k:v for k,v in mismatches.items() if v},
      'nir_only_cids':missing_joint,'optnir_only_cids':joint_only,
      'fit_flags':flags,
      'fixed_plot_cadence':cadence,
      'summaries':summaries,
      'missing_joint_reason':'Not encoded in the NIR FITRES row or SUBMIT.INFO; no per-CID failure cause inferred.'
    }
    (HERE/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    (HERE/'per-cohort-summary.json').write_text(json.dumps(summaries,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'counts':result['counts'],'nir_only_cids':missing_joint,
          'metadata_mismatch_counts':result['metadata_mismatch_counts'],
          'timing_all_common':{k:dict(mean=v['mean'],sd=v['sample_sd'],median=v['median'],q01=v['quantiles']['0.01'],q99=v['quantiles']['0.99']) for k,v in summaries['all_common'].items() if k.endswith('_day')}},indent=2))

if __name__=='__main__':main()
