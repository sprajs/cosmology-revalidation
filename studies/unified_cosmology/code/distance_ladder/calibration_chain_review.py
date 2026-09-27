"""Recover the official anchored flat-LCDM chain and inspect startup sensitivity.

The fixed row-prefix cuts are descriptive diagnostics, not a new convergence
qualification, posterior reweighting, or a claim about the paper's processing.
No chain points are removed from the downloaded original.
"""
import argparse
import configparser
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

import numpy as np

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/unified-cosmology/calibration-chain-review'
BASE='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon+_Data/5_COSMOLOGY/chains/'
PINS={'README':('README','a3d0e33a6bb4172263bf810982692033586b6099d17f795c5fc2bfeeeaea8e0f'),
      'Pantheon+SH0ES_FlatLambdaCDM.txt.gz':('Pantheon+SH0ES/Pantheon+SH0ES_FlatLambdaCDM.txt.gz','815125575137e3dc8b87a2ef60c6e8922501aae5a4d06ad583b1aef25cb3fe1b')}
CUTS=[0.,.1,.2,.3,.5]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def relative(path):return str(Path(path).resolve().relative_to(ROOT))
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def config_section(header,name):
    block=header.split('## START_OF_'+name+'_INI\n',1)[1].split('## END_OF_'+name+'_INI',1)[0]
    lines='\n'.join(line[3:] if line.startswith('## ') else '' for line in block.splitlines())
    config=configparser.ConfigParser(interpolation=None);config.read_string(lines)
    return {section:dict(config[section]) for section in config.sections()}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--acquire',action='store_true')
    ap.add_argument('--quadrature',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-chain-review.json')
    args=ap.parse_args();WORK.mkdir(parents=True,exist_ok=True)
    records={}
    for name,(upstream,expected) in PINS.items():
        path=WORK/name
        if not path.exists():
            if not args.acquire:raise FileNotFoundError(str(path)+'; rerun --acquire')
            with urllib.request.urlopen(BASE+upstream,timeout=120) as response:body=response.read()
            assert hashlib.sha256(body).hexdigest()==expected
            path.write_bytes(body)
        assert sha(path)==expected
        records[relative(path)]={'url':BASE+upstream,'sha256':expected,'bytes':path.stat().st_size}
    chain_path=WORK/'Pantheon+SH0ES_FlatLambdaCDM.txt.gz'
    header=[]
    with gzip.open(chain_path,'rt') as stream:
        for line in stream:
            if not line.startswith('#'):break
            header.append(line)
    header=''.join(header);header_path=WORK/'header.txt';header_path.write_text(header)
    values=config_section(header,'VALUES');params=config_section(header,'PARAMS')
    priors=config_section(header,'PRIORS');assert not priors
    columns=header.splitlines()[0][1:].split()
    assert columns==['cosmological_parameters--omega_m','cosmological_parameters--h0','supernova_params--m','prior','post']
    data=np.loadtxt(chain_path);assert data.shape==(160000,5) and np.isfinite(data).all()
    boxes={key:[float(v) for v in values[section][name].split()][::2]
           for key,section,name in [('Omega_m','cosmological_parameters','omega_m'),('h0','cosmological_parameters','h0'),('M','supernova_params','m')]}
    expected_prior=-sum(np.log(hi-lo) for lo,hi in boxes.values())
    prior_error=float(np.max(abs(data[:,3]-expected_prior)));assert prior_error<1e-12
    for j,key in enumerate(['Omega_m','h0','M']):assert np.all((data[:,j]>=boxes[key][0])&(data[:,j]<=boxes[key][1]))
    producer=json.loads(args.quadrature.read_text());assert producer['status']=='qualified_conditional_SN_only_LCDM_quadrature'
    summaries=[]
    for fraction in CUTS:
        removed=int(len(data)*fraction);part=data[removed:];physical=part[:,:3].copy();physical[:,1]*=100
        summary={key:{'mean':float(np.mean(physical[:,j])),'sd':float(np.std(physical[:,j],ddof=0)),
                      'quantiles':np.quantile(physical[:,j],[.025,.16,.5,.84,.975]).tolist()}
                 for j,key in enumerate(['Omega_m','H0_km_s_Mpc','M'])}
        summaries.append({'discarded_initial_fraction':fraction,'removed_rows':removed,'retained_rows':len(part),
                         'summary':summary,'H0_Omega_covariance_population':float(np.cov(physical[:,:2].T,bias=True)[0,1]),
                         'difference_from_quadrature':{key:{'mean':summary[key]['mean']-producer['posterior'][key]['mean'],
                                                            'sd':summary[key]['sd']-producer['posterior'][key]['sd']}
                                                      for key in ['Omega_m','H0_km_s_Mpc']}})
    result={'status':'descriptive_released_chain_and_embedded_config_review_no_new_inference',
            'official_commit':'c447f0fea703fcd0fff57de5000947b5ca81286b','chain_rows':len(data),'columns':columns,
            'header_params_ini':params,'header_values_ini':values,'header_priors_ini':priors,
            'parameter_boxes':boxes,'constant_log_prior':float(expected_prior),'constant_log_prior_max_error':prior_error,
            'row_prefix_sensitivity':summaries,'quantile_probabilities':[.025,.16,.5,.84,.975],
            'first16_rows_M_range':[float(data[:16,2].min()),float(data[:16,2].max())],
            'first16_rows_logpost_range':[float(data[:16,4].min()),float(data[:16,4].max())],
            'last_half_rows_logpost_quantiles':np.quantile(data[len(data)//2:,4],[.025,.5,.975]).tolist(),
            'comparison_limits':['The complete public chain includes an obvious initialization transient in M and log posterior. Fixed prefix sensitivities are reported; no preferred burn fraction was selected to force agreement.',
                                 'Neither the original paper burn/summary procedure nor exact historical data bytes, CAMB source and sampling qualification are certified by the embedded header.',
                                 'After10%-50%prefix removal the released-chain H0 standard deviations are compared with deterministic quadrature; this is descriptive numerical consistency, not independent observations.',
                                 'Historical priors are uniform Omega_m[0.1,0.9], h0[0.55,0.91], M[-20,-18], unlike the bounded late-time quadrature priors and improper flat-M measure.',
                                 'The header declares omnuh2=0.00083, massive_nu=3 and massless_nu=0.046. This does not certify the executed species configuration: the inspected pre-chain public CSL interface accepts num_massive_neutrinos rather than massive_nu and explicitly ignores massless_nu. The exact historical installed interface/CAMB identity is absent. The reviewed SN-only quadrature declares matter+Lambda without radiation; no exact historical theory-density replay is claimed.',
                                 'The reported prior is normalized in h0=H0/100. Its constant changes under physical H0 units and does not supply extra observational information.'],
            'input_sources':records,'source_sha256':{relative(Path(__file__)):sha(__file__)},
            'quadrature_path':relative(args.quadrature),'quadrature_sha256':sha(args.quadrature),
            'header_path':relative(header_path),'header_sha256':sha(header_path),
            'CMB_calls':0,'background_calls':0,'new_sampling_steps':0,'posterior_qualification':False}
    for p,record in records.items():assert sha(ROOT/p)==record['sha256']
    write(args.output,result)
    print(json.dumps({'status':result['status'],'row_prefix_sensitivity':summaries},indent=2))


if __name__=='__main__':main()
