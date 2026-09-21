#!/usr/bin/env python3
"""Propagation diagnostic using published weights, not an independent cosmology fit."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official'
ap=argparse.ArgumentParser();ap.add_argument('--sample',default='recovered');args=ap.parse_args()
d=pd.read_csv(OUT/f'results/{args.sample}_comparison.csv')
result={'sample':args.sample,'interpretation':'Fixed published alpha,beta,bias corrections and diagonal published MUERR_FINAL weights. Linear Omega_m response after profiling one intercept is a diagnostic, not rerun BBC/cosmology or new physical correction.','subsets':{}}
for label,subset in [('all',np.ones(len(d),bool)),('DES',d.IDSURVEY==10),('PIa_gt_0.999',d.PROB_SNNV19>.999),('BEAMS_Ia_gt_0.999',d.PROBCC_BEAMS<.001)]:
    q=d.loc[subset];w=1/q.MUERR_FINAL.to_numpy()**2
    r=q.delta_Tripp_fixed_coefficients.to_numpy();b=q.delta_mB.to_numpy();z=q.zHD_published.to_numpy()
    derivative=(FlatLambdaCDM(H0=70,Om0=.35195+.0001).distmod(z).value-FlatLambdaCDM(H0=70,Om0=.35195-.0001).distmod(z).value)/.0002
    mean=lambda x:float(np.sum(w*x)/np.sum(w))
    dc=derivative-mean(derivative);rc=r-mean(r)
    result['subsets'][label]={'n':len(q),'n_abs_dmB_gt_.01':int(np.sum(abs(b)>.01)),'max_abs_dmB':float(abs(b).max()),'weighted_mean_delta_Tripp_mag':mean(r),'weighted_rms_delta_Tripp_mag':float(np.sqrt(mean(r*r))),'weighted_centered_rms_delta_Tripp_mag':float(np.sqrt(mean(rc*rc))),'n_abs_delta_Tripp_gt_.01':int(np.sum(abs(r)>.01)),'n_abs_delta_Tripp_over_published_sigma_gt_.1':int(np.sum(abs(r)/q.MUERR_FINAL>.1)),'linear_Omega_m_shift_diagonal_weights_profiled_intercept':float(np.sum(w*dc*rc)/np.sum(w*dc*dc))}
path=OUT/f'diagnostics/{args.sample}_weighted_gate.json';path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
