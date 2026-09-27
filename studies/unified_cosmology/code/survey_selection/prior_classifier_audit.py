#!/usr/bin/env python3
"""Read-only audit that prior age-injection campaigns used measured native timing."""
import json,sys,re
import numpy as np
import pandas as pd
from astropy.io import fits
from common import ROOT,WORK,RESULTS,sha
sys.path.insert(0,str(ROOT));from lib.records import fitres

def main():
    prior=ROOT/'.work/survey-physics';paths=sorted((prior/'fits').glob('SP*/classifier.json'))+sorted((prior/'scaled-bbc/fits').glob('SPB_*/classifier.json'));checks=[];inputs={};total=0
    for p in paths:
        r=json.loads(p.read_text());folder=p.parent;prep=folder/'classifier-preprocessing.csv';nativepath=folder/'fit.SNANA.TEXT';fitpath=folder/'fit.FITRES.TEXT'
        assert sha(prep)==r['preprocessing_sha256'] and sha(nativepath)==r['native_peak_sha256'] and sha(fitpath)==r['native_fit_sha256']
        m=pd.read_csv(prep,dtype={'CID':str}).set_index('CID');d=fitres(nativepath);delta=float(np.max(abs(m.peak_mjd_measured-d.loc[m.index,'PKMJDINI'])))
        assert delta<1e-9
        txt=(folder/'fit.nml').read_text();option=int(re.search(r'OPT_SETPKMJD\s*=\s*(\d+)',txt).group(1));assert option==20 and option&2048==0
        home=prior/'scaled-bbc' if folder.parent.parent.name=='scaled-bbc' else prior
        hp=home/'simulations'/folder.name/(folder.name+'_HEAD.FITS')
        with fits.open(hp) as f:
            h=f[1].data;lookup={str(row['SNID']).strip():i for i,row in enumerate(h)};ix=[lookup[c] for c in m.index]
            simpeak=h['SIM_PEAKMJD'][ix].astype(float);zdiff=np.abs(h['REDSHIFT_FINAL'][ix]-h['SIM_REDSHIFT_CMB'][ix])
            peakdiff=np.abs(m.peak_mjd_measured.to_numpy()-simpeak)
        checks.append({'campaign':folder.name,'classified':len(m),'stored_peak_vs_native_PKMJDINI_max_days':delta,'peak_estimator_option':option,'uses_simulation_truth_peak_bit2048':False,'measured_peak_differs_from_true_by_gt_0_01day':int((peakdiff>.01).sum()),'median_absolute_measured_minus_true_peak_days':float(np.median(peakdiff)),'measured_redshift_differs_from_truth':int((zdiff>1e-7).sum())})
        total+=len(m)
        for q in [p,prep,nativepath,fitpath,folder/'fit.nml',hp]:inputs[str(q.relative_to(ROOT))]=sha(q)
    for name in ['model.pt','cli_args.json','data_norm.json']:
        old=prior/'SNDATA_ROOT/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19'/name
        new=WORK/'SNDATA_ROOT_2026-04-10/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19'/name
        assert sha(old)==sha(new);inputs[str(old.relative_to(ROOT))]=sha(old);inputs[str(new.relative_to(ROOT))]=sha(new)
    source=ROOT/'studies/host_ages/code/survey_physics/classify.py';scaled=source.parent/'scaled_classifier.py';native=prior/'SNANA/src/snana.F90'
    result={'status':'passed_no_prior_simulated_truth_peak_substitution','code_sha256':sha(__file__),'campaigns':len(checks),'classified_occurrences':total,'checks':checks,'inputs':inputs,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [source,scaled,native]},'interpretation':['Historical classifiers consume fit.SNANA.TEXT PKMJDINI, with option20=16(maxfluxclump)+4(noabort), not bit2048(SIM_PEAKMJD).','Source passes observed HEADREDSHIFT_FINAL and its error, not SIM_REDSHIFT_CMB.','Modelweights, modelarguments and normalization files match the newly closed observed-data model byteforbyte.','No old scientificresult needs replacement for the public-HEAD timing error discovered in the new diagnostic baseline.','This does not retroactively establish production retraining, CC mixture or population-selection validity.']}
    (RESULTS/'prior-classifier-audit.json').write_text(json.dumps(result,indent=2)+'\n');print('campaigns',len(checks),'occurrences',total)
if __name__=='__main__':main()
