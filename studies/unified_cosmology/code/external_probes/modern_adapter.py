"""Contemporary CMB comparison with explicit release/convention choices.

Run in .modern-venv after modern_acquire.py. This is an executable public-data
combination, not an exact reproduction of unavailable author configuration.
"""
import argparse
from adapter import external_info


def modern_info(model='cpl',bao=True,calibration='official_planck',act_lmax=6500):
    from candl.interface import CandlCobayaLikelihood
    assert calibration in ('official_planck','paper_literal')
    assert act_lmax in (6500,8500)
    info=external_info('lite',model,bao=bao,lensing=False)
    info['likelihood'].pop('planck_2018_highl_plik.TTTEEE_lite_native')
    info['likelihood'].update({
        'act_dr6_cmbonly.PlanckActCut':{},
        'act_dr6_cmbonly.ACTDR6CMBonly':{
            'input_file':'dr6_data_cmbonly.fits','version':'v1.0',
            'lmax_theory':9000,'ell_cuts':{p:[600,act_lmax] for p in ('TT','TE','EE')}},
        'SPT2025_TnE':{'external':CandlCobayaLikelihood,
            'data_set_file':'spt_candl_data.SPT3G_D1_TnE_lite',
            'clear_internal_priors':True,'feedback':False},
        'act_dr6_lenslike.ACTDR6LensLike':{
            'variant':'actplanck_baseline','lens_only':False,'lmax':4000,
            'version':'v1.2','apply_hartlap':True,'nsims_act':792.,'nsims_planck':400.,
            'no_like_corrections':False,'mock':False},
        'SPT2023_lensing':{'external':CandlCobayaLikelihood,
            'data_set_file':'candl_data.SPT3G_2018_Lens_and_CMB',
            'lensing':True,'clear_internal_priors':True,'feedback':False},
    })
    p=info['params']
    p['A_planck']={'prior':{'dist':'norm','loc':1.,'scale':.0025},'ref':1.,'proposal':.0005}
    p['P_act']={'prior':{'dist':'norm','loc':1.,'scale':.003},'ref':1.,'proposal':.0006}
    p['A_act']={'value':'lambda A_planck: A_planck' if calibration=='official_planck'
                else 'lambda P_act: P_act'}
    for name,sigma in [('Tcal',.0036),('Ecal',.0095)]:
        p[name]={'prior':{'dist':'norm','loc':1.,'scale':sigma},'ref':1.,'proposal':sigma/3}
    p['A_fg']={'prior':{'min':0.,'max':2.},'ref':.5,'proposal':.1}
    info['theory']['camb']['extra_args'].update(lmax=9000,lens_margin=1250,
        lens_potential_accuracy=4,AccuracyBoost=1,lSampleBoost=1,lAccuracyBoost=1,
        halofit_version='mead2016')
    return info


if __name__=='__main__':
    from cobaya.yaml import yaml_dump
    p=argparse.ArgumentParser()
    p.add_argument('--model',choices=['lcdm','cpl'],default='cpl')
    p.add_argument('--calibration',choices=['official_planck','paper_literal'],default='official_planck')
    p.add_argument('--act-lmax',type=int,choices=[6500,8500],default=6500)
    args=p.parse_args()
    print(yaml_dump(modern_info(args.model,calibration=args.calibration,act_lmax=args.act_lmax)))
