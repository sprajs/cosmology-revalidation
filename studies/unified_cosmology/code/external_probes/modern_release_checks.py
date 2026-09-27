"""Execute ACT authors' released CAMB reference cases with pinned dependencies."""
import json
import numpy as np
from cobaya.model import get_model
from adapter import PACKAGES
from acquire import sha,HERE,RESULTS

def main():
    out={}
    act={'packages_path':str(PACKAGES),'params':{'ombh2':.022,'omch2':.117,'ns':.96,'As':2e-9,'tau':.065,'cosmomc_theta':104.09e-4},
         'theory':{'camb':{'extra_args':{'lmax':9000,'lens_potential_accuracy':8,'min_l_logl_sampling':6000}}},
         'likelihood':{'act_dr6_cmbonly.ACTDR6CMBonly':{'params':{'A_act':1.,'P_act':1.}},
                       'act_dr6_cmbonly.PlanckActCut':{'params':{'A_planck':1.}}}}
    m=get_model(act);log=m.loglikes()[0]
    out['ACT_DR6_official_reference']={'actual_ACT_loglike':float(log[0]),'expected_ACT':-395.48,
        'actual_ACT_Planckcut_loglike':float(sum(log)),'expected_ACT_Planckcut':-962.27,
        'official_absolute_tolerance':.01,
        'official_numpy_isclose_relative_tolerance':1e-5,
        'passed':bool(np.isclose(log[0],-395.48,atol=.01,rtol=1e-5)
                      and np.isclose(sum(log),-962.27,atol=.01,rtol=1e-5))}
    lens={'packages_path':str(PACKAGES),
      'params':{'ombh2':.02219218,'omch2':.1203058,'tau':.06574325,'ns':.9625356,'H0':67.02393,'As':2.15086031154146e-9,'omnuh2':.00064},
      'theory':{'camb':{'extra_args':{'lmax':10000,'lens_margin':1250,'lens_potential_accuracy':4,
                  'AccuracyBoost':1,'lSampleBoost':1,'lAccuracyBoost':1,'halofit_version':'mead2016'}}},
      'likelihood':{'act_dr6_lenslike.ACTDR6LensLike':{'variant':'actplanck_baseline','lens_only':False}}}
    mm=get_model(lens);ll=float(mm.loglikes()[0][0]);chi=-2*ll
    out['ACT_Planck_lensing_official_reference']={'actual_chi2':chi,'expected_chi2':21.20,
        'official_decimal_places':1,'passed':bool(abs(chi-21.20)<.05)}
    out['status']='passed' if all(v['passed'] for v in out.values()) else 'released_reference_discrepancy'
    out['code_sha256']=sha(__file__)
    out['scope']='Original package reference points use their own fixed cosmological/accuracy settings; discrepancies are retained rather than retuning test values.'
    (RESULTS/'modern-release-checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
