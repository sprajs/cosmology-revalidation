"""Optional exact algebraic compression of fixed ACT/Planck lensing responses.

No bandpowers, covariance, response terms, theory requests or priors change.
Only matrix multiplication order changes. The active modern adapter is untouched.
"""
import numpy as np
from act_dr6_lenslike.act_dr6_lenslike import ACTDR6LensLike,standardize,pp_to_kk


class BinnedResponse:
    def __init__(self,data):
        assert data['include_planck'] and data['likelihood_corrections']
        self.data=data;self.blocks=[]
        fid=np.asarray(data['fiducial_cl_kk'])
        for suffix,key in [('', 'binmat_act'),('_planck','binmat_planck')]:
            B=np.asarray(data[key]);scale=-2*fid.copy()
            scale[2:]/=data['fAL'+suffix][2:]
            # Native norm response divides only ell>=2; preserve ell0/1 exactly.
            weighted=B*scale[None,:]
            ops={s:B@data['dN1_'+s+suffix]+weighted@data['dAL_dC'+suffix][i]
                 for i,s in enumerate(['tt','ee','bb','te'])}
            self.blocks.append((B,B@data['dN1_kk'+suffix],ops))

    def predict(self,kk,cl):
        d=self.data;delta_kk=kk-d['fiducial_cl_kk'];pieces=[]
        for B,K,ops in self.blocks:
            row=B@kk+K@delta_kk
            for s in ['tt','ee','bb','te']:row=row+ops[s]@(cl[s]-d['fiducial_cl_'+s])
            pieces.append(row)
        return np.concatenate(pieces)

    def loglike(self,kk,cl,return_theory=False):
        b=self.predict(kk,cl);res=self.data['data_binned_clkk']-b
        value=float(-.5*(res@(self.data['cinv']@res)))
        return (value,b) if return_theory else value


class FastACTDR6LensLike(ACTDR6LensLike):
    """Restricted to the already validated current modern configuration."""
    def initialize(self):
        forbidden=['lens_only','no_like_corrections','limber','varying_cmb_alens',
                   'act_cmb_rescale','act_calib','no_actlike_cmb_corrections','mock']
        for name in forbidden:
            if getattr(self,name):raise ValueError(f'Exact compressed wrapper does not support {name}=True')
        if self.variant!='actplanck_baseline':raise ValueError('Only validated actplanck_baseline supported')
        super().initialize();self.binned_response=BinnedResponse(self.data)

    def loglike(self,cl,**params_values):
        ell=cl['ell'];kk=standardize(ell,pp_to_kk(cl['pp'],ell),self.trim_lmax)
        cmb={s:standardize(ell,cl[s],self.trim_lmax) for s in ['tt','ee','bb','te']}
        return self.binned_response.loglike(kk,cmb)


def use_fast_lensing(info):
    """Explicit opt-in on a caller-owned config; never modify active adapter."""
    key='act_dr6_lenslike.ACTDR6LensLike'
    info['likelihood'][key]=dict(info['likelihood'][key],external=FastACTDR6LensLike)
    return info
