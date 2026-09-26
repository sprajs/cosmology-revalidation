"""Source-law support only: no observed fluxes or fitted posteriors."""
from pathlib import Path
import ast,json,hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar,brentq
P=Path(__file__).resolve().parent;R=P.parents[3]
S=R/'runs/research_2026_09_26/bayesn_distance_identification/official-code/bayesn'
source=S/'spline_utils.py';tree=ast.parse(source.read_text());wanted={'invKD_irr','spline_coeffs_irr'}
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in wanted]
assert len(nodes)==2;ns={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),ns)
xk=np.array([0,1e4/26500,1e4/12200,1e4/6000,1e4/5470,1e4/4670,1e4/4110,1e4/2700,1e4/2600])
wave=np.linspace(3000,18500,15501);X=ns['spline_coeffs_irr'](1e4/wave,xk,ns['invKD_irr'](xk),allow_extrap=False)
def knots(rv):
    c2=-.824+4.717/rv;c1=2.030-3.007*c2
    d=xk[7:]**2/((xk[7:]**2-4.596**2)**2+(.99*xk[7:])**2)
    return np.array([-rv,.26469*rv/3.1-rv,.82925*rv/3.1-rv,-.422809+1.00270*rv+2.13572e-4*rv**2-rv,-.051354+1.00216*rv-7.35778e-5*rv**2-rv,.700127+1.00184*rv-3.32598e-5*rv**2-rv,1.19456+1.01707*rv-5.46959e-3*rv**2+7.97809e-4*rv**3-4.45636e-5*rv**4-rv,*(c1+c2*xk[7:]+3.23*d)])
def min_law(rv):
    sp=CubicSpline(xk,knots(rv),bc_type='natural')
    # Exact stationary-point search for each cubic in inverse wavelength.
    candidates=np.r_[1e4/18500,1e4/3000,sp.derivative().roots(extrapolate=False)]
    candidates=candidates[(candidates>=1e4/18500)&(candidates<=1e4/3000)]
    a=1+sp(candidates)/rv;k=np.argmin(a)
    return float(a[k]),float(1e4/candidates[k])
rows=[];error=0
for rv in [.5,.6,.7,1,1.2,2.886,5]:
    y=knots(rv);a=1+X@y/rv;sp=CubicSpline(xk,y,bc_type='natural');error=max(error,float(np.max(abs(a-(1+sp(1e4/wave)/rv)))))
    roots=sp.solve(-rv,extrapolate=False);roots=roots[(roots>=1e4/18500)&(roots<=1e4/3000)]
    mn,at=min_law(rv);rows.append(dict(RV=rv,min_A_over_AV=mn,min_at_angstrom=at,negative_grid_min=float(wave[a<0].min()) if (a<0).any() else None,negative_grid_max=float(wave[a<0].max()) if (a<0).any() else None,zero_crossings_angstrom=sorted((1e4/roots).tolist())))
assert error<1e-12
critical=brentq(lambda rv:min_law(rv)[0],.6,.7,xtol=1e-14)
paths=[source,S/'bayesn_model.py',R/'runs/research_2026_09_26/bayesn_distance_identification/model-files/BAYESN.M20/BAYESN.YAML',Path(__file__)]
result=dict(status='Pure source-law support check; no observations or posterior masses evaluated',wavelength_domain_angstrom=[3000,18500],AV_convention='AV>0 passive extinction screen; attenuation10^(-.4 A_lambda)',source_operator='get_spectra F99 knots and natural-spline interpolation; source spline functions independently compared with SciPy CubicSpline',source_spline_vs_scipy_max_error=error,rows=rows,critical_RV_for_nonnegative_over_declared_domain=critical,scope='Not a posterior statement, empirical validity certificate, or universal physical lower bound on RV. It is a support limit for this specified empirical law over this wavelength domain. Do not clamp curves or retrospectively alter a paper-reproduction prior.',sha256={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(P/'bayesn-dust-support.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
