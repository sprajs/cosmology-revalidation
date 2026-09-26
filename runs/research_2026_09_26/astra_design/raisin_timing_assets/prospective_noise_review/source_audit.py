from pathlib import Path
import hashlib,json,numpy as np
ROOT=Path.cwd(); OUT=ROOT/'runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review'
S=ROOT/'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/source/src'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ranges={'snlc_sim.c':[(827,842),(878,890),(10883,11025),(14818,14840),(16455,16470),(22543,22684),(22741,22895),(23022,23065),(23333,23439),(24661,24703)],'snlc_fit.car':[(4383,4409),(4930,4978),(5294,5312),(5500,5522),(8000,8036),(8060,8072),(10948,11072)],'genmag_snoopy.c':[(108,140),(249,292)],'sntools_modelgrid_read.c':[(514,539)],'sntools_fluxErrModels.c':[(727,763),(853,902)]}
evidence={}
for fn,rr in ranges.items():
 p=S/fn; lines=p.read_text().splitlines(); evidence[str(p.relative_to(ROOT))]={'sha256':sha(p),'excerpts':[{'first':a,'last':b,'lines':lines[a-1:b]} for a,b in rr]}
p=ROOT/'runs/research_2026_09_26/raisin_sign_source/raisin_cosmo/genSimlib.py'
acq=json.loads((p.parent.parent/'simulation-acquisition.json').read_text())['files']; entry=next(x for x in acq if x['source_path']=='raisin_cosmo/genSimlib.py')
blob=hashlib.sha1(b'blob '+str(p.stat().st_size).encode()+b'\0'+p.read_bytes()).hexdigest()
assert sha(p)==entry['sha256'] and blob==entry['git_blob']
evidence[str(p.relative_to(ROOT))]={'sha256':sha(p),'git_blob':blob,'url':entry['url'],'excerpts':[{'first':196,'last':259,'lines':p.read_text().splitlines()[195:259]}]}
(OUT/'source-evidence.json').write_text(json.dumps({'scope':'Source-only review, no native calls or observed likelihoods','source_commit':'10ec91297e4482d593cb5d3d055b10d4aa915071','files':evidence},indent=2)+'\n')
# Independent analytic check of noise-dependent weights in the small-timing expansion.
h=np.array([.3,.75,1.,.85,.55,.2]); h1=np.array([-.07,-.05,.015,.055,.07,.04]); h2=np.array([.008,-.006,-.014,-.004,.003,.006]); A=2.; K=2.5/np.log(10)
C0=np.diag([.05,.04,.035,.04,.05,.065]); W0=np.linalg.inv(C0); n0=np.array([.11,-.12,.08,.09,-.13,.12]); p=np.array([50.,40.,60.,50.,35.,40.]); eta0=.6
S0=h@W0@h; a1=h1@W0@h/S0; b=h1@W0@h1/S0; c=h2@W0@h/S0
rows=[]
for scale in [1e-2,5e-3,2.5e-3,1.25e-3]:
 n=scale*n0; eta=scale*eta0; C=C0+np.diag(n/p); W=np.linalg.inv(C); y=A*h+n; ht=h+eta*h1+.5*eta**2*h2
 est=lambda t:(t@W@y)/(t@W@t)
 exact=-K*np.log(est(ht)/est(h)); dW=-W0@np.diag(n/p)@W0
 da1=(h1@dW@h-a1*(h@dW@h))/S0; e0=h@W0@n/(A*S0); e1=h1@W0@n/(A*S0)
 approx=K*(a1*eta+eta*da1-eta*e1+a1*eta*e0+(b+.5*c-1.5*a1*a1)*eta*eta)
 rows.append({'scale':scale,'exact_paired_D':exact,'second_order':approx,'absolute_remainder':abs(exact-approx),'remainder_over_scale_cubed':abs(exact-approx)/scale**3})
ratios=[rows[i]['absolute_remainder']/rows[i+1]['absolute_remainder'] for i in range(3)]
assert all(7.7<x<8.3 for x in ratios)
(OUT/'weight-timing-algebra.json').write_text(json.dumps({'scope':'Deterministic synthetic algebra only, no SN residual outcomes','rows':rows,'successive_remainder_ratio':ratios,'cubic_remainder_check':True},indent=2)+'\n')
print(json.dumps({'construction_blob_pass':True,'source_files':len(evidence),'synthetic_cubic_ratios':ratios},indent=2))
