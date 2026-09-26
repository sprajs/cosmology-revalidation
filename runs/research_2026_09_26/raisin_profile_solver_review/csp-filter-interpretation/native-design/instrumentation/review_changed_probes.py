from pathlib import Path
import json,hashlib,gzip
import numpy as np
from astropy.io import fits
from parse_native import parse,ROW
O=Path(__file__).resolve().parent;P=O/'changed-pilot-probes-v3';ROOT=Path.cwd();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();proto=json.loads((P/'protocol.json').read_text())
for name,h in proto['hashes'].items():assert sha(ROOT/name)==h,name
arms={a['name']:parse(P/'fits'/a['name']/'fit.log') for a in json.loads((P/'derived-fixed-probes.json').read_text())['arms']};last=lambda bs:{b['CID']:b for b in bs};nom=last(arms['nominal']);identitycols=[ROW.index(k) for k in ['MJD','dataF','data_error']]
def identity(a,b):return np.array_equal(a['array'][:,identitycols],b['array'][:,identitycols]) and [r['band'] for r in a['rows']]==[r['band'] for r in b['rows']]
def arithmetic(b):
 o=b['objective'];y=b['array'][:,4];f=b['array'][:,2];W=b['W'];res=y-f;Q=float(res@W@res);amp=float(f@W@y/(f@W@f));delta=float(-2.5*np.log10(amp)) if amp>0 else None
 return dict(CID=b['CID'],ITER=b['ITER'],n=b['NFITDATA'],D=o['D'],totalQ=o['totalQ'],dataQ=Q,priorQ=o['priorQ'],sigmaQ=o['sigmaQ'],Qclosure=Q+o['priorQ']+o['sigmaQ']-o['totalQ'],frozen_C_delta_D=delta,nonpositive_model_count=int(np.sum(f<=0)),W_symmetry_abs=float(abs(W-W.T).max()),C_min_eigenvalue=float(np.linalg.eigvalsh(b['C']).min()))
allq={a:[arithmetic(b) for b in bs] for a,bs in arms.items()};starts=[]
for name,shift in [('start_minus',-.2),('start_plus',.2)]:
 for cid,b in last(arms[name]).items():
  first=next(x for x in arms[name] if x['CID']==cid);basefirst=next(x for x in arms['nominal'] if x['CID']==cid);db=first['entry']['D_entry']-basefirst['entry']['D_entry'];diff=b['objective']['D']-nom[cid]['objective']['D'];r=dict(arm=name,CID=cid,baseline_entry_D=basefirst['entry']['D_entry'],actual_entry_D=first['entry']['D_entry'],actual_entry_shift=db,final_D_difference=diff,final_identity=identity(b,nom[cid]));r['pass']=abs(db-shift)<1e-12 and abs(diff)<=.001 and r['final_identity'];starts.append(r)
scaling=[]
for cid in nom:
 center=last(arms[f'{cid}_center'])[cid]
 for label in ['minus','plus']:
  b=last(arms[f'{cid}_{label}'])[cid];D=b['objective']['D']-center['objective']['D'];scale=10**(-.4*D);f=center['array'][:,2];g=b['array'][:,2];err=float(np.max(abs(g-scale*f))/max(np.max(abs(scale*f)),1));perrow=float(np.max(abs(g/(scale*f)-1)));r=dict(CID=cid,arm=label,actual_center_D=center['objective']['D'],actual_D=b['objective']['D'],actual_delta_D=D,scale=scale,normalized_max_error=err,max_perrow_relative_error=perrow,identity=identity(b,center),C_relative_Frobenius_change=float(np.linalg.norm(b['C']-center['C'])/np.linalg.norm(center['C'])));r['pass']=r['identity'] and perrow<=2e-8;scaling.append(r)
 # Float32 center versus original fully precise output coordinates.
 f=nom[cid]['array'][:,2];g=center['array'][:,2];scale=10**(-.4*(center['objective']['D']-nom[cid]['objective']['D']));r=dict(CID=cid,arm='center_vs_nominal',actual_center_D=nom[cid]['objective']['D'],actual_D=center['objective']['D'],actual_delta_D=center['objective']['D']-nom[cid]['objective']['D'],scale=scale,normalized_max_error=float(np.max(abs(g-scale*f))/max(np.max(abs(scale*f)),1)),max_perrow_relative_error=float(np.max(abs(g/(scale*f)-1))),identity=identity(center,nom[cid]),C_relative_Frobenius_change=float(np.linalg.norm(center['C']-nom[cid]['C'])/np.linalg.norm(nom[cid]['C'])));r['pass']=r['identity'] and r['max_perrow_relative_error']<=2e-8;scaling.append(r)
iterations=[]
for name in ['nominal','start_minus','start_plus']:
 for cid in nom:
  bs=[b for b in arms[name] if b['CID']==cid];a,b=bs[-2:];r=dict(arm=name,CID=cid,delta_D=b['objective']['D']-a['objective']['D'],identity=identity(a,b),C_relative_Frobenius_change=float(np.linalg.norm(b['C']-a['C'])/np.linalg.norm(a['C'])),W_relative_Frobenius_change=float(np.linalg.norm(b['W']-a['W'])/np.linalg.norm(a['W'])),final_frozen_C_delta_D=arithmetic(b)['frozen_C_delta_D']);r['pass']=r['identity'] and abs(r['delta_D'])<=.001 and abs(r['final_frozen_C_delta_D'])<=.001;iterations.append(r)
model=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/snoopy.B18/SNooPy_B18.fits';kcor=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits'
with fits.open(model) as h:phase=[float(h['TREST-GRID'].data['TREST'].min()),float(h['TREST-GRID'].data['TREST'].max())];shape=[float(h['LUMI-GRID'].data['STRETCH'].min()),float(h['LUMI-GRID'].data['STRETCH'].max())]
with fits.open(kcor) as h:
 t=h['FilterTrans'].data;wave=t['wavelength (A)'];bandbounds={band:[float(wave[t['CSP-'+band]!=0].min()),float(wave[t['CSP-'+band]!=0].max())] for band in sorted({r['band'] for bs in arms.values() for b in bs for r in b['rows']})}
support=[]
for name,bs in arms.items():
 for i,b in enumerate(bs):
  x=b['array'];z=x[0,6];bounds=np.array([bandbounds[r['band']] for r in b['rows']])/(1+z);lo=float(bounds.min());hi=float(bounds.max());r=dict(arm=name,callback=i,CID=b['CID'],ITER=b['ITER'],n=b['NFITDATA'],phase_min=float(x[:,1].min()),phase_max=float(x[:,1].max()),shape=b['objective']['shape'],rest_nonzero_throughput_min=lo,rest_nonzero_throughput_max=hi,rest_native_mean_min=float(x[:,8].min()),rest_native_mean_max=float(x[:,8].max()));r['phase_inside_source_grid']=r['phase_min']>=phase[0] and r['phase_max']<=phase[1];r['shape_inside_source_grid']=shape[0]<=r['shape']<=shape[1];r['all_nonzero_throughput_inside_native_lambda_screen']=bool(lo>=x[:,11].min() and hi<=x[:,12].max());support.append(r)
science=lambda p:[l for l in p.read_text().splitlines() if l.startswith(('SN:','VARNAMES:'))];orig=ROOT/'runs/research_2026_09_26/csp_native_filter_response/fits/nominal';default=P/'fits/nominal';control=lambda p:[l for l in p.read_text().splitlines() if l.split() and (l.startswith('VARNAMES:') or len(l.split())>1 and l.split()[1]=='2005hc')]
same=all(control(default/name)==control(orig/name) for name in ['fit.FITRES.TEXT','fit.LCPLOT.TEXT'])
result=dict(protocol_sha256=sha(P/'protocol.json'),unchanged_control_exact_science_and_LCPLOT=same,starts=starts,scaling=scaling,iterations=iterations,all_callback_arithmetic=allq,source_phase_grid=phase,source_shape_grid=shape,observer_nonzero_throughput_bounds=bandbounds,support=support,max_abs_Q_closure=max(abs(q['Qclosure']) for qs in allq.values() for q in qs),all_finite_positive_cov_and_mean=all(q['nonpositive_model_count']==0 and q['C_min_eigenvalue']>0 for qs in allq.values() for q in qs),scope='Native historical numerical gates only. Exact reproduction and numerical stability do not establish goodness of fit, physical dust or calibration-independent distances.')
result['numerical_gate_pass']=bool(same and all(r['pass'] for r in starts+scaling+iterations) and result['max_abs_Q_closure']<1e-8 and result['all_finite_positive_cov_and_mean'] and all(r['phase_inside_source_grid'] and r['shape_inside_source_grid'] and r['all_nonzero_throughput_inside_native_lambda_screen'] for r in support))
baseline=parse(O/'probes-v3/fits/nominal/fit.log'); physical=[]
for a,b in zip(baseline,arms['nominal'],strict=True):
 assert (a['CID'],a['ITER'])==(b['CID'],b['ITER'])
 old=sorted((r['MJD'],r['dataF'],r['data_error'],('j' if a['CID']=='2004ef' and r['band']=='J' else r['band'])) for r in a['rows'])
 new=sorted((r['MJD'],r['dataF'],r['data_error'],r['band']) for r in b['rows'])
 physical.append(dict(CID=a['CID'],ITER=a['ITER'],n_old=a['NFITDATA'],n_new=b['NFITDATA'],identity_with_declared_relabel=old==new))
result['physical_row_checks']=physical;result['numerical_gate_pass']=result['numerical_gate_pass'] and all(r['identity_with_declared_relabel'] for r in physical)
result['default_response']=[dict(CID=cid,delta_D=nom[cid]['objective']['D']-ob['objective']['D'],old_D=ob['objective']['D'],new_D=nom[cid]['objective']['D'],old_totalQ=ob['objective']['totalQ'],new_totalQ=nom[cid]['objective']['totalQ']) for cid,ob in last(baseline).items()]
(P/'result.json').write_text(json.dumps(result,indent=2)+'\n');(P/'manifest.json').write_text(json.dumps({str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file() and not 'SNANA-v11_04k-output' in p.parts and p.name!='manifest.json'},indent=2)+'\n');print(json.dumps({k:result[k] for k in ['numerical_gate_pass','unchanged_control_exact_science_and_LCPLOT','max_abs_Q_closure','starts','scaling','iterations','observer_nonzero_throughput_bounds']},indent=2))
