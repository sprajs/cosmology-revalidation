from pathlib import Path
import json,csv,hashlib,time
import numpy as np
from scipy.integrate import quad
from astropy.io import fits
P=Path(__file__).resolve().parent; H=P.parent.parent; O=P.parent; D=H/'dark_ramp_execution'
def save(n,x): (P/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvsave(n,rows):
 with (P/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
files=[P/'protocol.json',Path(__file__),D/'fixed-union-bad.npy',D/'fixed-mask-result.json']
names=['idbx41onq','idbx43p7q']; data=[]
for n in names:
 p=O/'work'/n/(n+'_flt.fits');files.append(p)
 with fits.open(p,memmap=False) as f:
  assert f['SCI',1].header['BUNIT']=='COUNTS/S'
  data.append({k:np.array(f[k,1].data) for k in ['SCI','ERR','DQ','SAMP','TIME']})
m=~np.load(D/'fixed-union-bad.npy')[5:1019,5:1019]; assert m.sum()==993750
z=data[1]['SCI'].astype(float)-data[0]['SCI'].astype(float);v=data[0]['ERR'].astype(float)**2+data[1]['ERR'].astype(float)**2
qmap=np.empty(m.shape,'U1');qmap[:507,:507]='B';qmap[:507,507:]='C';qmap[507:,:507]='A';qmap[507:,507:]='D'
def stats(s,label):
 a=z[s];b=v[s];ss=np.sum(a*a);vv=np.sum(b)
 return dict(region=label,pixels=int(a.size),signed_mean=float(a.mean()),sum_squared_difference=float(ss),sum_quoted_pair_variance=float(vv),repeat_to_quoted_ratio=float(ss/vv),summed_variance_excess=float(ss-vv),sum_squared_standardized_difference=float(np.sum(a*a/b)),coherent_mean_fraction=float(a.size*a.mean()**2/ss))
regions=[stats(m,'all_active_fixed_mask')]+[stats(m&(qmap==q),q) for q in 'BCAD'];blocks=[]
for by in range(16):
 for bx in range(16):
  s=np.zeros(m.shape,bool);ys=slice(max(5,by*64)-5,min(1019,(by+1)*64)-5);xs=slice(max(5,bx*64)-5,min(1019,(bx+1)*64)-5);s[ys,xs]=m[ys,xs];blocks.append(stats(s,f'{by},{bx}'))
# Circle/rectangle area: integrate vertical overlap, splitting all boundary intersections.
def area(x,y,r):
 lo=max(x-.5,-r);hi=min(x+.5,r)
 if lo>=hi:return 0.
 if max(abs(x-.5),abs(x+.5))**2+max(abs(y-.5),abs(y+.5))**2<=r*r:return 1.
 if max(abs(x)-.5,0)**2+max(abs(y)-.5,0)**2>=r*r:return 0.
 pts=[lo,hi]
 for yy in [y-.5,y+.5]:
  if abs(yy)<r:
   h=np.sqrt(r*r-yy*yy)
   pts += [u for u in [-h,h] if lo<u<hi]
 pts=sorted(set(pts))
 def f(xx):
  h=np.sqrt(max(0.,r*r-xx*xx));return max(0.,min(y+.5,h)-max(y-.5,-h))
 return sum(quad(f,a,b,epsabs=2e-13,epsrel=2e-13)[0] for a,b in zip(pts[:-1],pts[1:]))
r=[3.118908382066277,9.35672514619883,15.594541910331383];off=np.arange(-17,18)
a=np.array([[area(x,y,r[0]) for x in off] for y in off]);b=np.array([[area(x,y,r[2])-area(x,y,r[1]) for x in off] for y in off]);assert abs(a.sum()-np.pi*r[0]**2)<1e-9
sites=[];weights=[]
for t in json.loads((D/'fixed-mask-result.json').read_text())['coverage']:
 y,x=t['y']-5,t['x']-5;sl=(slice(y-17,y+18),slice(x-17,x+18));ma=a*m[sl];mb=b*m[sl];ok=ma.sum()/a.sum()>=.9 and mb.sum()/b.sum()>=.75;assert ok==t['eligible']
 if not ok:continue
 w=ma-mb*ma.sum()/mb.sum();assert abs(w.sum())<1e-10
 zz=float(np.sum(w*z[sl]));vv=float(np.sum(w*w*v[sl]));sites.append(dict(raw_x=t['x'],raw_y=t['y'],difference_dn_per_nominal_s=zz,quoted_diagonal_pair_variance=vv,squared_difference=zz*zz,aperture_coverage=float(ma.sum()/a.sum()),annulus_coverage=float(mb.sum()/b.sum())));weights.append((t['x'],t['y'],w))
assert len(sites)==247
checks={}
for name,rows in [('calibrated-quadrants.csv',regions),('calibrated-blocks.csv',blocks),('calibrated-apertures.csv',sites)]:
 ref=list(csv.DictReader((O/name).open()));assert len(ref)==len(rows);diff={}
 for k in rows[0]:
  if k=='region':assert [x[k] for x in rows]==[x[k] for x in ref];continue
  err=max(abs(float(x[k])-float(y[k])) for x,y in zip(rows,ref));diff[k]=err
  assert np.allclose([x[k] for x in rows],[float(x[k]) for x in ref],rtol=1e-10,atol=1e-10),(name,k,err)
 checks[name]=diff;csvsave('independent-'+name,rows)
# Fixed categories retain every original selected pixel.
cats=[];keys=np.stack([data[0]['SAMP'][m],data[1]['SAMP'][m],data[0]['DQ'][m].astype(np.uint16),data[1]['DQ'][m].astype(np.uint16)],1)
for q in ['all']+list('BCAD'):
 sel=np.ones(len(keys),bool) if q=='all' else qmap[m]==q
 u,idx=np.unique(keys[sel],axis=0,return_inverse=True);zz=z[m][sel];vv=v[m][sel]
 for i,row in enumerate(u):
  k=idx==i;cats.append(dict(region=q,samp_old=int(row[0]),samp_new=int(row[1]),dq_old=int(row[2]),dq_new=int(row[3]),pixels=int(k.sum()),sum_d2=float(np.sum(zz[k]**2)),sum_V=float(np.sum(vv[k])),sum_d=float(np.sum(zz[k]))))
csvsave('native-category-decomposition.csv',cats)
yy,xx=np.where(m);order=np.lexsort((xx,yy,-z[m]**2));chosen=list(order[:20]);quadorders={}
for q in 'BCAD':
 qo=order[qmap[yy[order],xx[order]]==q];quadorders[q]=qo;chosen+=list(qo[:5])
chosen=list(dict.fromkeys(chosen));rank=np.empty(len(order),int);rank[order]=np.arange(1,len(order)+1)
infl=[];lookup=[]
for i in chosen:
 y,x=int(yy[i]),int(xx[i]);rr=dict(rank=int(rank[i]),quadrant=str(qmap[y,x]),raw_x=x+5,raw_y=y+5,trim_x=x,trim_y=y,difference=float(z[y,x]),variance=float(v[y,x]),d2=float(z[y,x]**2),fraction_full_d2=float(z[y,x]**2/np.sum(z[m]**2)))
 for j,n in enumerate(['old','new']):
  rr.update({n+'_'+k:float(data[j][k][y,x]) for k in data[j]})
 hits=[]
 for sx,sy,w in weights:
  if abs(x+5-sx)<=17 and abs(y+5-sy)<=17:
   ww=float(w[y+5-sy+17,x+5-sx+17])
   if ww: hits.append(dict(raw_x=sx,raw_y=sy,weight=ww,weighted_difference=ww*z[y,x]))
 rr['aperture_hits']=len(hits);infl.append(rr);lookup.append(dict(rank=int(rank[i]),raw_x=x+5,raw_y=y+5,apertures=hits))
csvsave('influence-pixels.csv',infl);save('influence-aperture-weights.json',lookup)
fracs=[]
for q,o in [('all',order)]+list(quadorders.items()):
 den=np.sum(z[yy[o],xx[o]]**2)
 for k in [1,5,10,20]:fracs.append(dict(region=q,K=k,sum_d2=float(np.sum(z[yy[o[:k]],xx[o[:k]]]**2)),fraction=float(np.sum(z[yy[o[:k]],xx[o[:k]]]**2)/den)))
csvsave('influence-cumulative.csv',fracs)
# Selected read histories only; no reranking or changes to original estimator.
hist=[]
for n in names:
 for kind in ['raw','ima']:
  p=O/'work'/n/(n+'_'+kind+'.fits');files.append(p)
  with fits.open(p,memmap=False) as f:
   for ev in range(8,0,-1):
    hdr=f['SCI',ev].header
    arrays={k:f[k,ev].data for k in ['SCI','ERR','DQ','SAMP','TIME']}
    for item in infl:
     y,x=item['raw_y'],item['raw_x'];rec=dict(root=n,kind=kind,rank=item['rank'],raw_x=x,raw_y=y,extver=ev,SAMPNUM=hdr.get('SAMPNUM'),SAMPTIME=hdr.get('SAMPTIME'),BUNIT=hdr.get('BUNIT',''))
     for k,arr in arrays.items():rec[k]=float(arr[y,x]) if arr is not None else float(f[k,ev].header['PIXVALUE'])
     hist.append(rec)
csvsave('selected-read-histories.csv',hist)
save('result.json',dict(pass_arithmetic=True,checks=checks,regions=regions,apertures=dict(n=len(sites),sum_d2=sum(x['squared_difference'] for x in sites),sum_V=sum(x['quoted_diagonal_pair_variance'] for x in sites)),selected_influence_count=len(infl),scope='Independent extraction; posthoc influence only; no changed mask, correction or native execution'))
save('manifest.json',{str(p):sha(p) for p in files+list(P.glob('*.csv'))+[P/'result.json',P/'influence-aperture-weights.json']})
print(json.dumps({'pass':True,'checks':checks,'top5':infl[:5],'categories_nonzero_dq':[x for x in cats if x['region']=='all' and (x['dq_old'] or x['dq_new'])]},indent=2))
