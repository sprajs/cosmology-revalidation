from pathlib import Path
import urllib.request,json,hashlib
import numpy as np
P=Path(__file__).parent;R=Path.cwd();fd=R/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/filters/CSP_TAMU_20180316'
rows=[];arrays={};sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for label,commit in [('pre2017','5576598c4603c48fe195c4901f572671feed785f'),('current','3b97ef08e48b4d535f587d19a14eb2d0ea3a7f33')]:
 for name in ['V_tel_ccd_atm_ext_1.2.dat','V_LC3014_tel_ccd_atm_ext_1.2.dat','V_LC3009_tel_ccd_atm_ext_1.2.dat']:
  u=f'https://raw.githubusercontent.com/obscode/snpy/{commit}/snpy/filters/filters/LCO/Swope/{name}'
  with urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0'}),timeout=30) as h:b=h.read(200000)
  path=P/f'{label}-{name}';path.write_bytes(b);a=np.loadtxt(path);arrays[label,name]=a;rows.append({'url':u,'path':str(path.relative_to(R)),'sha256':sha(path),'bytes':len(b)})
(P/'registry-acquisition.json').write_text(json.dumps(rows,indent=2)+'\n')
checks=[]
for label in ['pre2017','current']:
 for old,new in [('V_tel_ccd_atm_ext_1.2.dat','V_LC9844_tel_ccd_atm_ext_1.2.dat'),('V_LC3014_tel_ccd_atm_ext_1.2.dat','V_LC3014_tel_ccd_atm_ext_1.2.dat'),('V_LC3009_tel_ccd_atm_ext_1.2.dat','V_LC3009_tel_ccd_atm_ext_1.2.dat')]:
  a=arrays[label,old];b=np.loadtxt(fd/new);x=np.unique(np.r_[a[:,0],b[:,0]]);y=np.interp(x,a[:,0],a[:,1],left=0,right=0);z=np.interp(x,b[:,0],b[:,1],left=0,right=0)
  yn=y/np.trapezoid(x*y,x);zn=z/np.trapezoid(x*z,x)
  checks.append({'registry':label,'registry_file':old,'comparison_2018_file':new,'same_numeric_array':a.shape==b.shape and bool(np.array_equal(a,b)),'same_grid':a.shape==b.shape and bool(np.array_equal(a[:,0],b[:,0])),'max_raw_transmission_difference':float(np.max(abs(y-z))),'normalized_photon_halfL1':float(.5*np.trapezoid(x*abs(yn-zn),x))})
(P/'registry-comparison.json').write_text(json.dumps({'checks':checks,'scope':'Filter source/grid comparison only; no SN or stellar synthetic photometry, no observed magnitude outcomes'},indent=2)+'\n');print(json.dumps(checks,indent=2))
