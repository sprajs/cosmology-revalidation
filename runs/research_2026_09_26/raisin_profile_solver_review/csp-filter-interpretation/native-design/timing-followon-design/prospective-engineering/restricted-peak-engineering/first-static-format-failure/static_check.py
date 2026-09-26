from pathlib import Path
import csv,json,subprocess,os,hashlib
import numpy as np
P=Path(__file__).resolve().parent
rows=list(csv.DictReader((P.parent/'inputs-v2/cadence-ledger.csv').open()))
mjd=[float(r['mjd']) for r in rows if r['keep']=='True'];assert len(mjd)==117
h='''#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define IPAR_GRIDGEN_TREST 0
#define IPAR_GRIDGEN_SHAPEPAR 1
struct {int NBIN[2];double VALUE[2][3];} SNGRID_SNOOPY={{2,2},{{0,-20,70},{0,.69999998807907104,1.2999999523162842}}};
#include "build/src/prosp_peak_domain.h"
int main(int argc,char **argv){
 char cid[]="synthetic";int iter=1,n=117,model=7,snoopy=7,filt=2,epoch=1;
 double mjd[117]={TIMES};
 double z=.453000009059906,step=0,off=0,initial=57707.80078125;
 double klo=-20,khi=85,zmax=.69999998807907104,lo=-1,hi=1;
 double phase=0,tobs=0,D=42,shape=1,AV=0,RV=1.5180000066757202;
 const char *a=argc>1?argv[1]:"normal";
 if(!strcmp(a,"disabled")){n=-100;z=NAN;klo=NAN;model=-1;
  prosp_peak_bounds__(cid,&iter,&n,mjd,&z,&step,&off,&initial,&klo,&khi,&zmax,&model,&snoopy,&lo,&hi);
  prosp_peak_fcn_guard__(cid,&iter,&initial,&off,&z);
  prosp_physical_mean_guard__(cid,&iter,&filt,&epoch,&phase,&tobs,&z,&D,&shape,&AV,&RV,&klo,&khi,&zmax,&model,&snoopy);
  if(lo!=-1||hi!=1)return 5;puts("DISABLED_EXACT_NO_WRITE");return 0;}
 if(!strcmp(a,"offset")){off=57000;initial-=off;}
 if(!strcmp(a,"narrow")){klo=-18;khi=65;}
 if(!strcmp(a,"count"))n=116;
 if(!strcmp(a,"free_z"))step=.01;
 prosp_peak_bounds__(cid,&iter,&n,mjd,&z,&step,&off,&initial,&klo,&khi,&zmax,&model,&snoopy,&lo,&hi);
 if(!strcmp(a,"changed"))mjd[0]+=1;
 iter=2;
 prosp_peak_bounds__(cid,&iter,&n,mjd,&z,&step,&off,&initial,&klo,&khi,&zmax,&model,&snoopy,&lo,&hi);
 if(!strcmp(a,"fcn_out"))initial=nextafter(hi,INFINITY);
 prosp_peak_fcn_guard__(cid,&iter,&initial,&off,&z);
 int k,i;for(k=0;k<2;k++)for(i=0;i<n;i++){
  double peak=k?hi:lo;tobs=(mjd[i]-off)-peak;phase=tobs/(1+z);
  prosp_physical_mean_guard__(cid,&iter,&filt,&epoch,&phase,&tobs,&z,&D,&shape,&AV,&RV,&klo,&khi,&zmax,&model,&snoopy);
 }
 if(!strcmp(a,"bad_phase"))phase=nextafter(-20,-INFINITY);
 if(!strcmp(a,"bad_shape"))shape=1.4;
 if(!strcmp(a,"nan"))D=NAN;
 prosp_physical_mean_guard__(cid,&iter,&filt,&epoch,&phase,&tobs,&z,&D,&shape,&AV,&RV,&klo,&khi,&zmax,&model,&snoopy);
 puts("SYNTHETIC_PASS");return 0;
}
'''.replace('TIMES',','.join(format(v,'.17g') for v in mjd))
(P/'synthetic_harness.c').write_text(h)
build=subprocess.run(['cc','-O1','-std=c99','-Wall','-Wextra','-Werror','-Wformat=2','synthetic_harness.c','-lm','-o','synthetic_harness'],cwd=P,text=True,capture_output=True)
(P/'synthetic-build.log').write_text(build.stdout+build.stderr)
assert build.returncode==0,build.stderr
checks=[]
for arg,mode,expected in [('disabled',None,0),('disabled','0',0),('normal','1',0),('offset','1',0),('narrow','1',0),('count','1',86),('free_z','1',86),('changed','1',86),('fcn_out','1',86),('bad_phase','1',86),('bad_shape','1',86),('nan','1',86),('normal','yes',86)]:
 env=os.environ.copy();env.pop('PROSP_HARD_PEAK_DOMAIN',None);env.pop('PROSP_FIT_SUPPORT',None)
 if mode is not None:env['PROSP_HARD_PEAK_DOMAIN']=mode
 p=subprocess.run([str(P/'synthetic_harness'),arg],env=env,text=True,capture_output=True)
 name=arg+'-'+str(mode)
 (P/(name+'.log')).write_text(p.stdout+p.stderr)
 checks.append(dict(case=name,returncode=p.returncode,expected=expected,passed=p.returncode==expected))
 assert p.returncode==expected,(name,p.stdout,p.stderr)
# independent NumPy arithmetic; source-driven raw endpoints, both execution precisions.
a=np.array(mjd);z=float(np.float32(.453));lo=float(np.max(a-(1+z)*70));hi=float(np.min(a+(1+z)*20))
e=[]
for x in [lo,hi]:
 r8=(a-x)/(1+z);r4=(a-x).astype(np.float32)/np.float32(np.float32(1)+np.float32(z))
 e.append(dict(peak=x,R8_phase=[float(r8.min()),float(r8.max())],R4_mask_phase=[float(r4.min()),float(r4.max())]))
 assert np.all(r8>=-20)&np.all(r8<=70)&np.all(r4>=-20)&np.all(r4<=70)
source=(P/'build/src/snlc_fit.car').read_text()
lines=source.splitlines();long=[(i+1,l)for i,l in enumerate(lines) if len(l)>72 and not l.startswith(('c','C','*','+','!')) and ('PROSP' in l or 'Trange_KCOR' in l)]
assert not long,long
res={'status':'PASS','native_fits':0,'native_photons':0,'tests':checks,'independent_endpoints':e,'no_inward_change_for_current_metadata':True,'source_sha256':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'build/src/snlc_fit.car',P/'build/src/genmag_snoopy.c',P/'build/src/prosp_peak_domain.h']}}
(P/'static-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
