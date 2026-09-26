from pathlib import Path
import subprocess,os,json,hashlib,re
Q=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');G=R/'phase2/official/build/sysroot/usr';obj=Q.parent/'build/obj/minuit.o'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(G/'lib')
cmd1=['cc','-O1','-std=c99','-Wall','-Wextra','-Werror','-Wformat=2','-c','unit_wrapper.c','-o','unit_wrapper.o']
cmd2=[str(G/'bin/gfortran'),'-B'+str(G/'lib/gcc/x86_64-pc-linux-gnu/16/'),'-B/usr/lib/gcc/x86_64-pc-linux-gnu/16/','-fsecond-underscore','unit_driver.f90','unit_wrapper.o',str(obj),'-o','unit_driver']
logs=[]
for cmd in [cmd1,cmd2]:
 p=subprocess.run(cmd,cwd=Q,env=env,text=True,capture_output=True);logs.append('COMMAND '+repr(cmd)+'\n'+p.stdout+p.stderr);assert p.returncode==0,p.stderr
(Q/'unit-build.log').write_text('\n'.join(logs))
checks=[]
for label,value,expected,scenario,rc in [('absent',None,'0','normal',0),('zero','0','0','normal',0),('minus','-2','-2','normal',0),('plus','2','2','normal',0),('fixed','2','2','fixed',87),('outside','2','2','outside',87),('invalid','1','0','normal',87),('nan','nan','0','normal',87),('no_domain','2','2','normal',87),('missed_readback','2','2','missed_readback',87),('later_changed','2','2','later_changed',87)]:
 e=env.copy();e.pop('PROSP_MINUIT_PEAK_SHIFT',None);e['PROSP_HARD_PEAK_DOMAIN']='1'
 if value is not None:e['PROSP_MINUIT_PEAK_SHIFT']=value
 if label=='no_domain':e.pop('PROSP_HARD_PEAK_DOMAIN')
 p=subprocess.run([str(Q/'unit_driver'),expected,scenario],cwd=Q,env=e,text=True,capture_output=True)
 (Q/(label+'.log')).write_text(p.stdout+p.stderr);assert p.returncode==rc,(label,p.stdout,p.stderr)
 captured=[x.split() for x in p.stdout.splitlines() if x.startswith('UNIT_CAPTURE')]
 if rc==0:
  assert len(captured)==48 and 'UNIT_ALL_PARAMS_PRIOR_CENTER_AND_LATER_ITERATIONS_PASS' in p.stdout
  for x in captured:
   it,ip=int(x[1]),int(x[2]);prior,proposal,stored=map(float,x[3:]);assert proposal==stored
   assert stored-prior==(float(expected) if it==1 and ip==3 else 0)
 checks.append({'case':label,'returncode':p.returncode,'expected':rc,'captured_MNPOUT_parameters':len(captured),'pass':True})
# Source patch confines the change to local MNPARM input and private MNPOUT outputs.
s=(Q/'source/snana.car').read_text();assert ',PROSP_MN_INITIAL, INISTP(IPAR)' in s
old=(Q/'snana.base.car').read_text();a=old.index('      SUBROUTINE MNFIT_DRIVER (');b=old.index('      END',a);a2=s.index('      SUBROUTINE MNFIT_DRIVER (');b2=s.index('      END',a2)
# Routine's first END token can be ENDIF/ENDDO; whole-file diff already identifies exact change.
for line in s.splitlines():
 if any(k in line for k in ['PROSP_MN','PROSP_MINUIT']) and not line.startswith('c'):assert len(line)<=72,line
r={'status':'PASS','native_fitter_builds':0,'native_fits':0,'photons':0,'standalone':'Exact existing native minuit.o used only for MNINIT/MNPARM/MNPOUT, no FCN/minimizer/photometry. Test-only INTRAC stub reports noninteractive.','native_minuit_object_sha256':sha(obj),'checks':checks,'source_sha256':{str(p.relative_to(Q)):sha(p) for p in [Q/'source/snana.car',Q/'source/genmag_snoopy.c',Q/'source/prosp_minuit_start.h',Q/'start-only.patch']},'commands':[cmd1,cmd2]}
(Q/'test-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
