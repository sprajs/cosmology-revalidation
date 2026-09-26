"""Source-defined inverse print/storage cells. Metadata only, no fitted outcomes."""
from pathlib import Path
import numpy as np,json,hashlib,decimal,csv,gzip,ctypes,re,collections
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');B=P.parent/'baseline_2021';S=P.parent/'snana_v11_04d/source/src';LIB=R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib'
D=decimal.Decimal;libc=ctypes.CDLL(None)
def fmt_c(x,fmt):
 b=ctypes.create_string_buffer(128);libc.snprintf(b,128,('%'+fmt).encode(),ctypes.c_double(float(x)));return b.value.decode()
def candidates(token,fmt):
 t=D(token);unit=D(1).scaleb(t.as_tuple().exponent);lo=t-unit/2;hi=t+unit/2
 x=np.float32(float(lo));a=[]
 if D(float(x))<lo:x=np.nextafter(x,np.float32(np.inf))
 while D(float(x))<=hi:
  if fmt_c(x,fmt)==token:a.append(float(x))
  x=np.nextafter(x,np.float32(np.inf))
  assert len(a)<100000
 assert a,(token,fmt,lo,hi)
 assert all(format(x,fmt)==fmt_c(x,fmt)==token for x in a)
 return a

def parse_simlib():
 out=collections.defaultdict(list);cid=None
 for line in LIB.read_text().splitlines():
  if line.startswith('LIBID:'):cid=line.split()[1]
  if line.startswith('S:'):
   t=line.split();out[cid].append({'MJD':t[1],'band':t[3]})
 return out

def main():
 assert not (P/'precision-ledger.json').exists()
 led=[];simlib=parse_simlib();byobj={}
 for cid in map(str,range(1,9)):
  path=B/'fits/baseline/data/DES_RAISIN_SIM'/f'{cid}.DAT';lines=path.read_text().splitlines();head={line.split(':',1)[0]:line.split(':',1)[1].strip().split()[0] for line in lines if ':' in line and line.split(':',1)[1].strip() and not line.startswith('OBS:')};epoch=0
  for line in lines:
   if not line.startswith('OBS:'):continue
   t=line.split();epoch+=1
   for key,index,fmt in [('MJD',1,'.3f'),('FLUXCAL',4,'.4e'),('FLUXCALERR',5,'.4e')]:
    token=t[index];a=candidates(token,fmt)
    if key=='MJD':
     # R8 native epoch -> R4 archive datum -> %.3f; take inverse R4 cell.
     low=(float(np.nextafter(np.float32(a[0]),np.float32(-np.inf)))+a[0])/2
     high=(float(np.nextafter(np.float32(a[-1]),np.float32(np.inf)))+a[-1])/2
     low=float(np.nextafter(low,np.inf));high=float(np.nextafter(high,-np.inf))
     # Tie-open endpoints guarantee the actual formatted archive token.
     assert fmt_c(np.float32(low),fmt)==token and fmt_c(np.float32(high),fmt)==token
     matches=[r for r in simlib[head['SIM_LIBID']] if r['band']==t[2] and low<=float(r['MJD'])<=high]
    else:
     low,high=a[0],a[-1];matches=[]
     if key=='FLUXCALERR':assert low>0
    led.append({'CID':cid,'epoch':epoch,'band':t[2],'field':key,'token':token,'printf_format':fmt,'float32_state_count':len(a),'float32_low':a[0],'float32_high':a[-1],'input_low':low,'input_high':high,'baseline_input':float(token),'SIMLIB_candidates':matches,'native_storage':'R8 time before archival R4; full precision interval' if key=='MJD' else 'R4 likelihood data'})
  for key,fmt in [('PEAKMJD','.4f'),('REDSHIFT_HELIO','.5f'),('REDSHIFT_FINAL','.5f'),('MWEBV','.5e')]:
   token=head[key];a=candidates(token,fmt);led.append({'CID':cid,'epoch':None,'band':None,'field':key,'token':token,'printf_format':fmt,'float32_state_count':len(a),'float32_low':a[0],'float32_high':a[-1],'input_low':a[0],'input_high':a[-1],'baseline_input':float(np.float32(float(token))),'SIMLIB_candidates':[],'native_storage':'R4 header state; no pre-storage subcell influences native likelihood'})
  byobj[cid]={'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'epochs':epoch,'source_SIM_LIBID':head['SIM_LIBID']}
 # Synthetic inverse-format checks, including negative and scientific values.
 tests=[]
 for val,fmt in [(57384.3125,'.3f'),(-1.2345,'.4e'),(0.23987,'.4e'),(0.46693,'.5f'),(0.00897533,'.5e')]:
  token=fmt_c(np.float32(val),fmt);a=candidates(token,fmt);assert float(np.float32(val)) in a;tests.append({'value':val,'token':token,'states':len(a)})
 time=[x for x in led if x['field']=='MJD'];heads=[x for x in led if x['epoch'] is None]
 out={'scope':'Metadata/source precision audit only. No new fits or outcome scores.','MJDOFF':0.0,'cohort':list(byobj),'objects':byobj,'coordinates':led,'source_gate':{'R4_epoch_before_writer':'snlc_fit.car:10335–10337 and4684,18118; VMJD REAL8 declared17996','text_writer':'sntools_output_text.c:1573 (MJD %.3f),1584–1593 (flux/error %.4le)','MJDOFF':'snana.car:352 PARAMETER MJDOFF=0.0','header_storage':'snana.car:815–886 REAL, sntools_dataformat_text.c:1770–1786 and1901–1903','header_print':'sntools_output_text.c:477–501; MJD4dec,z5dec,smallnoninteger5e'},'unit_tests':tests,'summary':{'epoch_coordinates':len(time)*3,'header_coordinates':len(heads),'unique_R4_time_states':all(x['float32_state_count']==1 for x in time),'time_cell_width_min':min(x['input_high']-x['input_low'] for x in time),'time_cell_width_max':max(x['input_high']-x['input_low'] for x in time),'source_SIMLIB_unique_times':sum(len(x['SIMLIB_candidates'])==1 for x in time),'total_times':len(time),'header_active_coordinates':sum(x['input_low']!=x['input_high'] for x in heads),'all_peak_unique_R4':all(x['float32_state_count']==1 for x in heads if x['field']=='PEAKMJD')}}
 (P/'precision-ledger.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
