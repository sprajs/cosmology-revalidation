from pathlib import Path
import json,hashlib,re,collections,csv
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=Path(__file__).resolve().parent;S=P.parent/'snana_v11_04d/source/src';sim=R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib';mp=R/'runs/research_2026_09_26/raisin_simulation_assets/author/sim/inputs/DES/DES3YR_SIM_ERRORFUDGES.DAT';inp=R/'runs/research_2026_09_26/raisin_sign_source/sim/inputs/DES/sim_DES_SNOOPY.input'
groups={};maps=[];m=None
for line in mp.read_text().splitlines():
 t=line.split('#')[0].split()
 if not t:continue
 if t[0]=='DEFINE_FIELDGROUP:':groups[t[1]]=t[2]
 if t[0]=='MAPNAME:':m={'name':t[1],'band':'ALL','field':'ALL'}
 if t[0]=='BAND:':m['band']=t[1]
 if 'FIELD:' in t:m['field']=t[t.index('FIELD:')+1]
 if t[0]=='ENDMAP:':m['expanded_field']=groups.get(m['field'],m['field']);maps.append(m);m=None
libs={};cid=None
for line in sim.read_text().splitlines():
 t=line.split('#')[0].split()
 if not t:continue
 if t[0]=='LIBID:':cid=t[1];libs[cid]={'LIBID':cid,'field':None,'bands':collections.Counter(),'NOBS_header':None}
 if cid is not None and 'NOBS:' in t:libs[cid]['NOBS_header']=int(t[t.index('NOBS:')+1])
 if t[0]=='FIELD:':libs[cid]['field']=t[1]
 if t[0]=='S:':libs[cid]['bands'][t[3]]+=1
 if t[0]=='END_LIBID:':cid=None
rows=[]
for l in libs.values():
 assert l['field'] is not None and sum(l['bands'].values())==l['NOBS_header']
 for band,n in sorted(l['bands'].items()):
  matches=[i for i,m in enumerate(maps) if (m['band']=='ALL' or band in m['band']) and (m['expanded_field']=='ALL' or l['field'] in m['expanded_field'])]
  rows.append({'LIBID':l['LIBID'],'actual_FIELD':l['field'],'band':band,'raw_epoch_count':n,'matched_map_indices':matches,'scale_fallback_true':1 if not matches else None,'scale_fallback_reported':1 if not matches else None})
assert all(not x['matched_map_indices'] for x in rows)
evidence={}
for f,ranges in {'sntools_fluxErrModels.c':[(166,181),(730,761),(853,900)],'snlc_sim.c':[(1042,1050),(5800,5824),(14818,14839),(15680,15700),(16450,16468),(22919,22947),(23022,23065)]}.items():
 ls=(S/f).read_text().splitlines();evidence[f]={'sha256':hashlib.sha256((S/f).read_bytes()).hexdigest(),'excerpts':[{'range':[a,b],'lines':[f'{j}: {ls[j-1]}' for j in range(a,b+1)]} for a,b in ranges]}
res={'scope':'Metadata/source matching audit only, no generator execution or photon outcomes.','source_commit':'10ec91297e4482d593cb5d3d055b10d4aa915071','input_FLUXERRMODEL_FILE':next(l.strip() for l in inp.read_text().splitlines() if l.startswith('FLUXERRMODEL_FILE:')),'declared_groups':groups,'maps':maps,'SIMLIB_entries':len(libs),'unique_actual_field_tokens':len({l['field'] for l in libs.values()}),'all_raw_epochs':sum(x['raw_epoch_count'] for x in rows),'tested_band_field_pairs':len(rows),'matching_pairs':0,'fallback':'get_FLUXERRMODEL initializes both returned errors to FLUXERR_IN; index<0 returns directly. Not a fatal mismatch. Fudge SQSIG_F=0 at scale1. No REDCOV declaration in map; defaultNREDCOV0. A configured FLUXERRMODEL_FILE suppresses legacy SIMLIB fluxerr fudge independently of whether individual fields match.','legacy_SIMLIB_fluxerr_keys':[l for l in sim.read_text().splitlines() if 'FLUXERR' in l or 'TEMPLATE' in l],'source_filename_not_execution_proof':True,'source_evidence':evidence,'inputs_sha256':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [sim,mp,inp,Path(__file__)]}}
(P/'map-result.json').write_text(json.dumps(res,indent=2)+'\n');(P/'map-matching.json').write_text(json.dumps(rows,indent=2)+'\n')
with (P/'map-matching.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({k:v for k,v in res.items() if k not in ['source_evidence','inputs_sha256','maps']},indent=2))
