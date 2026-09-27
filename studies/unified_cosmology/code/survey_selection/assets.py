#!/usr/bin/env python3
"""Recover native selection/classifier/CC assets and verify actual probability support."""
import argparse,gzip,json,re,tarfile,urllib.request,urllib.error
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from common import ROOT,WORK,RESULTS,sha
from acquire import get
DATA=WORK/'SNDATA_ROOT_2026-04-10'


def maps(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    blocks=[];mjd=None;field=None
    with opener(path,'rt') as f:
        for line in f:
            l=line.split('#')[0].strip()
            if l.startswith('PEAKMJD_RANGE:'):mjd=list(map(float,l.split()[1:]));field=None
            if l.startswith('FIELDLIST:'):field=l.split()[1]
            if l.startswith('VARNAMES:'):blocks.append({'mjd':mjd,'fields':field,'names':l.split()[1:],'rows':[]})
            if l.startswith('HOSTEFF:'):blocks[-1]['rows'].append(list(map(float,l.split()[1:])))
    return blocks


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plasticc',action='store_true');args=parser.parse_args()
    # Explicit further extraction, repeatable from the checksummed official bundle.
    prefixes=['models/classifiers/DES-SN5YR/','models/SALT3/SALT3.DOVEKIE/','models/NON1ASED/','models/lensing/','models/SIMSED/','models/SNIa/','models/population_pdf/','models/VPEC/','standards/']
    members=[];allnames=[]
    with tarfile.open(WORK/'SNDATA_ROOT_2026-04-10.tar.gz') as tar:
        for m in tar:
            if not m.isfile():continue
            allnames.append(m.name)
            if not any(m.name.startswith(p) for p in prefixes):continue
            q=DATA/m.name
            assert '..' not in Path(m.name).parts and not Path(m.name).is_absolute()
            if not q.exists():
                q.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(m) as src,q.open('wb') as dst:
                    while b:=src.read(2**20):dst.write(b)
            members.append({'member':m.name,'bytes':m.size,'sha256':sha(q)})
    eff=DATA/'models/searcheff/SEARCHEFF_zHOST_DES-SN5YR.DAT';old=ROOT/'.work/survey-physics/SNDATA_ROOT/models/searcheff/SEARCHEFF_zHOST_DES-SN5YR.DAT'
    newmaps=maps(eff);oldmaps=maps(old) if old.exists() else None
    changes=[];gridchecks=[];rng=np.random.default_rng(9273102)
    for i,b in enumerate(newmaps):
        a=np.array(b['rows']);assert a.shape[1]==3
        x=np.unique(a[:,0]);y=np.unique(a[:,1]);z=np.empty((len(x),len(y)))
        for row in a:z[np.where(x==row[0])[0][0],np.where(y==row[1])[0][0]]=row[2]
        assert len(a)==len(x)*len(y) and np.all((z>=0)&(z<=1))
        # Native OPT_EXTRAP=1 clamps input coordinates. This is coordinate support,
        # not silent clipping of the resulting probability or changing map nodes.
        q=np.column_stack([rng.uniform(10,35,10000),rng.uniform(-3,6,10000)])
        q=np.clip(q,[x.min(),y.min()],[x.max(),y.max()]);p=RegularGridInterpolator((x,y),z)(q)
        assert p.min()>=-1e-14 and p.max()<=1+1e-14
        gridchecks.append({'fields':b['fields'],'mjd':b['mjd'],'nodes':len(a),'grid_shape':z.shape,'probability_min':float(p.min()),'probability_max':float(p.max()),'tested_points':len(p)})
        if oldmaps:
            oa=np.array(oldmaps[i]['rows']);assert np.array_equal(a[:,:2],oa[:,:2])
            for row in np.flatnonzero(a[:,2]!=oa[:,2]):changes.append({'fields':b['fields'],'mjd':b['mjd'],'r_kron':float(a[row,0]),'g_minus_r_kron':float(a[row,1]),'old':float(oa[row,2]),'new':float(a[row,2])})
    ob=DATA/'models/searcheff/SEARCHEFF_zHOST_DES-SN5YR_OBS.DAT.gz';oldob=old.parent/ob.name
    obsvalues=[r[-1] for b in maps(ob) for r in b['rows']]
    configs=DATA/'sample_input_files/DES-SN5YR';refs=[]
    for f in sorted(configs.rglob('*')):
        if f.suffix not in ['.yml','.nml','.input','.README']:continue
        for n,line in enumerate(f.read_text(errors='replace').splitlines(),1):
            line=line.split('#')[0]
            for token in re.findall(r'\$SNDATA_ROOT/[A-Za-z0-9_./+\-]+',line):
                rel=token.removeprefix('$SNDATA_ROOT/').rstrip('.');path=DATA/rel
                exists=path.exists() or Path(str(path)+'.gz').exists() or rel in allnames or rel+'.gz' in allnames or any(v.startswith(rel.rstrip('/')+'/') for v in allnames)
                refs.append({'config':str(f.relative_to(DATA)),'line':n,'path':rel,'present_in_archive':exists})
    # Fresh one-object batch checks current LFS service; pinned unchanged Git tree
    # means the previous250 pointer objects are still the same omitted realizations.
    prior=ROOT/'studies/host_ages/results/survey_physics/acquisition.json';p=json.loads(prior.read_text());o=p['objects'][0]
    request={'operation':'download','transfers':['basic'],'objects':[{'oid':o['oid_sha256'],'size':o['bytes']} ]}
    req=urllib.request.Request('https://github.com/des-science/DES-SN5YR.git/info/lfs/objects/batch',data=json.dumps(request).encode(),headers={'Accept':'application/vnd.git-lfs+json','Content-Type':'application/vnd.git-lfs+json'})
    try:
        with urllib.request.urlopen(req,timeout=30) as response:status=response.status;body=response.read().decode()
    except urllib.error.HTTPError as e:status=e.code;body=e.read().decode()
    plasticc=[]
    if args.plasticc:
        zp=WORK/'plasticc-zenodo.json';get('https://zenodo.org/api/records/6672739',zp);z=json.loads(zp.read_text())
        for f in z['files']:
            if f['key'] not in ['SIMSED.SNIax.tar.gz','SIMSED.SNIa-91bg.tar.gz']:continue
            rec=get(f['links']['self'],WORK/f['key']);alg,digest=f['checksum'].split(':');assert sha(rec['path'],alg)==digest
            rec['release_checksum']=f['checksum'];info=[]
            with tarfile.open(rec['path']) as t:
                for member in t:
                    if member.isfile():
                        info.append({'member':member.name,'bytes':member.size})
                        if member.name.endswith(('SED.INFO','README','README.md')):
                            q=WORK/'plasticc-inspection'/member.name;assert '..' not in Path(member.name).parts and not Path(member.name).is_absolute();q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(t.extractfile(member).read())
            rec['members']=len(info);rec['unpacked_bytes']=sum(m['bytes'] for m in info);plasticc.append(rec)
    source=ROOT/'.work/survey-physics/SNANA/src/sntools_gridmap.c'
    result={'status':'passed_probability_support_but_exact_production_open','code_sha256':sha(__file__),'official_host_efficiency':{'source_sha256':sha(eff),'old_source_sha256':sha(old) if old.exists() else None,'changed_cells':changes,'grid_checks':gridchecks,'native_boundary_semantics':'OPT_EXTRAP>0 replaces out-of-rangecoordinatesbyboundary±range*1e−12; multilinearinterpolationthereforepreservesboundednodalprobabilities. Not linear continuation beyondthegrid.','native_source_sha256':sha(source),'native_source_lines':[575,578],'OBS_variant_sha256':sha(ob),'OBS_variant_unchanged_from_previous':sha(ob)==sha(oldob) if oldob.exists() else None,'OBS_variant_probability_range':[min(obsvalues),max(obsvalues)]},'new_assets':members,'required_SNDATA_ROOT_references':refs,'current_LFS_batch':{'http_status':status,'body':body,'object':o,'prior250pointer_audit_sha256':sha(prior)},'plasticc_assets':plasticc,'missing_exact_production_state':['250 originalDovekie LFS mock files remain inaccessible ifbatch403; newseedmock generation is not identical.','Pippin EXTERNAL bias/CCprior outputs and originalclassifiercompletepreprocessingstates not suppliedasonefrozenexecutableproduction artifact.','PLAsTiCC publicSED libraries acquired separately are versioned substitutes until equalitywith model_libs_updates usedbyauthors is established; binariescanberebuiltbutnotcalledoriginal.'],'probability_model_scope':'Releasedempiricalhost-redshift efficiency model, not independentlyestimatedcompleteage/dustdependentselection.'}
    (RESULTS/'native-assets.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'nodes_changed':len(changes),'all_nodes':sum(x['nodes'] for x in gridchecks),'refs':len(refs),'missing_refs':len([x for x in refs if not x['present_in_archive']]),'LFS_status':status,'plasticc_files':len(plasticc)}))
if __name__=='__main__':main()
