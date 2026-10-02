#!/usr/bin/env python3
"""Structural source transport and retained-SDK invocation; production science is native."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
FOLDER=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('retained_sdk_identity',ROOT/'experiments/released-ladder/controller.py')
sdk_checks=importlib.util.module_from_spec(spec);spec.loader.exec_module(sdk_checks)
load=sdk_checks.load;sha256=sdk_checks.sha256
CONFIG_SHA256='a95a36d99c304f99d23e800cf22fa5bc943fe05cd79b4c42bdbf8fabef9b09e3'

def sources(directory, config):
    result={}
    for name,expected in config['sources'].items():
        path=directory/name
        if path.stat().st_size!=expected['bytes'] or sha256(path)!=expected['sha256']:
            raise ValueError('source differs: '+name)
        result[name]=path
    return result

def transport(admitted, store):
    table=list(csv.DictReader(admitted['sn-table.dat'].read_text().splitlines(),delimiter=' ',skipinitialspace=True))
    if len(table)!=1701:raise ValueError('source row count')
    selected=[(i,r) for i,r in enumerate(table) if r['IDSURVEY']=='1']
    if len(selected)!=321 or len({r['CID'] for _,r in selected})!=321:raise ValueError('selected row identities')
    if selected[0][0]!=690 or selected[0][1]['CID']!='6057' or selected[0][1]['zHEL']!='0.06707':raise ValueError('selected photometry redshift identity')
    for _,r in selected:
        if not all(math.isfinite(float(r[k])) for k in ['zHD','zHEL']) or float(r['zHD'])<=.01:
            raise ValueError('selected redshift domain')
    (store/'rows.tsv').write_text(''.join(f"{i} {r['zHD']} {r['zHEL']}\n" for i,r in selected))
    (store/'row-identities.json').write_text(json.dumps([{'source_row_zero_based':i,'CID':r['CID'],'IDSURVEY':r['IDSURVEY']} for i,r in selected],indent=2)+'\n')
    raw=admitted['wfc3-filter-original.fits'].read_bytes()
    def header(offset):
        parsed={}
        for at in range(offset,offset+2880,80):
            card=raw[at:at+80].decode('ascii');key=card[:8].strip()
            if key=='END':return parsed
            if card[8:10]=='= ':parsed[key]=card[10:].split('/')[0].strip().strip("'").strip()
        raise ValueError('FITS header end missing')
    primary,extension=header(0),header(2880)
    if primary.get('SIMPLE')!='T' or primary.get('NAXIS')!='0':raise ValueError('primary FITS shape')
    expected={'XTENSION':'BINTABLE','BITPIX':'8','NAXIS':'2','NAXIS1':'8','NAXIS2':'910','TFIELDS':'2','TTYPE1':'WAVELENGTH','TTYPE2':'THROUGHPUT','TFORM1':'E','TFORM2':'E','TUNIT1':'ANGSTROMS','TUNIT2':'TRANSMISSION','PCOUNT':'0'}
    if any(extension.get(k)!=v for k,v in expected.items()):raise ValueError('typed FITS table differs')
    pairs=list(struct.iter_unpack('>ff',raw[5760:5760+910*8]))
    if any(not math.isfinite(x) or not math.isfinite(t) or x<=0 or not 0<=t<=1 for x,t in pairs):raise ValueError('response domain')
    if any(a[0]>=b[0] for a,b in zip(pairs,pairs[1:])) or pairs[0][1]!=0 or pairs[-1][1]!=0:raise ValueError('response support')
    (store/'passband-angstrom.tsv').write_text(''.join(f'{x:.17g} {t:.17g}\n' for x,t in pairs))

def execute(args):
    store=args.output.resolve();store.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();record={'status':'started','variant':'conditional-sdss-observer-historical-optical/v1','engine_identity':sdk_checks.SDK_IDENTITY,'gates':{'execution':'unassessed','numerical':'unassessed','observation':'blocked','inference':'not attempted','measured_calibration_law':'blocked'}}
    before=None;admitted=None
    try:
        for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
            os.environ[key]='1'
        record['reference_threads']={k:os.environ[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']}
        packet=load(FOLDER/'experiment.json');config=load(FOLDER/'config.json')
        if sha256(FOLDER/'config.json')!=CONFIG_SHA256:raise ValueError('reviewed configuration changed')
        (store/'config.json').write_bytes((FOLDER/'config.json').read_bytes())
        expected_candidates={'candidate-sn.json':packet['origin']['sha256'],'candidate-response.json':packet['limitations'][0].split('sha256=')[1]}
        for name,digest in expected_candidates.items():
            if sha256(FOLDER/name)!=digest:raise ValueError('candidate changed')
        record['source_candidate_revision']=packet['origin']['revision'];record['candidate_hashes']=expected_candidates
        record['packet_source_hashes']={p.name:sha256(p) for p in sorted(FOLDER.iterdir()) if p.is_file()}
        record['retained_sdk_checker_sha256']=sha256(Path(sdk_checks.__file__))
        record['consumer_revision']=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
        record['consumer_dirty']=bool(subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip())
        if record['consumer_dirty']:raise ValueError('consumer requires clean committed source')
        before=sdk_checks.fingerprint(args.engine_source,args.sdk);record['sdk_before']=before
        admitted=sources(args.sources,config);record['input_before']={n:sha256(p) for n,p in admitted.items()}
        transport(admitted,store)
        command=['/usr/bin/c++','-std=c++20','-O2','-fno-fast-math','-ffp-contract=off','-I',str(args.sdk/'include'),str(FOLDER/'consumer.cpp'),str(args.sdk/'lib/libirred_core.a'),'-o',str(store/'consumer')]
        record['compiler_command']=command
        sdk_checks.child(command,store,'compile',120)
        record['consumer_executable_sha256']=sha256(store/'consumer')
        native=sdk_checks.child([str(store/'consumer'),str(store/'rows.tsv'),str(store/'passband-angstrom.tsv'),str(store/'passband-metre.tsv')],store,'native',30)
        result=json.loads(native)
        if result['variant']!=config['variant'] or len(result['background'])!=645 or len(result['photometry'])!=3 or not all(result['invalid_controls'].values()):raise ValueError('native output admission')
        record['gates']['execution']='accepted'
        reference=json.loads(sdk_checks.child([str(args.reference_python),str(FOLDER/'reference.py'),str(store)],store,'reference',120))
        if reference['status']!='accepted':raise ValueError('reference comparisons failed')
        record['reference']=reference;record['gates']['numerical']='accepted';record['status']='accepted_conditional_deterministic_control'
    except Exception as error:
        record['status']='failed';record['error']=str(error)
    finally:
        errors=[]
        try:
            after=sdk_checks.fingerprint(args.engine_source,args.sdk);record['sdk_after']=after
            if before is not None and after!=before:errors.append('SDK changed')
        except Exception as error:errors.append('SDK final verification: '+str(error))
        if admitted:
            record['input_after']={n:sha256(p) for n,p in admitted.items()}
            if record['input_after']!=record['input_before']:errors.append('source changed')
        if 'packet_source_hashes' in record:
            after={p.name:sha256(p) for p in sorted(FOLDER.iterdir()) if p.is_file()}
            if after!=record['packet_source_hashes']:errors.append('packet changed')
            if sha256(Path(sdk_checks.__file__))!=record['retained_sdk_checker_sha256']:errors.append('SDK checker changed')
        if errors:record['status']='failed';record['integrity_errors']=errors
        record['elapsed_seconds']=time.monotonic()-started
        record['generated_hashes']={p.name:sha256(p) for p in store.iterdir() if p.is_file()}
        (store/'record.json').write_text(json.dumps(record,indent=2)+'\n')
        for path in store.iterdir():
            if path.is_file():path.chmod(0o444)
    print(json.dumps({'status':record['status'],'gates':record['gates'],'record_sha256':sha256(store/'record.json')}))
    return 0 if record['status']=='accepted_conditional_deterministic_control' else 1

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--sdk',type=Path,required=True);parser.add_argument('--engine-source',type=Path,required=True);parser.add_argument('--sources',type=Path,required=True);parser.add_argument('--reference-python',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    raise SystemExit(execute(parser.parse_args()))
