#!/usr/bin/env python3
"""Freeze source/config evidence and compare published DES/Dovekie scaffold."""
import hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/salt_dust_audit/snana_config'
OLD=ROOT/'sources/repos/des-science__DES-SN5YR@1.3'
NEW=ROOT/'sources/repos/des-science__DES-SN5YR'
SRC=ROOT/'sources/repos/RickKessler__SNANA/src'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for p in sorted((OLD/'7_PIPPIN_FILES').rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(OLD);q=NEW/rel
        rows.append({'relative_path':str(rel),'original_sha256':sha(p),'current_sha256':sha(q)if q.exists()else None,'identical':q.exists() and p.read_bytes()==q.read_bytes()})
    old_mask=int(re.search(r'^opt_biascor=(\d+)',(OLD/'7_PIPPIN_FILES/base_files/bbc/BBC_des5yr.input').read_text(),re.M)[1])
    defs=dict((name,int(val)) for name,val in re.findall(r'^#define (MASK_BIASCOR_\w+)\s+(\d+)\s', (SRC/'SALT2mu.c').read_text(),re.M))
    bits={k:v for k,v in defs.items()if v>0 and v&(v-1)==0 and old_mask&v}
    assert sum(bits.values())==old_mask
    fitopts=(NEW/'7_PIPPIN_FILES/base_files/lcfit/fitopts.yml').read_text()
    weights=[float(x)for x in re.findall(r'^\s*cal_\d+:\s*(\S+)',fitopts,re.M)]
    evidence=[]
    needles={
      SRC/'MWgaldust.h':['OPT_MWCOLORLAW_FITZ99_APPROX','OPT_MWCOLORLAW_FITZ99_EXACT','RVMIN_FITZ99'],
      SRC/'MWgaldust.c':['MWEBV_OUT = 0.86','MWEBV_ERR_OUT = 0.05','Only reliable for RV=3.1'],
      SRC/'genmag_SEDtools.c':['OPT_COLORLAW     = MWXT_SEDMODEL.OPT_COLORLAW','LAMREST    = LAMOBS/(1.0 + z)','XT_FRAC    = pow(TEN,arg)'],
      SRC/'sntools_genPDF.c':['get_VAL_RANGE_genPDF(IMAP','prob /= FUNMAX'],
      SRC/'SALT2mu.c':['muBiasErr = 0.0','muerrsq += muCOVadd','mu      -= muBias'],
      NEW/'7_PIPPIN_FILES/base_files/lcfit/lcfit_desSMP_5yr.nml':['OPT_MWCOLORLAW','OPT_MWEBV','RV_MWCOLORLAW','FITMODEL_NAME','DELCHI2_REJECT','RESTLAMBDA_FITRANGE'],
      NEW/'7_PIPPIN_FILES/base_files/bbc/BBC_des5yr.input':['opt_biascor','interp_biascor_logmass','CUTWIN x1ERR','p7=','p8='],
    }
    for p,finds in needles.items():
        lines=p.read_text().splitlines()
        evidence.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p),'matches':[{'line':j,'text':line}for j,line in enumerate(lines,1)if any(x in line for x in finds)]})
    result={'comparison':'Public tag1.3 PIPPIN files versus acquired current DES-Dovekie repository tree. These files are release scaffolding, not a claim they are the exact executed Dovekie jobs.', 'N_PIPPIN_files':len(rows),'N_byte_identical':sum(r['identical']for r in rows),'file_comparison':rows,'bbc_mask':{'value':old_mask,'enabled_bits':bits}, 'calibration_weights_in_current_scaffold':{'N':len(weights),'amplitude_weights':weights,'sum_squared':sum(x*x for x in weights),'equal_amplitude_for_unit_total':len(weights)**-.5},'findings':{'scaffold_matches_original':all(r['identical']for r in rows),'Dovekie_execution_provenance_missing':'Current scaffold still names SALT3.DES5YR, uses x1ERR<1.0 and 9 calibration amplitude weights .3. Paper reports DOVEKIE retraining, x1ERR<1.15 and weights 1/sqrt(9). Do not run this scaffold and call it exact Dovekie reproduction.'}, 'source_evidence':evidence}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    files=[Path(__file__).resolve(),ROOT/'catalog/repositories/RickKessler__SNANA.commit.json',ROOT/'catalog/repositories/des-science__DES-SN5YR.commit.json',ROOT/'catalog/repositories/des-science__DES-SN5YR@1.3.commit.json',SRC/'MWgaldust.c',SRC/'MWgaldust.h',SRC/'snlc_sim.c',SRC/'genmag_SEDtools.c',SRC/'genmag_SALT2.c',SRC/'sntools_genPDF.c',SRC/'SALT2mu.c',ROOT/'papers/text/2511.07517v3.txt',ROOT/'papers/text/2401.02945v2.txt',ROOT/'papers/text/2112.04456v2.txt',OUT/'results.json']
    (OUT/'manifest.json').write_text(json.dumps({'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in files]},indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k not in ['source_evidence','file_comparison']},indent=2))

if __name__=='__main__':main()
