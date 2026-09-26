#!/usr/bin/env python3
"""Audit full generated attempts without assuming CID is a unique attempt ID."""
from pathlib import Path
import argparse,json,re
import pandas as pd
from astropy.io import fits
R=Path(__file__).resolve().parents[3];O=R/'phase2/literature/simulations'
p=argparse.ArgumentParser();p.add_argument('version',nargs='?',default='PH2_pilot02_P21');a=p.parse_args();folder=O/'outputs'/a.version
lines=(folder/(a.version+'.DUMP')).read_text().splitlines();cols=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'));rows=[x.split()[1:] for x in lines if x.startswith('SN:')]
d=pd.DataFrame(rows,columns=cols).apply(pd.to_numeric);d.insert(0,'generated_attempt_index',range(1,len(d)+1));readme=(folder/(a.version+'.README')).read_text();head=fits.getdata(folder/(a.version+'_HEAD.FITS'),1)
generated=int(re.search(r'NGENLC_TOT:\s+(\d+)',readme).group(1));written=int(re.search(r'NGENLC_WRITE:\s+(\d+)',readme).group(1));passrows=d[(d.SIM_EFFMASK==5)&(d.CUTMASK==4095)]
headids={int(x) for x in head['SNID']};passedids=set(map(int,passrows.CID));assert generated==len(d);assert written==len(head)==len(passrows);assert len(passedids)==len(passrows);assert headids==passedids
record={'version':a.version,'generated':generated,'written':written,'all_generated_dump_rows':len(d),'duplicate_cid_count':int(d.CID.duplicated().sum()),'generated_key':'one-based generated_attempt_index, notCID','head_pass_unique_cid_match':True,'sim_effmask_counts':{str(k):int(v) for k,v in d.SIM_EFFMASK.value_counts().items()},'cutmask_counts':{str(k):int(v) for k,v in d.CUTMASK.value_counts().items()},'acceptance':written/generated,'limits':['Early rejected attempts reuse CID; join HEAD only toDUMPpassrows, neverallattempts byCID.','Columns notcomputed before anearly rejection have sentinel values; these are notphysicalobservations.','Mask5+cut4095 reproduces writtenmembership in thisconfiguration; measuredfluxfit selection not yet included.']}
(O/(a.version+'-denominator.json')).write_text(json.dumps(record,indent=2)+'\n');d.to_csv(O/(a.version+'-generated.csv.gz'),index=False,compression={'method':'gzip','mtime':0});print(json.dumps(record,indent=2))
