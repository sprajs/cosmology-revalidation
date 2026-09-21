#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
import jsonschema
ROOT=Path(__file__).resolve().parents[3];D=ROOT/'phase2/literature'
data=json.loads((D/'literature-map.json').read_text());schema=json.loads((D/'literature-map.schema.json').read_text())
jsonschema.Draft202012Validator.check_schema(schema);jsonschema.validate(data,schema)
assert len({e['id'] for e in data['entries']})==len(data['entries'])
testids={t['id'] for t in data['preregistered_tests']}; checked=0
for e in data['entries']:
 assert set(e['preregistered_tests'])<=testids
 assert e['methods_read'] or 'summary' in e['reading_status']
 for s in e['sources']:
  path=ROOT/s['path'];assert path.is_file(),path
  assert hashlib.sha256(path.read_bytes()).hexdigest()==s['sha256'],path;checked+=1
report={'schema':'Draft 2020-12','validator':'jsonschema 4.25.1','entries':len(data['entries']),'source_hash_checks':checked,'test_references_resolved':True,'reading_status_counts':data['reading_status_counts'],'result':'PASS','scientific_limit':'Schema/hash validation does not certify a full-method reanalysis of every entry; reading depth is explicit.'}
(D/'literature-map-validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
