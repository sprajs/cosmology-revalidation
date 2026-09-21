from pathlib import Path
import json,hashlib
from astropy.io.votable import parse_single_table
O=Path(__file__).resolve().parents[3]/'phase2/assumptions/image_access';out=[]
for p in sorted(O.glob('sia_des*.xml')):
 t=parse_single_table(str(p)).to_table();rows=[{k:str(r[k]) for k in t.colnames} for r in t];q=O/(p.stem+'.json');q.write_text(json.dumps(rows,indent=2)+'\n');out.append({'source':p.name,'source_sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'parsed_rows':len(t),'output':q.name,'output_sha256':hashlib.file_digest(q.open('rb'),'sha256').hexdigest()})
(O/'parsed_response_summary.json').write_text(json.dumps({'note':'Original access manifest parse_error was raised while serializing vector-valued masked table columns to CSV, not by the HTTP query or VOTable contents. This parser preserves every field as a JSON string.','responses':out},indent=2)+'\n')
