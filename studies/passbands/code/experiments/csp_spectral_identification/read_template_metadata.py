from pathlib import Path
import struct,zlib,json,hashlib
p=Path(__file__).parent;b=(p/'template-readme-local-range').read_bytes();a=struct.unpack('<4s5H3I2H',b[:30]);r=json.loads((p/'NIR_Ia_template.zip.members.json').read_text())[0]
nl,xl=a[-2:];data=zlib.decompress(b[30+nl+xl:30+nl+xl+r['compressed']],-15)
assert len(data)==r['size'] and zlib.crc32(data)==r['CRC']
(p/'readme_template.txt').write_bytes(data)
(p/'template-readme-acquisition.json').write_text(json.dumps({'url':'https://csp.obs.carnegiescience.edu/data/NIR_Ia_template.zip','range':'0-767','archive_downloaded':False,'member':r,'sha256':hashlib.sha256(data).hexdigest(),'original_attempt_failure':'Local header uses data-descriptor flag 8 and compressed size0; initial decompression used this zero size and failed before any output. Corrected extraction uses verified central-directory compressed size671 and CRC.'},indent=2)+'\n')
print(data.decode())
print((p/'spectra/readme_data.txt').read_text())
