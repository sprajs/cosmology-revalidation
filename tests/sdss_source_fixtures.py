"""Original tiny byte fixtures; no downloaded/source scientific fixture copies."""
import gzip,struct
def card(key,value):
    if isinstance(value,str):value="'"+value.replace("'","''")+"'"
    return (f'{key:<8}= {value}'.ljust(80)).encode('ascii')
def head(cards):
    body=b''.join(cards)+b'END'.ljust(80)
    return body+b' '*((-len(body))%2880)
def fits_one_row(extra=None):
    primary={'SIMPLE':'T','BITPIX':'8','NAXIS':'0'}
    extension={'XTENSION':'BINTABLE','BITPIX':'8','NAXIS':'2','NAXIS1':'20','NAXIS2':'1','PCOUNT':'0','GCOUNT':'1','TFIELDS':'2','TTYPE1':'SNID','TFORM1':'16A','TUNIT1':'','TTYPE2':'FLUXCAL','TFORM2':'1E','TUNIT2':''}
    if extra:extension.update(extra)
    p=head([card(k,v) for k,v in primary.items()])
    e=head([card(k,v) for k,v in extension.items()])
    raw=b'6057'+b'\x00'*12+struct.pack('>f',-3.25)
    return gzip.compress(p+e+raw+b'\0'*((-len(raw))%2880),mtime=0),[primary,extension],raw
