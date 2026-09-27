"""Pin primary calibration and signed-data availability evidence."""
import json,re,urllib.request
from datetime import datetime,timezone
from pathlib import Path
from acquire import ROOT,WORK,OUT,sha
from recover import fetch

def main():
    sources=[]
    urls={
      'ozdes/dr2.html':'https://docs.datacentral.org.au/ozdes/overview/dr2/',
      'ozdes/faq.html':'https://docs.datacentral.org.au/ozdes/overview/faq/',
      'yse-api-docs.html':'https://yse-pz.readthedocs.io/en/latest/api.html',
      'yse-pz-public/api.rst':'https://raw.githubusercontent.com/davecoulter/YSE_PZ/58f3e6a1622ec5755e5322aee2d00f3941510749/docs/api.rst',
      'yse-publications.html':'https://yse.ucsc.edu/publications/',
      'yse-prior-record.json':'https://zenodo.org/api/records/7186195',
      'yse-pz-public/README.md':'https://raw.githubusercontent.com/davecoulter/YSE_PZ/58f3e6a1622ec5755e5322aee2d00f3941510749/README.md',
      'yse-pz-public/api_views.py':'https://raw.githubusercontent.com/davecoulter/YSE_PZ/58f3e6a1622ec5755e5322aee2d00f3941510749/YSE_App/api_views.py'}
    for name,url in urls.items():
        path=WORK/name;path.parent.mkdir(exist_ok=True,parents=True)
        try:sources.append(fetch(url,path))
        except Exception as exc:sources.append(dict(url=url,status='fetch_failed',error=type(exc).__name__+': '+str(exc)))
    old=json.loads((WORK/'yse-prior-record.json').read_text())
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),sources=sources,
       yse_prior_release_files=[q['key']for q in old.get('files',[])],
       signed_YSE_recovery='The official survey publications page links DR1 Zenodo7317476. Its README explicitly discards negative flux. Prior record7186195 returns no downloadable files. Public YSE-PZ API documentation requires basic-auth credentials and pinned API model views use IsAuthenticated. The open software is not an unrestricted observational database. No credentials or author contact used; no signed upstream forced-flux table recovered through these sources. This bounded search does not establish that no other signed observations exist.',
       ozdes_calibration='Primary FAQ says spectra are not flux calibrated, including relative calibration; native units normalized counts per wavelength. Therefore count-band ratios are not physical Dn4000 without response information. DR2 warns dichroic continuum discontinuities from arm scale-factor failures.',
       magnitude_semantics='YSE README calls MAG AB magnitude and MAGERR magnitude error. Paper asinh/luptitude discussion concerns galaxy photo-z input, not an explicit SN epoch MAG conversion. Low-SNR metadata and possible sentinel/upper-limit conventions are not reverse-engineered into flux selection. Numerical MAG-vs-logFLUX discrepancies alone are not established photometric errors.')
    (OUT/'primary-source-audit.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
