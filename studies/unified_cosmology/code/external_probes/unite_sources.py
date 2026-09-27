"""Bounded, source-linked availability audit of the versioned Unite release.

Download papers and only Unite-specific repositories/deposits actually linked
in them. HTTP failures are evidence, never evidence of global unavailability.
No cosmological posterior or figure is converted into substitute observations.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
from html import unescape
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.error
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
SOURCES = {
    'abstract-v2.html': 'https://arxiv.org/abs/2609.05053v2',
    'paper-v2.html': 'https://arxiv.org/html/2609.05053v2',
    'paper-v2.pdf': 'https://arxiv.org/pdf/2609.05053v2',
    'paper-v1.html': 'https://arxiv.org/html/2609.05053v1',
    'paper-v1.pdf': 'https://arxiv.org/pdf/2609.05053v1',
    'source-v1.tar': 'https://arxiv.org/src/2609.05053v1',
    'source-v2.tar': 'https://arxiv.org/src/2609.05053v2',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


class TextLinks(HTMLParser):
    def __init__(self):
        super().__init__(); self.text = []; self.links = []
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'a' and attrs.get('href', '').startswith('http'):
            self.links.append(unescape(attrs['href']))
    def handle_data(self, text):
        self.text.append(text)


def acquire(cache, output):
    cache = cache.resolve()
    assert cache.is_relative_to(ROOT/'.work') and not cache.exists()
    assert not output.exists(), 'Use fresh paths; previous retrieval evidence is immutable.'
    cache.mkdir(parents=True)
    resources = {}; parsed = {}; links = {}; pdf_extraction = {}
    def fetch(item):
        name, url = item
        record = {'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()}
        try:
            request = urllib.request.Request(url, headers={'User-Agent':'cosmology-revalidation/source-audit'})
            with urllib.request.urlopen(request, timeout=55) as response:
                body = response.read(12_000_001)
                assert len(body) <= 12_000_000, 'Bounded acquisition limit exceeded.'
                record.update(status_code=response.status, resolved_url=response.url,
                              content_type=response.headers.get('Content-Type'))
        except urllib.error.HTTPError as error:
            body = error.read(1_000_000)
            record.update(status_code=error.code, resolved_url=error.url,
                          error=repr(error), content_type=error.headers.get('Content-Type'))
        except Exception as error:
            return name, dict(record, status='retrieval_failed', error=repr(error))
        path = cache/name; path.write_bytes(body)
        record.update(status='retrieved' if record['status_code']==200 else 'http_error',
                      path=rel(path), sha256=sha(path), bytes=len(body))
        return name, record
    with ThreadPoolExecutor(max_workers=3) as pool:
        resources.update(pool.map(fetch, SOURCES.items()))
    source_members = {}
    for name, record in list(resources.items()):
        if record['status'] != 'retrieved': continue
        path = ROOT/record['path']; text = ''
        if name.endswith('.html'):
            parser = TextLinks(); parser.feed(path.read_text(errors='replace'))
            text = ' '.join(parser.text); links[name] = sorted(set(parser.links))
        elif name.endswith('.pdf'):
            assert path.read_bytes().startswith(b'%PDF'), 'Successful PDF request returned non-PDF.'
            converted = subprocess.run(['pdftotext','-layout',str(path),'-'], capture_output=True, text=True, check=True)
            text = converted.stdout
            text_path=cache/(name+'.txt'); text_path.write_text(text)
            pdf_extraction[name]={'path':rel(text_path),'sha256':sha(text_path),
                                  'returncode':converted.returncode,'stderr':converted.stderr}
        elif name.endswith('.tar'):
            try:
                with tarfile.open(fileobj=io.BytesIO(path.read_bytes())) as archive:
                    source_members[name] = [{'name':m.name,'bytes':m.size} for m in archive.getmembers()]
                    contents=[]
                    for member in archive.getmembers():
                        if member.isfile() and member.name.endswith(('.tex','.bib','.txt')) and member.size<2_000_000:
                            contents.append(archive.extractfile(member).read().decode(errors='replace'))
                    text='\n'.join(contents)
            except tarfile.TarError as error:
                source_members[name]={'parse_error':repr(error)}
        parsed[name]=' '.join(text.split())
        links.setdefault(name, sorted(set(re.findall(r'https?://[^\s<>"{}\\]+', text))))
    all_links=sorted(set(link.rstrip('.,;)') for values in links.values() for link in values))
    relevant=[link for link in all_links if ('github.com' in link or 'zenodo.org' in link)
              and ('unite' in link.lower() or '2609.05053' in link)]
    # PDF text can wrap a URL at an underscore. Do not turn that fragment into
    # an independent missing repository when HTML supplies the full hyperlink.
    html_links={link.rstrip('.,;)') for name,values in links.items() if name.endswith('.html') for link in values}
    fragments=[link for link in relevant if link not in html_links and any(other.startswith(link) and other!=link for other in html_links)]
    relevant=[link for link in relevant if link not in fragments]
    repositories={}; deposits={}
    for index, link in enumerate(relevant):
        if 'github.com/' in link:
            parts=urllib.parse.urlparse(link).path.strip('/').split('/')
            if len(parts)<2: continue
            repository='/'.join(parts[:2]); key='repo-'+str(index)
            if any(repository.lower()==known.lower() for known in repositories): continue
            name, record=fetch((key+'.json','https://api.github.com/repos/'+repository))
            resources[name]=record
            entry={'linked_url':link,'repository':repository,'metadata_resource':name}
            if record['status']=='retrieved':
                metadata=json.loads((ROOT/record['path']).read_text()); branch=metadata['default_branch']
                for suffix, endpoint in [('tree',f'/git/trees/{urllib.parse.quote(branch,safe="")}?recursive=1'),('releases','/releases')]:
                    n,r=fetch((key+'-'+suffix+'.json','https://api.github.com/repos/'+repository+endpoint)); resources[n]=r
                    entry[suffix+'_resource']=n
                    if suffix=='tree' and r['status']=='retrieved':
                        tree=json.loads((ROOT/r['path']).read_text());entry['tree_commit']=tree['sha'];entry['tree_truncated']=tree['truncated']
                        entry['candidate_data_paths']=[v['path'] for v in tree['tree'] if any(s in v['path'].lower() for s in ['cov','hubble','likelihood','.dat','.fits'])]
            repositories[repository]=entry
        elif 'zenodo.org' in link:
            match=re.search(r'/records?/(\d+)',link)
            if match:
                name,record=fetch((f'zenodo-{match[1]}.json','https://zenodo.org/api/records/'+match[1]));resources[name]=record;deposits[match[1]]=name
    v2=parsed.get('paper-v2.html','');abstract=parsed.get('abstract-v2.html','')
    phrases=['publicly available upon acceptance','Removed inactive datalinks','2884']
    checks={phrase:phrase.lower() in (v2+' '+abstract).lower() for phrase in phrases}
    assert all(checks.values()), 'Versioned availability statement/sample/revision no longer matches; inspect before conclusion.'
    record={'status':'bounded_author_source_availability_audit_completed',
            'retrieved_utc':datetime.now(timezone.utc).isoformat(), 'resources':resources,
            'versioned_statement_checks':checks,'source_archive_members':source_members,
            'pdf_text_extraction':pdf_extraction,'excluded_PDF_wrapped_URL_fragments':fragments,
            'author_links_by_resource':links,'unite_specific_linked_interfaces':relevant,
            'repository_checks':repositories,'deposit_checks':deposits,
            'reported_likely_SN_count':2884,
            'usable_joint_distance_likelihood_acquired':False,
            'conclusion':'The versioned v2 article promises public data, likelihood and host measurements on acceptance. The checked author-linked interfaces did not supply an audited joint ordered distance vector, covariance and likelihood. HTTP failures or inactive links are scoped retrieval outcomes, not proof of global absence.',
            'integration_limit':'Unite overlaps Pantheon+ and DES-SN5YR/Dovekie and must be treated as a replacement SN analysis, not an additional independent factor. Published parameter posteriors and digitized figures are not substitutes for its distance likelihood.',
            'missing_for_target':['ordered SN identifiers and distances/redshifts','matching statistical+systematic covariance or exact joint likelihood','sample/normalization/anchoring semantics','release provenance and overlap metadata'],
            'model_calls':0,'source_sha256':{rel(__file__):sha(__file__)},
            'runtime':{'pdftotext':subprocess.run(['pdftotext','-v'],capture_output=True,text=True).stderr.splitlines()[0]}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); result=acquire(args.cache,args.output)
    print(json.dumps({'status':result['status'],'resources':len(result['resources']),
                      'linked_interfaces':result['unite_specific_linked_interfaces']}))
