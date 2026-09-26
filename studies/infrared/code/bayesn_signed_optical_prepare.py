"""Freeze signed optical inputs and a separately sealed NIR payload for the two-object BayeSN gate."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from ruamel.yaml import YAML

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / 'runs/research_2026_09_26/astra_design/bayesn_signed_pilot'
OLD = ROOT / 'runs/research_2026_09_26/bayesn_distance_identification'
OUT = ROOT / 'runs/research_2026_09_26/bayesn_signed_optical_pilot'
CIDS = ('DES16E2clk', 'DES16X3cry')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records(path):
    return list(csv.DictReader(path.open()))


def row_at(root, source, lineno):
    line = (root / source).read_text().splitlines()[int(lineno)-1]
    a = line.split()
    assert a[0] == 'OBS:'
    return dict(MJD=float(a[1]), band=a[2], flux=float(a[4]), error=float(a[5]), text=line)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    design_sha = sha(DESIGN / 'pilot-design-protocol.json')
    assert design_sha == '544131bf73dcd7de5069a1e25ca10f9cc3fa00e36ebabc0644a610429395d433'
    cohort = {r['CID']: r for r in records(DESIGN / 'frozen-cohort.csv')}
    mask = [r for r in records(DESIGN / 'frozen-mask-metadata.csv') if r['CID'] in CIDS and r['primary_keep'] == 'True']
    optical, nir = [], []
    inputs = {}
    for cid in CIDS:
        c = cohort[cid]
        inputs[c['raw_path']] = sha(ROOT / c['raw_path'])
        inputs[c['NIR_path']] = sha(ROOT / c['NIR_path'])
        for m in [x for x in mask if x['CID'] == cid]:
            obs = row_at(ROOT, m['source_path'], m['source_line'])
            assert obs['band'] == m['band'] and obs['MJD'] == float(m['MJD']) and obs['error'] > 0
            t = (obs['MJD']-float(c['optical_trigger']))/(1+float(c['zHEL']))
            assert abs(t-float(m['trigger_rest_time'])) < 1e-9
            assert 10 < t < 30
            row = dict(cid=cid, source_path=m['source_path'], source_line=int(m['source_line']),
                       MJD=obs['MJD'], band=obs['band'], flux=obs['flux'], error=obs['error'],
                       trigger_rest_time=t)
            (optical if m['role'] == 'optical' else nir).append(row)
    assert [sum(x['cid'] == cid for x in optical) for cid in CIDS] == [15,14]
    assert [sum(x['cid'] == cid for x in nir) for cid in CIDS] == [6,5]
    optical_payload = dict(scope='Signed optical only; forward runner must not read sealed NIR payload.',
                           design_sha256=design_sha, objects=[dict(cid=cid,zHEL=float(cohort[cid]['zHEL']),
                           zHD=float(cohort[cid]['zHD']),MWEBV=float(cohort[cid]['MWEBV']),
                           trigger=float(cohort[cid]['optical_trigger']),mu_LCDM=float(cohort[cid]['mu_LCDM']),
                           sigma_external=float(cohort[cid]['sigma_external'])) for cid in CIDS], rows=optical)
    (OUT/'optical-payload.json').write_text(json.dumps(optical_payload,indent=2)+'\n')
    (OUT/'nir-sealed.json').write_text(json.dumps(dict(scope='Sealed NIR outcomes; no forward/algebra run may read this file.',rows=nir),indent=2)+'\n')

    # Exact released DES J throughput and AB reference; old bridge files are immutable.
    kcor = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_DES_NIR.fits'
    inputs[str(kcor.relative_to(ROOT))] = sha(kcor)
    with fits.open(kcor) as h:
        ft = h['FilterTrans'].data
        zp = h['ZPoff'].data
        names = [str(x).strip() for x in zp['Filter Name']]
        index = names.index('WFC3_IR_F125W-J')
        z = zp[index]
        assert str(z['Primary Name']).strip() == 'AB'
        assert float(z['Primary Mag']) == 0 and float(z['ZPoff(Primary)']) == 0 and float(z['ZPoff(SNpot)']) == 0
        curve = np.column_stack([np.asarray(ft.field(0),float), np.asarray(ft['WFC3_IR_F125W-J'],float)])
    assert np.all(np.isfinite(curve)) and np.all(np.diff(curve[:,0])>0) and np.min(curve[:,1])>=0
    np.savetxt(OUT/'RAISIN_DES_J.dat',curve,fmt='%.12g')
    cfg = YAML(typ='safe').load((OLD/'release-filter-config.yaml').read_text())
    cfg['filters']['RAISIN_DES_J'] = dict(magsys='DES_AB',magzero=0.,path=str((OUT/'RAISIN_DES_J.dat').resolve()))
    YAML().dump(cfg,(OUT/'release-filter-config.yaml').open('w'))
    inputs[str((DESIGN/'pilot-design-protocol.json').relative_to(ROOT))] = design_sha
    inputs[str((DESIGN/'frozen-cohort.csv').relative_to(ROOT))] = sha(DESIGN/'frozen-cohort.csv')
    inputs[str((DESIGN/'frozen-mask-metadata.csv').relative_to(ROOT))] = sha(DESIGN/'frozen-mask-metadata.csv')
    inputs[str((OLD/'release-filter-config.yaml').relative_to(ROOT))] = sha(OLD/'release-filter-config.yaml')
    inputs[str((OLD/'official-code/bayesn/bayesn_model.py').relative_to(ROOT))] = sha(OLD/'official-code/bayesn/bayesn_model.py')
    result = dict(design_sha256=design_sha,inputs=inputs,outputs={p.name:sha(p) for p in
       [OUT/'optical-payload.json',OUT/'nir-sealed.json',OUT/'RAISIN_DES_J.dat',OUT/'release-filter-config.yaml']},
       optical_counts={cid:sum(x['cid']==cid for x in optical) for cid in CIDS},
       sealed_NIR_counts={cid:sum(x['cid']==cid for x in nir) for cid in CIDS},
       reference=dict(name='WFC3_IR_F125W-J',magsys='DES_AB',magzero=0.,positive_A=[float(curve[curve[:,1]>0,0].min()),float(curve[curve[:,1]>0,0].max())]))
    (OUT/'adapter-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(optical_counts=result['optical_counts'],NIR_sealed_counts=result['sealed_NIR_counts'],J_positive_A=result['reference']['positive_A'])))


if __name__ == '__main__':
    main()
