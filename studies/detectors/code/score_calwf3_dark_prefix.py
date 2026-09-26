"""Fixed-mask descriptive repeat slope/ERR comparison; no output DQ selection."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import time
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
H = ROOT/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
OUT = H/'calwf3_dark_prefix'
D = H/'dark_ramp_execution'
NAMES = ['idbx41onq', 'idbx43p7q']


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(name, value):
    with (OUT/name).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def csvsave(name, rows):
    with (OUT/name).open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    start = time.monotonic()
    release = json.loads((OUT/'score-release.json').read_text())
    assert release['approved'] and release['executor_sha256'] == sha(__file__)
    freeze = json.loads((OUT/'score-freeze.json').read_text())
    assert sha(OUT/'score-freeze.json') == release['freeze_sha256']
    for path, digest in freeze['inputs'].items():
        assert sha(path) == digest, path
    mask = ~np.load(D/'fixed-union-bad.npy')[5:1019, 5:1019]
    assert mask.shape == (1014, 1014) and int(mask.sum()) == 993750
    data, source, flags, reads = {}, [], [], []
    for root in NAMES:
        flt = OUT/'work'/root/(root+'_flt.fits')
        with fits.open(flt, memmap=False) as f:
            arrays = {key: f[key,1].data.copy() for key in ('SCI','ERR','DQ','SAMP','TIME')}
            assert all(x.shape == mask.shape for x in arrays.values())
            assert f['SCI',1].header['BUNIT'] == f['ERR',1].header['BUNIT'] == 'COUNTS/S'
            assert f[0].header['EXPFLAG'] == 'NORMAL'
            assert f[0].header['NSAMP'] == 8
            assert f['SCI',1].header['LTV1'] == f['SCI',1].header['LTV2'] == 0
            source.append(dict(root=root, units='DN per nominal second',
                               calibration_flags={k:v for k,v in f[0].header.items() if k.endswith('CORR')},
                               nonfinite={k:int(np.sum(~np.isfinite(arrays[k][mask]))) for k in ('SCI','ERR','TIME')},
                               negative_error=int(np.sum(arrays['ERR'][mask]<0)),
                               zero_error=int(np.sum(arrays['ERR'][mask]==0)),
                               sample_counts={str(k):int(v) for k,v in zip(*np.unique(arrays['SAMP'][mask],return_counts=True))},
                               DQ_counts={str(k):int(v) for k,v in zip(*np.unique(arrays['DQ'][mask],return_counts=True))}))
        data[root] = arrays
        for bit in [1<<i for i in range(16)]:
            flags.append(dict(root=root, stage='FLT', extver=1, bit=bit,
                              fixed_mask_pixels=int(np.sum((arrays['DQ'][mask].astype(np.uint16)&bit)!=0))))
        ima = OUT/'work'/root/(root+'_ima.fits')
        with fits.open(ima, memmap=False) as f:
            assert len(f) == 41
            for extver in range(8,0,-1):
                dq = f['DQ',extver].data
                # Native IMA is untrimmed; same fixed detector mask is translated once.
                assert dq.shape == (1024,1024)
                selected = dq[5:1019,5:1019][mask].astype(np.uint16)
                reads.append(dict(root=root, extver=extver,
                                  SAMPNUM=f['SCI',extver].header['SAMPNUM'],
                                  SAMPTIME=f['SCI',extver].header['SAMPTIME'],
                                  nonzero_DQ_pixels=int(np.sum(selected!=0))))
                for bit in [1<<i for i in range(16)]:
                    flags.append(dict(root=root, stage='IMA', extver=extver, bit=bit,
                                      fixed_mask_pixels=int(np.sum((selected&bit)!=0))))
    save('score-state-ledger.json', source)
    csvsave('native-DQ-bit-counts.csv', flags)
    csvsave('native-read-counts.csv', reads)
    assert all(all(n==0 for n in row['nonfinite'].values()) and row['negative_error']==0 for row in source)
    old, new = [data[k] for k in NAMES]
    difference = new['SCI'].astype(np.float64)-old['SCI'].astype(np.float64)
    variance = old['ERR'].astype(np.float64)**2 + new['ERR'].astype(np.float64)**2
    assert np.all(variance[mask]>0), 'Primary fixed region has nonpositive quoted pair variance'

    def summarize(selected, label):
        x, v = difference[selected], variance[selected]
        assert x.size > 0 and np.all(v>0)
        sq, vv = float(np.dot(x,x)), float(np.sum(v))
        return dict(region=label, pixels=int(x.size), signed_mean=float(np.mean(x)),
                    sum_squared_difference=sq, sum_quoted_pair_variance=vv,
                    repeat_to_quoted_ratio=sq/vv, summed_variance_excess=sq-vv,
                    sum_squared_standardized_difference=float(np.sum(x*x/v)),
                    coherent_mean_fraction=float(x.size*np.mean(x)**2/sq) if sq>0 else 0.)

    regions = [summarize(mask, 'all_active_fixed_mask')]
    quads = [('B',slice(0,507),slice(0,507)), ('C',slice(0,507),slice(507,1014)),
             ('A',slice(507,1014),slice(0,507)), ('D',slice(507,1014),slice(507,1014))]
    for name, yy, xx in quads:
        select = np.zeros(mask.shape,bool)
        select[yy,xx] = mask[yy,xx]
        regions.append(summarize(select,name))
    blocks=[]
    for by in range(16):
        yy=slice(max(5,64*by)-5,min(1019,64*(by+1))-5)
        for bx in range(16):
            xx=slice(max(5,64*bx)-5,min(1019,64*(bx+1))-5)
            select=np.zeros(mask.shape,bool)
            select[yy,xx]=mask[yy,xx]
            blocks.append(summarize(select,f'{by},{bx}'))
    spec=importlib.util.spec_from_file_location('fixed_geometry',D/'fixed_mask.py')
    geometry=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
    radius=int(np.ceil(geometry.R2+.5))
    offsets=range(-radius,radius+1)
    a=np.array([[geometry.pixel_circle_fraction(x,y,geometry.R0) for x in offsets] for y in offsets])
    b=np.array([[geometry.pixel_circle_fraction(x,y,geometry.R2)-geometry.pixel_circle_fraction(x,y,geometry.R1)
                 for x in offsets] for y in offsets])
    sites=[]
    coverage=json.loads((D/'fixed-mask-result.json').read_text())['coverage']
    for site in coverage:
        y,x=site['y']-5,site['x']-5
        yy,xx=slice(y-radius,y+radius+1),slice(x-radius,x+radius+1)
        good=mask[yy,xx]
        ma,mb=good*a,good*b
        eligible=ma.sum()/a.sum()>=.9 and mb.sum()/b.sum()>=.75
        assert eligible==site['eligible']
        if not eligible:
            continue
        w=ma-mb*ma.sum()/mb.sum()
        assert abs(w.sum())<1e-8
        d=float(np.sum(w*difference[yy,xx]))
        v=float(np.sum(w*w*variance[yy,xx]))
        sites.append(dict(raw_x=site['x'],raw_y=site['y'],difference_dn_per_nominal_s=d,
                          quoted_diagonal_pair_variance=v,squared_difference=d*d,
                          aperture_coverage=site['aperture_coverage'],annulus_coverage=site['annulus_coverage']))
    assert len(sites)==247
    csvsave('calibrated-quadrants.csv',regions)
    csvsave('calibrated-blocks.csv',blocks)
    csvsave('calibrated-apertures.csv',sites)
    sum_square=sum(x['squared_difference'] for x in sites)
    sum_var=sum(x['quoted_diagonal_pair_variance'] for x in sites)
    result=dict(arithmetic_pass=True,elapsed_seconds=time.monotonic()-start,
                pair=NAMES,regions=regions,apertures=dict(n=247,sum_squared_difference=sum_square,
                   sum_quoted_diagonal_variance=sum_var,ratio=sum_square/sum_var),
                scope='Single conditional NORMAL dark pair; native adaptive CR/slope operator with original dark switches, DN/nominal seconds',
                unresolved='Native internal segment/power/variance-term ledger not exported; physical clocks, detector covariance, full science recipe and historical SN photometry not validated',
                no_correction='No ERR rescaling, independent-pixel significance, or cosmological propagation')
    save('score-result.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
