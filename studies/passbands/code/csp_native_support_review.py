"""Independent, rowwise native-support audit using source-grid and KCOR bytes."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from astropy.io import fits
from review_csp_native_states import parse

ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log', type=Path)
    ap.add_argument('out', type=Path)
    a = ap.parse_args()
    a.out.mkdir(exist_ok=False, parents=True)
    model = RELEASE / 'model/snoopy.B18/SNooPy_B18.fits'
    kcor = RELEASE / 'kcor/kcor_CSPDR3_BD17.fits'
    with fits.open(model) as h:
        phase = np.asarray(h['TREST-GRID'].data['TREST'], float)
        shape = np.asarray(h['LUMI-GRID'].data['STRETCH'], float)
    records, _ = parse(a.log)
    bands = {r['band'] for s in records for r in s['rows']}
    with fits.open(kcor) as h:
        t = h['FilterTrans'].data
        wave = np.asarray(t.field(0), float)
        bounds = {b: (float(wave[t['CSP-'+b] != 0].min()),
                      float(wave[t['CSP-'+b] != 0].max())) for b in bands}
    results = []
    for i, s in enumerate(records):
        x = np.asarray([r['values'] for r in s['rows']])
        rest_bounds = np.asarray([bounds[r['band']] for r in s['rows']]) / (1+x[:, 6, None])
        ok_wave = (rest_bounds[:, 0] >= x[:, 11]) & (rest_bounds[:, 1] <= x[:, 12])
        ok_phase = (x[:, 1] >= phase.min()) & (x[:, 1] <= phase.max())
        ok_shape = shape.min() <= s['objective'][4] <= shape.max()
        results.append(dict(callback=i, CID=s['CID'], iteration=s['iteration'],
                            n=s['n'], phase_min=float(x[:, 1].min()),
                            phase_max=float(x[:, 1].max()), shape=s['objective'][4],
                            shape_supported=bool(ok_shape), phase_supported=bool(ok_phase.all()),
                            full_nonzero_throughput_supported=bool(ok_wave.all()),
                            positive_mean=bool((x[:, 2] > 0).all()),
                            unsupported_rows=[dict(index=s['rows'][j]['index'],
                                                   band=s['rows'][j]['band'],
                                                   phase=float(x[j, 1]),
                                                   rest_filter_bounds=rest_bounds[j].tolist(),
                                                   native_limits=x[j, 11:13].tolist())
                                              for j in range(s['n']) if not (ok_wave[j] and ok_phase[j])]))
    answer = dict(pass_all_support=all(r['shape_supported'] and r['phase_supported'] and
                                      r['full_nonzero_throughput_supported'] and r['positive_mean']
                                      for r in results),
                  scope='Finite interpolation and native wavelength-screen support only; does not validate empirical model or spectral training coverage.',
                  phase_bounds=[float(phase.min()), float(phase.max())],
                  shape_bounds=[float(shape.min()), float(shape.max())],
                  observer_filter_nonzero_bounds=bounds, callbacks=results)
    (a.out/'result.json').write_text(json.dumps(answer, indent=2)+'\n')
    (a.out/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    (a.out/'manifest.json').write_text(json.dumps(dict(inputs={str(p.resolve()): sha(p) for p in
        [a.log, model, kcor, Path(__file__).with_name('review_csp_native_states.py')]},
        outputs={p.name: sha(p) for p in a.out.iterdir() if p.is_file() and p.name!='manifest.json'}), indent=2)+'\n')
    print(json.dumps(dict(pass_all_support=answer['pass_all_support'], callbacks=len(results),
                         unsupported=[r for r in results if not r['phase_supported'] or not r['full_nonzero_throughput_supported'] or not r['shape_supported']]), indent=2))

if __name__ == '__main__': main()
