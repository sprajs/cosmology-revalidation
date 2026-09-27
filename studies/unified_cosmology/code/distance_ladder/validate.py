"""Verify current ladder identities, independent calculations and overlap geometry."""
import ast
import json
import re
import numpy as np
import pandas as pd
from scipy import linalg
from common import ROOT, HERE, WORK, OUT, sha, relative, write


def main():
    checked = {}
    def verify(path, digest):
        resolved = ROOT/path
        assert resolved.is_file(), path
        assert sha(resolved) == digest, path
        checked[path] = digest

    records = {}
    names = ['acquisition', 'reconstruction', 'cepheid-factor', 'overlap', 'independent-review', 'independent-factor-review', 'equation-visual-review']
    for name in names:
        records[name] = json.loads((OUT/(name+'.json')).read_text())
        record = records[name]
        assert record['status'].startswith('passed_') or record['status'] == 'visually_confirmed_printed_sign_inconsistency'
        for key, value in record.items():
            if key.endswith('sha256') and isinstance(value, dict):
                for path, digest in value.items():
                    verify(path, digest)
        if isinstance(record.get('source_sha256'), str):
            verify(relative(HERE/(name.replace('-', '_')+'.py')), record['source_sha256'])
    for entry in records['acquisition']['files']:
        verify(entry['path'], entry['sha256'])
        assert (ROOT/entry['path']).stat().st_size == entry['bytes']
    visual = records['equation-visual-review']
    verify(visual['source_pdf'], visual['source_pdf_sha256'])
    verify(relative(OUT/'independent-review.json'), visual['matrix_review_sha256'])
    baseline, review = records['reconstruction'], records['independent-review']
    assert abs(baseline['H0']['median_km_s_Mpc']-review['H0_median']) < 1e-8
    assert baseline['validation']['degrees_of_freedom'] == baseline['rows']-baseline['parameters'] == review['df']
    for row in baseline['released_cross_covariance_blocks']:
        if row['left'] == 'SN_calibrators':
            assert row['nonzero'] == 77*277 and row['maximum_absolute_correlation'] > .27
    factor, factor_review = records['cepheid-factor'], records['independent-factor-review']
    assert factor['removed_SN_rows'] == factor_review['removed_SN_rows'] == 354
    assert factor['data_rows'] == factor_review['retained_rows'] == 3138
    assert factor['host_order'] == [item['host'] for item in factor_review['host_labels']]
    assert [item['Cepheid_rows'] for item in factor['host_mapping']] == [item['Cepheid_rows'] for item in factor_review['host_labels']]
    covariance = np.array(factor['covariance_distance_modulus'])
    assert covariance.shape == (37, 37) and np.allclose(covariance, covariance.T, rtol=0, atol=1e-16)
    linalg.cholesky(covariance)
    assert np.count_nonzero(covariance-np.diag(np.diag(covariance))) > 0
    saved = np.load(WORK/'cepheid-host-factor.npz')
    assert np.array_equal(saved['covariance_mu'], covariance)
    assert np.array_equal(saved['mean_mu'], np.array(factor['mean_distance_modulus']))
    assert saved['host'].tolist() == factor['host_order']
    # Independent vector atan2 separation, compared to the producer's Astropy angle.
    o = pd.read_csv(WORK/'SN-overlap.csv', dtype={'pantheon_CID': str})
    p = pd.read_csv(ROOT/'.work/unified-cosmology/inference/pantheon/Pantheon+SH0ES.dat', sep=r'\s+', dtype={'CID': str})
    d = pd.read_csv(ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-ledger.csv').set_index('row')
    dr = d.loc[o.dovekie_row]; pr = p.iloc[o.pantheon_row]
    def unit(ra, dec):
        a, b = np.radians(np.asarray(ra)), np.radians(np.asarray(dec))
        return np.column_stack((np.cos(a)*np.cos(b), np.sin(a)*np.cos(b), np.sin(b)))
    uv, vv = unit(dr.RA, dr.DEC), unit(pr.RA, pr.DEC)
    separation = np.degrees(np.arctan2(np.linalg.norm(np.cross(uv, vv), axis=1), (uv*vv).sum(axis=1)))*3600
    angular_error = float(np.max(abs(separation-o.angular_separation_arcsec)))
    assert angular_error < 1e-7
    dz = abs(dr.zHEL.to_numpy()-pr.zHEL.to_numpy())
    assert np.max(abs(dz-o.absolute_delta_zHEL)) < 1e-14
    accepted = o.status == 'exact_alias_sky_z'
    assert np.all((separation[accepted] <= 1) & (dz[accepted] <= .001))
    assert np.array_equal(o.IS_CALIBRATOR.to_numpy(), pr.IS_CALIBRATOR.to_numpy())
    assert np.array_equal(o.USED_IN_SH0ES_HF.to_numpy(), pr.USED_IN_SH0ES_HF.to_numpy())
    matched = o[accepted]
    assert len(matched) == records['overlap']['accepted']['measurement_pairs'] == 301
    assert matched.dovekie_physical_ID.nunique() == 273
    hf = matched[matched.USED_IN_SH0ES_HF == 1]
    assert len(hf) == 95 and hf.dovekie_physical_ID.nunique() == 76
    assert int(matched.IS_CALIBRATOR.sum()) == 0
    parsed = []
    for source in sorted(HERE.glob('*.py')):
        ast.parse(source.read_text(), filename=str(source)); parsed.append(relative(source))
    docs = [HERE/'README.md', ROOT/'studies/unified_cosmology/notes/distance-ladder.md']
    links = []
    for doc in docs:
        for target in re.findall(r'\]\(([^)]+)\)', doc.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            path = (doc.parent/target.split('#')[0]).resolve()
            assert path.exists(), (doc, target)
            links.append(relative(path))
    source_files = sorted(HERE.glob('*'))
    result_files = [OUT/(name+'.json') for name in names]
    result = {'status': 'passed', 'acquired_files': len(records['acquisition']['files']),
              'verified_hash_paths': len(checked), 'AST_checked_scripts': parsed, 'local_links_checked': len(links),
              'independent_overlap_angle_max_difference_arcsec': angular_error,
              'retained_H0_km_s_Mpc': baseline['H0']['median_km_s_Mpc'],
              'H0_is_independent_cosmology_prior': False, 'SN_free_factor_hosts': 37,
              'overlap_events': 273, 'SH0ES_flagged_overlap_events': 76,
              'scientific_limits': ['Conditional compressed observations and original calibration/selection/anchor assumptions.',
                                    'The current Pantheon SN table is not numerically identical to frozen ladder ordinates.',
                                    'No cross-release calibration covariance or complete raw-SN integration has been inferred.',
                                    'No cosmological targets, likelihoods or chains were modified or evaluated.'],
              'dependencies_sha256': {relative(path): sha(path) for path in source_files+result_files+docs if path.is_file()},
              'verified_input_and_result_sha256': checked}
    write(OUT/'validation.json', result)
    print(json.dumps({k: result[k] for k in ['status','acquired_files','verified_hash_paths','local_links_checked','independent_overlap_angle_max_difference_arcsec']}))


if __name__ == '__main__':
    main()
