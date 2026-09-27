"""Conservative event overlap; event identity is not cross-release covariance."""
import re
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from common import ROOT, HERE, WORK, OUT, sha, relative, write

INPUTS = {
    '.work/unified-cosmology/survey-selection/normalized/dovekie-ledger.csv': '281b1da906b1f822643cda501c4e81ebd614e5e13f65f0c0d079ae06d26d46a2',
    '.work/unified-cosmology/inference/pantheon/Pantheon+SH0ES.dat': '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8',
}


def normalize(value):
    text = str(value).strip().lower()
    if text in {'', 'nan', 'unknown', 'null', '-9'}:
        return ''
    return text[2:] if re.match(r'^sn[12][0-9]{3}', text) else text


def counts(frame):
    return {'measurement_pairs': len(frame), 'Pantheon_measurement_rows': int(frame.pantheon_row.nunique()),
            'Dovekie_events': int(frame.dovekie_physical_ID.nunique()), 'Pantheon_CID_strings': int(frame.pantheon_CID.nunique())}


def main():
    for path, digest in INPUTS.items():
        assert sha(ROOT/path) == digest, path
    dp, pp = [ROOT/p for p in INPUTS]
    d = pd.read_csv(dp, dtype={'CID': str, 'IAU_alias': str, 'transient_alias': str}).fillna('')
    p = pd.read_csv(pp, sep=r'\s+', dtype={'CID': str})
    assert d.physical_ID.is_unique
    keys = p.CID.map(normalize).to_numpy()
    lookup = {}
    for i, key in enumerate(keys):
        if key:
            lookup.setdefault(key, []).append(i)
    valid = (p.RA >= 0) & (p.RA < 360) & (p.DEC >= -90) & (p.DEC <= 90)
    coords = SkyCoord(ra=p.RA.where(valid, 0).to_numpy()*u.deg, dec=p.DEC.where(valid, 0).to_numpy()*u.deg)
    rows = []
    for _, event in d.iterrows():
        aliases = {normalize(event[k]): k for k in ('CID', 'IAU_alias', 'transient_alias') if normalize(event[k])}
        valid_event = 0 <= float(event.RA) < 360 and -90 <= float(event.DEC) <= 90
        sep = SkyCoord(float(event.RA)*u.deg, float(event.DEC)*u.deg).separation(coords).arcsec if valid_event else np.full(len(p), np.inf)
        dz = abs(p.zHEL.to_numpy()-float(event.zHEL))
        candidates = {i for key in aliases for i in lookup.get(key, [])}
        candidates.update(np.flatnonzero(valid.to_numpy() & (sep <= 1) & (dz <= .001)).tolist())
        for i in sorted(candidates):
            row = p.iloc[i]
            key = keys[i]
            alias = key in aliases and (not key.isdigit() or int(row.IDSURVEY) == int(event.IDSURVEY))
            geometric = bool(valid.iloc[i] and valid_event and sep[i] <= 1 and dz[i] <= .001)
            if not alias and not geometric:
                continue
            status = ('exact_alias_sky_z' if geometric else 'exact_alias_failed_corroboration') if alias else 'sky_z_only_epoch_unverified'
            rows.append({'dovekie_row': int(event.row), 'dovekie_physical_ID': event.physical_ID,
                         'dovekie_CID': event.CID, 'dovekie_IDSURVEY': int(event.IDSURVEY),
                         'pantheon_row': i, 'pantheon_CID': row.CID, 'pantheon_IDSURVEY': int(row.IDSURVEY),
                         'alias_column': aliases.get(key, ''), 'status': status,
                         'angular_separation_arcsec': float(sep[i]), 'absolute_delta_zHEL': float(dz[i]),
                         'IS_CALIBRATOR': int(row.IS_CALIBRATOR), 'USED_IN_SH0ES_HF': int(row.USED_IN_SH0ES_HF),
                         'pantheon_PKMJD': float(row.PKMJD), 'pantheon_zHD_gt_0_01': bool(row.zHD > .01)})
    frame = pd.DataFrame(rows)
    accepted = frame[frame.status == 'exact_alias_sky_z']
    assert accepted.groupby('pantheon_row').dovekie_physical_ID.nunique().max() == 1
    assert accepted.groupby('dovekie_physical_ID').pantheon_CID.nunique().max() == 1
    output = WORK/'SN-overlap.csv'
    frame.to_csv(output, index=False)
    result = {
        'status': 'passed_conservative_overlap_inventory_not_joint_covariance',
        'gates': {'exact_alias': 'Case normalized; strip SN prefix only before four-digit year; purely numeric CID also requires identical survey ID.',
                  'separation_arcsec_max': 1, 'absolute_delta_zHEL_max': .001,
                  'unverified': 'Sky/redshift-only candidates are preserved but not accepted without epoch identity.'},
        'accepted': counts(accepted), 'accepted_calibrators': counts(accepted[accepted.IS_CALIBRATOR == 1]),
        'accepted_SH0ES_HF_flag': counts(accepted[accepted.USED_IN_SH0ES_HF == 1]),
        'alias_failed_corroboration': counts(frame[frame.status == 'exact_alias_failed_corroboration']),
        'sky_only_epoch_unverified': counts(frame[frame.status == 'sky_z_only_epoch_unverified']),
        'Pantheon_release_flag_counts': {key: {'measurement_rows': int((p[key] == 1).sum()), 'CID_strings': int(p.loc[p[key] == 1, 'CID'].nunique())}
                                        for key in ['IS_CALIBRATOR', 'USED_IN_SH0ES_HF']},
        'integration_limits': [
            'The flagged current Pantheon cohort is an event-overlap screen, not a proven one-to-one map of all rows in the earlier released ladder matrix.',
            'Corrected magnitudes differ between the products: e.g. current SN2011fe first two table values 9.74571/9.80286 versus ladder first two calibrator ordinates about9.752/9.808. No replacement or mixed covariance is performed.',
            'Strict position cuts may reject true aliases due to released coordinate rounding; all failed candidates remain available.',
            'The available full ladder covariance contains nonzero calibrator/Hubble-flow cross terms. Dovekie-to-ladder transformed calibration covariance or a shared latent calibration model is not supplied by these products.',
            'No calibrator event overlap under this gate does not establish independence of common surveys, standardization, calibration or population assumptions.',
        ],
        'input_sha256': INPUTS,
        'source_sha256': {relative(path): sha(path) for path in [HERE/'overlap.py', HERE/'common.py']},
        'output_sha256': {relative(output): sha(output)},
        'upstream_provenance': ['studies/unified_cosmology/results/survey_selection/acquisition.json',
                                'studies/unified_cosmology/results/inference/pantheon-interface.json'],
    }
    write(OUT/'overlap.json', result)
    print({key: result[key] for key in ['status', 'accepted', 'accepted_calibrators', 'accepted_SH0ES_HF_flag']})


if __name__ == '__main__':
    main()
