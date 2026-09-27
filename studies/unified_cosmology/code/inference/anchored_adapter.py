"""Separate anchored Pantheon+SH0ES target; no Dovekie/H0-prior combination."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import numpy as np
from cobaya.likelihood import Likelihood

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent/'distance_ladder'))
from calibration_interface import ReleasedCalibration, WORK as CALIBRATION_WORK, SOURCES, COMMIT
from modern_fast import configuration as original_configuration, identify as original_identify
from target_identity import canonical, digest

SAMPLE = 'pantheon_shoes_anchored'
DATA = CALIBRATION_WORK/'Pantheon+SH0ES.dat'
CALIBRATION_RECORD = ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-interface.json'
DESIGN = HERE/'anchored-design.json'


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


class AnchoredReleasedSN(Likelihood):
    """Cobaya distance provider for the independently audited release interface."""
    data_file: str = str(DATA)
    input_params = ['epsilon']
    type = 'SN'
    speed = 100.
    params = {'sn_chi2': {'derived': True}}

    def initialize(self):
        path = Path(self.data_file).resolve()
        assert path.name == DATA.name, 'Use the pinned raw release, not a distance-only NPZ.'
        self.release = ReleasedCalibration(path.parent)
        assert len(self.release.data) == 1657 and int(self.release.calibrator.sum()) == 77
        self.background_z, self.inverse_z = np.unique(self.release.z_hd_noncalibrator, return_inverse=True)
        assert len(self.inverse_z) == 1580 and np.all(self.background_z > .01)

    def get_requirements(self):
        return {'angular_diameter_distance': {'z': self.background_z}}

    def logp(self, epsilon=0., _derived=None):
        if epsilon != 0.:
            raise ValueError('This anchored target fixes luminosity evolution to zero.')
        unique_da = np.asarray(self.provider.get_angular_diameter_distance(self.background_z))
        assert unique_da.shape == self.background_z.shape
        result = self.release.evaluate(unique_da[self.inverse_z])
        if _derived is not None:
            _derived['sn_chi2'] = result['chi2']
        return result['loglike']


def calibration_lineage(record_path=CALIBRATION_RECORD):
    """Fail closed if the independent release audit or any input changed."""
    record_path = Path(record_path).resolve()
    record = json.loads(record_path.read_text())
    assert record['status'] == 'passed_released_Gaussian_interface_checks_no_cosmological_fit'
    assert record['released_target_usable'] is True
    assert record['independent_host_factor_replacement_certified'] is False
    assert record['official_commit'] == COMMIT
    expected_inputs = {relative(CALIBRATION_WORK/name): value[1] for name, value in SOURCES.items()}
    assert {path: item['sha256'] for path, item in record['input_sources'].items()} == expected_inputs
    producer = HERE.parent/'distance_ladder/calibration_interface.py'
    assert record['source_sha256'] == {relative(producer): digest(producer)}
    selection = record['selection']
    assert [selection[k] for k in ['selected_rows', 'calibrator_rows', 'noncalibrator_rows']] == [1657, 77, 1580]
    assert record['likelihood']['extra_SH0ES_or_Cepheid_factor'] is False
    assert record['likelihood']['normalized_evidence_claim'] is False
    inputs = {}
    for path, item in record['input_sources'].items():
        assert digest(ROOT/path) == item['sha256'], 'Anchored release input changed: '+path
        inputs[path] = item['sha256']
    for mapping in [record['source_sha256'], record['dependencies_sha256'], record['output_sha256'], record['validation']['source_sha256']]:
        for path, expected in mapping.items():
            assert digest(ROOT/path) == expected, 'Anchored audit evidence changed: '+path
            inputs[path] = expected
    inputs[relative(record_path)] = digest(record_path)
    return {'audit_record': relative(record_path), 'audit_sha256': digest(record_path),
            'official_commit': record['official_commit'], 'selection': selection,
            'likelihood': record['likelihood'], 'input_and_audit_sha256': inputs,
            'calibration_decomposition': 'Released covariance retained whole; embedded Cepheid decomposition is not certified and no extra calibration factor is added.'}


def configuration(model='lcdm', evolution='none', sample=SAMPLE, calibration='official_planck', surrogate=None, data_file=DATA):
    assert model in {'lcdm', 'cpl'} and evolution == 'none'
    assert sample == SAMPLE and calibration == 'official_planck'
    # Original factory constructs configuration only. No Dovekie data are read;
    # its placeholder SN block is completely replaced before model construction.
    info = original_configuration(model=model, evolution='none', sample='dovekie',
                                  calibration=calibration, surrogate=surrogate)
    info['likelihood']['released_sn'] = {'external': AnchoredReleasedSN, 'data_file': str(Path(data_file).resolve())}
    assert info['params']['epsilon'] == 0.
    return info


def identify(info, sample_file=DATA, model_file=None, calibration_record=CALIBRATION_RECORD):
    sample_file = Path(sample_file).resolve()
    lineage = calibration_lineage(calibration_record)
    assert relative(sample_file) in lineage['input_and_audit_sha256']
    config = canonical(info)
    assert config['likelihood']['released_sn'] == {'external': __name__+'.AnchoredReleasedSN', 'data_file': str(sample_file)}
    assert config['params']['epsilon'] == 0.
    assert set(config['likelihood']) == set(original_configuration()['likelihood'])
    record = original_identify(info, sample_file, model_file)
    for path in [Path(__file__), DESIGN, HERE.parent/'distance_ladder/calibration_interface.py']:
        record['source_sha256'][relative(path)] = digest(path)
    record['versions']['pandas'] = importlib.metadata.version('pandas')
    record['anchored_calibration'] = lineage
    record['sample_semantics'] = {'name': SAMPLE, 'replaces_other_SN_compilations': True,
        'extra_H0_or_Cepheid_factor': False, 'absolute_calibration_included': True,
        'flat_M_measure': 'Same improper flat dM convention; full Gaussian and intercept determinant retained, no evidence claim.'}
    record.pop('identity')
    record['identity'] = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return record
