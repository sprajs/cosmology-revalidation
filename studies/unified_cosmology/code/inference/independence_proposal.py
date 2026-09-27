"""Normalized frozen Gaussian/Student proposal; training changes no prior."""
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[4]
DESIGN = Path(__file__).with_name('independence-sampling-design.json')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


class FrozenMixture:
    def __init__(self, mean, covariance, names):
        self.mean = np.asarray(mean, dtype=float)
        self.covariance = np.asarray(covariance, dtype=float)
        self.names = list(names)
        self.d = len(self.mean)
        assert len(self.names) == len(set(self.names)) == self.d
        assert self.covariance.shape == (self.d, self.d)
        assert np.isfinite(self.mean).all() and np.isfinite(self.covariance).all()
        assert np.allclose(self.covariance, self.covariance.T, rtol=1e-12, atol=0)
        self.L = np.linalg.cholesky(self.covariance)
        self.logdet = 2 * np.log(np.diag(self.L)).sum()
        self.df = 5.
        self.gaussian_variance_factor = 1.05**2
        self.student_scale_factor = 4 * (self.df - 2) / self.df

    def logpdf(self, points):
        x = np.asarray(points, dtype=float)
        white = solve_triangular(self.L, (np.atleast_2d(x)-self.mean).T, lower=True).T
        radius2 = np.sum(white**2, axis=1)
        g = -.5 * (self.d*np.log(2*np.pi) + self.logdet
                   + self.d*np.log(self.gaussian_variance_factor)
                   + radius2/self.gaussian_variance_factor)
        t = (gammaln((self.df+self.d)/2) - gammaln(self.df/2)
             - .5*(self.d*np.log(self.df*np.pi) + self.logdet
                    + self.d*np.log(self.student_scale_factor))
             - (self.df+self.d)/2 * np.log1p(radius2/self.student_scale_factor/self.df))
        value = np.logaddexp(np.log(.9)+g, np.log(.1)+t)
        return float(value[0]) if x.ndim == 1 else value

    def draw(self, rng):
        student = bool(rng.random() >= .9)
        z = rng.normal(size=self.d)
        if student:
            z *= np.sqrt(self.student_scale_factor*self.df/rng.chisquare(self.df))
        else:
            z *= np.sqrt(self.gaussian_variance_factor)
        return self.mean + self.L@z, 'student_t5' if student else 'gaussian'

    @classmethod
    def read(cls, path, expected_sha256):
        assert digest(path) == expected_sha256, 'Frozen proposal changed.'
        with np.load(path, allow_pickle=False) as data:
            return cls(data['mean'], data['cov'], data['names'].tolist())


def freeze(folder, output, manifest_folder=None):
    """Copy complete source rows before learning; snapshots remain unqualified."""
    folder, output = Path(folder).resolve(), Path(output).resolve()
    assert not output.exists(), 'Refuse to overwrite a frozen proposal.'
    manifest_folder = Path(manifest_folder).resolve() if manifest_folder else folder
    manifests = [json.loads((manifest_folder/f'run-{rank}.json').read_text()) for rank in range(4)]
    target = manifests[0]['target_identity']
    assert all(m['target_identity'] == target for m in manifests)
    names = [k for k, v in target['configuration']['params'].items()
             if isinstance(v, dict) and 'prior' in v]
    assert not any(target['configuration']['params'][n].get('periodic', False) for n in names)
    # Actual scientific bytes are checked before consuming an active-chain snapshot.
    from measurement_summary import verify_current_target
    source_inputs = verify_current_target(manifests[0])
    output.mkdir(parents=True)
    copies = {}; training = []; sizes = []
    for rank in range(4):
        for name in [f'run-{rank}.json', f'chain.{rank+1}.txt']:
            source = (manifest_folder if name.startswith('run-') else folder)/name; content = source.read_bytes()
            if name.endswith('.txt'):
                content = content[:content.rfind(b'\n')+1]
            destination = output/name; destination.write_bytes(content)
            copies[name] = {'sha256': digest(destination), 'source': relative(source)}
        path = output/f'chain.{rank+1}.txt'
        headers = path.open().readline().lstrip('#').split()
        rows = np.loadtxt(path, ndmin=2)
        weights = rows[:, 0]
        assert (weights > 0).all() and np.array_equal(weights, weights.astype(int))
        expanded = np.repeat(rows[:, [headers.index(n) for n in names]], weights.astype(int), axis=0)
        start, stop = int(.5*len(expanded)), int(.9*len(expanded))
        assert stop-start >= 100
        training.append(expanded[start:stop])
        sizes.append({'rank': rank, 'frequency_draws': len(expanded), 'training_start': start,
                      'training_stop': stop, 'withheld_rows': len(expanded)-stop})
    data = np.concatenate(training)
    mean = data.mean(axis=0); centered = data-mean
    covariance = centered.T@centered/len(data)
    covariance += 1e-8*np.diag(np.diag(covariance))
    FrozenMixture(mean, covariance, names)
    path = output/'proposal.npz'
    np.savez(path, mean=mean, cov=covariance, names=names)
    record = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope': 'Frozen proposal training only; no posterior qualification.',
              'parent_target_identity': target['identity'], 'parent_settings': manifests[0]['arguments'],
              'names': names, 'sizes': sizes, 'snapshot_files': copies,
              'scientific_input_sha256': source_inputs,
              'proposal_sha256': digest(path), 'design_sha256': digest(DESIGN),
              'source_sha256': digest(__file__)}
    (output/'proposal.json').write_text(json.dumps(record, indent=2)+'\n')
    return record
