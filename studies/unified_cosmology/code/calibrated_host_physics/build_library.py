#!/usr/bin/env python3
"""Build the separately compiled, nebular-free C3K_HR SSP basis."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from datetime import datetime, timezone
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
CODE = Path(__file__).parent
WORK = ROOT / '.work/unified-cosmology/calibrated-host-physics'
RESULT = ROOT / 'studies/unified_cosmology/results/calibrated_host_physics'
os.environ['SPS_HOME'] = str(ROOT / '.work/physical-ages/fsps')
import fsps


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--refined', action='store_true')
    args = parser.parse_args()
    design = json.loads((CODE / 'design.json').read_text())
    ages = np.array(design['ages_Gyr'])
    metals = np.array(design['metallicities_logZ_solar'])
    if args.refined:
        ages = np.sort(np.r_[ages, np.sqrt(ages[:-1] * ages[1:])])
        metals = np.sort(np.r_[metals, (metals[:-1] + metals[1:]) / 2])
    start = time.monotonic()
    pop = fsps.StellarPopulation(zcontinuous=1, sfh=0, imf_type=1,
        add_dust_emission=False, add_neb_emission=False, dust1=0., dust2=0.)
    libraries = [s.decode() if isinstance(s, bytes) else s for s in pop.libraries]
    assert libraries[1] == 'c3k_hr', libraries
    spectra, masses = [], []
    for metal in metals:
        pop.params['logzsol'] = metal
        rows, mass = [], []
        for age in ages:
            wave, flux = pop.get_spectrum(tage=age, peraa=True)
            assert abs(pop.formed_mass - 1.) < 1e-12
            assert np.isfinite(flux).all() and np.min(flux) >= 0
            rows.append(flux.copy())
            mass.append(pop.stellar_mass)
        spectra.append(rows)
        masses.append(mass)
        print('Completed C3K_HR logZ/Zsolar', metal, flush=True)
    suffix = '-refined' if args.refined else ''
    output = WORK / ('ssp-library' + suffix + '.npz')
    np.savez_compressed(output, wavelength_A=wave,
        spectra_Lsun_per_A_per_formed_Msun=np.array(spectra),
        ages_Gyr=ages, metallicity_logZ_solar=metals,
        surviving_mass_fraction=np.array(masses),
        resolution_sigma_km_s=np.array(pop.resolutions),
        emission_line_wavelength_A=np.array(pop.emline_wavelengths))
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic() - start, libraries=libraries,
        ages_Gyr=ages.tolist(), metallicities_logZ_solar=metals.tolist(),
        wavelengths=len(wave), spectra=len(ages) * len(metals),
        semantics='Vacuum stellar-only SSP Llambda per formed solar mass; surviving mass is distinct. No nebular emission, host or foreground attenuation in the library.',
        dependencies_sha256={str(p.relative_to(ROOT)): sha(p) for p in
            [Path(__file__), CODE / 'design.json', CODE / 'implementation-design.json', RESULT / 'acquisition.json']},
        output_sha256={str(output.relative_to(ROOT)): sha(output)})
    (RESULT / ('stellar-library' + suffix + '.json')).write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
