"""FP64 GPU matrix evaluation for the unchanged spectral proposal.

The frozen SpectralSurrogate supplies every background, unit conversion,
envelope check and native fallback. Only its dense matrix multiplication is
dispatched to cuBLAS. Native-CAMB posterior correction remains mandatory.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from spectral_surrogate import SpectralSurrogate
from modern_fast import configuration as cpu_configuration, identify as cpu_identify
from target_identity import ROOT, digest


class GPUMatrix:
    """Right matrix operand with an explicit NumPy matmul dispatch only."""
    def __init__(self, coefficients):
        import cupy as cp
        assert coefficients.dtype == np.float64 and coefficients.ndim == 2
        self.matrix = cp.asarray(coefficients)
        self.shape = coefficients.shape
        cp.cuda.Stream.null.synchronize()

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if ufunc is not np.matmul or method != '__call__' or kwargs \
                or len(inputs) != 2 or inputs[1] is not self:
            return NotImplemented
        import cupy as cp
        features = np.asarray(inputs[0])
        assert features.dtype == np.float64 and features.shape == (self.shape[0],)
        # asnumpy synchronizes and returns an ordinary float64 array. No reduced
        # precision, asynchronous buffer reuse or change to polynomial features.
        return cp.asnumpy(cp.asarray(features) @ self.matrix)


class GPUSpectralSurrogate(SpectralSurrogate):
    def initialize(self):
        super().initialize()
        self.coefficients = GPUMatrix(self.coefficients)


def configuration(*args, **kwargs):
    info = cpu_configuration(*args, **kwargs)
    if 'spectral_surrogate' in info['theory']:
        info['theory']['spectral_surrogate']['external'] = GPUSpectralSurrogate
    return info


def identify(info, sample_file, model_file=None):
    import cupy as cp
    record = cpu_identify(info, sample_file, model_file)
    record['source_sha256'][str(Path(__file__).relative_to(ROOT))] = digest(__file__)
    for package in ['cupy-cuda12x', 'nvidia-cublas-cu12', 'nvidia-cuda-runtime-cu12',
                    'nvidia-cuda-nvrtc-cu12']:
        record['versions'][package] = importlib.metadata.version(package)
    record['gpu_numerics'] = {
        'dtype': 'float64', 'CUDA_runtime': cp.cuda.runtime.runtimeGetVersion(),
        'CUDA_driver': cp.cuda.runtime.driverGetVersion(),
        'device_name': cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
        'scope': 'Spectral proposal matrix multiplication only; native physics unchanged.'}
    record.pop('identity')
    record['identity'] = hashlib.sha256(json.dumps(record, sort_keys=True,
        separators=(',', ':')).encode()).hexdigest()
    return record
