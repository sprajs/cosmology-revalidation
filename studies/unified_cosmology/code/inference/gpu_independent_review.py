"""Independent concurrent FP64 GEMV audit; this is not a chain-speed estimate.

Run with the isolated modern likelihood environment. Each spawned worker uses
one CPU thread and identical predeclared synthetic features in the CPU and GPU
phases. Initialization and the five CPU reference checks are outside the timed
100-vector loop; GPU host/device transfers and synchronization are inside it.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import hashlib
import importlib.metadata
import json
import multiprocessing as mp
from pathlib import Path
import subprocess
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
MODEL = ROOT / '.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
EXPECTED_MODEL = 'dfff858fc0d468ab04902e0608f123a2fecf761901768c76b9bb39f6843b34dc'
WORKERS, VECTORS, SEED = 8, 100, 272901


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def memory():
    command = ['nvidia-smi', '--query-gpu=name,memory.total,memory.used,memory.free,driver_version',
               '--format=csv,noheader,nounits']
    output = subprocess.check_output(command, text=True).strip()
    name, total, used, free, driver = (v.strip() for v in output.split(','))
    return dict(device=name, total_MiB=int(total), used_MiB=int(used),
                free_MiB=int(free), driver=driver)


def worker(mode, index, queue, start):
    try:
        from spectral_surrogate import polynomial
        with np.load(MODEL, allow_pickle=False) as source:
            coefficients = source['coefficients']
            exponents = source['exponents']
        features = polynomial(np.random.default_rng(SEED + index).normal(size=(VECTORS, 8)), exponents)
        matrix = coefficients
        runtime = None
        if mode == 'gpu':
            from modern_gpu import GPUMatrix
            import cupy as cp
            matrix = GPUMatrix(coefficients)
            runtime = dict(runtime=cp.cuda.runtime.runtimeGetVersion(),
                           driver=cp.cuda.runtime.driverGetVersion(),
                           allocated_bytes=int(cp.get_default_memory_pool().used_bytes()))
        warmup_error = 0.
        for vector in features[:2]:
            actual = vector @ matrix
            reference = vector @ coefficients
            warmup_error = max(warmup_error, float(np.max(abs(actual-reference)/(1+abs(reference)))))
        assert warmup_error < 1e-11
        queue.put(dict(event='ready', mode=mode, worker=index, pid=os.getpid(),
                       warmup_error=warmup_error, runtime=runtime,
                       matrix_shape=list(coefficients.shape), matrix_bytes=coefficients.nbytes))
        assert start.wait(timeout=180), 'benchmark start timeout'
        checkpoints = []
        began = time.perf_counter()
        for i, vector in enumerate(features):
            value = vector @ matrix
            if i % 20 == 0:
                checkpoints.append((i, value))
        elapsed = time.perf_counter()-began
        scaled_error, absolute_error = 0., 0.
        for i, actual in checkpoints:
            reference = features[i] @ coefficients
            scaled_error = max(scaled_error, float(np.max(abs(actual-reference)/(1+abs(reference)))))
            absolute_error = max(absolute_error, float(np.max(abs(actual-reference))))
        assert scaled_error < 1e-11
        queue.put(dict(event='done', mode=mode, worker=index, seconds=elapsed,
                       comparison_indices=[i for i, _ in checkpoints],
                       maximum_scaled_error=scaled_error, maximum_absolute_error=absolute_error))
    except BaseException:
        queue.put(dict(event='error', mode=mode, worker=index, traceback=traceback.format_exc()))


def phase(context, mode):
    queue, start = context.Queue(), context.Event()
    processes = [context.Process(target=worker, args=(mode, i, queue, start)) for i in range(WORKERS)]
    began = time.perf_counter()
    for process in processes:
        process.start()
    ready = []
    try:
        while len(ready) < WORKERS:
            event = queue.get(timeout=180)
            assert event['event'] == 'ready', event
            ready.append(event)
        allocated = memory()
        # Preserve room for the display/server; refuse the timed phase if eight
        # actual contexts leave less than 1 GiB of the reported usable memory.
        assert allocated['free_MiB'] > 1024, allocated
        print(json.dumps(dict(event='all_workers_allocated', mode=mode,
                              memory=allocated, maximum_warmup_error=max(r['warmup_error'] for r in ready))), flush=True)
        timed = time.perf_counter()
        start.set()
        done = []
        while len(done) < WORKERS:
            event = queue.get(timeout=180)
            assert event['event'] == 'done', event
            done.append(event)
        timed_wall = time.perf_counter()-timed
        for process in processes:
            process.join(timeout=30)
            assert process.exitcode == 0, process.exitcode
        return dict(mode=mode, workers=WORKERS, vectors_per_worker=VECTORS,
                    allocation_memory=allocated, ready=sorted(ready, key=lambda v:v['worker']),
                    workers_results=sorted(done, key=lambda v:v['worker']),
                    maximum_worker_seconds=max(r['seconds'] for r in done),
                    timed_phase_wall_seconds=timed_wall,
                    whole_pool_wall_seconds=time.perf_counter()-began,
                    after_context_exit_memory=memory())
    finally:
        for process in processes:
            if process.is_alive():
                process.terminate()
            process.join(timeout=10)


def main():
    assert digest(MODEL) == EXPECTED_MODEL
    existing = ROOT / 'studies/unified_cosmology/results/inference/gpu-validation.json'
    full = json.loads(existing.read_text())
    assert full['status'] == 'passed_numerical_equivalence' and len(full['cases']) == 24
    for path, expected in {**full['source_sha256'], **full['input_sha256']}.items():
        assert digest(ROOT/path) == expected, path
    assert max(max(row['maximum_errors'].values()) for row in full['cases']) < 1e-6
    before = memory()
    context = mp.get_context('spawn')
    # GPU first gives the requested safe-allocation answer promptly. Processes
    # exit fully before the matched CPU phase, so GPU memory is not held idle.
    gpu = phase(context, 'gpu')
    cpu = phase(context, 'cpu')
    result = dict(status='passed_independent_concurrent_numerical_review',
        workers=WORKERS, vectors_per_worker=VECTORS, synthetic_coordinate_seed_base=SEED,
        sampling='Independent standard-normal eight-coordinate features per worker; identical CPU/GPU draws.',
        maximum_allowed_scaled_error=1e-11, memory_before=before, memory_after=memory(),
        gpu=gpu, cpu=cpu,
        component_throughput_speedup=cpu['maximum_worker_seconds']/gpu['maximum_worker_seconds'],
        full_likelihood_review=dict(cases=24, report_sha256=digest(existing),
            maxima={k:max(row['maximum_errors'][k] for row in full['cases'])
                    for k in ['logpost','loglikes','logpriors','derived']},
            scope='Reviewed existing root test/source; no duplicated native or joint likelihood calls.'),
        source_review=['The only substituted arithmetic is FP64 polynomial-feature matrix multiplication.',
            'asnumpy synchronizes and returns an ordinary CPU float64 vector.',
            'Polynomial features, scaling, positivity/envelope gates, background units and exact fallback are inherited unchanged.',
            'Configuration equality after replacing only the external class is checked for all three targets by the bound root report.',
            'No covariance, scientific prior, nuisance prior or likelihood expression is changed.'],
        scope='Concurrent component safety/equivalence only. Timings exclude initialization and reference checks but include GPU transfers. '
              'Not a full-likelihood/chain speed estimate, posterior convergence or interpolation accuracy certificate; exact posterior correction remains required.',
        versions={name:importlib.metadata.version(name) for name in ['numpy','cupy-cuda12x','nvidia-cublas-cu12','nvidia-cuda-runtime-cu12','nvidia-cuda-nvrtc-cu12']},
        dependencies_sha256={str(path.relative_to(ROOT)):digest(path) for path in
            [Path(__file__),Path(__file__).with_name('modern_gpu.py'),Path(__file__).with_name('spectral_surrogate.py'),
             Path(__file__).with_name('modern_fast.py'),Path(__file__).with_name('gpu_validate.py'),MODEL,existing]})
    output = ROOT / 'studies/unified_cosmology/results/inference/gpu-independent-review.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], speedup=result['component_throughput_speedup'],
        maximum_gpu_scaled_error=max(row['maximum_scaled_error'] for row in gpu['workers_results']),
        memory_after=result['memory_after']), indent=2), flush=True)


if __name__ == '__main__':
    main()
