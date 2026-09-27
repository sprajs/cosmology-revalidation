#!/usr/bin/env python3
"""Benchmark actual DES stellar-population forward matrices on CPU and CUDA.

Run --prepare with the science environment, then without it in the isolated
GPU environment. No device drivers or system CUDA installation are changed.
"""
import os,sys,sysconfig,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/physical-ages'
if '--prepare' not in sys.argv and not os.environ.get('COSMO_GPU_LIBRARIES_READY'):
    site=Path(sysconfig.get_paths()['purelib'])
    paths=[str(site/'nvidia'/name/'lib') for name in ['cuda_nvrtc','cuda_runtime','cublas']]
    env=dict(os.environ);env['LD_LIBRARY_PATH']=':'.join(paths+[env.get('LD_LIBRARY_PATH','')]);env['COSMO_GPU_LIBRARIES_READY']='1'
    os.execve(sys.executable,[sys.executable,*sys.argv],env)
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np

def main():
    path=WORK/'gpu-forward-fixture.npz'
    if '--prepare' in sys.argv:
        from fit import inputs
        from model import Grid,Library
        grid=Grid(Library(response_path=ROOT/'.work/host-transport/filters/des-deep-responses.npz'))
        items,_,_=inputs('des');kernels=np.concatenate([grid.kernels(x['z'],x['ebv']).T for x in items],axis=1)
        np.savez_compressed(path,spectra=grid.nu_flat,kernels=kernels,ids=[x['id'] for x in items],redshift=[x['z'] for x in items])
        dependencies=[Path(__file__).with_name(n) for n in ['model.py','fit.py','design.json']]+[WORK/'ssp-library.npz',ROOT/'.work/host-transport/filters/des-deep-responses.npz',ROOT/'.work/host-transport/des-deep-photometry.csv']
        (WORK/'gpu-forward-fixture.json').write_text(json.dumps({'fixture_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'dependencies_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}},indent=2)+'\n')
        print('Prepared actual 331-host, 68-nuisance-cell, 14-SSP forward calculation');return
    import cupy as cp
    fixture_record=json.loads((WORK/'gpu-forward-fixture.json').read_text())
    assert fixture_record['fixture_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    for name,expected in fixture_record['dependencies_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    data=np.load(path);a=np.ascontiguousarray(data['spectra']);b=np.ascontiguousarray(data['kernels'])
    def clock(operation,repeat=5):
        samples=[];result=None
        for _ in range(repeat):
            start=time.perf_counter();result=operation();samples.append(time.perf_counter()-start)
        return result,float(np.median(samples))
    cold=time.perf_counter();cp.cuda.Device(0).use();ga=cp.asarray(a);gb=cp.asarray(b);first=cp.asnumpy(ga@gb);cold=time.perf_counter()-cold
    reference,cpu=clock(lambda:a@b)
    def gpu_sync():
        result=ga@gb;cp.cuda.Stream.null.synchronize();return result
    _,gpu=clock(gpu_sync)
    gpu_output,end_to_end=clock(lambda:cp.asnumpy(cp.asarray(a)@cp.asarray(b)))
    rel=float(np.max(abs(gpu_output-reference)/np.maximum(abs(reference),1e-300)))
    assert rel<2e-12
    # Float32 is benchmarked but not admitted to the scientific fits.
    fa=cp.asarray(a.astype(np.float32));fb=cp.asarray(b.astype(np.float32))
    foutput,fseconds=clock(lambda:cp.asnumpy(fa@fb))
    frel=float(np.max(abs(foutput-reference)/np.maximum(abs(reference),1e-300)))
    single,cpusingle=clock(lambda:a@b[:,:8],20)
    def single_gpu():return cp.asnumpy(ga@gb[:,:8])
    _,gpusingle=clock(single_gpu,20)
    gpu_path=WORK/'des-gpu-forward-flux.npy';np.save(gpu_path,gpu_output)
    result=dict(device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),driver_version=cp.cuda.runtime.driverGetVersion(),runtime_version=cp.cuda.runtime.runtimeGetVersion(),cupy=cp.__version__,numpy=np.__version__,
        shapes={'stellar_spectra':list(a.shape),'combined_observer_kernels':list(b.shape)},objects=len(data['ids']),
        cold_GPU_seconds=cold,median_seconds={'CPU_float64_all_hosts':cpu,'GPU_float64_resident':gpu,'GPU_float64_with_host_transfers':end_to_end,'GPU_float32_resident_with_output_transfer':fseconds,'CPU_float64_one_host':cpusingle,'GPU_float64_one_host_resident_inputs_and_output_transfer':gpusingle},
        double_max_relative_difference=rel,float32_max_relative_difference=frel,CPU_to_GPU_batch_speedup_including_transfers=cpu/end_to_end,
        interpretation='Actual DES forward-model matrix products, not random benchmark arrays. GPU result calculated and retained locally; published bounds still use checked float64 CPU products. Small conic optimizations dominate the full analysis and remain compiled CPU work; matrix speedup is not whole-program speedup.',
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),dependencies_sha256=fixture_record['dependencies_sha256'],GPU_result_sha256=hashlib.sha256(gpu_path.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'studies/host_ages/results/physical_ages/compute-benchmark.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
