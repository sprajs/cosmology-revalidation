"""Persistent exact native USRFUN wrapper; each process owns one fixed SN state."""
from pathlib import Path
import subprocess,os,time,json,re
import numpy as np
import paired_refit as p
import native_profile_gate_v2 as ng
O=p.OUT;OUT=O/'fixed-c-profile';BUILD=OUT/'native-build'

class Oracle:
    def __init__(self,path,version='RSR_B',pars=None):
        self.path=path;self.counter=0;self.elapsed=0.;self.calls=[]
        if pars is None:pars=np.load(O/'native-profile-gate-v2/fixed/native.npz')['parameters']
        ng.config(path,version,pars)
        env=os.environ.copy();env.update(SNANA_DIR=str(BUILD),SNDATA_ROOT=str(p.ROOT/'phase2/official/inputs/SNDATA_ROOT'),
          LD_LIBRARY_PATH=str(p.ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RAISIN_PROBE_STREAM='1')
        self.log=(path/'fit.log').open('x')
        self.process=subprocess.Popen([str(BUILD/'bin/snlc_fit.exe'),str(path/'fit.nml')],cwd=path,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        self.started=time.monotonic()
        while True:
            line=self.process.stdout.readline()
            if not line:raise RuntimeError('Native oracle exited before READY; read '+str(path/'fit.log'))
            self.log.write(line)
            if line.startswith('RAISIN_READY:'):
                self.n=int(line.split()[1]);break
        self.initialization=time.monotonic()-self.started;self.log.flush()

    def mean(self,pars):
        pars=np.asarray(pars,float);assert pars.shape==(4,) and np.isfinite(pars).all()
        self.counter+=1;t=time.monotonic()
        self.process.stdin.write(str(self.counter)+' '+' '.join(format(x,'.17g') for x in pars)+'\n');self.process.stdin.flush()
        while True:
            line=self.process.stdout.readline()
            if not line:raise RuntimeError('Native oracle exited during mean query')
            self.log.write(line)
            if line.startswith('RAISIN_MEAN:'):
                q=line.split();assert int(q[1])==self.counter and int(q[2])==self.n
                a=np.array(list(map(float,q[3:])));assert np.array_equal(a[:4],pars)
                h=a[4:];assert h.shape==(self.n,) and np.isfinite(h).all()
                self.elapsed+=time.monotonic()-t
                return h

    def close(self):
        self.process.stdin.close()
        for line in self.process.stdout:self.log.write(line)
        self.process.wait(timeout=60);self.log.close()
        assert self.process.returncode==0
        return ng.export(self.path)

def gate():
    G=OUT/'oracle-gate';G.mkdir(exist_ok=False)
    prior=np.load(O/'native-profile-gate-v2/fixed/native.npz');base=prior['parameters']
    specs=[('fixed',base.copy()),('distance_plus',np.load(O/'native-profile-gate-v2/distance_plus/native.npz')['parameters']),
      ('distance_minus',np.load(O/'native-profile-gate-v2/distance_minus/native.npz')['parameters'])]
    for name,index,offset in [('shape_plus',1,.02),('shape_minus',1,-.02),('AV_plus',2,.1),('AV_minus',2,-.1)]:
        v=base.copy();v[index]+=offset;specs.append((name,v.astype(np.float32).astype(float)))
    protocol=dict(status='Frozen before oracle comparisons, benchmark or A/B profile outcomes',CID=ng.CID,reference='native-profile-gate-v2 fixed B',
      separate_native_probes=[dict(name=n,parameters=v.tolist()) for n,v in specs],
      benchmark='100 fixed-seed20260926 design vectors: shape .75–1.25, AV −.4–.8, D/peak fixed; repeat in deterministic random order and replay nominal after each batch',
      tolerances=dict(native_full_process_max_relative=2e-8,permuted_repeat_max_relative=2e-12,same_data_rows='exact native float32 values',same_C='exact array equality at returned reference state'),
      hashes={str(f.relative_to(p.ROOT)):p.sha(f) for f in [Path(__file__),OUT/'build-manifest.json',OUT/'mean-stream.patch',BUILD/'bin/snlc_fit.exe',O/'native-profile-gate-v2/manifest.json']})
    p.save(G/'protocol.json',protocol)
    standalone={}
    for n,v in specs:
        if n in ['fixed','distance_plus','distance_minus']:standalone[n]=np.load(O/'native-profile-gate-v2'/n/'native.npz')
        else:
            ng.config(G/n,'RSR_B',v);standalone[n]=ng.run(G/n)
    oracle=Oracle(G/'stream');dif={};saved={}
    for n,v in specs:
        h=oracle.mean(v);saved[n]=h
        e=float(np.max(abs(h-standalone[n]['model_flux']))/max(np.max(abs(h)),1e-30));assert e<2e-8,(n,e)
        dif[n]=e
    rng=np.random.default_rng(20260926);coords=np.tile(base,(100,1));coords[:,1]=rng.uniform(.75,1.25,100);coords[:,2]=rng.uniform(-.4,.8,100)
    t=time.monotonic();means=np.array([oracle.mean(v) for v in coords]);elapsed=time.monotonic()-t
    order=rng.permutation(100);repeat=np.empty_like(means)
    for i in order:repeat[i]=oracle.mean(coords[i])
    err=float(np.max(abs(repeat-means))/np.max(abs(means)));assert err<2e-12
    replay=oracle.mean(base);assert np.array_equal(replay,saved['fixed'])
    final=oracle.close()
    for k in ['band','MJD','data_flux','data_fluxerr','C','W']:
        assert np.array_equal(final[k],prior[k]),('reference state after oracle',k)
    np.savez_compressed(G/'design-replay.npz',coordinates=coords,means=means,permuted_replay=repeat,order=order)
    result=dict(status='PASS: exact native mean and state replay',initialization_seconds=oracle.initialization,N_native_vectors=oracle.counter,
      separate_process_relative_errors=dif,reordered_replay_relative_error=err,nominal_replay_exact=True,reference_data_and_C_exact=True,
      benchmark100_seconds=elapsed,seconds_per_vector=elapsed/100,projected100000_seconds=elapsed*1000,protocol_sha256=p.sha(G/'protocol.json'))
    p.save(G/'result.json',result);p.save(G/'manifest.json',dict(files_sha256={str(f.relative_to(G)):p.sha(f) for f in G.rglob('*') if f.is_file() and f.name!='manifest.json'}))
    print(json.dumps(result,indent=2))

if __name__=='__main__':gate()
