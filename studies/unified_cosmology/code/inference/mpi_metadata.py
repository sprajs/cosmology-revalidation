"""Keep common Cobaya metadata writes on rank zero, preserving collectives.

Cobaya 3.6.2 calls Output.check_and_dump_info on every MPI rank. Different
initial-reference strings can race when writing its common YAML/dill files.
Only those three output streams are suppressed on nonzero ranks; chain files,
likelihood inputs and the method's MPI communication are untouched.
"""
import builtins
import io
from pathlib import Path


def install(prefix):
    from cobaya import mpi, output
    shared={Path(str(prefix)+suffix).resolve() for suffix in
            ['.input.yaml','.updated.yaml','.updated.dill_pickle']}
    original=getattr(output,'open',builtins.open)

    def guarded_open(file,mode='r',*args,**kwargs):
        if not mpi.is_main_process() and mode.startswith('w'):
            try:
                suppress=Path(file).resolve() in shared
            except TypeError:
                suppress=False
            if suppress:
                return io.BytesIO() if 'b' in mode else io.StringIO()
        return original(file,mode,*args,**kwargs)

    output.open=guarded_open
    return {'common_metadata_writer_rank':0,
            'suppressed_paths':[str(p) for p in sorted(shared)]}


def validate():
    import hashlib
    import json
    from mpi4py import MPI
    from cobaya.output import get_output
    from cobaya.yaml import yaml_load_file
    root=Path(__file__).resolve().parents[4]
    folder=root/'.work/unified-cosmology/inference/metadata-validation'
    prefix=folder/'guarded'
    policy=install(prefix)
    with get_output(prefix=str(prefix),force=True) as out:
        for iteration in range(12):
            config={'params':{'x':{'prior':{'min':-1.,'max':1.},'ref':MPI.COMM_WORLD.rank*.001}},
                    'likelihood':{'placeholder':None},
                    'marker':str(MPI.COMM_WORLD.rank)*(101+iteration*37),
                    'output':str(prefix)}
            out.check_and_dump_info(config,config,check_compatible=False)
            MPI.COMM_WORLD.Barrier()
            if MPI.COMM_WORLD.rank==0:
                for suffix in ['input','updated']:
                    saved=yaml_load_file(str(prefix)+'.'+suffix+'.yaml')
                    assert saved['marker']=='0'*(101+iteration*37)
                    assert saved['params']['x']['ref']==0.
            MPI.COMM_WORLD.Barrier()
    if MPI.COMM_WORLD.rank==0:
        result={'status':'passed','MPI_ranks':MPI.COMM_WORLD.size,
                'rounds':12,'verified_common_yaml_writes':24,'policy':policy,
                'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'scope':'Metadata output only; no sampling density or chain collection changed.'}
        target=root/'studies/unified_cosmology/results/inference/metadata-validation.json'
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))


if __name__=='__main__':validate()
