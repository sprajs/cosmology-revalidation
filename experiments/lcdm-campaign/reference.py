"""Reference only: direct momentum GL panels, direct z and sqrt(a), scalar LDLT.

The numerical route shares ancestry with the original thermal peer reference.
SI constants, physical equations and explicit species conventions are shared
ancestry. The high-precision map starts with exact input binary64 coordinates;
production additionally rounds mapped fractions/temperatures to binary64.
No native result is an input. This admits only the frozen DESI control, and no
reference result is observational evidence or production physics.
"""
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import mpmath as mp
from mpmath.libmp.backend import BACKEND


REQUEST_SHA256 = "0e6fca28ec07187da6bcc91c0ae917d0b3caad7d16f91ceae9721eab2f937d15"
METHOD = "reference/direct-momentum-GL-direct-z-sqrt-a-scalar-LDLT/v1"
MAPPING = "physical-omega-h2-explicit-Kelvin-temperature-state-weight-SI/v1"
CONSTANTS = "SI2019-exact-h-c-kB-eV-IAU2012-AU-CODATA2018-G-fixed"
MODEL_ORDER = ["thermal-massless-supplied-drag", "thermal-massive-FD-supplied-drag"]
QUERY_ORDER = [f"DESI-DR2-row-{i}" for i in range(13)]
SETTINGS = {
    "coarse": dict(digits=60, momentum_nodes_per_panel=32, outer_nodes_per_panel=32, momentum_tail=128),
    "fine": dict(digits=90, momentum_nodes_per_panel=48, outer_nodes_per_panel=48, momentum_tail=160),
    # Cross refinements isolate each axis at common 90-digit precision.
    "momentum_refined": dict(digits=90, momentum_nodes_per_panel=48, outer_nodes_per_panel=32, momentum_tail=160),
    "outer_refined": dict(digits=90, momentum_nodes_per_panel=32, outer_nodes_per_panel=48, momentum_tail=128),
}
SOURCES = [
    {"path": "data/bao/desi_gaussian_bao_ALL_GCcomb_mean.txt", "bytes": 472,
     "sha256": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585",
     "source": "https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt"},
    {"path": "data/bao/desi_gaussian_bao_ALL_GCcomb_cov.txt", "bytes": 2547,
     "sha256": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509",
     "source": "https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt"},
]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _keys(value, expected, name):
    if type(value) is not dict or set(value) != set(expected):
        raise ValueError(name + " fields differ")


def _number(value, name, positive=False, nonnegative=False):
    if type(value) not in (int, float):
        raise ValueError(name + " requires a binary64 number")
    try:
        value = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(name + " is not binary64") from error
    if not math.isfinite(value) or (value != 0 and abs(value) < sys.float_info.min):
        raise ValueError(name + " is not finite normal binary64 or zero")
    if (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(name + " is outside the nonnegative reference domain")
    return value


def _integer(value, name, minimum=0, maximum=1_000_000_000):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(name + " requires a bounded integer")


def _duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON field " + key)
        result[key] = value
    return result


def load_json(path):
    path = Path(path)
    if path.stat().st_size > 65536:
        raise ValueError("reference input quota")
    blob = path.read_bytes()
    if len(blob) > 65536:
        raise ValueError("reference input quota")
    def nonfinite(value):
        raise ValueError("nonfinite JSON constant " + value)
    return json.loads(blob, object_pairs_hook=_duplicates, parse_constant=nonfinite)


def validate_model(model):
    _keys(model, {"id", "H0", "omega_b", "omega_c", "Tcmb", "omega_other", "species", "z_drag"}, "model")
    if type(model["id"]) is not str or model["id"] not in MODEL_ORDER:
        raise ValueError("unreviewed reference model identity")
    for key in ("H0", "Tcmb"):
        _number(model[key], key, positive=True)
    for key in ("omega_b", "omega_c", "omega_other", "z_drag"):
        _number(model[key], key, nonnegative=True)
    if type(model["species"]) is not list or len(model["species"]) > 1:
        raise ValueError("reference species quota")
    for species in model["species"]:
        _keys(species, {"mass_eV", "temperature_K", "g"}, "species")
        _number(species["mass_eV"], "mass_eV", nonnegative=True)
        _number(species["temperature_K"], "temperature_K", positive=True)
        _number(species["g"], "g", positive=True)


def validate_request(request, request_path):
    _keys(request, {"schema_version", "id", "operation", "engine", "models", "model_role", "budgets", "producer_policy", "resources"}, "request")
    _integer(request["schema_version"], "request schema", 1, 1)
    if request["id"] != "lcdm-campaign/thermal-desi-control-v1" or request["operation"] != "lcdm-campaign.thermal-control":
        raise ValueError("unreviewed reference request identity")
    _keys(request["engine"], {"revision", "build_id", "manifest_sha256", "archive_sha256", "cli_sha256", "route"}, "engine")
    if not all(type(value) is str for value in request["engine"].values()) or type(request["model_role"]) is not str:
        raise ValueError("reference request identity requires strings")
    if type(request["models"]) is not list or len(request["models"]) != 2:
        raise ValueError("exact two reference models required")
    for model in request["models"]:
        validate_model(model)
    if [model["id"] for model in request["models"]] != MODEL_ORDER:
        raise ValueError("reference model order differs")
    _keys(request["budgets"], {"ratio_absolute", "ratio_relative", "distance_absolute_mpc", "distance_relative", "density_absolute", "projection_absolute", "reference_fraction"}, "budgets")
    for name, value in request["budgets"].items():
        _number(value, name, positive=True)
    _keys(request["producer_policy"], {"absolute_mpc", "relative", "ratio_absolute", "ratio_relative", "moment_absolute", "moment_relative", "maximum_callbacks_per_point", "maximum_callbacks_per_model", "maximum_callbacks_total"}, "producer policy")
    for name, value in request["producer_policy"].items():
        if name.startswith("maximum_"):
            _integer(value, name, 1)
        else:
            _number(value, name, positive=True)
    _keys(request["resources"], {"jobs", "threads", "compile_seconds", "native_seconds", "reference_seconds", "memory_bytes", "output_bytes"}, "resources")
    for name, value in request["resources"].items():
        _integer(value, name, 1, 4_000_000_000)
    pinned_request = load_json(request_path)
    if sha256(request_path) != REQUEST_SHA256:
        raise ValueError("reference adjacent request byte identity differs")
    if canonical(request) != canonical(pinned_request):
        raise ValueError("reference supplied request differs from pinned snapshot")


def validate_axes(rows, covariance):
    if type(rows) is not list or len(rows) != 13:
        raise ValueError("full ordered13 reference queries required")
    for row in rows:
        _keys(row, {"z", "observed", "kind"}, "query")
        _number(row["z"], "query z", nonnegative=True)
        _number(row["observed"], "query observed", positive=True)
        _integer(row["kind"], "query kind", 0, 2)
        if row["z"] > 10:
            raise ValueError("query exceeds bounded DESI reference domain")
    if type(covariance) is not list or len(covariance) != 13 or any(type(row) is not list or len(row) != 13 for row in covariance):
        raise ValueError("full13 covariance dimensions differ")
    for row in covariance:
        for value in row:
            _number(value, "covariance")
    for i in range(13):
        if covariance[i][i] <= 0:
            raise ValueError("reference covariance diagonal not positive")
        for j in range(i):
            if float(covariance[i][j]).hex() != float(covariance[j][i]).hex():
                raise ValueError("reference covariance is not exactly symmetric binary64")


def validate_input(data, request_path=None, source_dir=None):
    _keys(data, {"schema_version", "request", "request_sha256", "rows", "covariance", "sources"}, "reference input")
    _integer(data["schema_version"], "reference schema", 1, 1)
    if data["request_sha256"] != REQUEST_SHA256:
        raise ValueError("reference request SHA differs")
    validate_request(data["request"], request_path or Path(__file__).with_name("request.json"))
    validate_axes(data["rows"], data["covariance"])
    if canonical(data["sources"]) != canonical(SOURCES):
        raise ValueError("exact ordered source identities differ")
    if source_dir is not None:
        blobs = []
        for i, source in enumerate(SOURCES):
            path = Path(source_dir) / f"original-input-{i}.txt"
            if path.stat().st_size != source["bytes"]:
                raise ValueError("reference original source bytes differ")
            blob = path.read_bytes()
            if len(blob) != source["bytes"] or hashlib.sha256(blob).hexdigest() != source["sha256"]:
                raise ValueError("reference original source bytes differ")
            blobs.append(blob)
        kinds = {"DM_over_rs": 0, "DH_over_rs": 1, "DV_over_rs": 2}
        rows = []
        for line in blobs[0].decode("utf-8").splitlines():
            if line.strip() and not line.lstrip().startswith("#"):
                z, observed, kind = line.split()
                rows.append(dict(z=float(z), observed=float(observed), kind=kinds[kind]))
        covariance = [list(map(float, line.split())) for line in blobs[1].decode("utf-8").splitlines()
                      if line.strip() and not line.lstrip().startswith("#")]
        if canonical(rows) != canonical(data["rows"]) or canonical(covariance) != canonical(data["covariance"]):
            raise ValueError("reference query/covariance order differs from original source bytes")


def runtime_fingerprint():
    if mp.__version__ != "1.3.0" or BACKEND != "python":
        raise ValueError("reference requires reviewed mpmath1.3.0 pure Python backend")
    root = Path(mp.__file__).resolve().parent
    paths = sorted(root.rglob("*.py"))
    if not paths or len(paths) > 256:
        raise ValueError("reference module inventory quota")
    inventory, total = [], 0
    for path in paths:
        if path.is_symlink() or not path.is_file():
            raise ValueError("reference module inventory requires regular sources")
        if total + path.stat().st_size > 8 * 1024 * 1024:
            raise ValueError("reference module source byte quota")
        blob = path.read_bytes()
        total += len(blob)
        if total > 8 * 1024 * 1024:
            raise ValueError("reference module source byte quota")
        inventory.append(dict(path=str(path), relative_path=str(path.relative_to(root)),
                              bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest()))
    executable = Path(sys.executable).absolute()
    return {
        "python": {"version": sys.version, "implementation": platform.python_implementation(),
                   "executable": {"path": str(executable), "resolved_path": str(executable.resolve()),
                                  "bytes": executable.stat().st_size, "sha256": sha256(executable)}},
        "mpmath": {"version": mp.__version__, "backend": BACKEND, "module_root": str(root),
                   "source_inventory": inventory, "source_inventory_sha256": hashlib.sha256(canonical(inventory)).hexdigest()},
    }


def binary64(value):
    """The exact binary64 input, with one high-precision rational division."""
    numerator, denominator = float(value).as_integer_ratio()
    return mp.mpf(numerator) / denominator


def legendre_rule(order, digits):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    tolerance = mp.mpf(10) ** (10 - digits)
    if len(nodes) != order or len(weights) != order:
        raise ValueError("reference Legendre dimensions differ")
    if any(not -1 < x < 1 for x in nodes) or any(w <= 0 for w in weights):
        raise ValueError("reference Legendre domain differs")
    if any(nodes[i] >= nodes[i + 1] for i in range(order - 1)):
        raise ValueError("reference Legendre root order differs")
    if abs(mp.fsum(weights) - 2) > tolerance:
        raise ValueError("reference Legendre weight normalization failed")
    for i, (x, weight) in enumerate(zip(nodes, weights)):
        if abs(mp.legendre(order, x)) > tolerance or abs(x + nodes[order - i - 1]) > tolerance or abs(weight - weights[order - i - 1]) > tolerance:
            raise ValueError("reference Legendre root/symmetry check failed")
    # Symmetry supplies odd moments; validate every nonconstant even moment
    # through the GL degree 2*n-1, including the high-degree weight constraint.
    for power in range(2, 2 * order, 2):
        if abs(mp.fsum(weight * node ** power for node, weight in zip(nodes, weights)) - mp.mpf(2) / (power + 1)) > tolerance:
            raise ValueError("reference Legendre polynomial exactness failed")
    return nodes, weights


def ldlt(covariance):
    if type(covariance) is not list or not 1 <= len(covariance) <= 13:
        raise ValueError("reference covariance dimensions differ")
    n = len(covariance)
    if any(type(row) is not list or len(row) != n for row in covariance):
        raise ValueError("reference covariance dimensions differ")
    for i in range(n):
        for j in range(n):
            _number(covariance[i][j], "covariance")
            if float(covariance[i][j]).hex() != float(covariance[j][i]).hex():
                raise ValueError("reference covariance is not exactly symmetric binary64")
    lower = [[mp.mpf(0)] * n for _ in range(n)]
    diagonal = []
    for j in range(n):
        lower[j][j] = mp.mpf(1)
        pivot = binary64(covariance[j][j]) - mp.fsum(lower[j][k] ** 2 * diagonal[k] for k in range(j))
        if not mp.isfinite(pivot) or pivot <= 0:
            raise ValueError("reference covariance not SPD")
        diagonal.append(pivot)
        for i in range(j + 1, n):
            lower[i][j] = (binary64(covariance[i][j]) - mp.fsum(lower[i][k] * lower[j][k] * diagonal[k] for k in range(j))) / pivot
    return lower, diagonal


def solve_ldlt(lower, diagonal, right):
    n = len(right)
    forward = []
    for i in range(n):
        forward.append(right[i] - mp.fsum(lower[i][k] * forward[k] for k in range(i)))
    solution = [mp.mpf(0)] * n
    for i in reversed(range(n)):
        solution[i] = forward[i] / diagonal[i] - mp.fsum(lower[k][i] * solution[k] for k in range(i + 1, n))
    return solution


def exponential_tail(order, endpoint):
    # Integral_endpoint^infinity q^order exp(-q) dq without subtraction.
    polynomial = mp.mpf(1)
    power = mp.mpf(1)
    for k in range(1, order + 1):
        power *= endpoint
        polynomial = power + k * polynomial
    return mp.exp(-endpoint) * polynomial


def evaluate(model, rows, covariance, digits, nq, no, tail):
    validate_model(model)
    validate_axes(rows, covariance)
    if (digits, nq, no, tail) not in {(60, 32, 32, 128), (90, 48, 48, 160), (90, 48, 32, 160), (90, 32, 48, 128)}:
        raise ValueError('unreviewed reference computational settings')
    with mp.workdps(digits):
        D=binary64
        H=D(model['H0']); h2=(H/100)**2
        c=mp.mpf(299792458); ev=mp.mpf('1.602176634e-19')
        kb=mp.mpf('1.380649e-23')/ev
        mpc=mp.mpf(648000000000)*149597870700/mp.pi
        hbar=mp.mpf('6.62607015e-34')/(2*mp.pi)
        critical=3*(H*1000/mpc)**2*c*c/(8*mp.pi*mp.mpf('6.67430e-11'))/ev*(hbar*c/ev)**3
        photons=mp.pi**2/15*(kb*D(model['Tcmb']))**4/critical
        radiation=D(model['omega_other'])/h2
        baryon=D(model['omega_b'])/h2; matter=baryon+D(model['omega_c'])/h2
        if not mp.isfinite(critical) or critical <= 0 or photons <= 0 or photons+radiation+matter > 1:
            raise ValueError('reference physical map outside flat nonnegative domain')
        xx,ww=legendre_rule(nq,digits)
        cuts=list(map(mp.mpf,['0','.01','.1','1','4']))+list(map(mp.mpf,range(8,tail+1,4)))
        if cuts[0] != 0 or cuts[-1] != tail:
            raise ValueError('reference momentum endpoint differs')
        nodes=[]
        for lo,hi in zip(cuts,cuts[1:]):
            for x,w in zip(xx,ww):
                q=(lo+hi)/2+(hi-lo)*x/2
                nodes.append((q*q,(hi-lo)*w*q*q/(2*(1+mp.exp(q)))))
        species=[(D(s['mass_eV'])/(kb*D(s['temperature_K'])),D(s['g'])*(kb*D(s['temperature_K']))**4/(2*mp.pi**2*critical)) for s in model['species']]
        def relic(a):
            return mp.fsum(pre*mp.fsum(w*mp.sqrt(q2+(ma*a)**2) for q2,w in nodes) for ma,pre in species)
        # For a<=1, sqrt(q^2+(ma)^2)<=q+m and occupation<=exp(-q).
        # Flat closure also changes with the tail: its scaled P correction is
        # nonnegative and bounded by this uniform envelope throughout [0,1].
        tail_scaled=mp.fsum(pre*(exponential_tail(3,tail)+ma*exponential_tail(2,tail)) for ma,pre in species)
        lam=1-photons-radiation-matter-relic(mp.mpf(1))
        if not mp.isfinite(lam) or lam < tail_scaled:
            raise ValueError('reference flat closure lacks nonnegative Lambda/tail margin')
        def P(a):
            if not 0 <= a <= 1:raise ValueError('reference scale factor outside [0,1]')
            value=photons+radiation+matter*a+relic(a)+lam*a**4
            if not mp.isfinite(value) or value <= 0:raise ValueError('reference scaled expansion not positive finite')
            return value
        def E(z):
            a=1/(1+z);return mp.sqrt(P(a))/a**2
        x,w=legendre_rule(no,digits)
        def integral(f,cuts):
            return mp.fsum((hi-lo)/2*mp.fsum(wi*f((lo+hi)/2+(hi-lo)*xi/2) for xi,wi in zip(x,w)) for lo,hi in zip(cuts,cuts[1:]))
        end=1/mp.sqrt(1+D(model['z_drag'])); loading=3*baryon/(4*photons)
        # P(0)>0: the sqrt(a) integrand tends to zero. Never evaluate E(a=0).
        ruler=c/1000/H/mp.sqrt(3)*integral(lambda u:2*u/mp.sqrt(P(u*u)*(1+loading*u*u)),[0,end/8,end/4,end/2,end])
        if not mp.isfinite(ruler) or ruler <= 0:raise ValueError('reference ruler not positive finite')
        predictions=[]; distances={}
        for row in rows:
            z=D(row['z']); key=float(row['z']).hex()
            if key not in distances:
                dm=c/1000/H*integral(lambda t:1/E(t),[0,z/4,z/2,z])
                dh=c/1000/(H*E(z));dv=mp.root(dm*dm*z*dh,3)
                distances[key]=(dm,dh,dv)
            predictions.append(distances[key][row['kind']]/ruler)
        # Scalar dense LDLT differs from production Cholesky. Full ordered
        # covariance is retained, with exact symmetry checked before factoring.
        n=len(rows);lower,diagonal=ldlt(covariance)
        residual=[D(r['observed'])-p for r,p in zip(rows,predictions)]
        inverse_residual=solve_ldlt(lower,diagonal,residual)
        quadratic=mp.fsum(r*v for r,v in zip(residual,inverse_residual));ld=mp.fsum(mp.log(d) for d in diagonal);norm=n*mp.log(2*mp.pi)
        # This bounds truncation only. GL error is measured by cross refinements.
        tail_relative=tail_scaled/(photons+radiation)
        ratio_relative=tail_relative/(mp.sqrt(1+tail_relative)+1)
        prediction_bounds=[abs(p)*ratio_relative for p in predictions]
        inverse_columns=[solve_ldlt(lower,diagonal,[mp.mpf(int(i==j)) for i in range(n)]) for j in range(n)]
        quadratic_bound=2*mp.fsum(abs(v)*b for v,b in zip(inverse_residual,prediction_bounds))
        quadratic_bound+=mp.fsum(abs(inverse_columns[j][i])*prediction_bounds[i]*prediction_bounds[j] for i in range(n) for j in range(n))
        strings=lambda vs:[mp.nstr(v,digits) for v in vs]
        if not all(mp.isfinite(v) for v in [ruler,*predictions,quadratic,ld,norm,tail_scaled,tail_relative,quadratic_bound]):
            raise ValueError('reference result not finite')
        return {'model':model['id'],'ruler_mpc':mp.nstr(ruler,digits),'predictions':strings(predictions),'density':dict(zip(['quadratic','log_determinant','normalization','log_density'],strings([quadratic,ld,norm,-(quadratic+ld+norm)/2]))),'digits':digits,'momentum_nodes_per_panel':nq,'outer_nodes_per_panel':no,'momentum_tail':tail,
                'tail_bounds':{'scaled_expansion':mp.nstr(tail_scaled,digits),'relative_scaled_expansion':mp.nstr(tail_relative,digits),'predictions':strings(prediction_bounds),'quadratic':mp.nstr(quadratic_bound,digits),'log_density':mp.nstr(quadratic_bound/2,digits)}}


def main(arguments):
    if arguments == ['--runtime-fingerprint']:
        print(json.dumps(runtime_fingerprint(),allow_nan=False))
        return
    if len(arguments) != 1:raise ValueError('one bounded reference input or --runtime-fingerprint required')
    p=Path(arguments[0]).resolve();d=load_json(p)
    validate_input(d,p.with_name('request.json'),p.parent)
    input_hash,script_hash=sha256(p),sha256(__file__);before=runtime_fingerprint()
    out={'schema_version':1,'method':METHOD,'mapping':MAPPING,'constants':CONSTANTS,
         'mpmath_version':mp.__version__,'runtime_before':before,'model_order':MODEL_ORDER,'models':d['request']['models'],
         'query_order':QUERY_ORDER,'queries':d['rows'],'sources':d['sources'],'request_sha256':REQUEST_SHA256,
         'request_canonical_sha256':hashlib.sha256(canonical(d['request'])).hexdigest(),
         'reference_input_sha256':input_hash,'reference_script_sha256':script_hash,'settings':SETTINGS}
    for key,settings in SETTINGS.items():
        out[key]=[evaluate(model,d['rows'],d['covariance'],settings['digits'],settings['momentum_nodes_per_panel'],settings['outer_nodes_per_panel'],settings['momentum_tail']) for model in d['request']['models']]
    out['runtime_after']=runtime_fingerprint()
    validate_input(load_json(p),p.with_name('request.json'),p.parent)
    if before != out['runtime_after'] or sha256(p) != input_hash or sha256(__file__) != script_hash:
        raise ValueError('reference source/runtime/input identity changed during execution')
    print(json.dumps(out,allow_nan=False))


if __name__=='__main__':
    main(sys.argv[1:])
