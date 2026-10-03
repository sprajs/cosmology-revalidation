"""SN-specific full-n fixtures/truth; import performs no numerical work."""
from fractions import Fraction as F
import transport as t

NATIVE_CONTROLS = ("analytic", "conditioning", "indefinite", "nonsymmetric", "nonfinite", "wrong-order", "missing-row", "original-B-refusal")
REFERENCE_CONTROLS = ("analytic", "conditioning")

def analytic_truth():
    v = [F(1 if i % 2 == 0 else -1, 4) for i in range(t.N)]
    r = [F(i % 17 - 8, 16) for i in range(t.N)]
    one = [F(1)] * t.N
    def inverse(x):
        vx = sum(a*b for a,b in zip(v,x))
        denom = 2 + sum(a*a for a in v)
        return [(a-b*vx/denom)/2 for a,b in zip(x,v)]
    wx = inverse(one)
    g = sum(wx)
    rows = []
    for shift in (F(0), F(1), F(-1,2), F(0)):
        residual = [x+shift for x in r]
        m = sum(a*b for a,b in zip(wx,residual))/g
        adjusted = [x-m for x in residual]
        solved = inverse(adjusted)
        q = sum(a*b for a,b in zip(adjusted,solved))
        rows.append({"M": m, "q": q, "score": -q/2, "adjusted": adjusted})
    return rows

def synthetic(name):
    t.need(name in REFERENCE_CONTROLS, "closed reference control")
    if name == "analytic":
        c = [[(2.0 if i==j else 0.0)+(0.0625 if i%2==j%2 else -0.0625) for j in range(t.N)] for i in range(t.N)]
    else:
        c = [[(1.0 if i==j else 0.0) for j in range(t.N)] for i in range(t.N)]
        c[-1][-1] = 2.0**-48
    r = [(i%17-8)/16 for i in range(t.N)]
    return c, [[x+shift for x in r] for shift in (0,1,-0.5,0)]

def analytic_admission(native, reference):
    from decimal import Decimal, localcontext
    t.need(native["status"] == "accepted" and reference["status"] == "accepted", "full-n analytic route refusal")
    truth = analytic_truth()
    with localcontext() as ctx:
        ctx.prec = 96
        def exact(f):
            return Decimal(f.numerator)/Decimal(f.denominator)
        for i, expected in enumerate(truth):
            # Every adjusted lane is checked, not just a selected coordinate.
            p, q = native["rows"][i]["payload"], reference["rows"][i]
            for value, wanted in ((p["M_coefficient"],expected["M"]),(q["M_coefficient"],expected["M"]),
                                  (p["relative_profile_score"],expected["score"]),(q["relative_profile_score"],expected["score"])):
                t.need(abs(Decimal(str(value))-exact(wanted)) <= Decimal("1e-6"), "full-n analytic scalar")
            t.need(len(p["adjusted_residuals"]) == len(q["adjusted_residuals"]) == t.N, "analytic no-drop")
            for a,b,wanted in zip(p["adjusted_residuals"],q["adjusted_residuals"],expected["adjusted"]):
                t.need(abs(Decimal(str(a))-exact(wanted)) <= Decimal("1e-6") and abs(Decimal(str(b))-exact(wanted)) <= Decimal("1e-6"), "analytic lane/shift")
        for rows, key in ((native["rows"],"payload"),(reference["rows"],None)):
            scalar = [(r[key] if key else r)["relative_profile_score"] for r in rows]
            t.need(all(abs(Decimal(str(v))-Decimal(str(scalar[0]))) <= Decimal("1e-6") for v in scalar), "constant-shift profile invariance")
    return {"full_n": t.N, "Sherman_Morrison": True, "all_adjusted_lanes_checked": True, "constant_shift": True,
            "scalar_and_lane_absolute_screen": "1e-6", "physical_qualification": False}

def asymmetry_algebra_control():
    # Exact two-lane adversary; no raw-B inversion or numerical score route.
    a,u,v = F(1,2),F(2),F(-1)
    qs = (u-v)**2/2
    qp = qs/(1+a*a)
    t.need(qp != qs and qp == qs/(1+(-a)**2), "skew-sign adversary")
    t.need(qs/(1+F(0)**2)==qs,"zero-skew limit")
    tiny_skew,tiny_lambda=F(1,2**40),F(1,2**80)
    weighted=tiny_skew/tiny_lambda
    qs_tiny=(u-v)**2/(2*tiny_lambda)
    qp_tiny=qs_tiny/(1+weighted*weighted)
    # S0=eI, A01=a: precision shrinks by1/(1+(a/e)^2).
    # This exact contrasting score remains sign-invariant and recovers zero skew.
    t.need(weighted>1 and qp_tiny<qs_tiny/2 and qp_tiny==qs_tiny/(1+(-weighted)**2) and qs_tiny/(1+F(0)**2)==qs_tiny,"tiny-skew exact profile contrast")
    return {"S0_q": str(qs), "effective_inverse_quadratic_q": str(qp), "different_targets": True,
            "reversed_skew_invariant": True, "zero_skew_limit":True,
            "near_singular_absolute_smallness_adversary":str(weighted),"near_singular_S0_relative_profile":str(-qs_tiny/2),
            "near_singular_effective_relative_profile":str(-qp_tiny/2),"raw_B_score_executed": False}

def projection_controls():
    import math
    import assemble
    lower=1.0;odd=math.nextafter(lower,math.inf);even=math.nextafter(odd,math.inf)
    t.need(assemble.rn_even((F.from_float(lower)+F.from_float(odd))/2).hex()==lower.hex(),"lower-even midpoint tie")
    t.need(assemble.rn_even((F.from_float(odd)+F.from_float(even))/2).hex()==even.hex(),"upper-even midpoint tie")
    a,b,exact=assemble.pair_projection(-0.0,0.0)
    t.need(a.hex()==(-0.0).hex() and b.hex()==(0.0).hex() and exact is None,"equal signed-zero pair retained")
    equal=[1.0,-0.0,0.0,2.0]
    t.need(list(assemble.projection_updates(equal,2))==[] and [x.hex() for x in equal]==[(1.0).hex(),(-0.0).hex(),(0.0).hex(),(2.0).hex()],"actual projection traversal preserves equal cells/diagonals")
    diagonal=(1.0).hex()
    source=[[1.0,3.0],[1.0,1.0]]
    updates=list(assemble.projection_updates([x for row in source for x in row],2))
    t.need(len(updates)==1 and updates[0][0:2]==(1,0),"one projected unordered pair only")
    _,_,_,_,a,exact=updates[0];b=a
    determinant=F.from_float(source[0][0])*F.from_float(source[1][1])-exact*exact
    t.need(a==b==2.0 and determinant==-3 and source[0][0].hex()==source[1][1].hex()==diagonal,"positive-diagonal projected non-SPD adversary")
    return {"both_tie_parities":True,"equal_bits_and_diagonal":True,"positive_diagonal_projection_determinant":str(determinant),"native_full_n_indefinite_control_uses_this_block":True}

def transport_fault_controls():
    # Bounded mocks exercise real read/failure serialization; no source/table IO.
    import hashlib
    import stat
    from types import SimpleNamespace
    from unittest import mock
    from pathlib import Path
    def metadata(size=4,stamp=1):
        return SimpleNamespace(st_dev=1,st_ino=2,st_size=size,st_mode=stat.S_IFREG|0o444,st_mtime_ns=stamp,st_ctime_ns=stamp)
    pin={"path":"/fixture/source","bytes":4,"sha256":hashlib.sha256(b"abcd").hexdigest()}
    cases=(
      ("changed-byte",[b"abce",b""],[metadata()]*3,[metadata()]*3,4,True,hashlib.sha256(b"abce").hexdigest()),
      ("stat-only",[b"abcd",b""],[metadata(),metadata(stamp=2),metadata(stamp=2)],[metadata(),metadata(stamp=2),metadata(stamp=2)],4,True,pin["sha256"]),
      ("unlink",[b"abcd",b""],[metadata()]*3,[metadata(),FileNotFoundError("fixture unlink"),FileNotFoundError("fixture unlink")],4,True,pin["sha256"]),
      ("partial-read",[b"ab",OSError("fixture read interruption")],[metadata()]*3,[metadata()]*3,2,False,hashlib.sha256(b"ab").hexdigest()),
      ("early-size",[],[metadata(size=5)]*3,[metadata(size=5)]*3,0,False,hashlib.sha256(b"").hexdigest()),
      ("overflow",[b"abcdx"],[metadata()]*3,[metadata()]*3,5,False,hashlib.sha256(b"abcdx").hexdigest())
    )
    results=[]
    for name,chunks,descriptor,linked,count,eof,digest in cases:
        with mock.patch.object(Path,"is_dir",return_value=True),mock.patch.object(Path,"is_symlink",return_value=False),mock.patch.object(Path,"lstat",side_effect=linked),mock.patch.object(t.os,"open",return_value=77),mock.patch.object(t.os,"close"),mock.patch.object(t.os,"fstat",side_effect=descriptor),mock.patch.object(t.os,"read",side_effect=chunks):
            try:t.read(pin);t.need(False,"expected source refusal")
            except t.SourceIdentityError as exc:
                record=t.document(t.encoded(t.failure(exc)))
                e=record["source_identity"]
                t.need(e["expected_authority"]==pin and e["bytes_consumed"]==count and e["eof_observed"] is eof and e["consumed_sha256"]==digest and e["accepted_complete_identity"] is False,"failure consumed identity/serialization")
                results.append({"id":name,"record":record})
    with mock.patch.object(Path,"is_dir",return_value=True),mock.patch.object(Path,"is_symlink",return_value=False),mock.patch.object(Path,"lstat",return_value=metadata()),mock.patch.object(t.os,"open",side_effect=OSError("fixture unopened")):
        try:t.read(pin);t.need(False,"expected unopened refusal")
        except t.SourceIdentityError as exc:
            e=exc.evidence;t.need(e["descriptor_opened"] is False and e["bytes_consumed"] is e["consumed_sha256"] is e["eof_observed"] is None,"unopened earns no stream")
            results.append({"id":"unopened","record":t.failure(exc)})
    return {"synthetic_nonqualification":True,"records":results}
