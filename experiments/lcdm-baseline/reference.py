"""Reference-only algorithms; never production physics or an inference engine.

60-digit Decimal Legendre P8 Newton roots/weights (derived at runtime), composite
GL8 directly in redshift and u=sqrt(a), dense scalar LDLT. Production instead
uses adaptive log-redshift/scale-factor integration and Cholesky. No external
reference library or copied constants. Inputs represent actual binary64 values.
"""
from decimal import Decimal as D, localcontext


def decimal(x):
    return D.from_float(float(x))


def legendre(x):
    previous, value = D(1), x
    for n in range(2, 9):
        previous, value = value, ((2*n-1)*x*value-(n-1)*previous)/n
    return value, 8*(x*value-previous)/(x*x-1)


def nodes():
    # Rational initial guesses isolate the four positive roots of P8.
    result = []
    for guess in ("0.18", "0.53", "0.80", "0.96"):
        x = D(guess)
        for _ in range(40):
            p, derivative = legendre(x)
            new = x-p/derivative
            if abs(new-x) < D("1e-55"):
                x = new
                break
            x = new
        p, derivative = legendre(x)
        if abs(p) > D("1e-52"):
            raise ValueError("Reference Legendre root failed")
        result.append((x, 2/((1-x*x)*derivative*derivative)))
    return result


def quadrature(function, end, panels, rule):
    h = end/(2*panels)
    return sum((h*w*(function((2*j+1)*h-h*x)+function((2*j+1)*h+h*x))
                for j in range(panels) for x, w in rule), D(0))


def pi():
    def atan(x):
        total, power = D(0), x
        for n in range(1000):
            term = power/(2*n+1)
            total += term if n % 2 == 0 else -term
            if abs(term) < D("1e-58"):
                return total
            power *= x*x
        raise ValueError("Reference pi failed")
    return 16*atan(D(1)/5)-4*atan(D(1)/239)


def ldlt(covariance, residual):
    n = len(residual)
    lower = [[D(0)]*n for _ in range(n)]
    diagonal = []
    for j in range(n):
        lower[j][j] = D(1)
        pivot = covariance[j][j]-sum((lower[j][k]**2*diagonal[k] for k in range(j)), D(0))
        if pivot <= 0:
            raise ValueError("Reference covariance not SPD")
        diagonal.append(pivot)
        for i in range(j+1, n):
            lower[i][j] = (covariance[i][j]-sum((lower[i][k]*lower[j][k]*diagonal[k]
                                               for k in range(j)), D(0)))/pivot
    transformed = []
    for i in range(n):
        transformed.append(residual[i]-sum((lower[i][k]*transformed[k] for k in range(i)), D(0)))
    quadratic = sum((x*x/p for x, p in zip(transformed, diagonal)), D(0))
    logdet = sum((p.ln() for p in diagonal), D(0))
    normalization = n*(2*pi()).ln()
    return dict(quadratic=quadratic, log_determinant=logdet, normalization=normalization,
                log_density=-(quadratic+logdet+normalization)/2)


def evaluate(request, observed, covariance, panels):
    with localcontext() as context:
        context.prec = 60
        m = {k: decimal(v) for k, v in request["model"].items() if k != "drag_origin"}
        matter, radiation = m["omega_m"], m["omega_r"]
        lam = 1-matter-radiation
        loading = 3*m["omega_b"]/(4*m["omega_gamma"])
        scale = D("299792.458")/m["h0_km_s_mpc"]
        rule = nodes()
        def expansion(z):
            y = 1+z
            return (radiation*y**4+matter*y**3+lam).sqrt()
        def distance(z):
            return scale*quadrature(lambda t: 1/expansion(t), z, panels, rule)
        def sound(u):
            a = u*u
            return 2*u/(radiation+matter*a+lam*a**4).sqrt()/(1+loading*a).sqrt()
        ruler = scale/D(3).sqrt()*quadrature(sound, 1/(1+m["z_drag"]).sqrt(), panels, rule)
        background = []
        for z in request["redshifts"]:
            z = decimal(z)
            dm = distance(z)
            background.append(dict(z=z, E=expansion(z), DM=dm, DL=(1+z)*dm))
        predictions = []
        for row in request["rows"]:
            z = decimal(row["z"])
            dm, dh = distance(z), scale/expansion(z)
            value = {"DM_over_rs": dm/ruler, "DH_over_rs": dh/ruler,
                     "DV_over_rs": (dm*dm*z*dh)**(D(1)/3)/ruler}[row["observable"]]
            predictions.append(value)
        residual = [decimal(y)-mu for y, mu in zip(observed, predictions)]
        density = ldlt([[decimal(v) for v in row] for row in covariance], residual)
        return dict(background=background, predictions=predictions, **density)


def serialize(value):
    if isinstance(value, D):
        return str(value)
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serialize(v) for v in value]
    return value


def deserialize(value):
    if isinstance(value, str):
        return D(value)
    if isinstance(value, dict):
        return {k: deserialize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [deserialize(v) for v in value]
    return value


if __name__ == "__main__":
    import json
    from pathlib import Path
    import sys
    source = Path(sys.argv[1])
    if source.stat().st_size > 65536:
        raise ValueError("Reference input quota exceeded")
    inputs = json.loads(source.read_text())
    if set(inputs) != {"request", "observed", "covariance"}:
        raise ValueError("Reference input fields differ")
    output = {name: serialize(evaluate(inputs["request"], inputs["observed"], inputs["covariance"], panels))
              for name, panels in (("coarse", 128), ("fine", 256))}
    print(json.dumps(output, allow_nan=False))
