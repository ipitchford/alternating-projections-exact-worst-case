"""Finite stress checks for the exact alternating-projection formula.

All polygon projections and input geometry use fractions.Fraction.  Comparison
with the scalar maximum uses 80-digit Decimal arithmetic, not interval
arithmetic.  This is a reproducible finite diagnostic, not a universal proof.
"""
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
from itertools import combinations
from pathlib import Path
import hashlib
import json
import random
import platform
import sys

# Explicit checks: unlike `assert`, these remain active under `python -O`.
CHECKS = {"evaluated": 0, "failed": 0}
FAILURES = []


def check(condition, label):
    CHECKS["evaluated"] += 1
    if not condition:
        CHECKS["failed"] += 1
        if len(FAILURES) < 50:
            FAILURES.append(label)
        print("CHECK FAILED:", label, flush=True)
    return bool(condition)


def manuscript_sha256():
    """sha256 of ../manuscript.md (bundle root) if present, else None."""
    path = Path(__file__).resolve().parents[1] / "manuscript.md"
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def output_dir():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "results"),
                        help="directory for the JSON report (default: <bundle>/results)")
    out = Path(parser.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    return out


def finish(results, out, name):
    """Record the check tally, write JSON, print it, and exit 1 on any failure."""
    results["checks_evaluated"] = CHECKS["evaluated"]
    results["checks_failed"] = CHECKS["failed"]
    results["first_failures"] = FAILURES
    results["status"] = "PASS" if CHECKS["failed"] == 0 else "FAIL"
    results["python_optimize_flag"] = sys.flags.optimize
    rendered = json.dumps(results, indent=2) + "\n"
    (out / name).write_text(rendered)
    print(rendered, end="")
    print(f"{name}: {CHECKS['evaluated']} explicit checks, {CHECKS['failed']} failed")
    if CHECKS["failed"]:
        raise SystemExit(1)


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))

def norm2(x):
    return dot(x, x)

def dist2(x, y):
    return norm2(tuple(a-b for a, b in zip(x, y)))

class Polyhedron:
    def __init__(self, faces):
        self.faces = faces
        check(all(norm2(a) == 1 and b >= 1 for a, b in faces), 'check_convex_polygons.py:75 all(norm2(a) == 1 and b >= 1 for a, b in faces)')
        self.vertices = []
        for (a, b), (c, d) in combinations(faces, 2):
            det = a[0]*c[1] - a[1]*c[0]
            if det:
                v = ((b*c[1]-a[1]*d)/det, (a[0]*d-b*c[0])/det)
                if self.contains(v):
                    self.vertices.append(v)

    def contains(self, x):
        return all(dot(a, x) <= b for a, b in self.faces)

    def project(self, x):
        if self.contains(x):
            return x
        candidates = list(self.vertices)
        for a, b in self.faces:
            excess = dot(a, x)-b
            p = tuple(xi-excess*ai for xi, ai in zip(x, a))
            if self.contains(p):
                candidates.append(p)
        check(candidates, 'check_convex_polygons.py:96 candidates')
        return min(candidates, key=lambda p: dist2(x, p))

def decimal(x):
    return D(x.numerator)/D(x.denominator)

def worst_error(R2, N):
    """Unique quartic maximiser; delta=1."""
    if R2 <= 2:
        return D(0)
    with localcontext() as ctx:
        ctx.prec = 80
        S = (decimal(R2)-1).sqrt()
        m = 2*N-1
        lo, hi = D(0), S-1
        for _ in range(300):
            d = (lo+hi)/2
            z = S-d
            polynomial = z**4-1-4*m*z*d
            if polynomial > 0:
                lo = d
            else:
                hi = d
        d = (lo+hi)/2
        z = S-d
        c = (z*z-1)/(z*z+1)
        return +(d*c**m)

def main():
    out = output_dir()
    rng = random.Random(20261101)
    upper = [(F(0), F(1)), (F(3,5), F(4,5)), (F(-3,5), F(4,5)),
             (F(5,13), F(12,13)), (F(-5,13), F(12,13))]
    lower = [(F(4,5), F(-3,5)), (F(3,5), F(-4,5)),
             (F(5,13), F(-12,13)), (F(7,25), F(-24,25)),
             (F(8,17), F(-15,17))]
    counts = {'instances': 0, 'full_cycle_comparisons': 0,
              'nonzero_final_errors': 0, 'zero_final_errors': 0,
              'threshold_instances': 0, 'multi_facet_instances': 0}
    best_ratio = D(0)
    best_case = None
    for case in range(160):
        # The first 40 are halfspace controls; the rest have multiple facets.
        size = 1 if case < 40 else rng.randint(2, 4)
        faces_a = [(upper[0], F(1))]
        faces_b = [(rng.choice(lower), F(1))]
        for normal in rng.sample(upper[1:], size-1):
            faces_a.append((normal, F(rng.randint(3, 10), 2)))
        for normal in rng.sample(lower, size-1):
            faces_b.append((normal, F(rng.randint(3, 10), 2)))
        A, B = Polyhedron(faces_a), Polyhedron(faces_b)
        intersection = Polyhedron(faces_a+faces_b)
        if case % 16 == 0:
            x = (F(1), F(rng.choice([-1, 0, 1])))
        else:
            x = (F(rng.randint(2, 12), 1), F(rng.randint(-4, 5), 2))
        R2 = norm2(x)
        counts['instances'] += 1
        counts['multi_facet_instances'] += int(size > 1)
        counts['threshold_instances'] += int(R2 <= 2)
        for N in range(1, 9):
            x = B.project(A.project(x))
            error2 = dist2(x, intersection.project(x))
            W = worst_error(R2, N)
            with localcontext() as ctx:
                ctx.prec = 80
                actual = decimal(error2).sqrt()
                check(actual <= W*(1+D('1e-55'))+D('1e-70'), "check_convex_polygons.py:163 actual <= W*(1+D('1e-55'))+D('1e-70')")
                if W and actual/W > best_ratio:
                    best_ratio = actual/W
                    best_case = {'case': case, 'N': N, 'R_squared': str(R2)}
            counts['full_cycle_comparisons'] += 1
            counts['nonzero_final_errors'] += int(error2 > 0)
            counts['zero_final_errors'] += int(error2 == 0)
            if R2 <= 2:
                check(error2 == 0, 'check_convex_polygons.py:171 error2 == 0')
    result = {'python': platform.python_version(),
              'seed': 20261101, 'decimal_precision': 80, 'counts': counts,
              'largest_observed_error_over_formula': str(best_ratio),
              'case_at_largest_ratio': best_case,
              'scope': 'Exact rational polygon projections; rounded high-precision scalar comparisons. Finite diagnostics, not a proof or an exhaustive extremal search.'}
    result['manuscript_sha256'] = manuscript_sha256()
    finish(result, out, 'convex_polygon_checks.json')

if __name__ == '__main__':
    main()
