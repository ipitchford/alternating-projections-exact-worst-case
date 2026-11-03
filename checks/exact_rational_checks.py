"""Independent exact-rational diagnostics for the finite-horizon proof.

No floating-point calculation is used.  These are finitely many checks of
identities and inequalities; they are not a formal proof of the theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
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


ANGLES = sorted({
    (F(0), F(1)), (F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)),
    (F(5, 13), F(12, 13)), (F(12, 13), F(5, 13)),
    (F(8, 17), F(15, 17)), (F(7, 25), F(24, 25)),
    (F(20, 29), F(21, 29)), (F(9, 41), F(40, 41)),
})


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def inverse_identity(r2, t, delta, c):
    # Verifies the positive square-root branch of the displayed inverse.
    root = r2 * c + delta * (delta + t)
    check(root >= 0, 'exact_rational_checks.py:74 root >= 0')
    check(root * root == (r2 - delta * delta) * (r2 - (delta + t) ** 2), 'exact_rational_checks.py:75 root * root == (r2 - delta * delta) * (r2 - (delta + t) ** 2)')
    check((root - delta * (delta + t)) / r2 == c, 'exact_rational_checks.py:76 (root - delta * (delta + t)) / r2 == c')
    check(r2 * c * c + 2 * delta * (delta + t) * c + (delta + t) ** 2 + delta * delta - r2 == 0, 'exact_rational_checks.py:77 r2 * c * c + 2 * delta * (delta + t) * c + (delta + t) ** 2 + delta * delta - r2 == 0')


def state_checks():
    counts = {"inverse_branches": 0, "half_map_invariance": 0,
              "transition_polynomials": 0, "value_candidates": 0,
              "fixed_residual_radius_monotonicity": 0}
    for delta in [F(1, 2), F(1), F(2)]:
        for c, s in ANGLES:
            check(c * c + s * s == 1, 'exact_rational_checks.py:86 c * c + s * s == 1')
            vertex = delta * (1 + c) / s
            for error in [F(1, 4), F(1), F(3), F(10)]:
                p = vertex + error
                r2 = p * p + delta * delta
                t = s * error
                check(0 <= t <= p - delta, 'exact_rational_checks.py:92 0 <= t <= p - delta')
                inverse_identity(r2, t, delta, c)
                counts["inverse_branches"] += 1

                rstar2 = r2 - t * t - 2 * delta * t
                pstar = c * p + delta * s
                tstar = c * t
                check(rstar2 == pstar * pstar + delta * delta, 'exact_rational_checks.py:99 rstar2 == pstar * pstar + delta * delta')
                check(tstar == pstar * s - delta * (1 + c), 'exact_rational_checks.py:100 tstar == pstar * s - delta * (1 + c)')
                check(pstar - vertex == c * error, 'exact_rational_checks.py:101 pstar - vertex == c * error')
                check(0 <= tstar <= pstar - delta, 'exact_rational_checks.py:102 0 <= tstar <= pstar - delta')
                inverse_identity(rstar2, tstar, delta, c)
                counts["inverse_branches"] += 1
                counts["half_map_invariance"] += 1

                for fraction in [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]:
                    q = c * fraction
                    polynomial = r2 * q * q + 2 * delta * (delta + t) * q + (delta + t) ** 2 + delta * delta - r2
                    check(polynomial <= 0, 'exact_rational_checks.py:110 polynomial <= 0')
                    counts["transition_polynomials"] += 1

                kappa_star = 1 - 2 * delta * delta / rstar2
                for cp, sp in ANGLES:
                    if c <= cp <= kappa_star:
                        b = delta * (1 + cp) / sp
                        dpre, dpost = p - b, pstar - b
                        check(dpost >= 0 and dpre >= 0, 'exact_rational_checks.py:118 dpost >= 0 and dpre >= 0')
                        check(pstar <= cp * p + delta * sp, 'exact_rational_checks.py:119 pstar <= cp * p + delta * sp')
                        for horizon in [1, 2, 3, 5, 10]:
                            check(dpost * cp ** (horizon - 1) <= dpre * cp ** horizon, 'exact_rational_checks.py:121 dpost * cp ** (horizon - 1) <= dpre * cp ** horizon')
                            counts["value_candidates"] += 1

        # At fixed residual cap, increasing the inverse angle constructs an
        # increasing-radius comparison and increases every value candidate.
        for t in [F(1, 10), F(1), F(5)]:
            for c0, s0 in ANGLES:
                p0 = (t + delta * (1 + c0)) / s0
                for c1, s1 in ANGLES:
                    if c1 >= c0:
                        p1 = (t + delta * (1 + c1)) / s1
                        check(p1 >= p0, 'exact_rational_checks.py:132 p1 >= p0')
                        for horizon in [0, 1, 2, 5]:
                            check((t / s1) * c1 ** horizon >= (t / s0) * c0 ** horizon, 'exact_rational_checks.py:134 (t / s1) * c1 ** horizon >= (t / s0) * c0 ** horizon')
                            counts["fixed_residual_radius_monotonicity"] += 1
    return counts


def metric_checks():
    counts = {"metric_normal_decompositions": 0,
              "metric_coordinate_inequalities": 0,
              "metric_sharp_equality_cases": 0}
    for delta in [F(1, 2), F(1), F(2)]:
        for c, s in ANGLES:
            if c == 0:
                continue  # The d=t degenerate case is handled separately.
            for sigma, a0 in ANGLES:
                if not sigma < s:
                    continue
                a, b = (-sigma, a0), (s, -c)
                K = s * a0 - c * sigma
                check(K > 0, 'exact_rational_checks.py:152 K > 0')
                z_l = delta * (s + sigma) / K
                z_e = (delta + c * z_l) / s
                check(dot(a, (z_e, z_l)) == delta, 'exact_rational_checks.py:155 dot(a, (z_e, z_l)) == delta')
                check(dot(b, (z_e, z_l)) == delta, 'exact_rational_checks.py:156 dot(b, (z_e, z_l)) == delta')
                for error in [F(1, 3), F(1), F(4)]:
                    U, V = error * c / K, error * a0 / K
                    check(U > 0 and V > 0, 'exact_rational_checks.py:159 U > 0 and V > 0')
                    check(tuple(U * x + V * y for x, y in zip(a, b)) == (error, F(0)), 'exact_rational_checks.py:160 tuple(U * x + V * y for x, y in zip(a, b)) == (error, F(0))')
                    counts["metric_normal_decompositions"] += 1
                    t = error * s
                    for dilation in [F(1), F(2)]:
                        for orthogonal in [F(0), delta]:
                            y = (dilation * z_e + error, dilation * z_l, orthogonal)
                            r2 = dot(y, y)
                            check(dilation * z_l >= delta, 'exact_rational_checks.py:167 dilation * z_l >= delta')
                            check(dilation * z_e >= delta * (1 + c) / s, 'exact_rational_checks.py:168 dilation * z_e >= delta * (1 + c) / s')
                            check(r2 >= delta * delta + (error + delta * (1 + c) / s) ** 2, 'exact_rational_checks.py:169 r2 >= delta * delta + (error + delta * (1 + c) / s) ** 2')
                            check(r2 > delta * delta + (delta + t) ** 2, 'exact_rational_checks.py:170 r2 > delta * delta + (delta + t) ** 2')
                            counts["metric_coordinate_inequalities"] += 1
                            if sigma == 0 and dilation == 1 and orthogonal == 0:
                                check(r2 == delta * delta + (error + delta * (1 + c) / s) ** 2, 'exact_rational_checks.py:173 r2 == delta * delta + (error + delta * (1 + c) / s) ** 2')
                                inverse_identity(r2, t, delta, c)
                                counts["metric_sharp_equality_cases"] += 1
    return counts


def main():
    out = output_dir()
    results = {
        "scope": "Finite exact-rational identities and inequalities, not a formal or exhaustive proof.",
        "arithmetic": "fractions.Fraction only; no floating point",
        "python": platform.python_version(),
        "counts": state_checks() | metric_checks(),
        "manuscript_sha256": manuscript_sha256(),
    }
    finish(results, out, "exact_rational_checks.json")


if __name__ == "__main__":
    main()
