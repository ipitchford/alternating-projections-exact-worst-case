#!/usr/bin/env python3
"""Independent exact-rational diagnostics for the final finite-horizon proof.

These finite checks supplement the analytic audit. They are not a proof of
the unrestricted theorem. No project code is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

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


def run():
    started = time.perf_counter()
    rational_parameters = sorted({F(a, b) for b in range(1, 9) for a in range(1, b + 1)})
    angles = [( (1-q*q)/(1+q*q), 2*q/(1+q*q) ) for q in rational_parameters]
    deltas = [F(2, 3), F(1), F(3, 2)]
    amplitudes = [F(0), F(1, 7), F(1), F(4)]
    counts = dict(state_identity_cases=0, endpoint_cases=0,
                  pointwise_value_comparisons=0, quadratic_reduction_cases=0,
                  normal_cone_metric_cases=0)

    for delta in deltas:
        for c, s in angles:
            check(c*c+s*s == 1 and s > 0, 'second_exact_checks.py:73 c*c+s*s == 1 and s > 0')
            w = delta*(1+c)/s
            for amplitude in amplitudes:
                p = w+amplitude
                r2 = p*p+delta*delta
                t = s*amplitude
                kappa = 1-2*delta*delta/r2
                check(F(0) <= c <= kappa, 'second_exact_checks.py:80 F(0) <= c <= kappa')
                check(F(0) <= t <= p-delta, 'second_exact_checks.py:81 F(0) <= t <= p-delta')
                polynomial = (r2*c*c+2*delta*(delta+t)*c
                              +(delta+t)**2+delta**2-r2)
                check(polynomial == 0, 'second_exact_checks.py:84 polynomial == 0')
                p_next = c*p+delta*s
                r2_next = r2-t*t-2*delta*t
                check(r2_next == p_next*p_next+delta*delta, 'second_exact_checks.py:87 r2_next == p_next*p_next+delta*delta')
                check(r2_next >= 2*delta*delta, 'second_exact_checks.py:88 r2_next >= 2*delta*delta')
                check(s*(p_next-w) == c*t, 'second_exact_checks.py:89 s*(p_next-w) == c*t')
                check(p_next-w == c*amplitude, 'second_exact_checks.py:90 p_next-w == c*amplitude')
                kappa_next = 1-2*delta*delta/r2_next
                check(c <= kappa_next <= kappa, 'second_exact_checks.py:92 c <= kappa_next <= kappa')
                check(c*t <= p_next-delta, 'second_exact_checks.py:93 c*t <= p_next-delta')
                counts['state_identity_cases'] += 1
                if not c or not t:
                    counts['endpoint_cases'] += 1
                for c_prime, s_prime in angles:
                    if c <= c_prime <= kappa_next:
                        w_prime = delta*(1+c_prime)/s_prime
                        next_amplitude = p_next-w_prime
                        start_amplitude = p-w_prime
                        check(next_amplitude >= 0 and start_amplitude >= 0, 'second_exact_checks.py:102 next_amplitude >= 0 and start_amplitude >= 0')
                        check(next_amplitude <= c_prime*start_amplitude, 'second_exact_checks.py:103 next_amplitude <= c_prime*start_amplitude')
                        for k in range(1, 7):
                            check(next_amplitude*c_prime**(k-1) <= start_amplitude*c_prime**k, 'second_exact_checks.py:105 next_amplitude*c_prime**(k-1) <= start_amplitude*c_prime**k')
                            counts['pointwise_value_comparisons'] += 1

    # Reconstruct Lemma 2 from the orthonormal-coordinate lower bound, rather
    # than evaluating c(r,t) with a numerical square root.
    for delta in deltas:
        for t in [F(1, 9), F(1), F(7, 3)]:
            for c, s in angles:
                if not c:
                    continue
                for ratio in [c/5, c/2, c]:
                    r2 = (delta+t)**2 + (delta*(1+c)+ratio*t)**2/(1-c*c)
                    reduced = (r2*ratio**2+2*delta*(delta+t)*ratio
                               +(delta+t)**2+delta**2-r2)
                    check(reduced <= 0, 'second_exact_checks.py:119 reduced <= 0')
                    check(r2 > delta**2+(delta+t)**2, 'second_exact_checks.py:120 r2 > delta**2+(delta+t)**2')
                    counts['quadratic_reduction_cases'] += 1

    # Construct active common-ball supporting normals a=(-sigma,alpha),
    # b=(s,-c), with y-z parallel to e=(1,0). Check the metric proof's
    # norm lower bound exactly, including its equality case sigma=0.
    for delta in deltas:
        for c, s in angles:
            if not c:
                continue  # This would be the separately handled d=t case.
            for sigma, alpha in angles:
                if sigma >= s:
                    continue
                K = s*alpha-c*sigma
                check(K > 0 and K <= s, 'second_exact_checks.py:134 K > 0 and K <= s')
                z_e = delta*(alpha+c)/K
                z_l = delta*(s+sigma)/K
                check(-sigma*z_e+alpha*z_l == delta, 'second_exact_checks.py:137 -sigma*z_e+alpha*z_l == delta')
                check(s*z_e-c*z_l == delta, 'second_exact_checks.py:138 s*z_e-c*z_l == delta')
                check(z_l >= delta and z_e >= delta*(1+c)/s, 'second_exact_checks.py:139 z_l >= delta and z_e >= delta*(1+c)/s')
                for distance in [F(1, 9), F(1), F(7, 3)]:
                    U, V = distance*c/K, distance*alpha/K
                    check(U > 0 and V > 0, 'second_exact_checks.py:142 U > 0 and V > 0')
                    check(-U*sigma+V*s == distance, 'second_exact_checks.py:143 -U*sigma+V*s == distance')
                    check(U*alpha-V*c == 0, 'second_exact_checks.py:144 U*alpha-V*c == 0')
                    residual = distance*s
                    r2 = (z_e+distance)**2+z_l*z_l
                    metric_lower_radius = delta**2+(distance+delta*(1+c)/s)**2
                    check(r2 >= metric_lower_radius, 'second_exact_checks.py:148 r2 >= metric_lower_radius')
                    check(r2 >= delta**2+(delta+residual)**2, 'second_exact_checks.py:149 r2 >= delta**2+(delta+residual)**2')
                    if sigma == 0:
                        check(r2 == metric_lower_radius, 'second_exact_checks.py:151 r2 == metric_lower_radius')
                    counts['normal_cone_metric_cases'] += 1

    return {
        'python': platform.python_version(),
        'arithmetic': 'fractions.Fraction; no floating-point comparisons',
        'angle_parameters': len(rational_parameters),
        'deltas': [str(x) for x in deltas],
        'amplitudes': [str(x) for x in amplitudes],
        'horizon_exponents_checked': [1, 2, 3, 4, 5, 6],
        'counts': counts,
        'elapsed_seconds': time.perf_counter()-started,
        'manuscript_sha256': manuscript_sha256(),
        'scope': ('Finite diagnostics for scalar identities, endpoint cases, '
                  'pointwise envelope comparison, quadratic reduction and '
                  'two-normal metric geometry. The universal claims rely '
                  'on the separately recorded analytic derivation.'),
    }


if __name__ == '__main__':
    out = output_dir()
    finish(run(), out, 'second_exact_checks.json')
