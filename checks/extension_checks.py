#!/usr/bin/env python3
"""Finite diagnostics for the two-radius extension (Theorems 7 and 8).

Setting: closed convex A containing the closed ball B(0,a) and B containing
B(0,b); alternating projections start with P_A (order AB) or with P_B (order
BA). For a point in X (radius rho_X) whose next projection is onto Y (radius
rho_Y) the scalar functions are p_X, T_XY, w_XY, D_XY, tau_XY, kappa_XY, c_XY,
E_XY and V^XY_k of the extension specification.

Arithmetic, stated per section:
  * fractions.Fraction (exact) for (i)-exact, (ii)-exact, (iii)-exact and (vi)-exact,
    using rational angles c=(1-q^2)/(1+q^2), s=2q/(1+q^2) and radii r with
    sqrt(r^2-rho_X^2) rational, so every identity is tested with no rounding;
  * decimal.Decimal with 60 significant digits (rounded, not interval arithmetic)
    for random irrational states, kappa/arccos comparisons, value-function
    maximisations, attainment at the optimal angle, the adversarial search,
    the Zamani-Glineur limit, the asymptotics and the tables.
Decimal comparisons use explicit tolerances recorded in the JSON report.

These are finitely many checks. They can expose transcription or
implementation errors; they do not prove the universally quantified theorems.
Standard library only. Usage: python3 -I -B checks/extension_checks.py [--out DIR]
"""
from decimal import Decimal as D, getcontext, localcontext
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import platform
import random
import sys
import time

PREC = 60
getcontext().prec = PREC
TOL = D('1e-45')            # absolute/relative tolerance for 60-digit Decimal identities
SEARCH_REL_TOL = D('1e-12')  # a search ratio above 1 + this is a failure (reported loudly)
SEED = 20261103

# Explicit checks: unlike `assert`, these remain active under `python -O`.
CHECKS = {"evaluated": 0, "failed": 0}
FAILURES = []
SECTION = {"name": None}
SECTION_COUNTS = {}


def check(condition, label):
    CHECKS["evaluated"] += 1
    SECTION_COUNTS[SECTION["name"]] = SECTION_COUNTS.get(SECTION["name"], 0) + 1
    if not condition:
        CHECKS["failed"] += 1
        if len(FAILURES) < 50:
            FAILURES.append(f"[{SECTION['name']}] {label}")
        print(f"CHECK FAILED [{SECTION['name']}]:", label, flush=True)
    return bool(condition)


def section(name):
    SECTION["name"] = name
    SECTION_COUNTS.setdefault(name, 0)
    print(f"\n=== {name} ===", flush=True)


def manuscript_sha256():
    path = Path(__file__).resolve().parents[1] / "manuscript.md"
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def output_dir():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "results"))
    out = Path(parser.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    return out


# ---------------------------------------------------------------- Decimal helpers
def dec(x):
    if isinstance(x, F):
        return D(x.numerator) / D(x.denominator)
    return D(x)


def dsqrt(x):
    if x < 0:
        if x > -TOL:
            return D(0)
        raise ValueError(f"square root of negative {x}")
    return x.sqrt()


def d_atan(x):
    with localcontext() as ctx:
        ctx.prec = PREC + 15
        k = 0
        while abs(x) > D("0.05"):
            x = x / (1 + (1 + x * x).sqrt())
            k += 1
        x2, term, total, n = x * x, x, x, 1
        while True:
            term *= -x2
            n += 2
            nt = term / n
            if abs(nt) < D(10) ** -(PREC + 12):
                break
            total += nt
        result = total * (2 ** k)
    return +result


def d_asin(x):
    return d_atan(x / dsqrt(1 - x * x))


def d_cos(x):
    with localcontext() as ctx:
        ctx.prec = PREC + 15
        x2, term, total, n = x * x, D(1), D(1), 0
        while True:
            n += 2
            term *= -x2 / (n * (n - 1))
            if abs(term) < D(10) ** -(PREC + 12):
                break
            total += term
    return +total


# ------------------------------------------------- scalar functions (generic order)
def p_of(r, rX):
    return dsqrt(r * r - rX * rX)


def T_of(r, rX, rY):
    return p_of(r, rX) - rY


def w_of(c, rX, rY):
    return (rY + c * rX) / dsqrt(1 - c * c)


def wprime(c, rX, rY):
    s = dsqrt(1 - c * c)
    return (rX + rY * c) / (s * s * s)


def D_of(r, c, rX, rY):
    return p_of(r, rX) - w_of(c, rX, rY)


def tau_of(r, c, rX, rY):
    return p_of(r, rX) * dsqrt(1 - c * c) - rY - c * rX


def kappa_of(r, rX, rY):
    r2 = r * r
    return (-rX * rY + dsqrt((r2 - rX * rX) * (r2 - rY * rY))) / r2


def c_of(r, t, rX, rY):
    r2 = r * r
    return (-rX * (rY + t) + dsqrt((r2 - rX * rX) * (r2 - (rY + t) ** 2))) / r2


def E_of(r, t, rX, rY):
    return D_of(r, c_of(r, t, rX, rY), rX, rY)


def bisect_decreasing(g, lo, hi, iters=210):
    """Root of a decreasing function with g(lo) >= 0 >= g(hi)."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def argmax_D_ck(r, k, rX, rY):
    """Maximiser over [0, kappa] of c^k D_XY(r,c) (log-concave; k >= 1)."""
    kap = kappa_of(r, rX, rY)
    return bisect_decreasing(lambda c: k * D_of(r, c, rX, rY) - c * wprime(c, rX, rY), D(0), kap)


def V_of(k, r, t, rX, rY):
    """V^XY_k(r,t) = max over c in [c_XY(r,t), kappa_XY(r)] of D_XY(r,c) c^k."""
    clo = c_of(r, t, rX, rY)
    if k == 0:
        return D_of(r, clo, rX, rY)
    c = max(argmax_D_ck(r, k, rX, rY), clo)
    return D_of(r, c, rX, rY) * c ** k


def W7(m, R, r1, r2):
    """Theorem 7: max_{0<=c<=kappa} c^m D(R,c); r1 = radius of the first-projected set."""
    if R * R <= r1 * r1 + r2 * r2:
        return D(0), None
    c = argmax_D_ck(R, m, r1, r2)
    return D_of(R, c, r1, r2) * c ** m, c


def W8(m, R, r1, r2):
    """Theorem 8: max_{0<=c<=kappa} c^m s(c) D(R,c) = max c^m tau(R,c)."""
    if R * R <= r1 * r1 + r2 * r2:
        return D(0), None
    p = p_of(R, r1)
    kap = kappa_of(R, r1, r2)

    def g(c):
        s = dsqrt(1 - c * c)
        return m * (p * s - r2 - c * r1) - c * (p * c / s + r1)
    c = bisect_decreasing(g, D(0), kap)
    return c ** m * tau_of(R, c, r1, r2), c


def close(x, y, scale=1):
    return abs(x - y) <= TOL * (1 + abs(dec(scale)))


# ------------------------------------------------------------- exact rational states
def pyth(q):
    return (1 - q * q) / (1 + q * q), 2 * q / (1 + q * q)


def rational_state(rng, rX, rY):
    """Exact state (c, s, p, r2, t) with c = c_XY(r,t) by construction."""
    q = F(rng.randint(1, 60), 60)
    c, s = pyth(q)
    amp = rng.choice([F(0), F(rng.randint(1, 40), rng.randint(1, 12))])
    p = (rY + c * rX) / s + amp
    return c, s, p, p * p + rX * rX, s * amp, amp


RADII = [(F(1), F(1)), (F(3, 10), F(1, 5)), (F(1, 5), F(3, 10)), (F(1), F(1, 2)),
         (F(1, 2), F(2)), (F(2), F(1, 3)), (F(1, 10), F(1))]


def orders(a, b):
    return [("AB", a, b), ("BA", b, a)]


# ------------------------------------------------------------------ sections
def sec_i_to_iv(rng, results):
    section("(i)-(iii) exact rational states")
    n_states = 0
    q_grid = [F(j, 24) for j in range(1, 25)]
    for a, b in RADII:
        for name, rX, rY in orders(a, b):
            for _ in range(40):
                c0, s0, p, r2, t, amp = rational_state(rng, rX, rY)
                n_states += 1
                # (i) the closed form (positive branch) equals c0, and c0 in [0, kappa]
                root = r2 * c0 + rX * (rY + t)
                check(root >= 0, f"{name} root>=0")
                check(rY + t <= p, f"{name} t <= T_XY(r)")
                check(root * root == (r2 - rX * rX) * (r2 - (rY + t) ** 2), f"{name} c_XY closed form")
                check(p * s0 - rY - c0 * rX == t, f"{name} tau_XY(r,c0) == t")
                check(c0 >= 0 and p - (rY + c0 * rX) / s0 >= 0, f"{name} 0 <= c0 <= kappa_XY(r)")
                # (ii) halfspace update identities
                rs2 = r2 - t * t - 2 * rY * t
                pY = c0 * p + rX * s0
                ts = c0 * t
                check(pY >= 0 and rs2 - rY * rY == pY * pY, f"{name} p_Y(r*) = c0 p_X + rho_X s0")
                root2 = rs2 * c0 + rY * (rX + ts)
                check(root2 >= 0 and rX + ts <= pY, f"{name} updated state admissible")
                check(root2 * root2 == (rs2 - rY * rY) * (rs2 - (rX + ts) ** 2), f"{name} c_YX(r*,t*) = c0")
                DYX = pY - (rX + c0 * rY) / s0
                DXY = p - (rY + c0 * rX) / s0
                check(DYX == c0 * DXY, f"{name} D_YX(r*,c0) = c0 D_XY(r,c0)")
                # (iii) key inequality on a grid of c' in [c0, kappa_XY(r)]
                h0 = c0 * p + rX * s0
                for qq in q_grid:
                    cp, sp = pyth(qq)
                    if cp < c0 or p - (rY + cp * rX) / sp < 0:
                        continue
                    lhs = pY - (rX + cp * rY) / sp
                    rhs = cp * (p - (rY + cp * rX) / sp)
                    hp = cp * p + rX * sp
                    check(lhs - rhs == h0 - hp, f"{name} D_YX - c'D_XY == h(c0)-h(c')")
                    check(hp >= h0, f"{name} h increasing on [c0,kappa]")
                    check(lhs <= rhs, f"{name} D_YX(r*,c') <= c' D_XY(r,c')")
                # kappa_YX(r*) <= kappa_XY(r) (Decimal; kappa is irrational in general)
                rs = dsqrt(dec(rs2))
                check(kappa_of(rs, dec(rY), dec(rX)) <= kappa_of(dsqrt(dec(r2)), dec(rX), dec(rY)) + TOL,
                      f"{name} kappa_YX(r*) <= kappa_XY(r)")
    results["exact_rational_states"] = n_states

    section("(i)-(ii) random irrational states (Decimal 60)")
    n_dec = 0
    for _ in range(150):
        a = dec(F(rng.randint(1, 30), rng.randint(1, 15)))
        b = dec(F(rng.randint(1, 30), rng.randint(1, 15)))
        for name, rX, rY in orders(a, b):
            rmin = dsqrt(rX * rX + rY * rY)
            r = rmin * (1 + dec(F(rng.randint(0, 400), 100)))
            t = T_of(r, rX, rY) * dec(F(rng.randint(0, 1000), 1000))
            c0 = c_of(r, t, rX, rY)
            kap = kappa_of(r, rX, rY)
            n_dec += 1
            check(close(tau_of(r, c0, rX, rY), t, r), f"{name} tau(r,c_XY(r,t)) == t")
            check(-TOL <= c0 <= kap + TOL, f"{name} c in [0,kappa]")
            check(close(E_of(r, t, rX, rY), t / dsqrt(1 - c0 * c0), r), f"{name} E = t/s(c)")
            s0 = dsqrt(1 - c0 * c0)
            rs = dsqrt(r * r - t * t - 2 * rY * t)
            check(close(p_of(rs, rY), c0 * p_of(r, rX) + rX * s0, r), f"{name} p_Y(r*) identity")
            check(close(c_of(rs, c0 * t, rY, rX), c0), f"{name} c_YX(r*,t*) = c0")
            check(close(D_of(rs, c0, rY, rX), c0 * D_of(r, c0, rX, rY), r), f"{name} D_YX = c0 D_XY")
            check(kappa_of(rs, rY, rX) <= kap + TOL, f"{name} kappa_YX(r*) <= kappa_XY(r)")
            check(close(kappa_of(r, rX, rY), kappa_of(r, rY, rX)), f"{name} kappa symmetric")
    results["decimal_random_states"] = n_dec

    section("(iv) kappa and c_XY closed forms vs arccos/arcsin formulas (Decimal 60)")
    worst = D(0)
    for _ in range(120):
        a = dec(F(rng.randint(1, 30), rng.randint(1, 15)))
        b = dec(F(rng.randint(1, 30), rng.randint(1, 15)))
        for name, rX, rY in orders(a, b):
            r = dsqrt(rX * rX + rY * rY) * (1 + dec(F(rng.randint(1, 500), 100)))
            kap = kappa_of(r, rX, rY)
            trig = d_cos(d_asin(rX / r) + d_asin(rY / r))
            worst = max(worst, abs(kap - trig))
            check(abs(kap - trig) <= TOL, f"{name} kappa closed form == cos(asin+asin)")
            check(abs(D_of(r, kap, rX, rY)) <= TOL * (1 + r), f"{name} D_XY(r,kappa) == 0")
            t = T_of(r, rX, rY) * dec(F(rng.randint(0, 100), 100))
            trig_c = d_cos(d_asin(rX / r) + d_asin((rY + t) / r))
            worst = max(worst, abs(c_of(r, t, rX, rY) - trig_c))
            check(abs(c_of(r, t, rX, rY) - trig_c) <= TOL, f"{name} c_XY closed form == cos(asin+asin)")
            # independent float sanity check against math.acos-free double evaluation
            kf = math.cos(math.asin(float(rX / r)) + math.asin(float(rY / r)))
            check(abs(float(kap) - kf) <= 1e-12, f"{name} kappa vs float trig")
    results["kappa_trig_max_abs_difference"] = format(worst, ".3e")


def sec_v(rng, results):
    section("(v) monotonicity of V^XY_k in r and t; Lemma 4' comparison (Decimal 60)")
    n = 0
    for a, b in [(F(1), F(1)), (F(3, 10), F(1, 5)), (F(1), F(1, 2))]:
        for name, rX, rY in orders(dec(a), dec(b)):
            rmin = dsqrt(rX * rX + rY * rY)
            rs_grid = [rmin * (1 + D(i) / 6) for i in range(0, 13)]
            for k in [0, 1, 2, 3, 5, 8]:
                for i in range(len(rs_grid) - 1):
                    r, r_next = rs_grid[i], rs_grid[i + 1]
                    T = T_of(r, rX, rY)
                    tvals = [T * D(j) / 8 for j in range(9)]
                    vals = [V_of(k, r, t, rX, rY) for t in tvals]
                    for j in range(8):
                        n += 1
                        check(vals[j] <= vals[j + 1] + TOL * (1 + r), f"{name} V_{k} nondecreasing in t")
                    for t, v in zip(tvals, vals):
                        n += 1
                        check(v <= V_of(k, r_next, t, rX, rY) + TOL * (1 + r), f"{name} V_{k} nondecreasing in r")
                        if k >= 1:
                            c0 = c_of(r, t, rX, rY)
                            rs = dsqrt(r * r - t * t - 2 * rY * t)
                            n += 1
                            check(V_of(k - 1, rs, c0 * t, rY, rX) <= v + TOL * (1 + r),
                                  f"{name} V^YX_{k-1}(r*,t*) <= V^XY_{k}(r,t)")
    results["value_function_comparisons"] = n


def halfspace_project(x, n, h):
    viol = x[0] * n[0] + x[1] * n[1] - h
    if viol <= 0:
        return x
    return (x[0] - viol * n[0], x[1] - viol * n[1])


def sec_vi(rng, results):
    section("(vi) attainment by the two-halfspace construction (exact, then Decimal at the optimum)")
    n_exact = 0
    for a, b in RADII:
        for q in [F(1, 7), F(1, 3), F(1, 2), F(2, 3), F(5, 6)]:
            c, s = pyth(q)
            for amp in [F(1, 9), F(1), F(7, 2)]:
                p = (b + a * c) / s + amp           # p = sqrt(R^2 - a^2), D_AB(R,c) = amp
                R2 = p * p + a * a
                nA, nB = (F(0), F(1)), (s, -c)
                z = ((b + c * a) / s, a)
                check(s * z[0] - c * z[1] == b and z[1] == a, "vertex on both boundary lines")
                x0 = (p, a)
                y = halfspace_project(x0, nA, a)
                check(y == x0, "x0 in A")
                for j in range(1, 9):
                    y = halfspace_project(y, nB if j % 2 else nA, b if j % 2 else a)
                    diff = (y[0] - z[0], y[1] - z[1])
                    beta = diff[0] / s
                    alpha = diff[1] + c * beta
                    check(alpha >= 0 and beta >= 0 and diff == (beta * s, alpha - beta * c),
                          "y_m - z in N_{A cap B}(z), so z = P_{A cap B}(y_m)")
                    check(diff[0] ** 2 + diff[1] ** 2 == (c ** j * amp) ** 2, f"dist(y_{j},A cap B) = c^m D_AB(R,c)")
                    other = (y[1] - a) if j % 2 else (s * y[0] - c * y[1] - b)
                    check(max(F(0), other) == c ** j * s * amp, f"residual at y_{j} = c^m s D_AB(R,c)")
                    n_exact += 1
                    if j <= 5 and amp == 1:
                        Rd = dsqrt(dec(R2))
                        check(dec(c ** j * amp) <= W7(j, Rd, dec(a), dec(b))[0] * (1 + TOL), "instance <= Theorem 7 value")
                        check(dec(c ** j * s * amp) <= W8(j, Rd, dec(a), dec(b))[0] * (1 + TOL), "instance <= Theorem 8 value")
    results["exact_attainment_trajectory_points"] = n_exact

    # Decimal trajectory at the optimal angle reproduces the maximum (equality).
    n_opt, worst = 0, D(0)
    for a, b in RADII:
        a, b = dec(a), dec(b)
        R = dsqrt(a * a + b * b) * D(3)
        for m in range(1, 9):
            for which in (7, 8):
                Wv, c = (W7 if which == 7 else W8)(m, R, a, b)
                s = dsqrt(1 - c * c)
                z = ((b + c * a) / s, a)
                y = (p_of(R, a), a)
                for j in range(1, m + 1):
                    y = halfspace_project(y, (s, -c) if j % 2 else (D(0), D(1)), b if j % 2 else a)
                dist = dsqrt((y[0] - z[0]) ** 2 + (y[1] - z[1]) ** 2)
                resid = max(D(0), (y[1] - a) if m % 2 else (s * y[0] - c * y[1] - b))
                val = dist if which == 7 else resid
                worst = max(worst, abs(val / Wv - 1))
                check(abs(val / Wv - 1) <= TOL, f"Theorem {which} attained at optimal c, m={m}")
                n_opt += 1
    results["optimal_angle_attainment_cases"] = n_opt
    results["optimal_angle_attainment_max_rel_error"] = format(worst, ".3e")

    # Cross-check the bisection maximiser against a dense grid.
    for a, b in [(D(1), D(1)), (D("0.3"), D("0.2")), (D("0.2"), D("0.3"))]:
        R = D(1) if a < 1 else D(3)
        for m in [1, 2, 5, 13]:
            for which, fn in ((7, W7), (8, W8)):
                Wv, _ = fn(m, R, a, b)
                kap = kappa_of(R, a, b)
                grid = []
                for i in range(1, 2000):
                    c = kap * i / 2000
                    grid.append(c ** m * (D_of(R, c, a, b) if which == 7 else tau_of(R, c, a, b)))
                g = max(grid)
                check(g <= Wv * (1 + TOL) and g >= Wv * (1 - D("1e-4")), f"bisection max vs grid, Thm {which}, m={m}")


# --------------------------------------------------------- adversarial search
def constraints_of(params, rho):
    cons = []
    for theta, extra in params["faces"]:
        n = (D(repr(math.cos(theta))), D(repr(math.sin(theta))))
        nn = dsqrt(n[0] ** 2 + n[1] ** 2)
        cons.append(("L", (n[0] / nn, n[1] / nn), rho * (1 + D(repr(extra)))))
    if params["disc"] is not None:
        phi, cfrac, extra = params["disc"]
        ctr = (rho * D(repr(cfrac * math.cos(phi))), rho * D(repr(cfrac * math.sin(phi))))
        rad = dsqrt(ctr[0] ** 2 + ctr[1] ** 2) + rho * (1 + D(repr(extra)))
        cons.append(("C", ctr, rad))
    return cons


def violation(x, con):
    if con[0] == "L":
        return x[0] * con[1][0] + x[1] * con[1][1] - con[2]
    return dsqrt((x[0] - con[1][0]) ** 2 + (x[1] - con[1][1]) ** 2) - con[2]


def proj_one(x, con):
    if con[0] == "L":
        return halfspace_project(x, con[1], con[2])
    ctr, rad = con[1], con[2]
    d = dsqrt((x[0] - ctr[0]) ** 2 + (x[1] - ctr[1]) ** 2)
    if d <= rad:
        return x
    return (ctr[0] + rad * (x[0] - ctr[0]) / d, ctr[1] + rad * (x[1] - ctr[1]) / d)


def boundary_meets(c1, c2):
    if c1[0] == "C" and c2[0] == "L":
        c1, c2 = c2, c1
    if c1[0] == "L" and c2[0] == "L":
        (a1, b1), h1 = c1[1], c1[2]
        (a2, b2), h2 = c2[1], c2[2]
        det = a1 * b2 - a2 * b1
        if abs(det) < D("1e-30"):
            return []
        return [((h1 * b2 - b1 * h2) / det, (a1 * h2 - a2 * h1) / det)]
    if c1[0] == "L":
        n, h = c1[1], c1[2]
        ctr, rad = c2[1], c2[2]
        off = h - (n[0] * ctr[0] + n[1] * ctr[1])
        if abs(off) > rad:
            return []
        half = dsqrt(rad * rad - off * off)
        p0 = (ctr[0] + off * n[0], ctr[1] + off * n[1])
        u = (-n[1], n[0])
        return [(p0[0] + half * u[0], p0[1] + half * u[1]), (p0[0] - half * u[0], p0[1] - half * u[1])]
    (x1, y1), r1 = c1[1], c1[2]
    (x2, y2), r2 = c2[1], c2[2]
    dx, dy = x2 - x1, y2 - y1
    d2 = dx * dx + dy * dy
    if d2 == 0:
        return []
    d = dsqrt(d2)
    if d > r1 + r2 or d < abs(r1 - r2):
        return []
    al = (r1 * r1 - r2 * r2 + d2) / (2 * d)
    hh = dsqrt(r1 * r1 - al * al)
    mx, my = x1 + al * dx / d, y1 + al * dy / d
    return [(mx + hh * dy / d, my - hh * dx / d), (mx - hh * dy / d, my + hh * dx / d)]


FEAS_TOL = D("1e-40")


def project(x, cons):
    """Euclidean projection onto a planar intersection of halfplanes and discs.

    At the projection at most two constraints need be active (planar KKT), so the
    projection is the nearest feasible point among: x itself, single-constraint
    projections, and pairwise boundary intersections.
    """
    if not cons:
        return x
    feasible = lambda p: all(violation(p, c) <= FEAS_TOL for c in cons)
    if feasible(x):
        return x
    cands = [proj_one(x, c) for c in cons]
    for i in range(len(cons)):
        for j in range(i + 1, len(cons)):
            cands.extend(boundary_meets(cons[i], cons[j]))
    cands = [p for p in cands if feasible(p)]
    if not check(bool(cands), "projection candidate set nonempty"):
        raise RuntimeError("no feasible projection candidate")
    return min(cands, key=lambda p: (x[0] - p[0]) ** 2 + (x[1] - p[1]) ** 2)


def dist(x, y):
    return dsqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)


def random_set(rng, kind):
    if kind == "poly":
        faces = [[rng.uniform(0, 2 * math.pi), 0.0 if rng.random() < 0.5 else rng.expovariate(1.0)]
                 for _ in range(rng.randint(1, 4))]
        return {"faces": faces, "disc": None}
    faces = [[rng.uniform(0, 2 * math.pi), 0.0 if rng.random() < 0.5 else rng.expovariate(1.0)]
             for _ in range(rng.randint(1, 2))]
    return {"faces": faces, "disc": [rng.uniform(0, 2 * math.pi), rng.uniform(0, 2), rng.uniform(0, 3)]}


def random_instance(rng):
    a, b = rng.choice(SEARCH_RADII)
    return {"a": a, "b": b, "A": random_set(rng, rng.choice(["poly", "poly", "disc"])),
            "B": random_set(rng, rng.choice(["poly", "poly", "disc"])),
            "psi": rng.uniform(0, 2 * math.pi), "Rfac": math.exp(rng.uniform(0, math.log(8)))}


def perturb(rng, inst, scale):
    new = json.loads(json.dumps(inst))
    for key in ("A", "B"):
        for face in new[key]["faces"]:
            face[0] += rng.gauss(0, 0.3 * scale)
            face[1] = 0.0 if rng.random() < 0.3 else max(0.0, face[1] + rng.gauss(0, 0.3 * scale))
        if new[key]["disc"] is not None:
            d = new[key]["disc"]
            d[0] += rng.gauss(0, 0.3 * scale)
            d[1] = max(0.0, d[1] + rng.gauss(0, 0.3 * scale))
            d[2] = max(0.0, d[2] + rng.gauss(0, 0.3 * scale))
    new["psi"] += rng.gauss(0, 0.3 * scale)
    new["Rfac"] = min(8.0, max(1.0, new["Rfac"] * math.exp(rng.gauss(0, 0.2 * scale))))
    return new


SEARCH_RADII = [(1.0, 1.0), (1.0, 0.5), (0.5, 1.0), (0.3, 0.2), (0.2, 0.3), (1.0, 0.1)]


def evaluate(inst, order, mmax=5):
    a, b = D(repr(inst["a"])), D(repr(inst["b"]))
    consA, consB = constraints_of(inst["A"], a), constraints_of(inst["B"], b)
    both = consA + consB
    first, second = (consA, consB) if order == "AB" else (consB, consA)
    r1, r2 = (a, b) if order == "AB" else (b, a)
    R = dsqrt(a * a + b * b) * D(repr(inst["Rfac"]))
    x0 = (R * D(repr(math.cos(inst["psi"]))), R * D(repr(math.sin(inst["psi"]))))
    R = dsqrt(x0[0] ** 2 + x0[1] ** 2)
    y = project(x0, first)
    out = []
    for m in range(1, mmax + 1):
        produced_by_second = (m % 2 == 1)
        y = project(y, second if produced_by_second else first)
        other = first if produced_by_second else second
        d_int = dist(y, project(y, both))
        d_res = dist(y, project(y, other))
        w7, _ = W7(m, R, r1, r2)
        w8, _ = W8(m, R, r1, r2)
        out.append((m, d_int, d_res, w7, w8))
    return out


def ratio(x, w):
    if w == 0:
        return D(0) if x <= D("1e-35") else D("Infinity")
    return x / w


def sec_vii(rng, results):
    section("(vii) randomized adversarial search on polygons and disc-halfspace sets (Decimal 60)")
    best = {}       # (order, kind, m) -> (ratio, instance)
    n_inst = 0

    def record(inst, order):
        nonlocal n_inst
        n_inst += 1
        vals = evaluate(inst, order)
        for m, d_int, d_res, w7, w8 in vals:
            for kind, x, w in (("dist", d_int, w7), ("resid", d_res, w8)):
                rr = ratio(x, w)
                if rr > 1 + SEARCH_REL_TOL:
                    print(f"!!! SEARCH VIOLATION: order={order} kind={kind} m={m} ratio={rr} instance={json.dumps(inst)}",
                          flush=True)
                check(rr <= 1 + SEARCH_REL_TOL, f"{order} {kind} m={m} ratio {rr} <= 1")
                key = (order, kind, m)
                if key not in best or rr > best[key][0]:
                    best[key] = (rr, inst)
        return vals

    pool = {"AB": [], "BA": []}
    for order in ("AB", "BA"):
        for _ in range(200):
            inst = random_instance(rng)
            record(inst, order)
            pool[order].append(inst)
    # Hill climbing: for each (order, kind, m) start from the best random instance.
    climbs = 0
    for order in ("AB", "BA"):
        for kind in ("dist", "resid"):
            for m in range(1, 6):
                cur_r, cur = best[(order, kind, m)]
                for step in range(16):
                    cand = perturb(rng, cur, 1.0 if step < 8 else 0.3)
                    vals = record(cand, order)
                    climbs += 1
                    _, d_int, d_res, w7, w8 = vals[m - 1]
                    rr = ratio(d_int, w7) if kind == "dist" else ratio(d_res, w8)
                    if rr > cur_r:
                        cur_r, cur = rr, cand
    table = {}
    overall = D(0)
    for (order, kind, m), (rr, inst) in sorted(best.items()):
        table[f"{order}/{kind}/m={m}"] = format(rr, ".12f")
        overall = max(overall, rr)
    results["search"] = {"instances_evaluated": n_inst, "hill_climb_steps": climbs,
                         "orders": ["AB", "BA"], "m_range": [1, 5],
                         "max_ratio_overall": format(overall, ".15f"),
                         "max_ratio_by_order_kind_m": table,
                         "failure_threshold": f"ratio > 1 + {SEARCH_REL_TOL}"}
    print(f"instances evaluated: {n_inst} (incl. {climbs} hill-climb steps)")
    print("max ratio (actual / theorem value), dist = Theorem 7, resid = Theorem 8:")
    for key, val in table.items():
        print(f"  {key:22s} {val}")
    print(f"overall max ratio: {format(overall, '.15f')}")


def sec_viii(results):
    section("(viii) delta -> 0 limit of Theorem 8 with m = 2N vs Zamani-Glineur")
    rows = []
    for N in range(1, 7):
        zg = dsqrt(D(2 * N) ** (2 * N) / D(2 * N + 1) ** (2 * N + 1))
        for e in (8, 16, 24):
            delta = D(10) ** -e
            val, _ = W8(2 * N, D(1), delta, delta)
            diff = abs(val - zg)
            check(diff <= 4 * delta, f"N={N} delta=1e-{e}: |W8 - ZG| <= 4 delta")
            rows.append({"N": N, "delta": f"1e-{e}", "theorem8": format(val, ".25e"),
                         "zamani_glineur": format(zg, ".25e"), "abs_diff": format(diff, ".3e")})
            if e == 24:
                print(f"N={N}: ZG={format(zg, '.20f')}  |Thm8(delta=1e-24) - ZG| = {format(diff, '.2e')}")
    results["zamani_glineur_limit"] = rows


def sec_ix(results):
    section("(ix) asymptotic ratio W_m e m / (w'(kappa) kappa^{m+1}) -> 1")
    rows = []
    e1 = D(1).exp()
    for a, b, R in [(D(1), D("0.5"), D(3)), (D("0.3"), D("0.2"), D(1))]:
        for name, r1, r2 in orders(a, b):
            kap = kappa_of(R, r1, r2)
            prev = None
            for m in (100, 1000, 10000):
                Wv, _ = W7(m, R, r1, r2)
                rr = Wv * e1 * m / (wprime(kap, r1, r2) * kap ** (m + 1))
                gap = abs(rr - 1)
                if prev is not None:
                    check(gap < prev, f"{name} a={a} b={b}: |ratio-1| decreases in m")
                prev = gap
                rows.append({"a": str(a), "b": str(b), "R": str(R), "order": name, "m": m,
                             "kappa": format(kap, ".15f"), "ratio": format(rr, ".10f")})
                print(f"a={a} b={b} R={R} order={name} m={m:>6}: ratio = {format(rr, '.8f')}")
            check(prev < D("0.01"), f"{name}: ratio within 1% at m=1e4")
    results["asymptotic_ratios"] = rows


def W_equal_over_R(ratio_, N):
    R = D(ratio_)
    return W7(2 * N - 1, R, D(1), D(1))[0] / R


def smallest_N(pred):
    hi = 1
    while not pred(hi):
        hi *= 2
    lo = hi // 2 + 1 if hi > 1 else 1
    while lo < hi:
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid + 1
    return hi


def sec_x(results):
    section("(x) practical cycle counts: smallest N with W_N(R,delta)/R <= eps")
    rows = []
    print(f"{'R/delta':>8} {'eps':>7} {'exact N':>9} {'(R^2/4d^2)ln(1/eps)':>21} {'A3 N':>8} {'classical N':>12}")
    for ratio_ in (10, 100):
        for eps in ("1e-3", "1e-6"):
            e = D(eps)
            n_exact = smallest_N(lambda N: W_equal_over_R(ratio_, N) <= e)
            check(n_exact == 1 or W_equal_over_R(ratio_, n_exact - 1) > e, "minimality of N")
            lead = D(ratio_) ** 2 / 4 * (-e.ln())
            q = 1 - 1 / D(ratio_ + 1) ** 2
            # smallest N with q^{(2N-1)/2} <= eps
            n_class = int(((2 * e.ln() / q.ln()) + 1) / 2)
            while dsqrt(q) ** (2 * n_class - 1) > e:
                n_class += 1
            while n_class > 1 and dsqrt(q) ** (2 * n_class - 3) <= e:
                n_class -= 1
            check(W_equal_over_R(ratio_, n_class) <= dsqrt(q) ** (2 * n_class - 1), "exact <= classical bound")
            check(n_exact <= n_class, "exact count <= classical count")
            # Corollary A3 bound with initial distance <= R: kappa^(2N-1)/sqrt(1-kappa^2) <= eps
            kap = 1 - 2 / D(ratio_) ** 2
            bound_a3 = lambda N: kap ** (2 * N - 1) / dsqrt(1 - kap * kap)
            n_a3 = smallest_N(lambda N: bound_a3(N) <= e)
            check(n_a3 == 1 or bound_a3(n_a3 - 1) > e, "minimality of Cor. A3 N")
            check(W_equal_over_R(ratio_, n_a3) <= bound_a3(n_a3), "exact <= Cor. A3 bound")
            check(n_exact <= n_a3 <= n_class, "exact <= Cor. A3 <= classical count")
            rows.append({"R_over_delta": ratio_, "eps": eps, "exact_N": n_exact,
                         "leading_order_count": format(lead, ".1f"), "corollary_A3_N": n_a3,
                         "classical_N": n_class})
            print(f"{ratio_:>8} {eps:>7} {n_exact:>9} {format(lead, '.1f'):>21} {n_a3:>8} {n_class:>12}")
    for ratio_ in (5, 10, 20):
        for N in (1, 10, 100, 1000):
            q = 1 - 1 / D(ratio_ + 1) ** 2
            check(W_equal_over_R(ratio_, N) <= dsqrt(q) ** (2 * N - 1), f"exact <= classical, R/delta={ratio_}, N={N}")
    results["cycle_counts"] = rows


def sec_xi(results):
    section("(xi) two-radius value vs Theorem 1 with delta = min(a,b)")
    rows = []
    print(f"{'a':>4} {'b':>4} {'R':>2} {'m':>2} {'two-radius W_m':>22} {'Thm 1 (delta=0.2)':>22} {'ratio':>9}")
    for a, b in ((D("0.3"), D("0.2")), (D("0.2"), D("0.3"))):
        R = D(1)
        for m in (1, 3, 5):
            two, _ = W7(m, R, a, b)
            eq, _ = W7(m, R, min(a, b), min(a, b))
            check(two <= eq, "two-radius value <= Theorem 1 with delta=min(a,b)")
            rows.append({"a": str(a), "b": str(b), "R": "1", "m": m,
                         "two_radius": format(two, ".12f"), "theorem1_delta_min": format(eq, ".12f"),
                         "ratio": format(two / eq, ".6f")})
            print(f"{a:>4} {b:>4} {R:>2} {m:>2} {format(two, '.15f'):>22} {format(eq, '.15f'):>22} {format(two / eq, '.6f'):>9}")
    results["two_radius_vs_equal"] = rows


def main():
    out = output_dir()
    started = time.perf_counter()
    rng = random.Random(SEED)
    results = {"python": platform.python_version(), "seed": SEED, "decimal_precision": PREC,
               "decimal_tolerance": str(TOL), "search_failure_threshold": f"1 + {SEARCH_REL_TOL}",
               "arithmetic": ("Fraction (exact) for rational-angle states in (i)-(iii) and the (vi) "
                              "trajectories; Decimal with 60 significant digits elsewhere (rounded, "
                              "not interval arithmetic)."),
               "scope": ("Finite diagnostics of the two-radius extension; not a proof of the "
                         "universally quantified Theorems 7 and 8.")}
    sec_i_to_iv(rng, results)
    sec_v(rng, results)
    sec_vi(rng, results)
    sec_vii(rng, results)
    sec_viii(results)
    sec_ix(results)
    sec_x(results)
    sec_xi(results)
    results["checks_by_section"] = SECTION_COUNTS
    results["checks_evaluated"] = CHECKS["evaluated"]
    results["checks_failed"] = CHECKS["failed"]
    results["first_failures"] = FAILURES
    results["status"] = "PASS" if CHECKS["failed"] == 0 else "FAIL"
    results["python_optimize_flag"] = sys.flags.optimize
    results["manuscript_sha256"] = manuscript_sha256()
    results["elapsed_seconds"] = round(time.perf_counter() - started, 2)
    (out / "extension_checks.json").write_text(json.dumps(results, indent=2) + "\n")
    print("\nchecks by section:")
    for k, v in SECTION_COUNTS.items():
        print(f"  {v:>7}  {k}")
    print(f"extension_checks.json: {CHECKS['evaluated']} explicit checks, {CHECKS['failed']} failed; "
          f"status {results['status']}; {results['elapsed_seconds']} s")
    if CHECKS["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
