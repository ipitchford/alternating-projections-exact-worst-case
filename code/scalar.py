"""High-precision evaluation of the exact equal-radius minimax value (Theorem 1).

Standard library only (decimal.Decimal). Used by plot_exact_results.py.
lower_scaled() is carried over unchanged in substance from the original
bundle's code/verify.py; its assertions are replaced by explicit exceptions so
the code behaves identically under `python -O`.
"""
from decimal import Decimal as D, localcontext


def lower_scaled(kappa, N, precision=70):
    """Return (W_N / asymptotic equivalent, maximising c, W_N) at delta=1.

    Here R^2 = 2/(1-kappa). The maximiser of c^m D(R,c), m=2N-1, is found by
    bisection on the critical quartic written in d = S - z (S=sqrt(R^2-1),
    z=sqrt((1+c)/(1-c))) to avoid loss of significance when z approaches S.
    """
    with localcontext() as ctx:
        ctx.prec = precision
        k = D(str(kappa))
        R2 = 2 / (1 - k)
        S = (R2 - 1).sqrt()
        m = 2 * N - 1

        def f(d):
            z = S - d
            return z**4 - 1 - 4 * m * z * d

        lo, hi = D(0), S - 1
        if not (f(lo) > 0 and f(hi) < 0):
            raise ArithmeticError('critical quartic is not bracketed')
        for _ in range(260):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        d = (lo + hi) / 2
        z = S - d
        c = (z * z - 1) / (z * z + 1)
        L = d * c**m
        C = R2 * R2 / (4 * S)
        asym = C * k**(2 * N) / (D(2) * D(1).exp() * N)
        residual = z**4 + 4 * m * z * z - 4 * m * S * z - 1
        if not abs(residual) < D('1e-55') * (1 + 4 * m * S * z):
            raise ArithmeticError('critical quartic residual too large')
        return +(L / asym), +c, +L


def exact_W_over_R(ratio, N, precision=70):
    """W_N(R,delta)/R for R/delta = ratio (> sqrt 2), via lower_scaled."""
    with localcontext() as ctx:
        ctx.prec = precision
        rho = D(str(ratio))
        if rho * rho <= 2:
            return D(0)
        kappa = 1 - 2 / (rho * rho)
        _, _, L = lower_scaled(str(kappa), N, precision)
        return L / rho  # delta = 1, so R = rho


def asymptotic_W_over_R(ratio, N, precision=70):
    """Sharp equivalent R^4 kappa^{2N} / (8 e delta^2 sqrt(R^2-delta^2) N), divided by R."""
    with localcontext() as ctx:
        ctx.prec = precision
        rho = D(str(ratio))
        kappa = 1 - 2 / (rho * rho)
        return rho**3 * kappa**(2 * N) / (8 * D(1).exp() * (rho * rho - 1).sqrt() * N)


def classical_W_over_R(ratio, N, precision=70):
    """Classical interior-point bound R (1 - delta^2/(R+delta)^2)^{(2N-1)/2}, divided by R."""
    with localcontext() as ctx:
        ctx.prec = precision
        rho = D(str(ratio))
        q = 1 - 1 / (rho + 1)**2
        return q.sqrt()**(2 * N - 1)
