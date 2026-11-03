"""Regenerate the three manuscript figures and their plotted data.

Figure 1 (geometry.png): the extremal two-halfspace construction.
Figure 2 (exact_envelopes.png): for R/delta in {5, 10, 20}, on a log scale,
  - the exact worst-case error W_N(R,delta)/R of Theorem 1 (solid),
  - its sharp asymptotic equivalent R^4 kappa^{2N}/(8 e delta^2 sqrt(R^2-delta^2) N), over R (dashed),
  - the classical interior-point bound (1 - delta^2/(R+delta)^2)^{(2N-1)/2}
    (Gubin-Polyak-Raik type; d_0 <= R and Fejer step d_{j+1}^2 <= d_j^2 - t_{j+1}^2) (dotted).
Figure 3 (critical_scaling.png): the limiting profile L(lambda) and exact values for N=5, 50.

All plotted values are evaluated from the scalar formulas with 70-digit
decimal arithmetic (code/scalar.py); rendering uses binary floating point.
Usage: python3 -I -B code/plot_exact_results.py   (needs numpy, matplotlib)
"""
from pathlib import Path
import os
import sys
import csv
import math
import tempfile
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # needed under `python -I`, which omits the script directory
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'ap-matplotlib-config'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scalar import lower_scaled, exact_W_over_R, asymptotic_W_over_R, classical_W_over_R
from plot_results import geometry

ROOT = HERE.parent
FIG = ROOT / 'figures'
PLUM, OCHRE, TEAL, INK = '#72466e', '#ad7628', '#267c78', '#24333d'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelcolor': INK, 'text.color': INK, 'axes.edgecolor': '#9aa3a8',
    'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': 'white',
    'axes.spines.top': False, 'axes.spines.right': False})


def profile(lam):
    if lam == 0:
        return 1.0, 0.0
    lo, hi = 0., 1.
    for _ in range(80):
        v = (lo + hi) / 2
        if v**3 > 8 * lam * (1 - v):
            hi = v
        else:
            lo = v
    v = (lo + hi) / 2
    return (1 - v) * math.exp(-4 * lam / v**2), v


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def figure2():
    ratios = [(5, PLUM), (10, OCHRE), (20, TEAL)]
    Ns = [1, 2, 3, 5] + list(range(10, 1001, 10))
    rows = []
    fig, ax = plt.subplots(figsize=(7.4, 4.0), constrained_layout=True)
    for ratio, col in ratios:
        exact, asym, classical = [], [], []
        for N in Ns:
            e, a, c = (exact_W_over_R(ratio, N), asymptotic_W_over_R(ratio, N),
                       classical_W_over_R(ratio, N))
            exact.append(float(e)); asym.append(float(a)); classical.append(float(c))
            rows.append({'R_over_delta': ratio, 'N': N, 'exact_W_over_R': format(e, '.20e'),
                         'asymptotic_equivalent_over_R': format(a, '.20e'),
                         'classical_bound_over_R': format(c, '.20e')})
        ax.semilogy(Ns, exact, '-', color=col, lw=2.0, label=rf'$R/\delta={ratio}$')
        ax.semilogy(Ns, asym, '--', color=col, lw=1.1, alpha=.9)
        ax.semilogy(Ns, classical, ':', color=col, lw=1.6)
    # Style keys for the line types.
    ax.plot([], [], '-', color=INK, lw=2.0, label=r'exact $W_N(R,\delta)$ (Theorem 1)')
    ax.plot([], [], '--', color=INK, lw=1.1,
            label=r'asymptotic $R^4\kappa^{2N}/(8e\delta^2\sqrt{R^2-\delta^2}\,N)$')
    ax.plot([], [], ':', color=INK, lw=1.6,
            label=r'classical interior-point bound (Gubin–Polyak–Raik type)')
    ax.set(xlabel='Full cycles N', ylabel=r'Error bound / $R$', ylim=(1e-15, 2), xlim=(0, 1000))
    ax.legend(frameon=False, loc='lower left', fontsize=7.6, ncol=2)
    ax.grid(axis='y', alpha=.15)
    fig.savefig(FIG / 'exact_envelopes.png', dpi=220)
    plt.close(fig)
    write_csv(FIG / 'exact_envelopes_data.csv', rows)


def figure3():
    lambdas = np.logspace(-5, math.log10(2), 100)
    rows = []
    fig, ax = plt.subplots(figsize=(7.4, 3.65), constrained_layout=True)
    ax.semilogx(lambdas, [profile(float(l))[0] for l in lambdas],
                color=PLUM, lw=2.1, label='Limiting profile')
    for N, col, ls in [(5, OCHRE, '--'), (50, TEAL, '-.')]:
        vals = []
        for lam in lambdas:
            eps = math.sqrt(float(lam) / N)
            if eps >= 1 / math.sqrt(2):
                val = 0.
            else:
                k = 1 - 2 * eps * eps
                _, _, L = lower_scaled(str(k), N)
                val = eps * float(L)
            limiting, v = profile(float(lam))
            vals.append(val)
            rows.append({'N': N, 'lambda': format(float(lam), '.16g'),
                         'epsilon': format(eps, '.16g'), 'exact_W_over_R': format(val, '.16g'),
                         'limiting_profile': format(limiting, '.16g'),
                         'optimal_profile_v': format(v, '.16g')})
        ax.semilogx(lambdas, vals, color=col, ls=ls, lw=1.5, label=f'Exact value, N={N}')
    ax.set(xlabel=r'$\lambda=N(\delta/R)^2$',
           ylabel=r'Worst-case error / $R$', ylim=(0, 1.03))
    ax.legend(frameon=False, loc='lower left')
    ax.grid(axis='y', alpha=.15)
    fig.savefig(FIG / 'critical_scaling.png', dpi=220)
    plt.close(fig)
    write_csv(FIG / 'critical_scaling_data.csv', rows)


def main():
    FIG.mkdir(exist_ok=True)
    geometry()
    figure2()
    figure3()
    print('Saved geometry.png, exact_envelopes.png, critical_scaling.png and the plotted data CSVs in figures/.')


if __name__ == '__main__':
    main()
