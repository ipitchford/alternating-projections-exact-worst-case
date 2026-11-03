"""Figure 1: the extremal two-halfspace geometry (equal radii, delta=1).

Plotting dependencies: numpy, matplotlib. Imported by plot_exact_results.py,
which regenerates all three figures.
"""
from pathlib import Path
import os
import csv
import math
import numpy as np
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ap-matplotlib-config')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
FIG.mkdir(exist_ok=True)
PLUM, OCHRE, TEAL, INK = '#72466e', '#ad7628', '#267c78', '#24333d'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.labelcolor': INK, 'text.color': INK,
                     'axes.edgecolor': '#9aa3a8', 'xtick.color': INK,
                     'ytick.color': INK, 'figure.facecolor': 'white',
                     'axes.spines.top': False, 'axes.spines.right': False})

def geometry():
    fig, ax = plt.subplots(figsize=(7.4, 4.25), constrained_layout=True)
    xx = np.linspace(-1.6, 3.7, 800)
    # A: y <= 1; B: .8*x - .6*y <= 1.
    lower_b = (.8 * xx - 1) / .6
    ax.fill_between(xx, lower_b, 1, where=(lower_b <= 1),
                    color=TEAL, alpha=.12, label=r'$A\cap B$')
    ax.plot(xx, np.ones_like(xx), color=PLUM, lw=1.7, label=r'$\partial A$')
    ax.plot(xx, lower_b, color=OCHRE, lw=1.7, label=r'$\partial B$')
    ax.add_patch(Circle((0, 0), 1, facecolor=TEAL, edgecolor=TEAL, alpha=.2))
    ax.plot(0, 0, 'o', color=TEAL, ms=4)
    ax.text(-.05, -.23, r'$q=0$', ha='center', fontsize=10)
    ax.annotate(r'$\delta=1$', xy=(0, 1), xytext=(-.72, .42),
                arrowprops={'arrowstyle': '-', 'color': TEAL}, color=TEAL)
    z = np.array([2., 1.])
    points = [np.array([3., 1.])]
    for j in range(1, 8):
        v = np.array([.6, .8]) if j % 2 else np.array([1., 0.])
        points.append(z + .6**j * v)
    pp = np.array(points)
    ax.plot(pp[:, 0], pp[:, 1], '-o', color=INK, lw=1.3, ms=3.3,
            label='Alternating projections', zorder=6)
    ax.scatter(*z, color=TEAL, marker='D', s=28, zorder=7)
    ax.annotate(r'$z=(2,1)$', xy=z, xytext=(1.48, .55),
                arrowprops={'arrowstyle': '-', 'color': TEAL}, color=TEAL)
    ax.annotate(r'$y_0=x_0$', xy=points[0], xytext=(3.03, .64))
    ax.annotate(r'$y_1$', xy=points[1], xytext=(2.46, 1.61))
    ax.annotate(r'$y_2$', xy=points[2], xytext=(2.41, .77))
    ax.set(xlim=(-1.4, 3.65), ylim=(-1.35, 2.05), xlabel=r'$x$', ylabel=r'$y$')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='upper left', frameon=False, fontsize=8.8)
    fig.savefig(FIG / 'geometry.png', dpi=220)
    plt.close(fig)

if __name__ == '__main__':
    geometry()
    print('Saved geometry.png')
