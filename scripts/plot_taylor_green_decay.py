"""Decay history of Taylor-Green vortex amplitude vs analytical solution.

Plots the peak |u| versus time step for the pure-LBM and the k-eps-augmented
runs alongside the analytical decay rate u(t) = U0 * exp(-2 nu k^2 t).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'taylor_green' / 'tg_history.csv'
KEPS_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'taylor_green_keps' / 'tg_keps_history.csv'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

pure = np.loadtxt(PURE_CSV, delimiter=',', skiprows=1)
keps = np.loadtxt(KEPS_CSV, delimiter=',', skiprows=1)

# columns: pure = step,u_max,u_max_theory ; keps = step,u_max,u_max_pure_theory,k_mean,eps_mean,nut_mean

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

LEGEND_KW = dict(
    loc='upper center', bbox_to_anchor=(0.5, -0.13),
    ncol=3, frameon=True, fontsize=11,
    handlelength=3.0, handletextpad=0.8, columnspacing=1.5,
)

# Left: log-scale decay
ax = axes[0]
ax.semilogy(pure[:, 0], pure[:, 2], 'k--', lw=2.0, label='解析解 $U_0 e^{-2\\nu k^2 t}$')
ax.semilogy(pure[:, 0], pure[:, 1], lw=2.2, color='C0', label='Pure LBM')
ax.semilogy(keps[:, 0], keps[:, 1], lw=2.2, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step')
ax.set_ylabel('peak $|u|$')
ax.set_title('Taylor-Green 渦の減衰（対数スケール）')
ax.grid(True, which='both', ls=':')
ax.legend(**LEGEND_KW)

# Right: ratio (numerical / analytical) — proves Pure LBM matches, k-eps decays faster
ax = axes[1]
ax.plot(pure[:, 0], pure[:, 1]/pure[:, 2], lw=2.2, color='C0', label='Pure LBM / 解析')
ax.plot(keps[:, 0], keps[:, 1]/keps[:, 2], lw=2.2, color='C3', label='($k$-$\\varepsilon$) / 解析')
ax.axhline(1.0, ls=':', color='k', lw=1.0)
ax.set_xlabel('step')
ax.set_ylabel('数値解 / 解析解 (Pure LBM 基準)')
ax.set_title('解析解との比（k-εで減衰加速）')
ax.grid(True, ls=':')
ax.legend(**dict(LEGEND_KW, ncol=2))

NX = 128  # match C source
nu0 = (1.0 - 0.5) / 3.0
plt.suptitle(f'Taylor-Green 渦 NX={NX} NY={NX} $\\tau=1.0$ $\\nu_0={nu0:.4f}$ $U_0=0.04$')
plt.tight_layout()

fname = 'taylor_green_decay.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved decay plot to {OUTPUT_DIR / fname}')
print(f'Saved decay plot to {PUBLISHED_ASSET_DIR / fname}')
print(f'Pure LBM final ratio: {pure[-1, 1]/pure[-1, 2]:.4f}')
print(f'k-eps   final ratio: {keps[-1, 1]/keps[-1, 2]:.4f}')
