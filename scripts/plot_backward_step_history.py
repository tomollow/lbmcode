"""Reattachment-length history of the backward-facing step.

Shows how the reattachment point x_R drifts downstream as the recirculation
zone develops. Compares pure LBM and k-eps overlay.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'backward_step' / 'step_history.csv'
KEPS_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'backward_step_keps' / 'step_keps_history.csv'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'

pure = np.loadtxt(PURE_CSV, delimiter=',', skiprows=1)
keps = np.loadtxt(KEPS_CSV, delimiter=',', skiprows=1)

STEP_LENGTH = 30
STEP_HEIGHT = 30
nu0 = (0.55 - 0.5) / 3.0

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

LEGEND_KW = dict(loc='upper center', bbox_to_anchor=(0.5, -0.13),
                 ncol=2, frameon=True, fontsize=11,
                 handlelength=3.0, handletextpad=0.8, columnspacing=1.5)

# Left: x_R/H vs time
ax = axes[0]
ax.plot(pure[:, 0], (pure[:, 2] - STEP_LENGTH) / STEP_HEIGHT, lw=2.2, color='C0', label='Pure LBM')
ax.plot(keps[:, 0], (keps[:, 2] - STEP_LENGTH) / STEP_HEIGHT, lw=2.2, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step'); ax.set_ylabel('$x_R / H$')
ax.set_title('再循環ゾーン長さの推移')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

# Right: u_max vs time
ax = axes[1]
ax.plot(pure[:, 0], pure[:, 1], lw=2.2, color='C0', label='Pure LBM')
ax.plot(keps[:, 0], keps[:, 1], lw=2.2, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step'); ax.set_ylabel('$|u|_{\\max}$')
ax.set_title('最大流速の推移')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

plt.suptitle('後方ステップ流れ — 再循環長と流速の発展（NX=240, NY=60, $\\tau$=0.55）')
plt.tight_layout()

fname = 'backward_step_history.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved history plot to {OUTPUT_DIR / fname}')
print(f'Saved history plot to {PUBLISHED_ASSET_DIR / fname}')
print(f'Pure final x_R/H  = {(pure[-1, 2] - STEP_LENGTH)/STEP_HEIGHT:.2f}')
print(f'k-eps final x_R/H = {(keps[-1, 2] - STEP_LENGTH)/STEP_HEIGHT:.2f}')
