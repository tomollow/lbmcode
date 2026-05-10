"""Growth-rate history of the Kelvin-Helmholtz instability.

Plots the peak |v| and the total turbulent kinetic energy proxy (u_rms) versus
time step for the pure-LBM and the k-eps-augmented runs.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'kelvin_helmholtz' / 'kh_history.csv'
KEPS_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'kelvin_helmholtz_keps' / 'kh_keps_history.csv'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

pure = np.loadtxt(PURE_CSV, delimiter=',', skiprows=1)  # step,v_max,u_rms
keps = np.loadtxt(KEPS_CSV, delimiter=',', skiprows=1)  # step,v_max,u_rms,k_mean,eps_mean,nut_mean

NX, NY = 256, 128
TAU = 0.52
nu0 = (TAU - 0.5) / 3.0

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

LEGEND_KW = dict(
    loc='upper center', bbox_to_anchor=(0.5, -0.13),
    ncol=2, frameon=True, fontsize=11,
    handlelength=3.0, handletextpad=0.8, columnspacing=1.5,
)

# Left: peak |v| (instability amplitude)
ax = axes[0]
ax.semilogy(pure[:, 0], pure[:, 1], lw=2.2, color='C0', label='Pure LBM')
ax.semilogy(keps[:, 0], keps[:, 1], lw=2.2, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step')
ax.set_ylabel('peak $|v|$')
ax.set_title('Kelvin-Helmholtz 摂動成長 (peak $|v|$)')
ax.grid(True, which='both', ls=':')
ax.legend(**LEGEND_KW)

# Right: u_rms decay (mean shear smoothing)
ax = axes[1]
ax.plot(pure[:, 0], pure[:, 2], lw=2.2, color='C0', label='Pure LBM')
ax.plot(keps[:, 0], keps[:, 2], lw=2.2, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step')
ax.set_ylabel('$u_{\\mathrm{rms}}$')
ax.set_title('平均せん断の平滑化 ($u_{\\mathrm{rms}}$)')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

plt.suptitle(f'Kelvin-Helmholtz 不安定性 NX={NX} NY={NY} $\\tau={TAU}$ $\\nu_0={nu0:.4f}$ $U_0=0.05$ $\\delta=4$')
plt.tight_layout()

fname = 'kelvin_helmholtz_growth.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved growth plot to {OUTPUT_DIR / fname}')
print(f'Saved growth plot to {PUBLISHED_ASSET_DIR / fname}')
print(f'Pure LBM: peak v_max = {pure[:, 1].max():.4e} at step {int(pure[pure[:, 1].argmax(), 0])}')
print(f'k-eps   : peak v_max = {keps[:, 1].max():.4e} at step {int(keps[keps[:, 1].argmax(), 0])}')
print(f'final v_max ratio (k-eps/pure) = {keps[-1, 1]/pure[-1, 1]:.4f}')
