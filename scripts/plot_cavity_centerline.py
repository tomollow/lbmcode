"""Centerline u(y) and v(x) profiles for cavity flow with Ghia (1982) reference.

Ghia data is hardcoded from Ghia, Ghia, Shin (1982) Table I/II for Re=400.
We compare against Re=400 even though our simulation is at Re~384 (the
parameters in the C source give Re = U_LID*NX/nu_0 = 384). The agreement is
qualitative — exact match would require U_LID=0.0521 (Re=400 exactly) and
longer convergence.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_DIR = ROOT_DIR / 'outputs' / 'sec4' / 'cavity'
KEPS_DIR = ROOT_DIR / 'outputs' / 'sec4' / 'cavity_keps'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

# Ghia, Ghia, Shin (1982) Re=400 reference (u along vertical centerline, v along horizontal)
GHIA_RE400_U = np.array([
    # y/L,  u/U_lid
    (0.0000,  0.00000),
    (0.0547, -0.08186),
    (0.0625, -0.09266),
    (0.0703, -0.10338),
    (0.1016, -0.14612),
    (0.1719, -0.24299),
    (0.2813, -0.32726),
    (0.4531, -0.17119),
    (0.5000, -0.11477),
    (0.6172, +0.02135),
    (0.7344, +0.16256),
    (0.8516, +0.29093),
    (0.9531, +0.55892),
    (0.9609, +0.61756),
    (0.9688, +0.68439),
    (0.9766, +0.75837),
    (1.0000, +1.00000),
])

GHIA_RE400_V = np.array([
    # x/L,  v/U_lid
    (0.0000, +0.00000),
    (0.0625, +0.18360),
    (0.0703, +0.19713),
    (0.0781, +0.20920),
    (0.0938, +0.22965),
    (0.1563, +0.28124),
    (0.2266, +0.30203),
    (0.2344, +0.30174),
    (0.5000, +0.05186),
    (0.8047, -0.38598),
    (0.8594, -0.44993),
    (0.9063, -0.33827),
    (0.9453, -0.22847),
    (0.9531, -0.19254),
    (0.9609, -0.15663),
    (0.9688, -0.12146),
    (1.0000, +0.00000),
])


def last_snapshot(directory, prefix):
    paths = sorted(directory.glob(f'{prefix}*.csv'),
                   key=lambda p: int(p.stem.replace(prefix, '')))
    return paths[-1]


def centerline_profiles(path, NX_expected=128, NY_expected=128):
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    # u along vertical centerline (x = NX/2)
    cu = df[df['x'] == NX // 2].sort_values('y')
    # v along horizontal centerline (y = NY/2)
    cv = df[df['y'] == NY // 2].sort_values('x')
    return NX, NY, cu['y'].values, cu['u'].values, cv['x'].values, cv['v'].values


U_LID = 0.05

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
LEGEND_KW = dict(loc='upper center', bbox_to_anchor=(0.5, -0.13),
                 ncol=3, frameon=True, fontsize=11,
                 handlelength=3.0, handletextpad=0.8, columnspacing=1.5)

# Left: u along vertical centerline
ax = axes[0]
NX, NY, y_pure, u_pure, _, _ = centerline_profiles(last_snapshot(PURE_DIR, 'cavity_snapshot_'))
_, _, y_keps, u_keps, _, _ = centerline_profiles(last_snapshot(KEPS_DIR, 'cavity_keps_snapshot_'))
ax.plot(u_pure / U_LID, y_pure / (NY-1), lw=2.0, color='C0', label='Pure LBM')
ax.plot(u_keps / U_LID, y_keps / (NY-1), '--', lw=2.0, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.plot(GHIA_RE400_U[:, 1], GHIA_RE400_U[:, 0], 'o', ms=5, color='k', label='Ghia (1982) Re=400')
ax.set_xlabel('$u / U_{\\mathrm{lid}}$')
ax.set_ylabel('$y / L$')
ax.set_title('垂直中心線 ($x = NX/2$) の $u$ プロファイル')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

# Right: v along horizontal centerline
ax = axes[1]
NX, NY, _, _, x_pure, v_pure = centerline_profiles(last_snapshot(PURE_DIR, 'cavity_snapshot_'))
_, _, _, _, x_keps, v_keps = centerline_profiles(last_snapshot(KEPS_DIR, 'cavity_keps_snapshot_'))
ax.plot(x_pure / (NX-1), v_pure / U_LID, lw=2.0, color='C0', label='Pure LBM')
ax.plot(x_keps / (NX-1), v_keps / U_LID, '--', lw=2.0, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.plot(GHIA_RE400_V[:, 0], GHIA_RE400_V[:, 1], 'o', ms=5, color='k', label='Ghia (1982) Re=400')
ax.set_xlabel('$x / L$')
ax.set_ylabel('$v / U_{\\mathrm{lid}}$')
ax.set_title('水平中心線 ($y = NY/2$) の $v$ プロファイル')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

plt.suptitle(f'蓋駆動キャビティ Ghia 比較 (NX={NX}, NY={NY}, $\\tau=0.55$, $Re \\approx 384$)')
plt.tight_layout()

fname = 'cavity_centerline.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved centerline plot to {OUTPUT_DIR / fname}')
print(f'Saved centerline plot to {PUBLISHED_ASSET_DIR / fname}')
