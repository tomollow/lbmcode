"""Vorticity snapshots of the Kelvin-Helmholtz instability over time.

Usage: python plot_kelvin_helmholtz_snapshots.py [pure|keps|les]   (default: pure)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

variant = sys.argv[1] if len(sys.argv) > 1 else 'pure'
if variant not in ('pure', 'keps', 'les'):
    raise SystemExit(f'unknown variant: {variant}')

if variant == 'pure':
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'kelvin_helmholtz'
    prefix = 'kh_snapshot_'
    suptitle_label = 'Pure LBM'
    out_fname = 'kelvin_helmholtz_snapshots.png'
elif variant == 'keps':
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'kelvin_helmholtz_keps'
    prefix = 'kh_keps_snapshot_'
    suptitle_label = 'LBM + $k$-$\\varepsilon$'
    out_fname = 'kelvin_helmholtz_snapshots_keps.png'
else:  # les
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'kelvin_helmholtz_les'
    prefix = 'kh_les_snapshot_'
    suptitle_label = 'LBM + Smagorinsky LES'
    out_fname = 'kelvin_helmholtz_snapshots_les.png'

snapshots = sorted(src_dir.glob(f'{prefix}*.csv'))
if not snapshots:
    raise SystemExit(f'no snapshots found in {src_dir}')

def step_of(p):
    return int(p.stem.replace(prefix, ''))

snapshots = sorted(snapshots, key=step_of)
if len(snapshots) > 6:
    idx = np.linspace(0, len(snapshots) - 1, 6).astype(int)
    snapshots = [snapshots[i] for i in idx]

arrays = []
for path in snapshots:
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    omega = df.pivot(index='y', columns='x', values='vorticity').reindex(
        index=range(NY), columns=range(NX)
    ).values
    arrays.append((step_of(path), omega))

# For KH, anchor color scale on the late-stage frame (mid-table) where the
# vortices are saturated — using the initial sharp shear layer would put the
# whole roll-up phase off-scale.
mid_frame_omega = arrays[len(arrays) // 2][1]
vmax = float(np.percentile(np.abs(mid_frame_omega), 99))

n = len(arrays)
ncols = 2
nrows = (n + ncols - 1) // ncols
fig, axs = plt.subplots(nrows, ncols, figsize=(13, 3.5 * nrows))
axs = np.array(axs).reshape(-1)

for ax, (step, omega) in zip(axs, arrays):
    im = ax.imshow(omega, origin='lower', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax, aspect='equal')
    ax.set_title(f'step = {step}, $|\\omega|_\\max$ = {np.abs(omega).max():.3e}')
    ax.set_xlabel('x'); ax.set_ylabel('y')
for ax in axs[len(arrays):]:
    ax.axis('off')

cbar = fig.colorbar(im, ax=axs.tolist(), orientation='horizontal',
                    fraction=0.04, pad=0.08, label='vorticity $\\omega = \\partial_x v - \\partial_y u$')

plt.suptitle(f'Kelvin-Helmholtz 不安定 渦度フィールド（{suptitle_label}）— カラースケールは中盤フレーム99%タイルで固定')

fname = out_fname
plt.savefig(OUTPUT_DIR / fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200, bbox_inches='tight')
print(f'Saved snapshots ({variant}) to {OUTPUT_DIR / fname}')
print(f'Saved snapshots ({variant}) to {PUBLISHED_ASSET_DIR / fname}')
