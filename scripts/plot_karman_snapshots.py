"""Vorticity snapshots of the Karman vortex street.

Renders 6 panels showing the wake evolution from start to vortex shedding.
Usage: python plot_karman_snapshots.py [pure|keps]   (default: pure)
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
if variant not in ('pure', 'keps'):
    raise SystemExit(f'unknown variant: {variant}')

if variant == 'pure':
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'karman'
    prefix = 'karman_snapshot_'
    suptitle_label = 'Pure LBM'
    out_fname = 'karman_snapshots.png'
else:
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'karman_keps'
    prefix = 'karman_keps_snapshot_'
    suptitle_label = 'LBM + $k$-$\\varepsilon$'
    out_fname = 'karman_snapshots_keps.png'

snapshots = sorted(src_dir.glob(f'{prefix}*.csv'),
                   key=lambda p: int(p.stem.replace(prefix, '')))
if not snapshots:
    raise SystemExit(f'no snapshots in {src_dir}')

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
    solid = df.pivot(index='y', columns='x', values='solid').reindex(
        index=range(NY), columns=range(NX)
    ).values
    omega_masked = np.ma.array(omega, mask=(solid > 0))
    step = int(path.stem.replace(prefix, ''))
    arrays.append((step, omega_masked))

# Anchor color scale on the last frame where shedding is fully developed
last_omega = arrays[-1][1]
vmax = float(np.percentile(np.abs(last_omega.compressed()), 99))

n = len(arrays)
fig, axs = plt.subplots(n, 1, figsize=(13, 3.0 * n))
axs = np.array(axs).reshape(-1)

for ax, (step, omega) in zip(axs, arrays):
    im = ax.imshow(omega, origin='lower', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax, aspect='equal')
    ax.set_title(f'step = {step}, $|\\omega|_\\max$ = {np.abs(omega).max():.3e}',
                 pad=8, fontsize=11)
    ax.set_xlabel('x', fontsize=10)
    ax.set_ylabel('y', fontsize=10)

# Adjust spacing so panel titles don't collide with the panel above's x-tick labels.
fig.subplots_adjust(hspace=0.6, top=0.95, bottom=0.07)
cbar = fig.colorbar(im, ax=axs.tolist(), orientation='horizontal',
                    fraction=0.03, pad=0.06,
                    label='vorticity $\\omega = \\partial_x v - \\partial_y u$')
fig.suptitle(f'Kármán 渦列 渦度フィールド（{suptitle_label}）', y=0.99, fontsize=12)

plt.savefig(OUTPUT_DIR / out_fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / out_fname, dpi=200, bbox_inches='tight')
print(f'Saved snapshots ({variant}) to {OUTPUT_DIR / out_fname}')
print(f'Saved snapshots ({variant}) to {PUBLISHED_ASSET_DIR / out_fname}')
