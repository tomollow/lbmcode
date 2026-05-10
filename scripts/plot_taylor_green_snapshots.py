"""Render vorticity snapshots from Taylor-Green CSV files as a 2x3 panel.

Usage: python plot_taylor_green_snapshots.py [pure|keps]   (default: pure)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

variant = sys.argv[1] if len(sys.argv) > 1 else 'pure'
if variant not in ('pure', 'keps'):
    raise SystemExit(f'unknown variant: {variant}')

if variant == 'pure':
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'taylor_green'
    prefix = 'tg_snapshot_'
    suptitle_label = 'Pure LBM'
    out_fname = 'taylor_green_snapshots.png'
else:
    src_dir = ROOT_DIR / 'outputs' / 'sec4' / 'taylor_green_keps'
    prefix = 'tg_keps_snapshot_'
    suptitle_label = 'LBM + $k$-$\\varepsilon$'
    out_fname = 'taylor_green_snapshots_keps.png'

snapshots = sorted(src_dir.glob(f'{prefix}*.csv'))
if not snapshots:
    raise SystemExit(f'no snapshots found in {src_dir}')

# Parse step from filename
def step_of(p):
    return int(p.stem.replace(prefix, ''))

snapshots = sorted(snapshots, key=step_of)
# Pick up to 6 panels evenly distributed
if len(snapshots) > 6:
    idx = np.linspace(0, len(snapshots) - 1, 6).astype(int)
    snapshots = [snapshots[i] for i in idx]

# Determine common color range across all snapshots for the vorticity field
all_vmax = 0.0
arrays = []
for path in snapshots:
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    omega = df.pivot(index='y', columns='x', values='vorticity').reindex(
        index=range(NY), columns=range(NX)
    ).values
    arrays.append((step_of(path), omega))
    all_vmax = max(all_vmax, float(np.abs(omega).max()))

# Symmetric colormap (RdBu) from -vmax to +vmax based on first frame for stable colors
vmax_init = float(np.abs(arrays[0][1]).max())

n = len(arrays)
ncols = 3
nrows = (n + ncols - 1) // ncols
fig, axs = plt.subplots(nrows, ncols, figsize=(13, 4.2 * nrows))
axs = np.array(axs).reshape(-1)

for ax, (step, omega) in zip(axs, arrays):
    im = ax.imshow(omega, origin='lower', cmap='RdBu_r',
                   vmin=-vmax_init, vmax=vmax_init, aspect='equal')
    ax.set_title(f'step = {step}, $|\\omega|_\\max$={np.abs(omega).max():.3e}')
    ax.set_xlabel('x'); ax.set_ylabel('y')
for ax in axs[len(arrays):]:
    ax.axis('off')

# Single horizontal colorbar at bottom
cbar = fig.colorbar(im, ax=axs.tolist(), orientation='horizontal',
                    fraction=0.04, pad=0.08, label='vorticity $\\omega = \\partial_x v - \\partial_y u$')

plt.suptitle(f'Taylor-Green 渦 渦度フィールド（{suptitle_label}）— カラースケールは初期 $|\\omega|_\\max$ で固定')

fname = out_fname
plt.savefig(OUTPUT_DIR / fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200, bbox_inches='tight')
print(f'Saved snapshots ({variant}) to {OUTPUT_DIR / fname}')
print(f'Saved snapshots ({variant}) to {PUBLISHED_ASSET_DIR / fname}')
