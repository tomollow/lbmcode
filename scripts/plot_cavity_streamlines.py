"""Streamlines and streamfunction contours of the lid-driven cavity at the final step.

Loads the last snapshot from outputs/sec4/cavity{,_keps}/ and renders a 2-panel
comparison of pure LBM vs k-eps with overlaid streamlines and color-mapped
streamfunction.
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


def last_snapshot(directory, prefix):
    paths = sorted(directory.glob(f'{prefix}*.csv'),
                   key=lambda p: int(p.stem.replace(prefix, '')))
    if not paths:
        raise SystemExit(f'no snapshots in {directory}')
    return paths[-1]


def to_grid(df, column, NX, NY):
    return df.pivot(index='y', columns='x', values=column).reindex(
        index=range(NY), columns=range(NX)
    ).values


def render(ax, path, title):
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    psi = to_grid(df, 'psi', NX, NY)
    u = to_grid(df, 'u', NX, NY)
    vy = to_grid(df, 'v', NX, NY)

    psi_min = float(psi.min())
    extent = max(abs(psi_min), 1e-12)
    im = ax.contourf(psi, levels=20, cmap='RdBu_r', vmin=-extent, vmax=extent)
    # streamlines
    ys, xs = np.mgrid[0:NY, 0:NX]
    ax.streamplot(xs, ys, u, vy, density=1.4, color='k', linewidth=0.6, arrowsize=0.7)
    ax.set_xlim(0, NX-1)
    ax.set_ylim(0, NY-1)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'{title}\n($\\psi_\\min = {psi_min:.4f}$, normalized: ${psi_min/(0.05*NX):.4f}$)')
    return im


fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
im_pure = render(axes[0], last_snapshot(PURE_DIR, 'cavity_snapshot_'), 'Pure LBM')
im_keps = render(axes[1], last_snapshot(KEPS_DIR, 'cavity_keps_snapshot_'), 'LBM + $k$-$\\varepsilon$')

cbar = fig.colorbar(im_pure, ax=axes.ravel().tolist(), orientation='horizontal',
                    fraction=0.04, pad=0.10, label='streamfunction $\\psi$ (lattice units)')
plt.suptitle('蓋駆動キャビティ流れ — 流線と流線関数（最終スナップショット, $Re \\approx 384$）')

fname = 'cavity_streamlines.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200, bbox_inches='tight')
print(f'Saved streamline plot to {OUTPUT_DIR / fname}')
print(f'Saved streamline plot to {PUBLISHED_ASSET_DIR / fname}')
