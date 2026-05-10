"""Streamlines for the backward-facing step at the final time step.

Renders pure LBM and k-eps side by side with streamfunction contours and
streamlines, with the solid step block masked out.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_DIR = ROOT_DIR / 'outputs' / 'sec4' / 'backward_step'
KEPS_DIR = ROOT_DIR / 'outputs' / 'sec4' / 'backward_step_keps'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)


def last_snapshot(directory, prefix):
    paths = sorted(directory.glob(f'{prefix}*.csv'),
                   key=lambda p: int(p.stem.replace(prefix, '')))
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
    u_arr = to_grid(df, 'u', NX, NY)
    v_arr = to_grid(df, 'v', NX, NY)
    solid = to_grid(df, 'solid', NX, NY)
    psi_masked = np.ma.array(psi, mask=(solid > 0))
    u_masked = np.ma.array(u_arr, mask=(solid > 0))
    v_masked = np.ma.array(v_arr, mask=(solid > 0))

    extent_psi = float(np.abs(psi_masked).max() + 1e-12)
    im = ax.contourf(psi_masked, levels=20, cmap='RdBu_r',
                     vmin=-extent_psi, vmax=extent_psi)
    ys, xs = np.mgrid[0:NY, 0:NX]
    ax.streamplot(xs, ys, u_masked, v_masked, density=1.4, color='k', linewidth=0.5,
                  arrowsize=0.6)
    # overlay step block
    ax.fill_between([0, 30], 0, 30, color='gray', alpha=0.85, zorder=5)
    ax.set_xlim(0, NX-1); ax.set_ylim(0, NY-1)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')

    # Find reattachment
    u_bottom = u_arr[1, :]   # y=1
    reattach = -1
    for x in range(30, NX):
        if u_bottom[x] > 0:
            reattach = x; break
    xR_over_H = (reattach - 30) / 30.0 if reattach > 0 else float('nan')
    ax.axvline(reattach, ls='--', color='lime', lw=1.5, label=f'reattach $x_R={reattach}$')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_title(f'{title} — $x_R/H$ = {xR_over_H:.2f}')
    return im


fig, axes = plt.subplots(2, 1, figsize=(13, 7))
im_pure = render(axes[0], last_snapshot(PURE_DIR, 'step_snapshot_'), 'Pure LBM')
im_keps = render(axes[1], last_snapshot(KEPS_DIR, 'step_keps_snapshot_'), 'LBM + $k$-$\\varepsilon$')

cbar = fig.colorbar(im_pure, ax=axes.ravel().tolist(), orientation='horizontal',
                    fraction=0.04, pad=0.10, label='streamfunction $\\psi$')
plt.suptitle('後方ステップ流れ — 流線と再循環ゾーン（最終スナップショット）')

fname = 'backward_step_streamlines.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200, bbox_inches='tight')
print(f'Saved streamlines plot to {OUTPUT_DIR / fname}')
print(f'Saved streamlines plot to {PUBLISHED_ASSET_DIR / fname}')
