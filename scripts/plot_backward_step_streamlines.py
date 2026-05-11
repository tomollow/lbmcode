"""Streamlines and streamfunction for the backward-facing step.

Without args, renders pure LBM and k-eps side by side (top/bottom) as
`backward_step_streamlines.png`. With a variant arg `[pure|keps|les]`, renders
a single-variant panel (streamfunction + reattach line, plus a $\\nu_t$ field
panel for variants that compute it) as `backward_step_streamlines_<variant>.png`.

Usage:
    python plot_backward_step_streamlines.py            # 2-panel compare
    python plot_backward_step_streamlines.py keps       # single keps
    python plot_backward_step_streamlines.py les        # single les
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

VARIANTS = {
    'pure': dict(src='backward_step',      prefix='step_snapshot_',      label='Pure LBM'),
    'keps': dict(src='backward_step_keps', prefix='step_keps_snapshot_', label='LBM + $k$-$\\varepsilon$'),
    'les':  dict(src='backward_step_les',  prefix='step_les_snapshot_',  label='LBM + Smagorinsky LES'),
}


def last_snapshot(name, prefix):
    directory = ROOT_DIR / 'outputs' / 'sec4' / name
    paths = sorted(directory.glob(f'{prefix}*.csv'),
                   key=lambda p: int(p.stem.replace(prefix, '')))
    if not paths:
        raise SystemExit(f'no snapshots in {directory}')
    return paths[-1]


def to_grid(df, column, NX, NY):
    return df.pivot(index='y', columns='x', values=column).reindex(
        index=range(NY), columns=range(NX)
    ).values


def find_reattach(u_arr):
    NX = u_arr.shape[1]
    u_bottom = u_arr[1, :]  # y=1
    for x in range(30, NX):
        if u_bottom[x] > 0:
            return x
    return -1


def render_streamlines(ax, path, title):
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
    ax.fill_between([0, 30], 0, 30, color='gray', alpha=0.85, zorder=5)
    reattach = find_reattach(u_arr)
    xR_over_H = (reattach - 30) / 30.0 if reattach > 0 else float('nan')
    ax.axvline(reattach, ls='--', color='lime', lw=1.5, label=f'reattach $x_R={reattach}$')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(0, NX-1); ax.set_ylim(0, NY-1); ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'{title} — $x_R/H$ = {xR_over_H:.2f}')
    return im, df, NX, NY, solid


variant = sys.argv[1] if len(sys.argv) > 1 else None

if variant is None:
    # Default: 2-panel pure-vs-keps comparison
    fig, axes = plt.subplots(2, 1, figsize=(13, 7))
    im, _, _, _, _ = render_streamlines(axes[0], last_snapshot(VARIANTS['pure']['src'], VARIANTS['pure']['prefix']), VARIANTS['pure']['label'])
    render_streamlines(axes[1], last_snapshot(VARIANTS['keps']['src'], VARIANTS['keps']['prefix']), VARIANTS['keps']['label'])
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
                 fraction=0.04, pad=0.10, label='streamfunction $\\psi$')
    plt.suptitle('後方ステップ流れ — 流線と再循環ゾーン（最終スナップショット）')
    out_fname = 'backward_step_streamlines.png'
elif variant in VARIANTS:
    cfg = VARIANTS[variant]
    path = last_snapshot(cfg['src'], cfg['prefix'])
    df_peek = pd.read_csv(path, nrows=1)
    has_nut = 'nut' in df_peek.columns

    fig, axes = plt.subplots(2 if has_nut else 1, 1, figsize=(13, 6.5 if has_nut else 3.5))
    axes = np.atleast_1d(axes)
    im_psi, df, NX, NY, solid = render_streamlines(axes[0], path, cfg['label'])
    plt.colorbar(im_psi, ax=axes[0], orientation='vertical',
                 fraction=0.03, pad=0.02, label='streamfunction $\\psi$')

    if has_nut:
        nut = to_grid(df, 'nut', NX, NY)
        nut_masked = np.ma.array(np.maximum(nut, 1e-12), mask=(solid > 0))
        floor = max(float(nut_masked.min()), 1e-8)
        ceil  = max(float(nut_masked.max()), floor * 10)
        im_nut = axes[1].imshow(nut_masked, origin='lower', cmap='viridis',
                                norm=matplotlib.colors.LogNorm(vmin=floor, vmax=ceil))
        axes[1].fill_between([0, 30], 0, 30, color='gray', alpha=0.85, zorder=5)
        axes[1].set_aspect('equal'); axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
        axes[1].set_title(f"{cfg['label']} — $\\nu_t$ field (log scale, max = {nut.max():.2e})")
        plt.colorbar(im_nut, ax=axes[1], orientation='vertical',
                     fraction=0.03, pad=0.02, label='$\\nu_t$ (lattice units)')

    plt.suptitle(f"後方ステップ流れ — 流線と $\\nu_t$ 場（{cfg['label']}, $Re_H \\approx 56$）")
    out_fname = f'backward_step_streamlines_{variant}.png'
else:
    raise SystemExit(f'unknown variant: {variant}')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / out_fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / out_fname, dpi=200, bbox_inches='tight')
print(f'Saved {out_fname} to {OUTPUT_DIR} and {PUBLISHED_ASSET_DIR}')
