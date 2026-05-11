"""Streamlines and streamfunction contours for the lid-driven cavity.

Without args, renders a side-by-side comparison of pure LBM vs k-eps as the
default `cavity_streamlines.png`. With a variant arg `[pure|keps|les]`, renders
a single-variant panel (streamfunction + streamlines, plus a $\\nu_t$ field
panel for variants that compute it) as `cavity_streamlines_<variant>.png`.

Usage:
    python plot_cavity_streamlines.py            # 2-panel compare
    python plot_cavity_streamlines.py keps       # single keps variant
    python plot_cavity_streamlines.py les        # single les variant
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
    'pure': dict(src='cavity',      prefix='cavity_snapshot_',      label='Pure LBM'),
    'keps': dict(src='cavity_keps', prefix='cavity_keps_snapshot_', label='LBM + $k$-$\\varepsilon$'),
    'les':  dict(src='cavity_les',  prefix='cavity_les_snapshot_',  label='LBM + Smagorinsky LES'),
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


def render_streamlines(ax, path, title):
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    psi = to_grid(df, 'psi', NX, NY)
    u = to_grid(df, 'u', NX, NY)
    vy = to_grid(df, 'v', NX, NY)
    psi_min = float(psi.min())
    extent = max(abs(psi_min), 1e-12)
    im = ax.contourf(psi, levels=20, cmap='RdBu_r', vmin=-extent, vmax=extent)
    ys, xs = np.mgrid[0:NY, 0:NX]
    ax.streamplot(xs, ys, u, vy, density=1.4, color='k', linewidth=0.6, arrowsize=0.7)
    ax.set_xlim(0, NX-1); ax.set_ylim(0, NY-1); ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'{title}\n($\\psi_\\min = {psi_min:.4f}$, normalized: ${psi_min/(0.05*NX):.4f}$)')
    return im, df, NX, NY


variant = sys.argv[1] if len(sys.argv) > 1 else None

if variant is None:
    # Default: 2-panel pure-vs-keps comparison (back-compat)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    im, _, _, _ = render_streamlines(axes[0], last_snapshot(VARIANTS['pure']['src'], VARIANTS['pure']['prefix']), VARIANTS['pure']['label'])
    render_streamlines(axes[1], last_snapshot(VARIANTS['keps']['src'], VARIANTS['keps']['prefix']), VARIANTS['keps']['label'])
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal',
                 fraction=0.04, pad=0.10, label='streamfunction $\\psi$ (lattice units)')
    plt.suptitle('蓋駆動キャビティ流れ — 流線と流線関数（最終スナップショット, $Re \\approx 384$）')
    out_fname = 'cavity_streamlines.png'
elif variant in VARIANTS:
    cfg = VARIANTS[variant]
    path = last_snapshot(cfg['src'], cfg['prefix'])
    df_peek = pd.read_csv(path, nrows=1)
    has_nut = 'nut' in df_peek.columns

    fig, axes = plt.subplots(1, 2 if has_nut else 1, figsize=(13 if has_nut else 7, 5.8))
    axes = np.atleast_1d(axes)
    im_psi, df, NX, NY = render_streamlines(axes[0], path, cfg['label'])
    plt.colorbar(im_psi, ax=axes[0], orientation='horizontal',
                 fraction=0.05, pad=0.12, label='streamfunction $\\psi$')

    if has_nut:
        nut = to_grid(df, 'nut', NX, NY)
        nut_safe = np.maximum(nut, 1e-12)
        floor = max(float(nut_safe.min()), 1e-8)
        ceil  = max(float(nut_safe.max()), floor * 10)
        im_nut = axes[1].imshow(nut_safe, origin='lower', cmap='viridis',
                                norm=matplotlib.colors.LogNorm(vmin=floor, vmax=ceil))
        axes[1].set_aspect('equal'); axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
        axes[1].set_title(f"{cfg['label']}\n$\\nu_t$ field (log scale, max = {nut.max():.2e})")
        plt.colorbar(im_nut, ax=axes[1], orientation='horizontal',
                     fraction=0.05, pad=0.12, label='$\\nu_t$ (lattice units)')

    plt.suptitle(f"蓋駆動キャビティ流れ — 流線と $\\nu_t$ 場（{cfg['label']}, $Re \\approx 384$）")
    out_fname = f'cavity_streamlines_{variant}.png'
else:
    raise SystemExit(f'unknown variant: {variant}')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / out_fname, dpi=200, bbox_inches='tight')
plt.savefig(PUBLISHED_ASSET_DIR / out_fname, dpi=200, bbox_inches='tight')
print(f'Saved {out_fname} to {OUTPUT_DIR} and {PUBLISHED_ASSET_DIR}')
