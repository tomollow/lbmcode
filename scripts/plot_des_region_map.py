"""SA-DES length-scale switch visualization for sec4 DES cases.

For each DES case (cavity, karman, backward_step), reads the final snapshot
CSV and renders a panel showing:
  - d_wall field (left): wall distance in lattice units
  - RANS/LES region map (right): cells where d_wall < C_DES*Delta=0.65 are
    in SA-RANS mode (red); cells where d_wall >= 0.65 are in the LES branch
    (blue). This is the DES97 hybrid switch geometry.

Usage:
    python plot_des_region_map.py            # render all 3 cases
    python plot_des_region_map.py cavity     # single case
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

C_DES = 0.65
DELTA_LES = 1.0
SWITCH = C_DES * DELTA_LES

CASES = {
    'cavity': dict(src='cavity_des', prefix='cavity_des_snapshot_',
                   label='Cavity ($Re \\approx 384$)', has_solid=False,
                   figsize=(11, 4.5)),
    'karman': dict(src='karman_des', prefix='karman_des_snapshot_',
                   label='Karman ($Re_D \\approx 127$)', has_solid=True,
                   figsize=(13, 3.5)),
    'step':   dict(src='backward_step_des', prefix='step_des_snapshot_',
                   label='Backward Step ($Re_H \\approx 56$)', has_solid=True,
                   figsize=(13, 3.5)),
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


def render_case(name, cfg):
    path = last_snapshot(cfg['src'], cfg['prefix'])
    df = pd.read_csv(path)
    NX = df['x'].max() + 1
    NY = df['y'].max() + 1
    d_wall = to_grid(df, 'd_wall', NX, NY)
    if cfg['has_solid']:
        solid = to_grid(df, 'solid', NX, NY).astype(bool)
        d_wall_m = np.ma.array(d_wall, mask=solid)
    else:
        solid = None
        d_wall_m = d_wall

    fig, axes = plt.subplots(1, 2, figsize=cfg['figsize'])

    # Left: d_wall field (log-ish to highlight RANS layer)
    vmax_d = float(np.percentile(d_wall_m.compressed() if hasattr(d_wall_m, 'compressed')
                                  else d_wall_m, 99))
    im_d = axes[0].imshow(d_wall_m, origin='lower', cmap='viridis',
                          vmin=0.5, vmax=max(vmax_d, 2.0))
    axes[0].set_title(f"{cfg['label']} — $d_\\text{{wall}}$ (LU)")
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    axes[0].set_aspect('equal')
    plt.colorbar(im_d, ax=axes[0], orientation='vertical',
                 fraction=0.04, pad=0.02, label='$d_\\text{wall}$')

    # Right: RANS/LES map. RANS=1 where d_wall < SWITCH, LES=0 otherwise.
    rans_mask = (d_wall < SWITCH).astype(float)
    if solid is not None:
        rans_mask = np.ma.array(rans_mask, mask=solid)
    cmap = mcolors.ListedColormap(['#3b6dd6', '#e85c5c'])
    im_r = axes[1].imshow(rans_mask, origin='lower', cmap=cmap, vmin=0, vmax=1)
    n_fluid = (~solid).sum() if solid is not None else NX*NY
    n_les = int(((d_wall >= SWITCH) & (~solid if solid is not None else np.ones_like(d_wall, dtype=bool))).sum())
    les_frac = n_les / n_fluid
    axes[1].set_title(f"RANS / LES branch (LES fraction = {les_frac:.3f})")
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    axes[1].set_aspect('equal')
    cbar = plt.colorbar(im_r, ax=axes[1], orientation='vertical',
                        fraction=0.04, pad=0.02, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(['LES', 'RANS'])

    plt.suptitle(f"SA-DES length-scale switch ($C_{{DES}}\\Delta = {SWITCH}$ LU)",
                 fontsize=12)
    plt.tight_layout()
    out_fname = f'des_region_map_{name}.png'
    plt.savefig(OUTPUT_DIR / out_fname, dpi=180, bbox_inches='tight')
    plt.savefig(PUBLISHED_ASSET_DIR / out_fname, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_fname} (LES fraction = {les_frac:.3f})')


target = sys.argv[1] if len(sys.argv) > 1 else None
if target is None:
    for name, cfg in CASES.items():
        render_case(name, cfg)
elif target in CASES:
    render_case(target, CASES[target])
else:
    raise SystemExit(f'unknown case: {target} (choose from {list(CASES)})')
