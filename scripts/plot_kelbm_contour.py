"""2D contour plots for kelbm output.

Variant `keps` (default) reads outputs/sec4/kelbm/kelbm_output.csv and renders
u, k, epsilon. Variant `les` reads outputs/sec4/kelbm_les/kelbm_les_output.csv
and renders u, nut. Output PNGs are suffixed with the variant.

Usage: python plot_kelbm_contour.py [keps|les]   (default: keps)
"""

import sys
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

variant = sys.argv[1] if len(sys.argv) > 1 else 'keps'
if variant not in ('keps', 'les'):
    raise SystemExit(f'unknown variant: {variant}')

if variant == 'keps':
    csv_path = OUTPUT_DIR / 'kelbm' / 'kelbm_output.csv'
    fields = [('u', 'jet'), ('k', 'plasma'), ('epsilon', 'plasma')]
    title_suffix = '(k-ε)'
else:
    csv_path = OUTPUT_DIR / 'kelbm_les' / 'kelbm_les_output.csv'
    fields = [('u', 'jet'), ('nut', 'viridis')]
    title_suffix = '(LES)'

df = pd.read_csv(csv_path)
NX = df['x'].max() + 1
NY = df['y'].max() + 1


def to_grid(var):
    return df.pivot(index='y', columns='x', values=var).reindex(
        index=range(NY), columns=range(NX)
    ).values


for name, cmap in fields:
    arr = to_grid(name)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(f'{name} の 2 次元分布コンター {title_suffix}')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    im = ax.imshow(arr, origin='lower', cmap=cmap)
    fig.colorbar(im, ax=ax, orientation='horizontal', label=name, pad=0.15, fraction=0.05)
    fig.tight_layout()
    fname = f'kelbm_contour_{name}_{variant}.png'
    plt.savefig(OUTPUT_DIR / fname, dpi=200)
    plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
    plt.close()

print(f'Saved {variant} contour plots to {OUTPUT_DIR} and {PUBLISHED_ASSET_DIR}')
