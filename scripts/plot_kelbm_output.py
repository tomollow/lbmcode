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

csv_path = ROOT_DIR / 'outputs' / 'sec4' / 'kelbm' / 'kelbm_output.csv'
df = pd.read_csv(csv_path)
NX = df['x'].max() + 1
NY = df['y'].max() + 1

# Vectorized: pivot the long-format CSV into 2D grids (much faster than iterrows)
def to_grid(var):
    return df.pivot(index='y', columns='x', values=var).reindex(
        index=range(NY), columns=range(NX)
    ).values

u = to_grid('u')
v = to_grid('v')
k = to_grid('k')
eps = to_grid('epsilon')

fig, axs = plt.subplots(2, 2, figsize=(10, 9))
im0 = axs[0, 0].imshow(u, origin='lower', cmap='jet')
axs[0, 0].set_title('u (x方向速度)')
plt.colorbar(im0, ax=axs[0, 0])

im1 = axs[0, 1].imshow(v, origin='lower', cmap='jet')
axs[0, 1].set_title('v (y方向速度)')
plt.colorbar(im1, ax=axs[0, 1])

im2 = axs[1, 0].imshow(k, origin='lower', cmap='plasma')
axs[1, 0].set_title('k (乱流エネルギー)')
plt.colorbar(im2, ax=axs[1, 0])

im3 = axs[1, 1].imshow(eps, origin='lower', cmap='plasma')
axs[1, 1].set_title(r'$\varepsilon$ (散逸率)')
plt.colorbar(im3, ax=axs[1, 1])

plt.tight_layout()
fname = 'kelbm_output_overview.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved overview plot to {OUTPUT_DIR / fname}')
print(f'Saved overview plot to {PUBLISHED_ASSET_DIR / fname}')
