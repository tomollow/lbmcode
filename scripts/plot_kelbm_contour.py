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

def to_grid(var):
    arr = np.zeros((NY, NX))
    for _, row in df.iterrows():
        arr[int(row['y']), int(row['x'])] = row[var]
    return arr

u = to_grid('u')
k = to_grid('k')
eps = to_grid('epsilon')

for arr, name, cmap in zip([u, k, eps], ['u', 'k', 'epsilon'], ['jet', 'plasma', 'plasma']):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(f'{name} の2次元分布コンター')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    im = ax.imshow(arr, origin='lower', cmap=cmap)
    fig.colorbar(im, ax=ax, orientation='horizontal', label=name, pad=0.15, fraction=0.05)
    fig.tight_layout()
    fname = f'kelbm_contour_{name}.png'
    plt.savefig(OUTPUT_DIR / fname, dpi=200)
    plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
    plt.close()

print(f'Saved contour plots to {OUTPUT_DIR} and {PUBLISHED_ASSET_DIR}')
