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

df_c = df[df['x'] == NX // 2].sort_values('y')
y = df_c['y'].values
u = df_c['u'].values
k = df_c['k'].values
eps = df_c['epsilon'].values

fig, axs = plt.subplots(3, 1, figsize=(7, 10))
axs[0].plot(u, y)
axs[0].set_xlabel('u')
axs[0].set_ylabel('y')
axs[0].set_title('中心断面 x=NX/2 の u 分布')
axs[0].grid()

axs[1].plot(k, y, color='orange')
axs[1].set_xlabel('k')
axs[1].set_ylabel('y')
axs[1].set_title('中心断面 x=NX/2 の k 分布')
axs[1].grid()

axs[2].plot(eps, y, color='green')
axs[2].set_xlabel(r'$\varepsilon$')
axs[2].set_ylabel('y')
axs[2].set_title(r'中心断面 x=NX/2 の $\varepsilon$ 分布')
axs[2].grid()

plt.tight_layout()
fname = 'kelbm_centerline.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved centerline plot to {OUTPUT_DIR / fname}')
print(f'Saved centerline plot to {PUBLISHED_ASSET_DIR / fname}')
