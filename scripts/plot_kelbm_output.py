
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

# ファイルパス
csv_path = os.path.join('outputs', 'sec4', 'kelbm', 'kelbm_output.csv')

# CSV 読み込み
df = pd.read_csv(csv_path)
NX = df['x'].max() + 1
NY = df['y'].max() + 1

# 格子状に変換
def to_grid(var):
    arr = np.zeros((NY, NX))
    for _, row in df.iterrows():
        arr[int(row['y']), int(row['x'])] = row[var]
    return arr

u = to_grid('u')
v = to_grid('v')
k = to_grid('k')
eps = to_grid('epsilon')

fig, axs = plt.subplots(2, 2, figsize=(10, 9))

im0 = axs[0,0].imshow(u, origin='lower', cmap='jet')
axs[0,0].set_title('u (x方向速度)')
plt.colorbar(im0, ax=axs[0,0])

im1 = axs[0,1].imshow(v, origin='lower', cmap='jet')
axs[0,1].set_title('v (y方向速度)')
plt.colorbar(im1, ax=axs[0,1])

im2 = axs[1,0].imshow(k, origin='lower', cmap='plasma')
axs[1,0].set_title('k (乱流エネルギー)')
plt.colorbar(im2, ax=axs[1,0])

im3 = axs[1,1].imshow(eps, origin='lower', cmap='plasma')
axs[1,1].set_title('epsilon (散逸率)')
plt.colorbar(im3, ax=axs[1,1])

plt.tight_layout()
plt.show()
