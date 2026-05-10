
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

csv_path = os.path.join('outputs', 'sec4', 'kelbm', 'kelbm_output.csv')
df = pd.read_csv(csv_path)
NX = df['x'].max() + 1
NY = df['y'].max() + 1

# y方向中心断面（x=NX//2）を抽出
df_c = df[df['x'] == NX//2].sort_values('y')
y = df_c['y'].values
u = df_c['u'].values
v = df_c['v'].values
k = df_c['k'].values
eps = df_c['epsilon'].values

fig, axs = plt.subplots(3, 1, figsize=(7, 10))
axs[0].plot(u, y, label='u')
axs[0].set_xlabel('u')
axs[0].set_ylabel('y')
axs[0].set_title('中心断面 x=NX/2 の u 分布')
axs[0].grid()

axs[1].plot(k, y, label='k', color='orange')
axs[1].set_xlabel('k')
axs[1].set_ylabel('y')
axs[1].set_title('中心断面 x=NX/2 の k 分布')
axs[1].grid()

axs[2].plot(eps, y, label='epsilon', color='green')
axs[2].set_xlabel('epsilon')
axs[2].set_ylabel('y')
axs[2].set_title('中心断面 x=NX/2 の epsilon 分布')
axs[2].grid()

plt.tight_layout()
plt.show()
