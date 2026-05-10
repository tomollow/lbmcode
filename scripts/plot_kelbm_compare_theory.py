
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

# kelbm.c のパラメータ（src/sec4/kelbm.c と同期）
TAU = 0.55
FORCE_X = 5e-6
nu0 = (TAU - 0.5) / 3.0
Cmu = 0.09

csv_path = ROOT_DIR / 'outputs' / 'sec4' / 'kelbm' / 'kelbm_output.csv'
df = pd.read_csv(csv_path)
NX = df['x'].max() + 1
NY = df['y'].max() + 1

df_c = df[df['x'] == NX // 2].sort_values('y')
y = df_c['y'].values.astype(float)
u = df_c['u'].values
k = df_c['k'].values
eps = df_c['epsilon'].values
nut = Cmu * k**2 / (eps + 1e-12)

yc = (NY - 1) / 2.0
H_half = NY / 2.0

# 層流 Poiseuille 予測（k-ε が無効ならこれ）
u_max_lam = FORCE_X * NY**2 / (8.0 * nu0)
u_lam = u_max_lam * (1.0 - ((y - yc) / H_half)**2)

# 摩擦速度（力学平衡から）
u_tau = (FORCE_X * H_half)**0.5
Re_tau = u_tau * H_half / nu0

fig, axes = plt.subplots(1, 2, figsize=(13, 6))

# (左) 速度プロファイル：LBM vs 層流予測
ax = axes[0]
ax.plot(u, y, lw=2, label='LBM ($k$-$\\varepsilon$ 有効)')
ax.plot(u_lam, y, '--', lw=1.5, label='層流予測（パラボラ型）')
ax.set_xlabel('u')
ax.set_ylabel('y')
ax.set_title(f'速度プロファイル（$Re_\\tau$≈{Re_tau:.0f}, $Re_{{max}}$≈{u.max()*NY/nu0:.0f}）')
ax.grid(True)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=True)

# (右) 渦粘性プロファイル：ν_t / ν_0
ax = axes[1]
ax.plot(nut / nu0, y, lw=2, color='C2')
ax.axvline(1.0, ls=':', color='k', lw=0.8, label='$\\nu_t = \\nu_0$')
ax.set_xlabel('$\\nu_t / \\nu_0$')
ax.set_ylabel('y')
ax.set_title(f'渦粘性プロファイル（壁: $\\nu_t/\\nu_0$={nut[0]/nu0:.2f}, 中心: {nut[len(nut)//2]/nu0:.2f}）')
ax.grid(True)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=1, frameon=True)

plt.suptitle(f'チャンネル流れ k-$\\varepsilon$ 検証 (NX={NX}, NY={NY}, $\\tau$={TAU}, $\\nu_0$={nu0:.4f})')
plt.tight_layout()
output_path = OUTPUT_DIR / 'kelbm_compare_theory.png'
published_path = PUBLISHED_ASSET_DIR / 'kelbm_compare_theory.png'
plt.savefig(output_path, dpi=200)
plt.savefig(published_path, dpi=200)
print(f'Saved plot to {output_path}')
print(f'Saved plot to {published_path}')
print(f'  u_max (LBM)       = {u.max():.5f}')
print(f'  u_max (laminar)   = {u_max_lam:.5f}')
print(f'  ratio LBM/laminar = {u.max()/u_max_lam:.3f}')
print(f'  Re_tau            = {Re_tau:.1f}')
print(f'  nu_t/nu_0 (max)   = {nut.max()/nu0:.3f}')
plt.show()
