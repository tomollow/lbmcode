"""Reproduce Fig. 4.5 (Sec. 4.6, Geier 2006 cascaded LBM) qualitatively.

Top panel:    speed magnitude |U| in the wake behind a rectangular obstacle
              (similar to Fig. 4.5(a)).
Bottom panel: radially-averaged 1D energy spectrum E(k) versus k built from
              wake-window snapshots, with a Kolmogorov -5/3 reference line
              (similar to Fig. 4.5(b)).

The book's figure is from Re = 1.4e6; we run a Re ~ O(10^3) demonstration that
the central-moment collision is stable in the tau -> 0.5 regime and produces a
broadband wake whose spectrum is power-law-like over a finite range.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / 'outputs' / 'sec4' / 'lbmcm_karman'
ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'
ASSET_DIR.mkdir(parents=True, exist_ok=True)


def load_snapshot(path: Path):
    df = pd.read_csv(path)
    nx = int(df['x'].max() + 1)
    ny = int(df['y'].max() + 1)

    def pivot(col):
        return df.pivot(index='y', columns='x', values=col).reindex(
            index=range(ny), columns=range(nx)).values

    return {
        'nx': nx, 'ny': ny,
        'u': pivot('u'), 'v': pivot('v'),
        'speed': pivot('speed'),
        'solid': pivot('solid'),
    }


def load_wake_frame(path: Path):
    df = pd.read_csv(path)
    xs = sorted(df['x'].unique())
    ys = sorted(df['y'].unique())
    nx = len(xs)
    ny = len(ys)
    u = df.pivot(index='y', columns='x', values='u').values
    v = df.pivot(index='y', columns='x', values='v').values
    return u, v, nx, ny


# ---- (a) speed-magnitude panel: last full snapshot --------------------------

snap_files = sorted(SRC_DIR.glob('lbmcm_karman_snapshot_*.csv'),
                    key=lambda p: int(re.search(r'(\d+)', p.stem).group(1)))
if not snap_files:
    raise SystemExit(f'no snapshots in {SRC_DIR}')

snap = load_snapshot(snap_files[-1])
speed = np.ma.array(snap['speed'], mask=(snap['solid'] > 0))

# ---- (b) energy spectrum from wake-window frames ---------------------------

wake_files = sorted(SRC_DIR.glob('lbmcm_karman_wake_*.csv'),
                    key=lambda p: int(re.search(r'(\d+)', p.stem).group(1)))
if not wake_files:
    raise SystemExit(f'no wake frames in {SRC_DIR}')

# Determine grid size from the first frame
u0, v0, NX_W, NY_W = load_wake_frame(wake_files[0])

# Time-mean velocity for detrending (subtract mean flow so the spectrum reflects
# the turbulent fluctuations, not the body-force-driven uniform component)
u_mean = np.zeros((NY_W, NX_W))
v_mean = np.zeros((NY_W, NX_W))
for path in wake_files:
    u, v, _, _ = load_wake_frame(path)
    u_mean += u
    v_mean += v
u_mean /= len(wake_files)
v_mean /= len(wake_files)

# 2D Hanning window to reduce spectral leakage from non-periodic wake field
wx = np.hanning(NX_W)
wy = np.hanning(NY_W)
window = np.outer(wy, wx)
window_norm = (window ** 2).sum()

# Radial wavenumber grid
kx = np.fft.fftfreq(NX_W) * NX_W   # cycles over the window, in lattice units
ky = np.fft.fftfreq(NY_W) * NY_W
KX, KY = np.meshgrid(kx, ky)
K = np.sqrt(KX ** 2 + KY ** 2)

# Bin into integer-k shells; max wavenumber limited by Nyquist of shorter side
k_max = int(min(NX_W, NY_W) // 2)
k_bins = np.arange(0, k_max + 1)
nbins = len(k_bins) - 1

E_k = np.zeros(nbins)
n_frames = 0
for path in wake_files:
    u, v, _, _ = load_wake_frame(path)
    up = (u - u_mean) * window
    vp = (v - v_mean) * window
    U_hat = np.fft.fft2(up)
    V_hat = np.fft.fft2(vp)
    E2d = 0.5 * (np.abs(U_hat) ** 2 + np.abs(V_hat) ** 2) / window_norm
    # Radial average -> times 2*pi*k_shell to get 1D E(k) (azimuthal area)
    for b in range(nbins):
        mask = (K >= k_bins[b]) & (K < k_bins[b + 1])
        if mask.any():
            E_k[b] += E2d[mask].mean() * 2 * np.pi * (b + 0.5)
    n_frames += 1
E_k /= n_frames

k_centers = k_bins[:-1] + 0.5
# Drop k=0 (DC component) for log-log plot
mask_pos = (k_centers >= 1) & (E_k > 0)
k_plot = k_centers[mask_pos]
E_plot = E_k[mask_pos]

# Reference -5/3 line, anchored near the inertial-range center
if len(k_plot) > 0:
    k_ref = k_plot[len(k_plot) // 3]
    E_ref = E_plot[len(k_plot) // 3]
    ref_line = E_ref * (k_plot / k_ref) ** (-5.0 / 3.0)
else:
    ref_line = None

# ---- assemble figure -------------------------------------------------------

fig = plt.figure(figsize=(11, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.7], hspace=0.32)

# (a) speed magnitude — match Fig. 4.5(a) layout (wide, grayscale-friendly)
axa = fig.add_subplot(gs[0])
vmax = float(np.percentile(speed.compressed(), 99.5))
im = axa.imshow(speed, origin='lower', cmap='magma', vmin=0, vmax=vmax,
                aspect='equal')
axa.set_title('(a) 長方形障害物周りの速さ $|\\mathbf{U}|$（中心モーメント LBM）',
              fontsize=13)
axa.set_xlabel('$x$ [lu]')
axa.set_ylabel('$y$ [lu]')
cb = fig.colorbar(im, ax=axa, orientation='horizontal',
                  shrink=0.55, pad=0.18, aspect=30)
cb.set_label('$|\\mathbf{U}|$ [lu]')

# (b) energy spectrum
axb = fig.add_subplot(gs[1])
axb.loglog(k_plot, E_plot, 'o', ms=5, mfc='C0', mec='C0',
           label='LBM (中心モーメント)')
if ref_line is not None:
    axb.loglog(k_plot, ref_line, '-', color='k', lw=1.4,
               label="Kolmogorov $k^{-5/3}$ 参照線")
axb.set_xlabel('$k$ [cycles per wake window]')
axb.set_ylabel('$E(k)$')
axb.set_title('(b) 乱流エネルギー $E(k)$ と波数 $k$', fontsize=13)
axb.grid(True, which='both', ls=':', alpha=0.6)
axb.legend(loc='lower left', frameon=True)

_final_step = int(re.search(r'(\d+)', snap_files[-1].stem).group(1))
fig.suptitle(f'カスケード／中心モーメント LBM による乱流計算  '
             f'($n_\\mathrm{{wake\\,frames}} = {len(wake_files)}$, '
             f'最終ステップ {_final_step})',
             fontsize=12)

out_png = ASSET_DIR / 'lbmcm_karman_spectrum.png'
fig.savefig(out_png, dpi=150, bbox_inches='tight')
fig.savefig(SRC_DIR / 'lbmcm_karman_spectrum.png', dpi=150, bbox_inches='tight')
print(f'wrote {out_png}')
