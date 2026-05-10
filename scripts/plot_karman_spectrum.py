"""Time series and frequency spectrum of the Karman vortex shedding.

Loads probe time series for both pure LBM and k-eps runs, plots:
- v_probe(t) showing the onset of shedding
- power spectrum highlighting the Strouhal frequency
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.rcParams['font.family'] = 'Yu Gothic'

ROOT_DIR = Path(__file__).resolve().parents[1]
PURE_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'karman' / 'karman_probe.csv'
KEPS_CSV = ROOT_DIR / 'outputs' / 'sec4' / 'karman_keps' / 'karman_keps_probe.csv'
OUTPUT_DIR = ROOT_DIR / 'outputs' / 'sec4'
PUBLISHED_ASSET_DIR = ROOT_DIR / 'docs' / 'assets' / 'sec4'

pure = np.loadtxt(PURE_CSV, delimiter=',', skiprows=1)
keps = np.loadtxt(KEPS_CSV, delimiter=',', skiprows=1)

D = 20
SAMPLE_INTERVAL = 5

LEGEND_KW = dict(loc='upper center', bbox_to_anchor=(0.5, -0.13),
                 ncol=2, frameon=True, fontsize=11,
                 handlelength=3.0, handletextpad=0.8, columnspacing=1.5)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Left: v_probe time series — see when shedding starts
ax = axes[0]
ax.plot(pure[:, 0], pure[:, 3], lw=1.0, color='C0', label='Pure LBM')
ax.plot(keps[:, 0], keps[:, 3], lw=1.0, color='C3', label='LBM + $k$-$\\varepsilon$')
ax.set_xlabel('step')
ax.set_ylabel('$v$ at probe')
ax.set_title('プローブ点 $v(t)$ — 渦放出の発達')
ax.grid(True, ls=':')
ax.legend(**LEGEND_KW)

# Right: power spectrum (late half)
ax = axes[1]


def spectrum(probe):
    N = len(probe)
    late = probe[3 * N // 4:, 3]
    late = late - late.mean()
    N2 = len(late)
    fft = np.fft.rfft(late)
    freqs = np.fft.rfftfreq(N2, SAMPLE_INTERVAL)
    return freqs, np.abs(fft)


fp, mp = spectrum(pure)
fk, mk = spectrum(keps)
# Peak detection
def peak_freq(freqs, mag):
    if mag[1:].size == 0:
        return 0.0
    idx = np.argmax(mag[1:]) + 1
    return freqs[idx]
fp_peak = peak_freq(fp, mp)
fk_peak = peak_freq(fk, mk)
u_pure = pure[-1, 1]
u_keps = keps[-1, 1]
St_pure = fp_peak * D / u_pure if u_pure > 0 else 0
St_keps = fk_peak * D / u_keps if u_keps > 0 else 0

ax.semilogy(fp[1:], mp[1:], lw=1.5, color='C0',
            label=f'Pure LBM (St = {St_pure:.3f})')
ax.semilogy(fk[1:], mk[1:], lw=1.5, color='C3',
            label=f'LBM + $k$-$\\varepsilon$ (St = {St_keps:.3f})')
ax.axvline(fp_peak, ls=':', color='C0', alpha=0.5)
ax.axvline(fk_peak, ls=':', color='C3', alpha=0.5)
ax.set_xlabel('frequency (per step)')
ax.set_ylabel('|FFT($v_{\\rm probe}$)|')
ax.set_title('パワースペクトル — Strouhal 周波数の検出')
ax.set_xlim(0, 0.005)
ax.grid(True, ls=':', which='both')
ax.legend(**LEGEND_KW)

plt.suptitle(f'Kármán 渦列 — 時系列とスペクトル（NX=360, NY=80, $D=20$, $Re_D \\approx {u_pure*D/((0.55-0.5)/3.0):.0f}$）')
plt.tight_layout()

fname = 'karman_spectrum.png'
plt.savefig(OUTPUT_DIR / fname, dpi=200)
plt.savefig(PUBLISHED_ASSET_DIR / fname, dpi=200)
print(f'Saved spectrum plot to {OUTPUT_DIR / fname}')
print(f'Saved spectrum plot to {PUBLISHED_ASSET_DIR / fname}')
print(f'Pure  Strouhal = {St_pure:.4f}, freq = {fp_peak:.5f}/step')
print(f'k-eps Strouhal = {St_keps:.4f}, freq = {fk_peak:.5f}/step')
print(f'(literature for Re~100 cylinder in unbounded flow: ~0.165)')
