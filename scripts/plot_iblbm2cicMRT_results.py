"""Result figure for src/sec6/iblbm2cicMRT.c (cylindrical Couette, IB-LBM).

Reads the converged interior velocity fields produced by

    cmd /c scripts\\run_one.cmd src\\sec6\\iblbm2cicMRT.c

from outputs/sec6/iblbm2cicMRT/ and builds a 4-panel figure:

  (a) |u| magnitude with the two immersed cylinders + velocity vectors
  (b) clockwise tangential velocity u_theta vs radius, with the analytic
      circular-Couette profile overlaid
  (c) azimuthally-averaged u_theta(r) vs analytic profile (benchmark)
  (d) absolute error |u_theta - u_analytic| vs radius inside the annulus

Also writes the radial benchmark table to
    docs/sec6/generated/iblbm2cicMRT_couette_profile.csv
and prints the relative L2 error that the C code reports as `err`.

This is the MRT + implicit-velocity-correction sibling of
plot_iblbm2cdfSRT_results.py; the field layout and analytic profile are
identical, only the run directory / output names differ.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec6" / "iblbm2cicMRT"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec6"
CSV_PATH = ROOT_DIR / "docs" / "sec6" / "generated" / "iblbm2cicMRT_couette_profile.csv"

# Source parameters (iblbm2cicMRT.c).
NX = NY = 50
RP_OUT = 70.0 / 200.0 * NX   # rp[0] = 17.5 (outer, stationary)
RP_IN = 45.0 / 200.0 * NX    # rp[1] = 11.25 (inner, rotating)
U0 = 0.01
CENTER = NX // 2             # 25 (integer center used by the source)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def load_grid(path: Path) -> np.ndarray:
    """Load the interior field; one lattice row (fixed j) per text line.

    The source writes  for j=1..ny-1 { for i=1..nx-1 { u[i][j] } }  so the
    array is data[row=j-1, col=i-1] over the interior 1..49.
    """
    rows = []
    for line in path.read_text().splitlines():
        vals = [float(v) for v in line.split()]
        if vals:
            rows.append(vals)
    return np.array(rows)


def analytic_utheta(r: np.ndarray) -> np.ndarray:
    """Stokes circular-Couette tangential velocity (clockwise positive)."""
    return U0 * (r / RP_OUT - RP_OUT / r) / (RP_IN / RP_OUT - RP_OUT / RP_IN)


def main() -> None:
    u = load_grid(RUN_DIR / "datau")   # (49, 49) rows=j-1, cols=i-1
    v = load_grid(RUN_DIR / "datav")
    ny_i, nx_i = u.shape

    # Physical coordinates of each interior sample (i = col+1, j = row+1).
    ci = np.arange(nx_i) + 1
    rj = np.arange(ny_i) + 1
    ii, jj = np.meshgrid(ci, rj)          # ii[row,col]=i, jj[row,col]=j
    dx = ii - CENTER
    dy = jj - CENTER
    r = np.sqrt(dx ** 2 + dy ** 2)
    r_safe = np.where(r == 0, 1.0, r)

    speed = np.sqrt(u ** 2 + v ** 2)
    # Clockwise tangential velocity, matching the source's ut definition.
    utheta = (u * dy - v * dx) / r_safe

    annulus = (r >= RP_IN) & (r <= RP_OUT)

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 9.6), constrained_layout=True)

    # (a) speed map + cylinders + velocity vectors.
    ax = axes[0, 0]
    im = ax.imshow(speed, origin="lower", cmap="viridis",
                   extent=[1, nx_i, 1, ny_i])
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(CENTER + RP_OUT * np.cos(th), CENTER + RP_OUT * np.sin(th),
            color="white", lw=1.6, label="外円筒 (静止)")
    ax.plot(CENTER + RP_IN * np.cos(th), CENTER + RP_IN * np.sin(th),
            color="red", lw=1.6, label="内円筒 (回転)")
    step = 3
    ax.quiver(ii[::step, ::step], jj[::step, ::step],
              u[::step, ::step], v[::step, ::step],
              color="white", scale_units="xy", angles="xy", width=0.004)
    ax.set_title(rf"(a) 速度の大きさ $|\mathbf{{u}}|$ (最大 {speed.max():.3e})")
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="upper right")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$|\mathbf{u}|$")

    # (b) scattered tangential velocity vs radius.
    ax = axes[0, 1]
    ax.scatter(r.ravel(), utheta.ravel(), s=6, color="#1f4e79", alpha=0.35,
               label="LBM 全格子点")
    rr = np.linspace(RP_IN, RP_OUT, 200)
    ax.plot(rr, analytic_utheta(rr), "-", color="#c0392b", lw=2.0,
            label="解析解 (annulus)")
    ax.axvline(RP_IN, color="0.5", ls=":", lw=1.0)
    ax.axvline(RP_OUT, color="0.5", ls=":", lw=1.0)
    ax.axhline(U0, color="0.7", ls="--", lw=0.8)
    ax.set_title(r"(b) 時計回り接線速度 $u_\theta(r)$ の散布")
    ax.set_xlabel("$r$ (格子単位)")
    ax.set_ylabel(r"$u_\theta$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    # (c) azimuthally-averaged profile vs analytic (the benchmark).
    nbins = 24
    edges = np.linspace(r.min(), r.max(), nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    rflat, uflat = r.ravel(), utheta.ravel()
    mean_prof = np.full(nbins, np.nan)
    std_prof = np.zeros(nbins)
    for b in range(nbins):
        sel = (rflat >= edges[b]) & (rflat < edges[b + 1])
        if sel.any():
            mean_prof[b] = uflat[sel].mean()
            std_prof[b] = uflat[sel].std()

    ax = axes[1, 0]
    ax.errorbar(centers, mean_prof, yerr=std_prof, fmt="o", color="#1f4e79",
                ms=4, capsize=2, label="LBM 方位平均")
    ax.plot(rr, analytic_utheta(rr), "-", color="#c0392b", lw=2.0,
            label="解析解")
    ax.axvspan(RP_IN, RP_OUT, color="#f6d6cf", alpha=0.4, label="比較領域 (annulus)")
    ax.set_title(r"(c) 方位平均 $\langle u_\theta\rangle(r)$ とベンチマーク")
    ax.set_xlabel("$r$ (格子単位)")
    ax.set_ylabel(r"$u_\theta$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    # (d) absolute error vs radius inside the annulus.
    ax = axes[1, 1]
    ra = r[annulus]
    err_pts = np.abs(utheta[annulus] - analytic_utheta(ra))
    ax.scatter(ra, err_pts, s=8, color="#7a1f12", alpha=0.5)
    ax.set_title(r"(d) annulus 内の絶対誤差 $|u_\theta-u_{\theta,\mathrm{exact}}|$")
    ax.set_xlabel("$r$ (格子単位)")
    ax.set_ylabel(r"$|\Delta u_\theta|$")
    ax.grid(alpha=0.3)

    # Relative L2 error over the annulus (the C code's `err`).
    num = np.sum((utheta[annulus] - analytic_utheta(r[annulus])) ** 2)
    den = np.sum(analytic_utheta(r[annulus]) ** 2)
    rel_l2 = float(np.sqrt(num / den))

    fig.suptitle(
        f"図 6.3 円筒 Couette 流の IB-LBM 解 (陰的補正・MRT, 相対 L2 誤差 err = {rel_l2:.4f})",
        fontsize=13)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "iblbm2cicMRT_results.png",
                 PUBLISHED_ASSET_DIR / "iblbm2cicMRT_results.png"):
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {path}")

    # Benchmark CSV: azimuthally-averaged profile vs analytic in the annulus.
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["r", "utheta_lbm_mean", "utheta_std",
                         "utheta_analytic", "abs_error", "rel_error"])
        for b in range(nbins):
            if np.isnan(mean_prof[b]):
                continue
            if not (RP_IN <= centers[b] <= RP_OUT):
                continue
            ana = float(analytic_utheta(np.array([centers[b]]))[0])
            abs_e = abs(mean_prof[b] - ana)
            rel_e = abs_e / abs(ana) if ana != 0 else float("nan")
            writer.writerow([f"{centers[b]:.4f}", f"{mean_prof[b]:.6e}",
                             f"{std_prof[b]:.6e}", f"{ana:.6e}",
                             f"{abs_e:.6e}", f"{rel_e:.6e}"])
    print(f"Saved CSV to {CSV_PATH}")
    print(f"Relative L2 error in annulus (matches C `err`) = {rel_l2:.6f}")
    print(f"max |u| = {speed.max():.6e}  (u0 = {U0})")


if __name__ == "__main__":
    main()
