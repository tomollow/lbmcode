"""Result figure for src/sec5/lbmlap.c (Laplace's law, static droplet).

Reads the converged fields produced by

    cmd /c scripts\\run_one.cmd src\\sec5\\lbmlap.c

from outputs/sec5/lbmlap/ and builds a 4-panel figure:

  (a) order-parameter field phi with the phi=0 interface contour
  (b) centerline phi(x) against the analytic tanh profile
  (c) spurious-current vector field over |u| magnitude
  (d) density rho along the horizontal centerline (the Laplace pressure jump)

Also writes the interface-profile benchmark table to
    docs/sec5/generated/lbmlap_interface_profile.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmlap"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"
CSV_PATH = ROOT_DIR / "docs" / "sec5" / "generated" / "lbmlap_interface_profile.csv"

# Source parameters (lbmlap.c).
NX = NY = 50
PHI0 = 1.0
WID = 5.0
RADIUS = NX * 0.25          # 12.5
CENTER = NX * 0.5           # 25.0

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def load_grid(path: Path) -> np.ndarray:
    """Load a whitespace-separated field, one lattice row (j) per text line."""
    rows = []
    for line in path.read_text().splitlines():
        vals = [float(v) for v in line.split()]
        if vals:
            rows.append(vals)
    return np.array(rows)


def main() -> None:
    phi = load_grid(RUN_DIR / "dataphi")        # (51, 51) rows=j, cols=i
    rho = load_grid(RUN_DIR / "datarho")         # (51, 51)
    u_in = load_grid(RUN_DIR / "datau")          # (49, 49) interior i,j = 1..49
    v_full = load_grid(RUN_DIR / "datav")        # (51, 51)
    phi_line = np.array([float(v) for v in
                         (RUN_DIR / "dataphi2D").read_text().split()])  # (51,)

    # Align interior velocity blocks: datav interior is rows/cols 1..49.
    v_in = v_full[1:NY, 1:NX]
    speed = np.sqrt(u_in ** 2 + v_in ** 2)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.2), constrained_layout=True)

    # (a) order-parameter field + interface contour.
    ax = axes[0, 0]
    im = ax.imshow(phi, origin="lower", cmap="coolwarm", vmin=-PHI0, vmax=PHI0)
    ax.contour(phi, levels=[0.0], colors="k", linewidths=1.2)
    ax.set_title(r"(a) 秩序変数場 $\phi$ と界面 ($\phi=0$)")
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$\phi$")

    # (b) centerline phi vs analytic tanh.
    ax = axes[0, 1]
    xx = np.arange(phi_line.size)
    r = np.abs(xx - CENTER)
    phi_analytic = -PHI0 * np.tanh(2.0 * (r - RADIUS) / WID)
    ax.plot(xx, phi_line, "o", color="#1f4e79", ms=4, label="LBM")
    ax.plot(xx, phi_analytic, "-", color="#c0392b",
            label=r"解析解 $\phi_0\tanh[2(R-|x-x_0|)/W]$")
    ax.set_title(r"(b) 中心線上の界面プロファイル")
    ax.set_xlabel("$i$ (中心線 $j=n_y/2$)")
    ax.set_ylabel(r"$\phi$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    # (c) spurious-current field.
    ax = axes[1, 0]
    im = ax.imshow(speed, origin="lower", cmap="viridis")
    step = 3
    ys, xs = np.mgrid[0:speed.shape[0]:step, 0:speed.shape[1]:step]
    ax.quiver(xs, ys, u_in[::step, ::step], v_in[::step, ::step],
              color="white", scale_units="xy", angles="xy")
    ax.set_title(rf"(c) 寄生流 $|\mathbf{{u}}|$ (最大 {speed.max():.2e})")
    ax.set_xlabel("$i-1$")
    ax.set_ylabel("$j-1$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$|\mathbf{u}|$")

    # (d) density along horizontal centerline (pressure jump).
    ax = axes[1, 1]
    rho_line = rho[NY // 2, :]
    dp = (rho[NY // 2, NX // 2] - rho[0, 0]) / 3.0
    ax.plot(np.arange(rho_line.size), rho_line, "-o", color="#1f4e79", ms=3)
    ax.axvspan(CENTER - RADIUS, CENTER + RADIUS, color="#f6d6cf", alpha=0.5,
               label="液滴内部")
    ax.set_title(rf"(d) 中心線密度 $\rho$ ($\Delta p={dp:.3e}$)")
    ax.set_xlabel("$i$ (中心線 $j=n_y/2$)")
    ax.set_ylabel(r"$\rho$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    fig.suptitle("図 5.1 静止液滴の自由エネルギー型 LBM 解 (既定: $\\sigma=10^{-4}$, $R=12.5$, $W=5$)",
                 fontsize=13)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "lbmlap_results.png",
                 PUBLISHED_ASSET_DIR / "lbmlap_results.png"):
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {path}")

    # Interface-profile benchmark CSV (a few representative points around the edge).
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["i", "r_from_center", "phi_lbm", "phi_analytic", "abs_error"])
        for i in range(phi_line.size):
            writer.writerow([i, f"{r[i]:.3f}", f"{phi_line[i]:.6e}",
                             f"{phi_analytic[i]:.6e}",
                             f"{abs(phi_line[i]-phi_analytic[i]):.6e}"])
    rms = float(np.sqrt(np.mean((phi_line - phi_analytic) ** 2)))
    print(f"Saved CSV to {CSV_PATH}")
    print(f"Interface profile RMS error vs tanh = {rms:.6e}")
    print(f"max |u| (spurious current) = {speed.max():.6e}")


if __name__ == "__main__":
    main()
