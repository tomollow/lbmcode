"""Result figure for src/sec5/lbmzalesak.c (Zalesak's rotating slotted disk).

Reads the final order-parameter field produced by

    cmd /c scripts\\run_one.cmd src\\sec5\\lbmzalesak.c

from outputs/sec5/lbmzalesak/dataphi (the field after exactly one full
revolution, 2500 steps) and compares it against the analytically reconstructed
initial slotted disk. Builds a 4-panel figure:

  (a) initial order-parameter field phi_0 with the phi=0 interface
  (b) final phi after one revolution with the phi=0 interface
  (c) overlay of the initial vs final phi=0 contours (shape preservation)
  (d) phi along the vertical slot column (i=24): initial vs final

Also writes the benchmark table to
    docs/sec5/generated/lbmzalesak_advection.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmzalesak"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"
CSV_PATH = ROOT_DIR / "docs" / "sec5" / "generated" / "lbmzalesak_advection.csv"

# Source parameters (lbmzalesak.c).
NX = NY = 50
RADIUS = NX * 0.4           # 20.0
SLOT_I = (23, 26)           # 93*nx//200 .. 107*nx//200  (integer division)
SLOT_J = (2, NY // 2)       # 2 .. 25
SLOT_COL = 24               # vertical cut inside the slot

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def load_grid(path: Path) -> np.ndarray:
    """Load dataphi: one lattice row j per text line, columns i = 0..nx."""
    rows = []
    for line in path.read_text().splitlines():
        vals = [float(v) for v in line.split()]
        if vals:
            rows.append(vals)
    return np.array(rows)            # shape (ny+1, nx+1), [j, i]


def initial_field() -> np.ndarray:
    """Reconstruct the exact initial slotted disk (matches the C integer math).

    phi = -1 everywhere; phi = +1 inside the disk of radius 0.4*nx centred at
    (nx/2, ny/2); then phi = -1 again inside the slot i in [23,26], j in [2,25].
    Returned with the same [j, i] layout as dataphi.
    """
    phi = -np.ones((NY + 1, NX + 1))
    for j in range(NY + 1):
        for i in range(NX + 1):
            dist = np.hypot(i - NX * 0.5, j - NY * 0.5)
            if dist <= NX * 0.4:
                phi[j, i] = 1.0
    for i in range(SLOT_I[0], SLOT_I[1] + 1):
        for j in range(SLOT_J[0], SLOT_J[1] + 1):
            phi[j, i] = -1.0
    return phi


def main() -> None:
    phi_f = load_grid(RUN_DIR / "dataphi")       # (51, 51) [j, i]
    phi_0 = initial_field()

    # Indicator (color) functions for the geometric error norm.
    h0 = (phi_0 > 0.0).astype(float)
    hf = (phi_f > 0.0).astype(float)
    area0 = float(h0.sum())
    areaf = float(hf.sum())
    # Zalesak / Rudman geometric error: E1 = sum|Hf - H0| / sum H0.
    e1 = float(np.abs(hf - h0).sum() / h0.sum())

    mass0 = float(phi_0.sum())
    massf = float(phi_f.sum())

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.4), constrained_layout=True)

    # (a) initial field.
    ax = axes[0, 0]
    im = ax.imshow(phi_0, origin="lower", cmap="coolwarm", vmin=-1, vmax=1)
    ax.contour(phi_0, levels=[0.0], colors="k", linewidths=1.2)
    ax.set_title(r"(a) 初期 $\phi_0$（切欠き円板）")
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$\phi$")

    # (b) final field after one revolution.
    ax = axes[0, 1]
    im = ax.imshow(phi_f, origin="lower", cmap="coolwarm", vmin=-1, vmax=1)
    ax.contour(phi_f, levels=[0.0], colors="k", linewidths=1.2)
    ax.set_title(r"(b) 1 周回転後の $\phi$（2500 歩）")
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$\phi$")

    # (c) contour overlay.
    ax = axes[1, 0]
    ax.add_patch(plt.Rectangle((0, 0), NX, NY, facecolor="#eef2f6",
                               edgecolor="none"))
    c0 = ax.contour(phi_0, levels=[0.0], colors="#1f4e79", linewidths=1.6,
                    linestyles="dashed")
    cf = ax.contour(phi_f, levels=[0.0], colors="#c0392b", linewidths=1.6)
    ax.set_aspect("equal")
    ax.set_xlim(0, NX)
    ax.set_ylim(0, NY)
    ax.set_title(rf"(c) 界面 $\phi=0$ の比較（幾何誤差 $E_1={e1*100:.1f}\%$）")
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    ax.plot([], [], color="#1f4e79", ls="dashed", label="初期")
    ax.plot([], [], color="#c0392b", label="1 周後")
    ax.legend(fontsize=9, loc="upper right")

    # (d) phi along the slot column i=24.
    ax = axes[1, 1]
    jj = np.arange(NY + 1)
    ax.plot(jj, phi_0[:, SLOT_COL], "o--", color="#1f4e79", ms=3, label="初期")
    ax.plot(jj, phi_f[:, SLOT_COL], "-", color="#c0392b", label="1 周後")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.axvspan(SLOT_J[0], SLOT_J[1], color="#f6d6cf", alpha=0.5, label="切欠き $j$ 範囲")
    ax.set_title(rf"(d) 切欠き列 $i={SLOT_COL}$ 上の $\phi(j)$")
    ax.set_xlabel("$j$")
    ax.set_ylabel(r"$\phi$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    fig.suptitle(
        "図 5.1 Zalesak 円板の 1 回転移流テスト（既定: $\\mathrm{Pe}=400$, $W=2$, $\\tau=0.75$）",
        fontsize=13)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "lbmzalesak_results.png",
                 PUBLISHED_ASSET_DIR / "lbmzalesak_results.png"):
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {path}")

    # Benchmark CSV.
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["quantity", "initial", "after_one_rev", "rel_change_pct"])
        writer.writerow(["filled_area_cells", f"{area0:.1f}", f"{areaf:.1f}",
                         f"{(areaf - area0) / area0 * 100:.3f}"])
        writer.writerow(["total_mass_sum_phi", f"{mass0:.6e}", f"{massf:.6e}",
                         f"{(massf - mass0) / abs(mass0) * 100:.3f}"])
        writer.writerow(["geometric_error_E1", "0.0", f"{e1:.6f}", ""])
        writer.writerow(["phi_max", f"{phi_0.max():.6f}", f"{phi_f.max():.6f}", ""])
        writer.writerow(["phi_min", f"{phi_0.min():.6f}", f"{phi_f.min():.6f}", ""])
    print(f"Saved CSV to {CSV_PATH}")
    print(f"area:  init={area0:.0f}  final={areaf:.0f}  "
          f"change={(areaf - area0) / area0 * 100:+.2f}%")
    print(f"mass:  init={mass0:.4e}  final={massf:.4e}  "
          f"change={(massf - mass0) / abs(mass0) * 100:+.3f}%")
    print(f"geometric error E1 = {e1:.4f}  ({e1*100:.2f}%)")
    print(f"phi range: init [{phi_0.min():.3f}, {phi_0.max():.3f}]  "
          f"final [{phi_f.min():.3f}, {phi_f.max():.3f}]")


if __name__ == "__main__":
    main()
