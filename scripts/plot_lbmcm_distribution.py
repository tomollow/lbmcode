"""Plot u, v, stream function distributions for src/sec4/lbmcm.c.

Reads the three output files (dataCMu, dataCMv, dataCMs) produced by the
CM @ Re=5000 run and produces:

- A 2x2 figure with: u/U_lid heatmap, v/U_lid heatmap, stream-function contours,
  and centerline profiles (u(L/2, y) and v(x, L/2)).
- Saved to docs/assets/sec4/lbmcm_distribution.png and the run directory.

The run directory is outputs/sec4/lbmcm/cm_re5000/ when populated by
scripts/run_lbmcm_compare.ps1; otherwise outputs/sec4/lbmcm/ is used as
a fallback (single-run layout).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_BASE = ROOT_DIR / "outputs" / "sec4" / "lbmcm"
RUN_DIR = (
    RUN_BASE / "cm_re5000"
    if (RUN_BASE / "cm_re5000" / "dataCMu").exists()
    else RUN_BASE
)
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec4"


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic",
    "Meiryo",
    "MS Gothic",
    "Noto Sans CJK JP",
    "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                rows.append([float(value) for value in stripped.split()])
    return np.array(rows, dtype=float)


def main() -> None:
    data_u = read_matrix(RUN_DIR / "dataCMu")
    data_v = read_matrix(RUN_DIR / "dataCMv")
    data_psi = read_matrix(RUN_DIR / "dataCMs")

    ny, nx = data_psi.shape
    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(x, y)

    mid_i = nx // 2
    mid_j = ny // 2

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.4), constrained_layout=True)

    u_ax, v_ax = axes[0]
    psi_ax, line_ax = axes[1]

    u_max = float(np.max(np.abs(data_u)))
    im_u = u_ax.imshow(
        data_u,
        origin="lower",
        extent=(0.0, 1.0, 0.0, 1.0),
        cmap="RdBu_r",
        vmin=-u_max,
        vmax=u_max,
        aspect="equal",
    )
    u_ax.set_title(r"水平速度 $u / U_{\rm lid}$")
    u_ax.set_xlabel(r"$x / L$")
    u_ax.set_ylabel(r"$y / L$")
    figure.colorbar(im_u, ax=u_ax, fraction=0.046, pad=0.04)

    v_max = float(np.max(np.abs(data_v)))
    im_v = v_ax.imshow(
        data_v,
        origin="lower",
        extent=(0.0, 1.0, 0.0, 1.0),
        cmap="RdBu_r",
        vmin=-v_max,
        vmax=v_max,
        aspect="equal",
    )
    v_ax.set_title(r"鉛直速度 $v / U_{\rm lid}$")
    v_ax.set_xlabel(r"$x / L$")
    v_ax.set_ylabel(r"$y / L$")
    figure.colorbar(im_v, ax=v_ax, fraction=0.046, pad=0.04)

    psi_min = float(data_psi.min())
    psi_max = float(data_psi.max())
    main_levels = np.linspace(psi_min, 0.0, 11)[:-1]
    sec_levels = np.linspace(0.0, psi_max, 8)
    levels = np.unique(np.concatenate([main_levels, sec_levels]))
    cs = psi_ax.contour(xx, yy, data_psi, levels=levels, colors="black", linewidths=0.9)
    psi_ax.clabel(cs, inline=True, fontsize=7, fmt="%.3f")
    psi_ax.set_aspect("equal")
    psi_ax.set_xlabel(r"$x / L$")
    psi_ax.set_ylabel(r"$y / L$")
    psi_ax.set_title(
        rf"流線関数 $\psi/(U_{{\rm lid}}\,L)$" "\n"
        rf"$\psi_{{\min}}={psi_min:.4f},\ \psi_{{\max}}={psi_max:.4f}$"
    )
    psi_ax.plot([0.0, 1.0], [1.0, 1.0], color="red", linewidth=1.8, solid_capstyle="butt")
    psi_ax.text(0.5, 1.02, r"$U_{\rm lid}\to$", ha="center", va="bottom", color="red")

    line_ax.plot(
        data_u[:, mid_i],
        y,
        color="C0",
        linewidth=1.6,
        marker="o",
        markersize=3.5,
        markerfacecolor="white",
        label=r"$u(L/2,\,y)/U_{\rm lid}$",
    )
    line_ax.plot(
        x,
        data_v[mid_j, :],
        color="C3",
        linewidth=1.6,
        marker="s",
        markersize=3.2,
        markerfacecolor="white",
        label=r"$v(x,\,L/2)/U_{\rm lid}$",
    )
    line_ax.axhline(0.0, color="0.7", linewidth=0.7)
    line_ax.axvline(0.0, color="0.7", linewidth=0.7)
    line_ax.grid(True, linestyle=":", alpha=0.6)
    line_ax.set_xlabel(r"$u/U_{\rm lid}$ または $x/L$")
    line_ax.set_ylabel(r"$y/L$ または $v/U_{\rm lid}$")
    line_ax.set_title("中心線プロファイル")
    line_ax.legend(loc="best", fontsize=9)

    figure.suptitle(
        "lbmcm.c: D2Q9 + 中心モーメント衝突, Re = 5000, lid-driven cavity (51×51)",
        fontsize=13,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmcm_distribution.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcm_distribution.png"
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    figure.savefig(published_path, dpi=200, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")

    print(f"u/U_lid range  : [{data_u.min():.4f}, {data_u.max():.4f}]")
    print(f"v/U_lid range  : [{data_v.min():.4f}, {data_v.max():.4f}]")
    print(f"psi range      : [{psi_min:.6f}, {psi_max:.6f}]")
    print(f"u(L/2, y) max  : {data_u[:, mid_i].max():.4f} at y/L={y[np.argmax(data_u[:, mid_i])]:.3f}")
    print(f"u(L/2, y) min  : {data_u[:, mid_i].min():.4f} at y/L={y[np.argmin(data_u[:, mid_i])]:.3f}")
    print(f"v(x, L/2) max  : {data_v[mid_j, :].max():.4f} at x/L={x[np.argmax(data_v[mid_j, :])]:.3f}")
    print(f"v(x, L/2) min  : {data_v[mid_j, :].min():.4f} at x/L={x[np.argmin(data_v[mid_j, :])]:.3f}")


if __name__ == "__main__":
    main()
