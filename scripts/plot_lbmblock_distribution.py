"""Plot 2D distribution and mesh layout for src/sec4/lbmblock.c.

Reads u3000_coarse / u3000_fine / v3000_coarse / v3000_fine from
outputs/sec4/lbmblock/ and produces:

- lbmblock_distribution.png: 2D heatmaps of u/u_t on coarse and fine grids,
  plus center-line profile vs analytical, at t=3000.
- lbmblock_mesh.png: schematic of the multi-block grid layout showing
  coarse grid, fine grid, overlap region, interface, walls and periodic
  boundaries.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec4" / "lbmblock"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec4"

NX = 20
NY = 32
M = 2
NXF = M * NX  # 40
NYF = NY  # 32
NXC = NX  # 20
NYC = NY // M + 1  # 17
UT = 0.01

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


def plot_distribution() -> None:
    uc = read_matrix(RUN_DIR / "u3000_coarse")  # shape (nyc+1, nxc+1)
    uf = read_matrix(RUN_DIR / "u3000_fine")    # shape (nyf+1, nxf+1)

    assert uc.shape == (NYC + 1, NXC + 1)
    assert uf.shape == (NYF + 1, NXF + 1)

    # Physical y of fine: y_f / 2  in [0, ny/2 = 16]
    # Physical y of coarse: y_c + ny/2 - 1 in [ny/2 - 1, ny] = [15, 32]
    xf = np.arange(NXF + 1) * 0.5  # 0..20 in physical units
    yf = np.arange(NYF + 1) * 0.5  # 0..16
    xc = np.arange(NXC + 1) * 1.0  # 0..20
    yc = np.arange(NYC + 1) * 1.0 + (NY / 2 - 1)  # 15..32

    fig = plt.figure(figsize=(13.5, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.05, 0.95])

    ax_2d = fig.add_subplot(gs[0, 0:2])
    vmin, vmax = 0.0, 1.05

    pcm_f = ax_2d.pcolormesh(
        xf,
        yf,
        uf / UT,
        cmap="viridis",
        shading="auto",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )
    pcm_c = ax_2d.pcolormesh(
        xc,
        yc,
        uc / UT,
        cmap="viridis",
        shading="auto",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    # mark interface (physical y = ny/2 = 16) and overlap band (15..16)
    ax_2d.axhline(NY / 2, color="white", linestyle="--", linewidth=0.9, alpha=0.85)
    ax_2d.axhline(NY / 2 - 1, color="white", linestyle=":", linewidth=0.7, alpha=0.7)
    ax_2d.text(0.5, NY / 2 + 0.4, "interface (粗 $y_c=1$ / 細 $y_f=n_{y,f}$)",
               color="white", fontsize=8.5)
    ax_2d.text(0.5, NY / 2 - 0.9, "overlap 下端 (粗 $y_c=0$ / 細 $y_f=n_{y,f}-2$)",
               color="white", fontsize=8.5)

    ax_2d.set_xlim(0, NX)
    ax_2d.set_ylim(0, NY)
    ax_2d.set_aspect("equal")
    ax_2d.set_xlabel("$x$（物理単位）")
    ax_2d.set_ylabel("$y$（物理単位）")
    ax_2d.set_title("$u / u_t$ 分布 (t=3000)  ─ 細格子 (下) + 粗格子 (上)")
    cb = fig.colorbar(pcm_f, ax=ax_2d, shrink=0.85, pad=0.02)
    cb.set_label("$u / u_t$")

    # right: centerline profile
    ax_p = fig.add_subplot(gs[0, 2])
    j_mid_f = NXF // 2
    j_mid_c = NXC // 2
    ax_p.plot(uf[:, j_mid_f] / UT, yf / NY, color="C0", linewidth=1.6,
              marker="o", markersize=3.0, markerfacecolor="white",
              label="細格子 ($x=n_{x,f}/2$)")
    ax_p.plot(uc[:, j_mid_c] / UT, yc / NY, color="C3", linewidth=1.6,
              marker="s", markersize=3.0, markerfacecolor="white",
              label="粗格子 ($x=n_{x,c}/2$)")
    y_exact = np.linspace(0, NY, 401)
    ax_p.plot(y_exact / NY, y_exact / NY, color="black", linewidth=1.2,
              linestyle="--", label="解析解 $u_t y/H$")
    ax_p.axhline(0.5, color="gray", linewidth=0.6, linestyle=":")
    ax_p.set_xlim(-0.05, 1.05)
    ax_p.set_ylim(0.0, 1.0)
    ax_p.set_xlabel("$u / u_t$")
    ax_p.set_ylabel("$y / H$")
    ax_p.set_title("中心鉛直線プロファイル")
    ax_p.grid(True, linestyle=":", alpha=0.5)
    ax_p.legend(loc="lower right", fontsize=8.5)

    fig.suptitle("lbmblock.c: Multi-block LBM 分布結果 (Couette, t=3000)", fontsize=13)

    out = RUN_DIR / "lbmblock_distribution.png"
    pub = PUBLISHED_ASSET_DIR / "lbmblock_distribution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(pub, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")
    print(f"Saved {pub}")


def plot_mesh() -> None:
    fig, ax = plt.subplots(figsize=(13.0, 8.5), constrained_layout=True)

    # Fine grid nodes: physical y in [0, 16], step 0.5
    xf = np.arange(NXF + 1) * 0.5
    yf = np.arange(NYF + 1) * 0.5
    Xf, Yf = np.meshgrid(xf, yf)

    # Coarse grid nodes: physical y in [15, 32], step 1.0
    xc = np.arange(NXC + 1) * 1.0
    yc = np.arange(NYC + 1) * 1.0 + (NY / 2 - 1)
    Xc, Yc = np.meshgrid(xc, yc)

    # Overlap band y in [15, 16]
    overlap = patches.Rectangle(
        (0, NY / 2 - 1), NX, 1.0,
        linewidth=0, facecolor="#fff2b3", alpha=0.7, zorder=0,
    )
    ax.add_patch(overlap)

    # Fine grid lines
    for x in xf:
        ax.plot([x, x], [yf.min(), yf.max()], color="#1f6feb",
                linewidth=0.4, alpha=0.5, zorder=1)
    for y in yf:
        ax.plot([xf.min(), xf.max()], [y, y], color="#1f6feb",
                linewidth=0.4, alpha=0.5, zorder=1)
    ax.scatter(Xf.flatten(), Yf.flatten(), s=6, color="#1f6feb",
               zorder=3, label=f"細格子ノード ({NXF + 1}×{NYF + 1}, $\\delta x = 0.5$)")

    # Coarse grid lines
    for x in xc:
        ax.plot([x, x], [yc.min(), yc.max()], color="#d12f2f",
                linewidth=0.7, alpha=0.6, zorder=2)
    for y in yc:
        ax.plot([xc.min(), xc.max()], [y, y], color="#d12f2f",
                linewidth=0.7, alpha=0.6, zorder=2)
    ax.scatter(Xc.flatten(), Yc.flatten(), s=22, color="#d12f2f",
               marker="s", zorder=4,
               label=f"粗格子ノード ({NXC + 1}×{NYC + 1}, $\\delta x = 1$)")

    # Interface lines
    ax.axhline(NY / 2, color="#222", linestyle="--", linewidth=1.2, zorder=5)
    ax.axhline(NY / 2 - 1, color="#222", linestyle=":", linewidth=1.1, zorder=5)

    # Walls
    ax.plot([0, NX], [0, 0], color="black", linewidth=3.5, zorder=6)
    ax.plot([0, NX], [NY, NY], color="black", linewidth=3.5, zorder=6)
    # Lid arrow + label (above the wall, outside plot area on top)
    ax.annotate("", xy=(NX * 0.85, NY + 1.4), xytext=(NX * 0.15, NY + 1.4),
                arrowprops=dict(arrowstyle="->", color="black", linewidth=2.0))
    ax.text(NX * 0.5, NY + 2.5, "上壁: $u = u_t$ (動壁)", ha="center",
            fontsize=12, fontweight="bold")
    ax.text(NX * 0.5, -2.5, "下壁: $u = 0$ (固定壁)", ha="center",
            fontsize=12, fontweight="bold")

    # Periodic boundary annotation — moved further out
    for x_anchor in (-2.5, NX + 2.5):
        ax.annotate("", xy=(x_anchor, 2.0), xytext=(x_anchor, NY - 2.0),
                    arrowprops=dict(arrowstyle="<->", color="#0a7", linewidth=1.4))
    ax.text(-3.4, NY * 0.5, "x 周期境界", color="#0a7",
            fontsize=11, rotation=90, ha="center", va="center", fontweight="bold")
    ax.text(NX + 3.4, NY * 0.5, "x 周期境界", color="#0a7",
            fontsize=11, rotation=90, ha="center", va="center", fontweight="bold")

    # Labels for interface lines — to the right with arrows pointing in
    ax.annotate("interface\n粗 $y_c = 1$\n細 $y_f = n_{y,f}$",
                xy=(NX, NY / 2), xytext=(NX + 5.5, NY / 2 + 1.5),
                fontsize=10, ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#222", linewidth=0.9),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#222", lw=0.8))
    ax.annotate("overlap 下端\n粗 $y_c = 0$\n細 $y_f = n_{y,f} - 2$",
                xy=(NX, NY / 2 - 1), xytext=(NX + 5.5, NY / 2 - 2.5),
                fontsize=10, ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#222", linewidth=0.9),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#222", lw=0.8))

    # Region labels — placed away from interface clutter
    ax.text(NX * 0.5, NY * 0.28, "細格子領域 (下半分)",
            fontsize=13, ha="center", color="#1f6feb", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#1f6feb", alpha=0.85, lw=1.0))
    ax.text(NX * 0.5, NY * 0.78, "粗格子領域 (上半分)",
            fontsize=13, ha="center", color="#d12f2f", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#d12f2f", alpha=0.85, lw=1.0))
    # Overlap label — placed to the LEFT to avoid central node clutter
    ax.annotate("オーバーラップ層\n(両格子で重複)",
                xy=(NX * 0.4, NY / 2 - 0.5), xytext=(-9.0, NY / 2 - 0.5),
                fontsize=10, ha="left", va="center", color="#7a5b00", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#7a5b00", linewidth=0.9),
                bbox=dict(boxstyle="round,pad=0.3", fc="#fff2b3", ec="#7a5b00", lw=0.8))

    ax.set_xlim(-10.0, NX + 11.5)
    ax.set_ylim(-4.5, NY + 5.0)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$（物理単位）")
    ax.set_ylabel("$y$（物理単位）")
    ax.set_title(f"lbmblock.c: 多重格子 (Multi-block) 配置  "
                 f"$n_x \\times n_y = {NX}\\times{NY}$, refinement ratio $m={M}$",
                 fontsize=13)
    # Legend below the plot to avoid overlapping with top-wall label
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              fontsize=10, framealpha=0.95)
    ax.grid(False)

    out = RUN_DIR / "lbmblock_mesh.png"
    pub = PUBLISHED_ASSET_DIR / "lbmblock_mesh.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(pub, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")
    print(f"Saved {pub}")


def main() -> None:
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    plot_mesh()
    plot_distribution()


if __name__ == "__main__":
    main()
