"""Schematic for the lbmtherm.c advection-diffusion benchmark.

Three-panel layout:
  (a) Channel geometry: periodic in x, walls at y = 0 and y = h, uniform
      horizontal carrier velocity u_0, cosine-modulated wall temperature
      T_wall(x) = cos(k x).
  (b) Sub-grid wall offset q: zoom showing the lattice column near the
      bottom wall (ghost row j = 0, first interior j = 1, second j = 2)
      and the wall a distance q above the ghost row.
  (c) Analytical temperature mode T_a(x, y) at Pe = 20 for orientation.

Run from the repo root:
    python scripts/plot_lbmtherm_schematic.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec3" / "lbmtherm"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


# Parameters used purely for illustration (match lbmtherm.c defaults).
NX = 64
NY = 64
PE = 20.0
TAUG = 0.56
Q = 0.7


def draw_layout(ax: plt.Axes) -> None:
    """Panel (a): channel geometry, periodic BC, wall temperature."""
    ax.set_xlim(-0.20, 1.25)
    ax.set_ylim(-0.20, 1.30)
    ax.set_aspect("equal")
    ax.axis("off")

    flow_color = "0.20"

    channel_l, channel_r = 0.08, 0.92
    channel_b, channel_t = 0.18, 0.82
    wall_thickness = 0.045

    # Cosine-modulated wall temperature: paint with diverging colormap along x.
    n_band = 200
    xs = np.linspace(channel_l, channel_r, n_band)
    cos_vals = np.cos(2 * np.pi * (xs - channel_l) / (channel_r - channel_l))
    cmap = plt.get_cmap("coolwarm")

    for k in range(n_band - 1):
        c = cmap(0.5 + 0.5 * cos_vals[k])
        ax.add_patch(Rectangle((xs[k], channel_t),
                               xs[k + 1] - xs[k], wall_thickness,
                               facecolor=c, edgecolor="none"))
        ax.add_patch(Rectangle((xs[k], channel_b - wall_thickness),
                               xs[k + 1] - xs[k], wall_thickness,
                               facecolor=c, edgecolor="none"))

    for y_edge in (channel_t, channel_t + wall_thickness,
                   channel_b, channel_b - wall_thickness):
        ax.plot([channel_l, channel_r], [y_edge, y_edge],
                color="black", linewidth=0.6)

    # Periodic side edges (dashed).
    ax.plot([channel_l, channel_l],
            [channel_b - wall_thickness, channel_t + wall_thickness],
            color="black", linewidth=0.8, linestyle="--")
    ax.plot([channel_r, channel_r],
            [channel_b - wall_thickness, channel_t + wall_thickness],
            color="black", linewidth=0.8, linestyle="--")

    # Uniform advection arrows (carrier velocity u_0).
    y_arrows = np.linspace(channel_b + 0.10, channel_t - 0.10, 4)
    for y_arr in y_arrows:
        ax.annotate(
            "", xy=(channel_l + 0.55, y_arr), xytext=(channel_l + 0.10, y_arr),
            arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": flow_color},
        )
    ax.text((channel_l + channel_r) / 2, channel_b + 0.02,
            r"$\mathbf{u} = (u_0,\; 0)$ (given)",
            ha="center", va="bottom", fontsize=10, color=flow_color)

    # Periodic BC indication.
    ax.text(channel_l - 0.05, (channel_b + channel_t) / 2,
            "周期", ha="right", va="center", fontsize=10, color="0.30",
            rotation=90)
    ax.text(channel_r + 0.05, (channel_b + channel_t) / 2,
            "周期", ha="left", va="center", fontsize=10, color="0.30",
            rotation=90)

    # Wall labels.
    ax.text((channel_l + channel_r) / 2, channel_t + wall_thickness + 0.04,
            r"上壁: $T_{\mathrm{wall}}(x) = \cos(k x)$",
            ha="center", va="bottom", fontsize=11)
    ax.text((channel_l + channel_r) / 2, channel_b - wall_thickness - 0.04,
            r"下壁: $T_{\mathrm{wall}}(x) = \cos(k x)$",
            ha="center", va="top", fontsize=11)

    # Channel height.
    ax.annotate(
        "", xy=(channel_l - 0.03, channel_t), xytext=(channel_l - 0.03, channel_b),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.20"},
    )
    ax.text(channel_l - 0.05, (channel_b + channel_t) / 2,
            r"$h$", ha="right", va="center", fontsize=11)

    # Wavelength label.
    ax.annotate(
        "", xy=(channel_r, channel_t + 0.21), xytext=(channel_l, channel_t + 0.21),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.20"},
    )
    ax.text((channel_l + channel_r) / 2, channel_t + 0.22,
            r"$L_x = 2\pi / k$",
            ha="center", va="bottom", fontsize=11)

    # Coordinate axes.
    ox, oy = -0.15, -0.10
    ax.annotate("", xy=(ox + 0.13, oy), xytext=(ox, oy),
                arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    ax.annotate("", xy=(ox, oy + 0.13), xytext=(ox, oy),
                arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    ax.text(ox + 0.14, oy, r"$x$", ha="left", va="center", fontsize=11)
    ax.text(ox, oy + 0.14, r"$y$", ha="center", va="bottom", fontsize=11)

    ax.set_title("(a) チャネルと境界条件", fontsize=11)


def draw_subgrid(ax: plt.Axes) -> None:
    """Panel (b): zoom showing sub-grid wall offset q."""
    ax.set_xlim(-0.10, 1.10)
    ax.set_ylim(-0.10, 1.10)
    ax.set_aspect("equal")
    ax.axis("off")

    # Lattice columns 0..4 along x.
    node_xs = np.linspace(0.10, 0.90, 5)
    y_ghost = 0.18
    y_first = 0.40
    y_second = 0.55
    y_third = 0.70
    y_wall = y_ghost + Q * (y_first - y_ghost)

    # Faint lattice grid.
    for nx_ in node_xs:
        ax.plot([nx_, nx_], [y_ghost, y_third + 0.10],
                color="0.85", linewidth=0.5)
    for y_ in (y_ghost, y_first, y_second, y_third):
        ax.plot([node_xs[0] - 0.04, node_xs[-1] + 0.04], [y_, y_],
                color="0.85", linewidth=0.5)

    # Nodes.
    for nx_ in node_xs:
        ax.plot(nx_, y_ghost, "o", markersize=7,
                markerfacecolor="white", markeredgecolor="0.40")
        ax.plot(nx_, y_first, "o", markersize=7,
                markerfacecolor="0.30", markeredgecolor="0.30")
        ax.plot(nx_, y_second, "o", markersize=7,
                markerfacecolor="0.30", markeredgecolor="0.30")
        ax.plot(nx_, y_third, "o", markersize=7,
                markerfacecolor="0.30", markeredgecolor="0.30")

    ax.text(node_xs[-1] + 0.06, y_ghost, r"$j = 0$ (ghost)",
            ha="left", va="center", fontsize=10, color="0.40")
    ax.text(node_xs[-1] + 0.06, y_first, r"$j = 1$",
            ha="left", va="center", fontsize=10)
    ax.text(node_xs[-1] + 0.06, y_second, r"$j = 2$",
            ha="left", va="center", fontsize=10)
    ax.text(node_xs[-1] + 0.06, y_third, r"$j = 3$",
            ha="left", va="center", fontsize=10)

    # The wall.
    ax.plot([node_xs[0] - 0.05, node_xs[-1] + 0.05], [y_wall, y_wall],
            color="#c0392b", linewidth=2.4)
    ax.text(node_xs[0] - 0.05, y_wall + 0.025, r"壁面 ($T_{\mathrm{wall}}$)",
            ha="left", va="bottom", fontsize=10, color="#c0392b")

    # q distance arrow.
    qx = node_xs[1]
    ax.annotate(
        "", xy=(qx, y_wall), xytext=(qx, y_ghost),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.1, "color": "0.20"},
    )
    ax.text(qx - 0.04, 0.5 * (y_wall + y_ghost),
            r"$q$", ha="right", va="center", fontsize=12)

    # Distance from wall to first interior node = 1 - q.
    qx2 = node_xs[2]
    ax.annotate(
        "", xy=(qx2, y_first), xytext=(qx2, y_wall),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.1, "color": "0.20"},
    )
    ax.text(qx2 + 0.04, 0.5 * (y_wall + y_first),
            r"$1 - q$", ha="left", va="center", fontsize=10)

    # Lattice spacing.
    lx0 = node_xs[3]; lx1 = node_xs[4]
    ax.annotate(
        "", xy=(lx1, y_third + 0.05), xytext=(lx0, y_third + 0.05),
        arrowprops={"arrowstyle": "<|-|>", "lw": 0.8, "color": "0.30"},
    )
    ax.text(0.5 * (lx0 + lx1), y_third + 0.06,
            r"$\Delta x = 1$", ha="center", va="bottom", fontsize=9,
            color="0.30")

    ax.text(0.5, 1.0,
            r"$q$: 壁面と最隣接 ghost 行との距離 (コード既定 $q = 0.7$)",
            ha="center", va="top", fontsize=10, color="0.20")
    ax.text(0.5, 0.04,
            r"半過程 BB は $q = 0.5$ で 2 次精度。$q \neq 0.5$ では 3 セル補間",
            ha="center", va="top", fontsize=9, color="0.20")

    ax.set_title("(b) サブグリッド壁オフセット $q$", fontsize=11)


def draw_analytical(ax: plt.Axes) -> None:
    """Panel (c): analytical temperature mode T_a(x, y)."""
    chi = (TAUG - 0.5) / 3.0
    h = (NY - 2) + 2.0 * Q
    u0 = PE * chi / h
    k = 2.0 * np.pi / NX
    beta = k * np.sqrt(1.0 + 1j * u0 / (chi * k))

    x = np.linspace(0.0, NX, 200)
    y = np.linspace(0.0, h, 200)
    xx, yy = np.meshgrid(x, y)
    phase = np.exp(1j * k * xx)
    num = np.sinh(beta * yy) + np.sinh(beta * (h - yy))
    den = np.sinh(beta * h)
    t_field = np.real(phase * num / den)

    cs = ax.contourf(xx / NX, yy / h, t_field,
                     levels=np.linspace(-1.0, 1.0, 21), cmap="coolwarm")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x / L_x$")
    ax.set_ylabel(r"$y / h$")
    ax.set_title(r"(c) 解析解 $T_a$ ($Pe = 20$)", fontsize=11)
    plt.colorbar(cs, ax=ax, shrink=0.85, label=r"$T_a$")


def main() -> None:
    fig = plt.figure(figsize=(15.5, 5.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.05])

    ax_a = fig.add_subplot(gs[0, 0])
    draw_layout(ax_a)

    ax_b = fig.add_subplot(gs[0, 1])
    draw_subgrid(ax_b)

    ax_c = fig.add_subplot(gs[0, 2])
    draw_analytical(ax_c)

    fig.suptitle(
        r"lbmtherm.c: 周期チャネル advection-diffusion ($Pe = u_0 h / \chi$, $k = 2\pi/L_x$)",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    out_local = RUN_DIR / "lbmtherm_schematic.png"
    out_pub = PUBLISHED_ASSET_DIR / "lbmtherm_schematic.png"
    fig.savefig(out_local, dpi=220, bbox_inches="tight")
    fig.savefig(out_pub, dpi=220, bbox_inches="tight")
    print(f"Saved {out_local}")
    print(f"Saved {out_pub}")


if __name__ == "__main__":
    main()
