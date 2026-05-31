"""Analytical-model schematic for src/sec5/lbmzalesak.c (Zalesak's rotating disk).

Draws the periodic square domain, the slotted disk initialised at the centre,
the rigid-body (solid-body) rotation velocity field, the centre of rotation,
and the key parameters (radius R, slot, Peclet number, one-revolution period).

Outputs:
    outputs/sec5/lbmzalesak/lbmzalesak_schematic.png
    docs/assets/sec5/lbmzalesak_schematic.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Wedge


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmzalesak"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"

# Source parameters (lbmzalesak.c).
NX = NY = 50
RADIUS = NX * 0.4            # 20.0
SLOT_I = (23, 26)           # 93*nx/200 .. 107*nx/200 with integer division
SLOT_J = (2, NY // 2)       # 2 .. 25


def slotted_disk_patch(cx: float, cy: float, scale: float):
    """Return matplotlib path vertices (normalised) for the slotted disk."""
    # Build a filled indicator on a fine grid, then take its 0.5 contour as a
    # polygon for a crisp schematic fill.
    n = 400
    xs = np.linspace(0, NX, n)
    ys = np.linspace(0, NY, n)
    gx, gy = np.meshgrid(xs, ys)
    inside = np.sqrt((gx - NX * 0.5) ** 2 + (gy - NY * 0.5) ** 2) <= RADIUS
    slot = (gx >= SLOT_I[0]) & (gx <= SLOT_I[1]) & (gy >= SLOT_J[0]) & (gy <= SLOT_J[1])
    return (gx, gy, inside & ~slot)


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(-0.22, 1.34)
    axis.set_ylim(-0.20, 1.16)
    axis.set_aspect("equal")
    axis.axis("off")

    disk_color = "#d6604d"
    outside_color = "#4393c3"
    periodic_color = "#7d7d7d"
    rot_color = "#2c7a3f"

    dom_l, dom_b, dom_r, dom_t = 0.0, 0.0, 1.0, 1.0

    # Outside phase fill + periodic boundary (dashed = periodic).
    axis.add_patch(plt.Rectangle((dom_l, dom_b), 1.0, 1.0,
                                 facecolor=outside_color, alpha=0.14,
                                 edgecolor="none"))
    axis.add_patch(plt.Rectangle((dom_l, dom_b), 1.0, 1.0,
                                 facecolor="none", edgecolor=periodic_color,
                                 linewidth=1.6, linestyle=(0, (6, 4))))

    # Slotted disk (normalised coordinates: divide lattice index by NX).
    gx, gy, mask = slotted_disk_patch(0.5, 0.5, 1.0)
    axis.contourf(gx / NX, gy / NY, mask.astype(float), levels=[0.5, 1.5],
                  colors=[disk_color], alpha=0.55)
    axis.contour(gx / NX, gy / NY, mask.astype(float), levels=[0.5],
                 colors=[disk_color], linewidths=1.8)

    # Phase labels.
    axis.text(0.5, 0.70, r"円盤 $\phi=+1$", ha="center", va="center",
              fontsize=12, color="#7a1f12")
    axis.text(0.12, 0.90, r"外部 $\phi=-1$", ha="left", va="center",
              fontsize=11, color="#1f5d8c")
    axis.text(0.5, 0.255, r"切欠き (slot)", ha="center", va="center",
              fontsize=10, color="#7a1f12")

    # Centre of rotation.
    axis.plot([0.5], [0.5], "o", color="black", ms=4)
    axis.text(0.515, 0.515, r"回転中心", ha="left", va="bottom",
              fontsize=9, color="0.2")

    # Radius arrow (toward upper-left so it avoids the slot).
    ang = np.deg2rad(125.0)
    axis.annotate("", xy=(0.5 + 0.4 * np.cos(ang), 0.5 + 0.4 * np.sin(ang)),
                  xytext=(0.5, 0.5),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.5, "color": "black"})
    axis.text(0.5 + 0.24 * np.cos(ang) - 0.02, 0.5 + 0.24 * np.sin(ang),
              r"$R=0.4\,n_x$", ha="right", va="center", fontsize=12)

    # Solid-body rotation velocity field (sampled arrows, CCW).
    u0 = 0.04
    pts = []
    for r in (0.30, 0.46):
        for k in range(8):
            th = 2 * np.pi * k / 8
            pts.append((0.5 + r * np.cos(th), 0.5 + r * np.sin(th)))
    for (px, py) in pts:
        # u = -omega (y-yc), v = omega (x-xc)  ->  CCW
        ux, vy = -(py - 0.5), (px - 0.5)
        mag = np.hypot(ux, vy)
        ux, vy = ux / mag * 0.075, vy / mag * 0.075
        axis.add_patch(FancyArrowPatch((px, py), (px + ux, py + vy),
                       arrowstyle="-|>", mutation_scale=10,
                       color=rot_color, lw=1.1, alpha=0.85))

    # CCW indicator arc (top-right corner).
    axis.add_patch(Wedge((1.16, 0.92), 0.10, 20, 300, width=0.001,
                         edgecolor=rot_color, facecolor="none", lw=1.6))
    axis.add_patch(FancyArrowPatch((1.16 + 0.10 * np.cos(np.deg2rad(300)),
                                    0.92 + 0.10 * np.sin(np.deg2rad(300))),
                                   (1.16 + 0.10 * np.cos(np.deg2rad(285)),
                                    0.92 + 0.10 * np.sin(np.deg2rad(285))),
                   arrowstyle="-|>", mutation_scale=12, color=rot_color, lw=1.6))
    axis.text(1.16, 0.92, "CCW", ha="center", va="center", fontsize=9,
              color=rot_color)

    # Velocity field annotation.
    axis.text(1.04, 0.60, r"剛体回転場", ha="left", va="center",
              fontsize=12, color=rot_color)
    axis.text(1.04, 0.50,
              r"$u = -u_0\pi\!\left(\dfrac{j}{n_y}-\dfrac{1}{2}\right)$",
              ha="left", va="center", fontsize=11, color=rot_color)
    axis.text(1.04, 0.40,
              r"$v = +u_0\pi\!\left(\dfrac{i}{n_y}-\dfrac{1}{2}\right)$",
              ha="left", va="center", fontsize=11, color=rot_color)
    axis.text(1.04, 0.30,
              r"$\omega = u_0\pi/n_y$", ha="left", va="center",
              fontsize=11, color="0.3")
    axis.text(1.04, 0.21,
              r"1 周 $=2\pi/\omega=2n_y/u_0=2500$ 歩",
              ha="left", va="center", fontsize=10, color="0.3")
    axis.text(1.04, 0.11, r"$\mathrm{Pe}=400$", ha="left", va="center",
              fontsize=11, color="0.3")

    # Periodic BC labels.
    axis.text(0.5, dom_t + 0.04, "周期境界", ha="center", va="bottom",
              fontsize=10, color=periodic_color)
    axis.text(dom_l - 0.03, 0.5, "周期境界", ha="right", va="center",
              fontsize=10, color=periodic_color, rotation=90)

    # Coordinate axes.
    ox, oy = -0.16, -0.10
    axis.annotate("", xy=(ox + 0.13, oy), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.annotate("", xy=(ox, oy + 0.13), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.text(ox + 0.14, oy, r"$x\,(i)$", ha="left", va="center", fontsize=11)
    axis.text(ox, oy + 0.14, r"$y\,(j)$", ha="center", va="bottom", fontsize=11)

    # Domain size label.
    axis.annotate("", xy=(dom_r, dom_b - 0.07), xytext=(dom_l, dom_b - 0.07),
                  arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.2"})
    axis.text(0.5, dom_b - 0.10, r"$n_x = n_y = 50$", ha="center", va="top",
              fontsize=10)


def main() -> None:
    figure, axis = plt.subplots(figsize=(8.0, 6.2), constrained_layout=True)
    draw_schematic(axis)
    axis.set_title(
        r"図 5.0 Zalesak の円盤と剛体回転場（相場 D2Q9 LBM, 移流テスト）",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "lbmzalesak_schematic.png",
                 PUBLISHED_ASSET_DIR / "lbmzalesak_schematic.png"):
        figure.savefig(path, dpi=220, bbox_inches="tight")
        print(f"Saved schematic to {path}")


if __name__ == "__main__":
    main()
