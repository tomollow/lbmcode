"""Analytical-model schematic for src/sec6/iblbm2cicTRT.c.

Cylindrical (circular) Couette flow solved with an *implicit velocity
correction* immersed-boundary lattice Boltzmann method (IB-LBM) using the
TRT (two-relaxation-time) collision operator and the Guo (2002) forcing
term split consistently into even/odd parts.

Draws the periodic square fluid domain, the two concentric immersed
cylinders (outer stationary, inner rotating clockwise), the Lagrangian
marker points on each ring (ne factor 0.2, same as the MRT sibling), the
annular gap, the rotation arrows and the analytical tangential-velocity
profile.

Outputs:
    outputs/sec6/iblbm2cicTRT/iblbm2cicTRT_schematic.png
    docs/assets/sec6/iblbm2cicTRT_schematic.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec6" / "iblbm2cicTRT"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec6"

# Source parameters (iblbm2cicTRT.c).
NX = NY = 50
RP_OUT = 70.0 / 200.0 * NX   # rp[0] = 17.5 (outer, stationary)
RP_IN = 45.0 / 200.0 * NX    # rp[1] = 11.25 (inner, rotating)
U0 = 0.01
CENTER = NX / 2.0            # 25.0
NE_FACTOR = 0.2              # ne = (int)(2*pi*R*0.2)  -> 21 / 14 points

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(-8, 58)
    axis.set_ylim(-10, 58)
    axis.set_aspect("equal")
    axis.axis("off")

    fluid_color = "#cfe3f2"
    outer_color = "#555555"
    inner_color = "#c0392b"
    periodic_color = "#7d7d7d"

    # Fluid domain fill + periodic boundary (dashed).
    axis.add_patch(plt.Rectangle((0, 0), NX, NY, facecolor=fluid_color,
                                 alpha=0.5, edgecolor="none", zorder=0))
    axis.add_patch(plt.Rectangle((0, 0), NX, NY, facecolor="none",
                                 edgecolor=periodic_color, linewidth=1.6,
                                 linestyle=(0, (6, 4)), zorder=1))

    # Outer (stationary) cylinder ring.
    axis.add_patch(plt.Circle((CENTER, CENTER), RP_OUT, facecolor="none",
                              edgecolor=outer_color, linewidth=2.4, zorder=3))
    # Inner (rotating) cylinder ring + filled core.
    axis.add_patch(plt.Circle((CENTER, CENTER), RP_IN, facecolor=inner_color,
                              alpha=0.18, edgecolor=inner_color, linewidth=2.4,
                              zorder=3))

    # Lagrangian marker points (use the same counts as the source).
    ne_out = int(2.0 * np.pi * RP_OUT * NE_FACTOR)
    ne_in = int(2.0 * np.pi * RP_IN * NE_FACTOR)
    th = np.linspace(0, 2 * np.pi, ne_out, endpoint=False)
    axis.scatter(CENTER + RP_OUT * np.cos(th), CENTER + RP_OUT * np.sin(th),
                 s=22, color=outer_color, zorder=4)
    th = np.linspace(0, 2 * np.pi, ne_in, endpoint=False)
    axis.scatter(CENTER + RP_IN * np.cos(th), CENTER + RP_IN * np.sin(th),
                 s=22, color=inner_color, zorder=4)

    # Clockwise rotation arrows on the inner ring.
    for ang in np.deg2rad([45, 135, 225, 315]):
        x = CENTER + RP_IN * np.cos(ang)
        y = CENTER + RP_IN * np.sin(ang)
        # Clockwise tangent direction: (sin, -cos).
        tx, ty = np.sin(ang), -np.cos(ang)
        axis.annotate("", xy=(x + 3.0 * tx, y + 3.0 * ty), xytext=(x, y),
                      arrowprops={"arrowstyle": "-|>", "lw": 1.6,
                                  "color": inner_color}, zorder=5)

    # Radii arrows.
    axis.annotate("", xy=(CENTER + RP_IN * np.cos(np.deg2rad(20)),
                          CENTER + RP_IN * np.sin(np.deg2rad(20))),
                  xytext=(CENTER, CENTER),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.3, "color": "black"})
    axis.text(CENTER + 0.55 * RP_IN, CENTER + 2.0, r"$R_i$", fontsize=13)
    axis.annotate("", xy=(CENTER + RP_OUT * np.cos(np.deg2rad(-25)),
                          CENTER + RP_OUT * np.sin(np.deg2rad(-25))),
                  xytext=(CENTER, CENTER),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.3, "color": "black"})
    axis.text(CENTER + 0.62 * RP_OUT, CENTER - 5.5, r"$R_o$", fontsize=13)

    # Labels.
    axis.text(CENTER, CENTER, r"$\omega$ (時計回り)", ha="center", va="center",
              fontsize=10, color=inner_color)
    axis.text(CENTER, CENTER + RP_OUT + 2.2, "外円筒 (静止)", ha="center",
              va="bottom", fontsize=11, color=outer_color)
    axis.text(CENTER, CENTER + RP_IN - 3.2, "内円筒\n(回転 $u_0$)", ha="center",
              va="center", fontsize=10, color=inner_color)

    # Periodic BC labels.
    axis.text(NX / 2, NY + 1.0, "周期境界", ha="center", va="bottom",
              fontsize=10, color=periodic_color)
    axis.text(-1.0, NY / 2, "周期境界", ha="right", va="center",
              fontsize=10, color=periodic_color, rotation=90)

    # Analytical-solution annotation box.
    axis.text(NX + 2.5, NY * 0.82,
              "回転 Couette 流\n(Stokes 厳密解)", ha="left", va="center",
              fontsize=11)
    axis.text(NX + 2.5, NY * 0.59,
              r"$u_\theta(r)=u_0\dfrac{r/R_o-R_o/r}{R_i/R_o-R_o/R_i}$",
              ha="left", va="center", fontsize=12, color=inner_color)
    axis.text(NX + 2.5, NY * 0.40,
              r"$R_i \leq r \leq R_o$", ha="left", va="center", fontsize=10,
              color="0.3")
    axis.text(NX + 2.5, NY * 0.28,
              fr"$R_i={RP_IN:.2f},\ R_o={RP_OUT:.2f}$", ha="left", va="center",
              fontsize=10, color="0.3")
    axis.text(NX + 2.5, NY * 0.16,
              r"$u_0=0.01,\ \tau_+=10,\ \Lambda=9/8$ (TRT)", ha="left",
              va="center", fontsize=10, color="0.3")

    # Coordinate axes.
    ox, oy = -6.0, -7.0
    axis.annotate("", xy=(ox + 6, oy), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.annotate("", xy=(ox, oy + 6), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.text(ox + 6.5, oy, r"$x$", ha="left", va="center", fontsize=11)
    axis.text(ox, oy + 6.5, r"$y$", ha="center", va="bottom", fontsize=11)

    # Domain-size label.
    axis.annotate("", xy=(NX, -4.0), xytext=(0, -4.0),
                  arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.2"})
    axis.text(NX / 2, -5.5, r"$n_x = n_y = 50$", ha="center", va="top",
              fontsize=10)


def main() -> None:
    figure, axis = plt.subplots(figsize=(8.4, 6.4), constrained_layout=True)
    draw_schematic(axis)
    axis.set_title(
        r"図 6.4 円筒 Couette 流の解析モデル（陰的速度補正 IB-LBM, TRT）",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "iblbm2cicTRT_schematic.png",
                 PUBLISHED_ASSET_DIR / "iblbm2cicTRT_schematic.png"):
        figure.savefig(path, dpi=220, bbox_inches="tight")
        print(f"Saved schematic to {path}")


if __name__ == "__main__":
    main()
