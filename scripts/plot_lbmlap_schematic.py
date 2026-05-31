"""Analytical-model schematic for src/sec5/lbmlap.c (Laplace's law).

Draws a periodic square domain with a circular droplet of radius R, the
diffuse interface of width W, the inside/outside order-parameter values, the
pressure jump Delta p = sigma / R, and the coordinate system.

Outputs:
    outputs/sec5/lbmlap/lbmlap_schematic.png
    docs/assets/sec5/lbmlap_schematic.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmlap"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(-0.20, 1.30)
    axis.set_ylim(-0.22, 1.18)
    axis.set_aspect("equal")
    axis.axis("off")

    inside_color = "#d6604d"
    outside_color = "#4393c3"
    periodic_color = "#7d7d7d"

    dom_l, dom_b, dom_r, dom_t = 0.05, 0.05, 0.95, 0.95

    # Outside phase fill.
    axis.add_patch(plt.Rectangle((dom_l, dom_b), dom_r - dom_l, dom_t - dom_b,
                                 facecolor=outside_color, alpha=0.18,
                                 edgecolor="none"))
    # Periodic domain boundary (dashed = periodic).
    axis.add_patch(plt.Rectangle((dom_l, dom_b), dom_r - dom_l, dom_t - dom_b,
                                 facecolor="none", edgecolor=periodic_color,
                                 linewidth=1.6, linestyle=(0, (6, 4))))

    cx, cy = 0.5, 0.5
    radius = 0.27
    width = 0.06

    # Diffuse interface band (annulus, dotted guides at R +/- W/2).
    axis.add_patch(plt.Circle((cx, cy), radius + width / 2, facecolor="none",
                              edgecolor="0.5", linewidth=0.8, linestyle=":"))
    axis.add_patch(plt.Circle((cx, cy), radius - width / 2, facecolor="none",
                              edgecolor="0.5", linewidth=0.8, linestyle=":"))
    # Droplet interior.
    axis.add_patch(plt.Circle((cx, cy), radius, facecolor=inside_color,
                              alpha=0.45, edgecolor=inside_color, linewidth=1.8))

    # Phase labels.
    axis.text(cx, cy + 0.08, r"液滴内部", ha="center", va="center",
              fontsize=12, color="#7a1f12")
    axis.text(cx, cy - 0.02, r"$\phi = +\phi_0$", ha="center", va="center",
              fontsize=12, color="#7a1f12")
    axis.text(cx, cy - 0.10, r"$p_{\mathrm{in}}$", ha="center", va="center",
              fontsize=12, color="#7a1f12")
    axis.text(dom_l + 0.10, dom_t - 0.08, r"外部 $\phi=-\phi_0$,  $p_{\mathrm{out}}$",
              ha="left", va="center", fontsize=11, color="#1f5d8c")

    # Radius arrow.
    ang = np.deg2rad(35.0)
    axis.annotate("", xy=(cx + radius * np.cos(ang), cy + radius * np.sin(ang)),
                  xytext=(cx, cy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.5, "color": "black"})
    axis.text(cx + 0.45 * radius * np.cos(ang) + 0.01,
              cy + 0.45 * radius * np.sin(ang) + 0.03,
              r"$R$", ha="left", va="bottom", fontsize=13)

    # Interface width callout.
    axis.annotate("", xy=(cx + radius + width / 2, cy - 0.32),
                  xytext=(cx + radius - width / 2, cy - 0.32),
                  arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.25"})
    axis.text(cx + radius, cy - 0.30, r"界面厚さ $W$", ha="center", va="bottom",
              fontsize=10, color="0.25")
    axis.annotate("", xy=(cx + radius, cy - 0.30),
                  xytext=(cx + radius, cy - 0.06 - 0.0),
                  arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "0.5",
                              "linestyle": "dotted"})

    # Laplace law annotation.
    axis.text(1.02, 0.62, r"Laplace の法則", ha="left", va="center",
              fontsize=12, color="black")
    axis.text(1.02, 0.52, r"$\Delta p = p_{\mathrm{in}} - p_{\mathrm{out}} = \dfrac{\sigma}{R}$",
              ha="left", va="center", fontsize=13, color="#7a1f12")
    axis.text(1.02, 0.40, r"$p = \rho c_s^2 = \rho/3$", ha="left", va="center",
              fontsize=11, color="0.3")
    axis.text(1.02, 0.32, r"$\sigma$ : 表面張力", ha="left", va="center",
              fontsize=10, color="0.3")

    # Periodic BC labels.
    axis.text(0.5, dom_t + 0.04, "周期境界", ha="center", va="bottom",
              fontsize=10, color=periodic_color)
    axis.text(dom_l - 0.03, 0.5, "周期境界", ha="right", va="center",
              fontsize=10, color=periodic_color, rotation=90)

    # Coordinate axes.
    ox, oy = -0.14, -0.12
    axis.annotate("", xy=(ox + 0.13, oy), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.annotate("", xy=(ox, oy + 0.13), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.text(ox + 0.14, oy, r"$x$", ha="left", va="center", fontsize=11)
    axis.text(ox, oy + 0.14, r"$y$", ha="center", va="bottom", fontsize=11)

    # Domain size label.
    axis.annotate("", xy=(dom_r, dom_b - 0.06), xytext=(dom_l, dom_b - 0.06),
                  arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.2"})
    axis.text(0.5, dom_b - 0.085, r"$n_x = n_y = 50$", ha="center", va="top",
              fontsize=10)


def main() -> None:
    figure, axis = plt.subplots(figsize=(7.6, 6.2), constrained_layout=True)
    draw_schematic(axis)
    axis.set_title(
        r"図 5.0 静止液滴と Laplace の法則（自由エネルギー型 D2Q9 二相 LBM）",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    for path in (RUN_DIR / "lbmlap_schematic.png",
                 PUBLISHED_ASSET_DIR / "lbmlap_schematic.png"):
        figure.savefig(path, dpi=220, bbox_inches="tight")
        print(f"Saved schematic to {path}")


if __name__ == "__main__":
    main()
