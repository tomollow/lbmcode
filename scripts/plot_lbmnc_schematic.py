from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec3" / "lbmnc"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"


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


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(-0.18, 1.25)
    axis.set_ylim(-0.28, 1.15)
    axis.set_aspect("equal")
    axis.axis("off")

    hot_color = "#d6604d"
    cold_color = "#4393c3"
    adiabatic_color = "#b5b5b5"
    flow_color = "0.20"

    cavity_l, cavity_b = 0.10, 0.10
    cavity_r, cavity_t = 0.90, 0.90
    wall_thickness = 0.04

    # Hot left wall.
    axis.add_patch(plt.Rectangle((cavity_l - wall_thickness, cavity_b - wall_thickness),
                                 wall_thickness,
                                 (cavity_t - cavity_b) + 2 * wall_thickness,
                                 facecolor=hot_color, edgecolor="black", linewidth=0.7))
    # Cold right wall.
    axis.add_patch(plt.Rectangle((cavity_r, cavity_b - wall_thickness),
                                 wall_thickness,
                                 (cavity_t - cavity_b) + 2 * wall_thickness,
                                 facecolor=cold_color, edgecolor="black", linewidth=0.7))
    # Adiabatic top.
    axis.add_patch(plt.Rectangle((cavity_l - wall_thickness, cavity_t),
                                 (cavity_r - cavity_l) + 2 * wall_thickness,
                                 wall_thickness,
                                 facecolor=adiabatic_color, edgecolor="black", linewidth=0.7,
                                 hatch="//"))
    # Adiabatic bottom.
    axis.add_patch(plt.Rectangle((cavity_l - wall_thickness, cavity_b - wall_thickness),
                                 (cavity_r - cavity_l) + 2 * wall_thickness,
                                 wall_thickness,
                                 facecolor=adiabatic_color, edgecolor="black", linewidth=0.7,
                                 hatch="//"))

    # Cavity interior outline.
    axis.add_patch(plt.Rectangle((cavity_l, cavity_b),
                                 cavity_r - cavity_l, cavity_t - cavity_b,
                                 facecolor="white", edgecolor="black", linewidth=1.2))

    # Recirculation arrow (counter-clockwise: up on hot side, down on cold side).
    cx, cy = 0.50, 0.50
    radius = 0.22
    theta = np.linspace(np.deg2rad(190.0), np.deg2rad(-100.0), 220)
    x_curve = cx + radius * np.cos(theta)
    y_curve = cy + radius * np.sin(theta)
    axis.plot(x_curve, y_curve, color=flow_color, linewidth=1.6)
    axis.annotate(
        "",
        xy=(x_curve[-1], y_curve[-1]),
        xytext=(x_curve[-10], y_curve[-10]),
        arrowprops={"arrowstyle": "-|>", "lw": 1.6, "color": flow_color},
    )
    axis.text(cx, cy, "主循環", ha="center", va="center", fontsize=11, color=flow_color)

    # Up/down arrows near the side walls indicating buoyant rise and sink.
    axis.annotate(
        "",
        xy=(0.20, 0.78), xytext=(0.20, 0.22),
        arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": hot_color},
    )
    axis.text(0.13, 0.50, "上昇", ha="right", va="center", fontsize=10, color=hot_color, rotation=90)

    axis.annotate(
        "",
        xy=(0.80, 0.22), xytext=(0.80, 0.78),
        arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": cold_color},
    )
    axis.text(0.87, 0.50, "下降", ha="left", va="center", fontsize=10, color=cold_color, rotation=90)

    # Wall labels.
    axis.text(cavity_l - wall_thickness - 0.02, cavity_t + 0.03,
              r"高温壁 $T = 1$", ha="right", va="bottom", fontsize=11, color=hot_color)
    axis.text(cavity_r + wall_thickness + 0.02, cavity_t + 0.03,
              r"低温壁 $T = 0$", ha="left", va="bottom", fontsize=11, color=cold_color)
    axis.text(0.50, cavity_t + wall_thickness + 0.02,
              r"断熱壁 $\partial T/\partial y = 0$",
              ha="center", va="bottom", fontsize=10, color="0.30")
    axis.text(0.50, cavity_b - wall_thickness - 0.02,
              r"断熱壁 $\partial T/\partial y = 0$",
              ha="center", va="top", fontsize=10, color="0.30")

    # No-slip annotation (placed below cavity-length label).
    axis.text(0.50, -0.24, r"全壁面で $\mathbf{u} = 0$（half-way bounce-back）",
              ha="center", va="bottom", fontsize=9, color="0.30")

    # Gravity arrow.
    gx, gy_top = 1.07, 0.78
    axis.annotate(
        "",
        xy=(gx, gy_top - 0.30), xytext=(gx, gy_top),
        arrowprops={"arrowstyle": "-|>", "lw": 1.6, "color": "black"},
    )
    axis.text(gx + 0.02, gy_top - 0.15, r"$\mathbf{g}$", ha="left", va="center", fontsize=12)

    # Coordinate axes (placed outside the cavity, lower-left).
    ox, oy = -0.12, -0.16
    axis.annotate("", xy=(ox + 0.13, oy), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.annotate("", xy=(ox, oy + 0.13), xytext=(ox, oy),
                  arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"})
    axis.text(ox + 0.14, oy, r"$x$", ha="left", va="center", fontsize=11)
    axis.text(ox, oy + 0.14, r"$y$", ha="center", va="bottom", fontsize=11)

    # Cavity-side length label (placed below the bottom adiabatic wall).
    axis.annotate(
        "",
        xy=(cavity_r, -0.13), xytext=(cavity_l, -0.13),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.20"},
    )
    axis.text(0.50, -0.15, r"$L = n_x - 1$", ha="center", va="top", fontsize=10)
    axis.annotate(
        "",
        xy=(cavity_l - 0.10, cavity_t), xytext=(cavity_l - 0.10, cavity_b),
        arrowprops={"arrowstyle": "<|-|>", "lw": 1.0, "color": "0.20"},
    )
    axis.text(cavity_l - 0.12, 0.50, r"$L$", ha="right", va="center", fontsize=10)


def main() -> None:
    figure, axis = plt.subplots(figsize=(7.0, 6.4), constrained_layout=True)
    draw_schematic(axis)
    axis.set_title(
        r"側面加熱正方キャビティの自然対流（$Ra = 10^{4}$, $Pr = 0.71$）",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmnc_schematic.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmnc_schematic.png"
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved schematic to {output_path}")
    print(f"Saved schematic to {published_path}")


if __name__ == "__main__":
    main()
