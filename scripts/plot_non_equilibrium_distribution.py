from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrow, Rectangle


ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"


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


def add_arrow(axis: plt.Axes, start: tuple[float, float], delta: tuple[float, float], color: str, width: float = 0.035) -> None:
    axis.add_patch(
        FancyArrow(
            start[0],
            start[1],
            delta[0],
            delta[1],
            width=width,
            head_width=0.18,
            head_length=0.22,
            length_includes_head=True,
            color=color,
        )
    )


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 10.0)
    axis.set_ylim(0.0, 6.0)

    wall = Rectangle((1.0, 4.55), 8.0, 0.9, facecolor="#e8e8e8", edgecolor="none")
    axis.add_patch(wall)
    axis.plot([1.0, 9.0], [4.55, 4.55], color="black", linewidth=1.6)
    axis.text(5.0, 5.0, "壁面", ha="center", va="center", fontsize=15)

    fluid_node = (5.0, 3.0)
    wall_point = (5.0, 4.05)
    axis.add_patch(Circle(fluid_node, radius=0.14, facecolor="white", edgecolor="black", linewidth=1.4))
    axis.text(fluid_node[0], fluid_node[1] - 0.36, "流体節点", ha="center", va="top", fontsize=12)

    axis.add_patch(Circle(wall_point, radius=0.07, facecolor="black", edgecolor="black"))
    axis.text(wall_point[0] + 0.35, wall_point[1] + 0.02, "壁位置", ha="left", va="center", fontsize=12)
    axis.plot([fluid_node[0], wall_point[0]], [fluid_node[1], wall_point[1]], color="0.45", linestyle="--", linewidth=1.1)

    add_arrow(axis, (5.0, 3.25), (0.0, 0.95), "black")
    add_arrow(axis, (4.85, 3.15), (-1.0, 0.95), "black")
    add_arrow(axis, (5.15, 3.15), (1.0, 0.95), "black")
    axis.text(5.18, 4.0, r"既知: $f_4$", ha="left", va="center", fontsize=12)
    axis.text(3.55, 4.0, r"既知: $f_7$", ha="right", va="center", fontsize=12)
    axis.text(6.45, 4.0, r"既知: $f_8$", ha="left", va="center", fontsize=12)

    add_arrow(axis, (5.0, 4.15), (0.0, -0.95), "tab:blue")
    add_arrow(axis, (3.95, 4.1), (0.95, -0.9), "tab:blue")
    add_arrow(axis, (6.05, 4.1), (-0.95, -0.9), "tab:blue")
    axis.text(5.22, 3.52, r"未知: $f_2$", ha="left", va="center", fontsize=12, color="tab:blue")
    axis.text(4.15, 3.55, r"未知: $f_5$", ha="right", va="center", fontsize=12, color="tab:blue")
    axis.text(5.85, 3.55, r"未知: $f_6$", ha="left", va="center", fontsize=12, color="tab:blue")

    axis.text(1.3, 3.75, "Non-equilibrium\nbounce-back", ha="left", va="center", fontsize=15, weight="bold")
    axis.text(1.3, 3.0, r"$f_i = f_i^{eq} + f_i^{neq}$", ha="left", va="center", fontsize=16)
    axis.text(1.3, 2.45, r"$f_i^{neq}(x_w) = f_{\bar{i}}^{neq}(x_f)$", ha="left", va="center", fontsize=16, color="tab:blue")
    axis.text(1.3, 1.8, r"$f_i = f_i^{eq}(\rho_w, \mathbf{u}_w)$", ha="left", va="center", fontsize=15)
    axis.text(1.3, 1.35, r"$\qquad +\; f_{\bar{i}} - f_{\bar{i}}^{eq}(\rho_w, \mathbf{u}_w)$", ha="left", va="center", fontsize=15)

    axis.annotate(
        "反対方向の非平衡成分を\nそのまま未知分布へ写す",
        xy=(7.0, 3.9),
        xytext=(7.55, 2.25),
        fontsize=12,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "lw": 1.2, "color": "0.25"},
    )

    axis.text(8.7, 3.0, r"$\bar{i}$: 反対方向", ha="right", va="center", fontsize=11, color="0.35")
    axis.axis("off")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(10.8, 5.8))
    draw_schematic(axis)
    figure.tight_layout()

    output_path = OUTPUT_DIR / "non_equilibrium_distribution.png"
    published_path = PUBLISHED_ASSET_DIR / "non_equilibrium_distribution.png"
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()