from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec2" / "lbmcavi"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec2" / "generated"

GHIA_U_RE100 = [
    (0.9766, 0.84123),
    (0.7344, 0.00332),
    (0.5000, -0.20581),
    (0.2813, -0.15662),
    (0.1016, -0.06434),
]

GHIA_V_RE100 = [
    (0.9688, -0.05906),
    (0.8047, -0.24533),
    (0.5000, 0.05454),
    (0.2266, 0.17507),
    (0.0625, 0.09233),
]


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


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)

    wall_color = "0.84"
    line_color = "0.15"

    axis.add_patch(plt.Rectangle((0.10, 0.12), 0.80, 0.76, facecolor="white", edgecolor="none"))
    axis.add_patch(plt.Rectangle((0.06, 0.08), 0.10, 0.80, facecolor=wall_color, edgecolor="none"))
    axis.add_patch(plt.Rectangle((0.84, 0.08), 0.10, 0.80, facecolor=wall_color, edgecolor="none"))
    axis.add_patch(plt.Rectangle((0.06, 0.00), 0.88, 0.12, facecolor=wall_color, edgecolor="none"))

    axis.plot([0.16, 0.16], [0.12, 0.76], color=line_color, linewidth=1.2)
    axis.plot([0.84, 0.84], [0.12, 0.76], color=line_color, linewidth=1.2)
    axis.plot([0.16, 0.84], [0.12, 0.12], color=line_color, linewidth=1.2)
    axis.plot([0.16, 0.84], [0.76, 0.76], color=line_color, linewidth=1.2)

    axis.annotate(
        "",
        xy=(0.67, 0.83),
        xytext=(0.29, 0.83),
        arrowprops={"arrowstyle": "-|>", "lw": 1.3, "color": line_color},
    )
    axis.text(0.48, 0.86, r"$u^w$", ha="center", va="bottom", fontsize=12)

    theta = np.linspace(np.deg2rad(220.0), np.deg2rad(-35.0), 240)
    x_curve = 0.50 + 0.18 * np.cos(theta)
    y_curve = 0.44 + 0.18 * np.sin(theta)
    axis.plot(x_curve, y_curve, color="0.25", linewidth=1.4)
    axis.annotate(
        "",
        xy=(x_curve[-1], y_curve[-1]),
        xytext=(x_curve[-10], y_curve[-10]),
        arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": "0.25"},
    )
    axis.text(0.49, 0.53, "流れ", ha="center", va="center", fontsize=13, color="0.25")
    axis.text(0.50, 0.04, "壁", ha="center", va="bottom", fontsize=12, color="0.30")

    axis.text(0.02, -0.08, "(a) 模式図", transform=axis.transAxes, ha="left", va="top", fontsize=12)
    axis.axis("off")


def add_panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(0.02, 1.02, label, transform=axis.transAxes, ha="left", va="bottom", fontsize=12)


def build_ghia_rows(data_u: np.ndarray, data_v: np.ndarray) -> list[list[str]]:
    ny, nx = data_u.shape
    y_nodes = np.arange(1, ny + 1, dtype=float) / 50.0
    x_nodes = np.arange(1, nx + 1, dtype=float) / 50.0
    u_center = data_u[:, nx // 2]
    v_center = data_v[ny // 2, :]

    rows: list[list[str]] = []
    for coord, ghia_value in GHIA_U_RE100:
        code_value = float(np.interp(coord, y_nodes, u_center))
        rows.append([
            "u_vertical_centerline",
            f"{coord:.4f}",
            f"{code_value:.8f}",
            f"{ghia_value:.5f}",
            f"{abs(code_value - ghia_value):.8f}",
        ])

    for coord, ghia_value in GHIA_V_RE100:
        code_value = float(np.interp(coord, x_nodes, v_center))
        rows.append([
            "v_horizontal_centerline",
            f"{coord:.4f}",
            f"{code_value:.8f}",
            f"{ghia_value:.5f}",
            f"{abs(code_value - ghia_value):.8f}",
        ])

    return rows


def write_ghia_csv(rows: list[list[str]]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmcavi_ghia_re100_comparison.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["quantity", "coordinate", "lbmcavi", "ghia1982", "abs_error"])
        writer.writerows(rows)
    return csv_path


def main() -> None:
    data_u = read_matrix(RUN_DIR / "datau")
    data_v = read_matrix(RUN_DIR / "datav")
    data_psi = read_matrix(RUN_DIR / "datas")
    ghia_rows = build_ghia_rows(data_u, data_v)

    ny, nx = data_psi.shape
    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(x, y)

    mid_i = nx // 2
    mid_j = ny // 2

    figure = plt.figure(figsize=(12.0, 8.6), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.15, 1.0], wspace=0.18, hspace=0.16)
    schematic_axis = figure.add_subplot(grid[0, 0])
    cavity_axis = figure.add_subplot(grid[0, 1])
    u_axis = figure.add_subplot(grid[1, 0])
    v_axis = figure.add_subplot(grid[1, 1])

    draw_schematic(schematic_axis)

    levels = np.linspace(float(data_psi.min()), float(data_psi.max()), 15)
    cavity_axis.contour(xx, yy, data_psi, levels=levels, colors="black", linewidths=0.95)
    cavity_axis.set_aspect("equal")
    cavity_axis.set_xlabel(r"$x/L$")
    cavity_axis.set_ylabel(r"$y/L$")
    cavity_axis.set_title("流れ関数の等値線")
    cavity_axis.plot([0.0, 1.0], [1.0, 1.0], color="0.20", linewidth=1.6, solid_capstyle="butt")
    cavity_axis.text(0.50, 1.04, r"$u^w = 0.1$", ha="center", va="bottom")
    add_panel_label(cavity_axis, "(b) 流れ関数分布")

    u_axis.plot(data_u[:, mid_i], y, color="black", linewidth=1.5, marker="o", markersize=3.8, markerfacecolor="white")
    u_axis.set_xlim(-0.6, 1.05)
    u_axis.set_ylim(0.0, 1.0)
    u_axis.grid(True, linestyle=":", alpha=0.6)
    u_axis.set_xlabel(r"$u(x=L/2,y)/u_t$")
    u_axis.set_ylabel(r"$y/L$")
    u_axis.set_title("鉛直中心線速度")
    add_panel_label(u_axis, "(c)")

    v_axis.plot(x, data_v[mid_j, :], color="black", linewidth=1.5, marker="s", markersize=3.6, markerfacecolor="white")
    v_axis.set_xlim(0.0, 1.0)
    v_axis.grid(True, linestyle=":", alpha=0.6)
    v_axis.set_xlabel(r"$x/L$")
    v_axis.set_ylabel(r"$v(x,y=L/2)/u_t$")
    v_axis.set_title("水平中心線速度")
    add_panel_label(v_axis, "(d)")

    figure.suptitle("Lid-driven cavity flow by D2Q9 LBM, Re = 100", fontsize=14)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmcavi_streamfunction.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcavi_streamfunction.png"
    csv_path = write_ghia_csv(ghia_rows)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")
    print(f"Saved CSV to {csv_path}")


if __name__ == "__main__":
    main()