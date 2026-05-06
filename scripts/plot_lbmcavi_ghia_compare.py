from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec2" / "lbmcavi"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

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


def main() -> None:
    data_u = read_matrix(RUN_DIR / "datau")
    data_v = read_matrix(RUN_DIR / "datav")

    ny, nx = data_u.shape
    y_nodes = np.arange(1, ny + 1, dtype=float) / 50.0
    x_nodes = np.arange(1, nx + 1, dtype=float) / 50.0
    u_center = data_u[:, nx // 2]
    v_center = data_v[ny // 2, :]

    ghia_u_y = np.array([row[0] for row in GHIA_U_RE100], dtype=float)
    ghia_u = np.array([row[1] for row in GHIA_U_RE100], dtype=float)
    ghia_v_x = np.array([row[0] for row in GHIA_V_RE100], dtype=float)
    ghia_v = np.array([row[1] for row in GHIA_V_RE100], dtype=float)

    figure, (u_axis, v_axis) = plt.subplots(1, 2, figsize=(11.6, 4.8), constrained_layout=True)

    u_axis.plot(u_center, y_nodes, color="black", linewidth=1.7, label="LBMcode")
    u_axis.plot(
        ghia_u,
        ghia_u_y,
        linestyle="None",
        marker="o",
        markersize=6,
        markerfacecolor="white",
        markeredgecolor="black",
        label="Ghia (1982)",
    )
    u_axis.set_xlim(-0.28, 1.02)
    u_axis.set_ylim(0.0, 1.0)
    u_axis.set_xlabel(r"$u(x=L/2,y)/u_t$")
    u_axis.set_ylabel(r"$y/L$")
    u_axis.set_title("鉛直中心線上の x 方向速度")
    u_axis.grid(True, linestyle=":", alpha=0.6)
    u_axis.legend(loc="best")

    v_axis.plot(x_nodes, v_center, color="black", linewidth=1.7, label="LBMcode")
    v_axis.plot(
        ghia_v_x,
        ghia_v,
        linestyle="None",
        marker="s",
        markersize=5.8,
        markerfacecolor="white",
        markeredgecolor="black",
        label="Ghia (1982)",
    )
    v_axis.set_xlim(0.0, 1.0)
    v_axis.set_ylim(-0.30, 0.22)
    v_axis.set_xlabel(r"$x/L$")
    v_axis.set_ylabel(r"$v(x,y=L/2)/u_t$")
    v_axis.set_title("水平中心線上の y 方向速度")
    v_axis.grid(True, linestyle=":", alpha=0.6)
    v_axis.legend(loc="best")

    figure.suptitle("Lid-driven cavity flow: comparison with Ghia (1982), Re = 100", fontsize=13.5)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmcavi_ghia_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcavi_ghia_compare.png"
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()