from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

NX = 20
NY = 20
GX = 1.0e-5
TAU_VALUES = [0.56, 3.0]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic",
    "Meiryo",
    "MS Gothic",
    "Noto Sans CJK JP",
    "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False


def build_executable() -> Path:
    source = ROOT_DIR / "src" / "sec2" / "lbmpoi.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "lbmpoi.exe"


def run_case(executable: Path, tau: float) -> Path:
    output_dir = OUTPUT_DIR / f"lbmpoi_tau_{str(tau).replace('.', '_')}"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "run.log").open("w", encoding="utf-8") as log_file:
        subprocess.run(
            [str(executable), str(tau)],
            cwd=output_dir,
            check=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
    return output_dir / "data"


def read_profile(file_path: Path) -> np.ndarray:
    values: list[float] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped))
    return np.array(values, dtype=float)


def viscosity(tau: float) -> float:
    return (tau - 0.5) / 3.0


def analytical_umax(tau: float) -> float:
    return GX * NY * NY / (8.0 * viscosity(tau))


def analytical_profile(y_over_h: np.ndarray, tau: float) -> np.ndarray:
    return analytical_umax(tau) * 4.0 * y_over_h * (1.0 - y_over_h)


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 10.0)
    axis.set_ylim(0.0, 6.0)
    axis.add_patch(plt.Rectangle((1.0, 4.7), 6.8, 0.8, facecolor="0.9", edgecolor="none"))
    axis.add_patch(plt.Rectangle((1.0, 0.5), 6.8, 0.8, facecolor="0.9", edgecolor="none"))
    axis.plot([1.0, 7.8], [4.7, 4.7], color="black", linewidth=1.5)
    axis.plot([1.0, 7.8], [1.3, 1.3], color="black", linewidth=1.5)

    y_positions = np.linspace(1.5, 4.5, 9)
    y_center = 3.0
    half_height = 1.7
    for y_pos in y_positions:
        shape = max(0.0, 1.0 - ((y_pos - y_center) / half_height) ** 2)
        length = 1.0 + 3.6 * shape
        axis.arrow(
            1.4,
            y_pos,
            length,
            0.0,
            width=0.018,
            head_width=0.14,
            head_length=0.22,
            length_includes_head=True,
            color="black",
            alpha=0.9,
        )

    axis.arrow(
        5.2,
        3.0,
        1.1,
        0.0,
        width=0.035,
        head_width=0.22,
        head_length=0.28,
        length_includes_head=True,
        color="black",
    )
    axis.text(5.55, 3.35, r"$G_x$", ha="center", va="bottom")

    axis.annotate(
        "",
        xy=(8.2, 1.3),
        xytext=(8.2, 4.7),
        arrowprops={"arrowstyle": "<->", "lw": 1.2},
    )
    axis.text(8.45, 3.0, r"$h$", ha="left", va="center", fontsize=12)
    axis.text(4.4, 5.1, "壁", ha="center", va="center", fontsize=12)
    axis.text(4.4, 0.9, "壁", ha="center", va="center", fontsize=12)
    axis.text(0.7, 3.0, "周期的境界", ha="center", va="center", rotation=90)
    axis.text(8.95, 3.0, "周期的境界", ha="center", va="center", rotation=90)
    axis.axis("off")


def main() -> None:
    executable = build_executable()
    y_over_h = np.linspace(0.0, 1.0, NY + 1)

    figure = plt.figure(figsize=(11.0, 5.8))
    grid = figure.add_gridspec(1, 2, width_ratios=[1.15, 1.35])
    schematic_axis = figure.add_subplot(grid[0, 0])
    profile_axis = figure.add_subplot(grid[0, 1])

    draw_schematic(schematic_axis)

    marker_map = {
        0.56: "o",
        3.0: "x",
    }

    for tau in TAU_VALUES:
        data_path = run_case(executable, tau)
        normalized_by_case = read_profile(data_path)
        profile_axis.plot(
            normalized_by_case,
            y_over_h,
            linestyle="None",
            marker=marker_map[tau],
            markersize=7,
            markerfacecolor="none" if marker_map[tau] == "o" else None,
            markeredgewidth=1.2,
            color="black",
            label=f"tau={tau}",
        )

    analytical_curve = 4.0 * y_over_h * (1.0 - y_over_h)
    profile_axis.plot(
        analytical_curve,
        y_over_h,
        linestyle="-",
        color="black",
        linewidth=1.6,
        label="解析解",
    )

    profile_axis.set_xlim(0.0, 1.5)
    profile_axis.set_ylim(0.0, 1.0)
    profile_axis.set_xlabel(r"$u_x/u_c$")
    profile_axis.set_ylabel(r"$y/h$")
    profile_axis.set_title("Poiseuille flow の定常速度分布比較")
    profile_axis.grid(True, linestyle=":", alpha=0.6)
    profile_axis.legend(loc="best")

    figure.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "lbmpoi_tau_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmpoi_tau_compare.png"
    figure.savefig(output_path, dpi=200)
    figure.savefig(published_path, dpi=200)
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()