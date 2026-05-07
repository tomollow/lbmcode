from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2" / "fdlbm"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

NY = 22
GX = 1.0e-5
TAU = 0.06
H = NY - 2

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


def build_executable() -> Path:
    source = ROOT_DIR / "src" / "sec2" / "fdlbm.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "fdlbm.exe"


def run_case(executable: Path) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_DIR / "run.log").open("w", encoding="utf-8") as log_file:
        subprocess.run(
            [str(executable)],
            cwd=OUTPUT_DIR,
            check=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
    return OUTPUT_DIR / "data"


def read_profile(file_path: Path) -> np.ndarray:
    values: list[float] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped))
    return np.array(values, dtype=float)


def analytical_profile(point_count: int) -> tuple[np.ndarray, np.ndarray]:
    y = np.arange(point_count, dtype=float) - 1.0
    nu = TAU / 3.0
    u_max = GX * H * H / (8.0 * nu)
    u = GX * y * (H - y) / (2.0 * nu)
    return y, u, u_max


def main() -> None:
    executable = build_executable()
    data_path = run_case(executable)
    normalized_profile = read_profile(data_path)

    y_nodes = np.arange(len(normalized_profile), dtype=float) - 1.0
    analytical_y, analytical_u, u_max = analytical_profile(len(normalized_profile))
    profile = normalized_profile * u_max

    figure, axis = plt.subplots(figsize=(6.8, 5.2), constrained_layout=True)
    axis.plot(
        analytical_u,
        analytical_y,
        color="black",
        linewidth=1.8,
        label="解析解",
    )
    axis.plot(
        profile,
        y_nodes,
        linestyle="None",
        marker="o",
        markersize=6.5,
        markerfacecolor="white",
        markeredgecolor="black",
        label="FDLBM",
    )

    axis.text(
        0.06,
        0.95,
        "D2Q9, FTCS, gx = 1e-5, tau = 0.06, h = 20",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.75"},
    )
    axis.set_xlim(0.0, max(float(np.max(analytical_u)) * 1.05, float(np.max(profile)) * 1.05))
    axis.set_ylim(-1.0, H + 1.0)
    axis.set_xlabel(r"$u_x$")
    axis.set_ylabel(r"$y$")
    axis.set_title("FDLBM による Poiseuille flow の定常速度分布")
    axis.set_yticks(np.arange(0.0, H + 0.1, 5.0))
    axis.grid(True, linestyle=":", alpha=0.6)
    axis.legend(loc="best")

    top_axis = axis.secondary_xaxis("top", functions=(lambda value: value / u_max, lambda value: value * u_max))
    top_axis.set_xlabel(r"$u_x/u_{\max}$")
    right_axis = axis.secondary_yaxis("right", functions=(lambda value: value / H, lambda value: value * H))
    right_axis.set_ylabel(r"$y/h$")

    output_path = OUTPUT_DIR / "fdlbm_profile.png"
    published_path = PUBLISHED_ASSET_DIR / "fdlbm_profile.png"
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()