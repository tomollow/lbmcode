from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

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
    return output_dir / "fneq5"


def read_profile(file_path: Path) -> np.ndarray:
    values: list[float] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped))
    return np.array(values, dtype=float)


def normalization_scale(tau: float) -> float:
    return tau * GX * NY / (8.0 * (tau - 0.5))


def theoretical_profile(y_over_h: np.ndarray) -> np.ndarray:
    return 2.0 * y_over_h - 1.0


def main() -> None:
    executable = build_executable()
    y_over_h = np.linspace(0.0, 1.0, NY + 1)
    marker_map = {
        0.56: "o",
        3.0: "x",
    }

    figure, axis = plt.subplots(figsize=(6.3, 5.8))

    for tau in TAU_VALUES:
        file_path = run_case(executable, tau)
        fneq5 = read_profile(file_path)
        normalized = fneq5 / normalization_scale(tau)
        axis.plot(
            normalized,
            y_over_h,
            linestyle="None",
            marker=marker_map[tau],
            markersize=7,
            markerfacecolor="none" if marker_map[tau] == "o" else None,
            markeredgewidth=1.2,
            color="black",
            label=fr"tau={tau}",
        )

    axis.plot(
        theoretical_profile(y_over_h),
        y_over_h,
        linestyle="-",
        color="black",
        linewidth=1.6,
        label="理論値",
    )

    axis.set_xlim(-1.1, 1.1)
    axis.set_ylim(0.0, 1.0)
    axis.set_xlabel(r"$f_5^{neq}/f_{5,	au}^{neq}$")
    axis.set_ylabel(r"$y/h$")
    axis.set_title("非平衡分布関数の比較")
    axis.grid(True, linestyle=":", alpha=0.6)
    axis.legend(loc="best")

    figure.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "lbmpoi_fneq5_tau_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmpoi_fneq5_tau_compare.png"
    figure.savefig(output_path, dpi=220)
    figure.savefig(published_path, dpi=220)
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()