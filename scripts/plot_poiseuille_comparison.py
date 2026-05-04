from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs"
SECTION_OUTPUT_DIR = OUTPUT_DIR / "sec2"

CASES = [
    {
        "name": "LBM",
        "source": ROOT_DIR / "src" / "sec2" / "lbmpoi.c",
        "exe": BUILD_DIR / "lbmpoi.exe",
        "output_dir": SECTION_OUTPUT_DIR / "lbmpoi",
        "output_file": "data",
        "marker": "o",
    },
    {
        "name": "FDLBM",
        "source": ROOT_DIR / "src" / "sec2" / "fdlbm.c",
        "exe": BUILD_DIR / "fdlbm.exe",
        "output_dir": SECTION_OUTPUT_DIR / "fdlbm",
        "output_file": "data",
        "marker": "x",
    },
]


def build_case(source: Path) -> None:
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )


def run_case(exe: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "run.log").open("w", encoding="utf-8") as log_file:
        subprocess.run([str(exe)], cwd=output_dir, check=True, stdout=log_file, stderr=subprocess.STDOUT)


def read_profile(file_path: Path) -> np.ndarray:
    values = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped))
    return np.array(values, dtype=float)


def normalize_profile(profile: np.ndarray) -> np.ndarray:
    return profile / np.max(profile)


def analytical_profile(point_count: int) -> tuple[np.ndarray, np.ndarray]:
    y = np.linspace(0.0, 1.0, point_count)
    u = 4.0 * y * (1.0 - y)
    return y, u


def draw_schematic(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 10.0)
    axis.set_ylim(0.0, 4.0)
    axis.plot([1.0, 9.0], [3.0, 3.0], color="black", linewidth=2.0)
    axis.plot([1.0, 9.0], [1.0, 1.0], color="black", linewidth=2.0)
    axis.arrow(2.0, 2.0, 5.0, 0.0, width=0.06, head_width=0.25, head_length=0.35,
               length_includes_head=True, color="tab:blue")
    axis.arrow(2.0, 2.5, 4.0, 0.0, width=0.03, head_width=0.18, head_length=0.25,
               length_includes_head=True, color="tab:gray", alpha=0.7)
    axis.arrow(2.0, 1.5, 4.0, 0.0, width=0.03, head_width=0.18, head_length=0.25,
               length_includes_head=True, color="tab:gray", alpha=0.7)
    axis.text(5.0, 3.2, "upper plate", ha="center", va="bottom")
    axis.text(5.0, 0.8, "lower plate", ha="center", va="top")
    axis.text(7.4, 2.2, "flow direction", color="tab:blue")
    axis.text(1.2, 3.4, "2D parallel-plate channel", fontsize=12, weight="bold")
    axis.annotate("channel width H", xy=(8.8, 1.0), xytext=(8.8, 3.0),
                  arrowprops={"arrowstyle": "<->", "lw": 1.2},
                  ha="left", va="center")
    axis.axis("off")


def main() -> None:
    profiles = []
    analytical_y = None
    analytical_u = None

    for case in CASES:
      build_case(case["source"])
      run_case(case["exe"], case["output_dir"])
      profile = read_profile(case["output_dir"] / case["output_file"])
      y = np.linspace(0.0, 1.0, len(profile))
      profiles.append((case, y, normalize_profile(profile)))
      if analytical_y is None:
          analytical_y, analytical_u = analytical_profile(len(profile))

    figure = plt.figure(figsize=(8, 10))
    grid = figure.add_gridspec(2, 1, height_ratios=[1, 2])
    schematic_axis = figure.add_subplot(grid[0, 0])
    profile_axis = figure.add_subplot(grid[1, 0])

    draw_schematic(schematic_axis)

    for case, y, profile in profiles:
        profile_axis.plot(
            profile,
            y,
            linestyle="None",
            marker=case["marker"],
            markersize=7,
            markerfacecolor="none" if case["marker"] == "o" else None,
            label=case["name"],
        )

    profile_axis.plot(analytical_u, analytical_y, linestyle="-", color="black", linewidth=1.8,
                      label="Analytical")
    profile_axis.set_xlabel("Normalized velocity u / u_max")
    profile_axis.set_ylabel("Channel width direction y / H")
    profile_axis.set_title("Steady velocity profile in a 2D parallel-plate channel")
    profile_axis.grid(True, linestyle=":", alpha=0.6)
    profile_axis.legend(loc="best")
    profile_axis.set_xlim(left=0.0)
    profile_axis.set_ylim(0.0, 1.0)

    figure.tight_layout()
    SECTION_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = SECTION_OUTPUT_DIR / "poiseuille_profile_comparison.png"
    figure.savefig(output_path, dpi=200)
    print(f"Saved plot to {output_path}")


if __name__ == "__main__":
    main()