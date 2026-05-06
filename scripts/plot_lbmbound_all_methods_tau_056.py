from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2" / "lbmbound_tau_0_56_summary"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

TAU = 0.56
NY_STANDARD = 20
NY_INTERPOLATED = 21
Q = 0.25
H_INTERPOLATED = NY_INTERPOLATED - 2.0 + 2.0 * Q

CASES = [
    (1, "equilibrium", "Equilibrium"),
    (2, "on_grid_bounce_back", "On-grid BB"),
    (3, "inamuro_no_slip", "Inamuro"),
    (4, "zou_non_equilibrium", "Zou-He"),
    (5, "half_way_bounce_back", "Half-way BB"),
    (6, "interpolated_linear", "Interpolated BB (L)"),
    (7, "interpolated_quadratic", "Interpolated BB (Q)"),
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


def build_executable() -> Path:
    source = ROOT_DIR / "src" / "sec2" / "lbmbound.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "lbmbound.exe"


def run_case(executable: Path, flag: int, output_name: str) -> Path:
    output_dir = OUTPUT_DIR / output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    run_result = subprocess.run(
        [str(executable), str(TAU)],
        cwd=output_dir,
        input=f"{flag}\n",
        text=True,
        capture_output=True,
        check=True,
    )
    (output_dir / "run.log").write_text(run_result.stdout + run_result.stderr, encoding="utf-8")
    return output_dir / "data"


def read_profile(file_path: Path) -> np.ndarray:
    values: list[float] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped))
    return np.array(values, dtype=float)


def get_y_coordinates(flag: int) -> np.ndarray:
    if flag <= 4:
        return np.linspace(0.0, 1.0, NY_STANDARD + 1)
    if flag == 5:
        return (np.arange(1, NY_STANDARD + 1, dtype=float) - 0.5) / NY_STANDARD
    return (np.arange(1, NY_INTERPOLATED, dtype=float) - 1.0 + Q) / H_INTERPOLATED


def main() -> None:
    executable = build_executable()

    figure, profile_axis = plt.subplots(figsize=(7.2, 5.0))

    style_map = {
        1: ("o", "black"),
        2: ("s", "#444444"),
        3: ("^", "#111111"),
        4: ("D", "#666666"),
        5: ("x", "#222222"),
        6: ("+", "#555555"),
        7: ("v", "#777777"),
    }

    for flag, output_name, label in CASES:
        data_path = run_case(executable, flag, output_name)
        normalized_profile = read_profile(data_path)
        y_over_h = get_y_coordinates(flag)
        marker, color = style_map[flag]
        profile_axis.plot(
            normalized_profile,
            y_over_h,
            linestyle="None",
            marker=marker,
            markersize=5.8,
            markerfacecolor="none" if marker not in {"x", "+"} else None,
            markeredgewidth=1.0,
            color=color,
            label=label,
        )

    y_reference = np.linspace(0.0, 1.0, 200)
    analytical_curve = 4.0 * y_reference * (1.0 - y_reference)
    profile_axis.plot(
        analytical_curve,
        y_reference,
        linestyle="-",
        color="black",
        linewidth=1.4,
        label="解析解",
    )

    profile_axis.set_xlim(0.0, 1.5)
    profile_axis.set_xticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5])
    profile_axis.set_ylim(0.0, 1.0)
    profile_axis.set_xlabel(r"$u_x/u_c$")
    profile_axis.set_ylabel(r"$y/h$")
    profile_axis.set_title("7 手法の定常速度分布比較 ($\\tau=0.56$)")
    profile_axis.grid(True, linestyle=":", alpha=0.5)
    profile_axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)

    figure.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "lbmbound_all_methods_tau_056.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmbound_all_methods_tau_056.png"
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    figure.savefig(published_path, dpi=200, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()