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
TAU_VALUES = [0.56, 3.0]
FLAG = 2


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


def run_case(executable: Path, tau: float) -> Path:
    output_dir = OUTPUT_DIR / f"lbmbound_ongrid_tau_{str(tau).replace('.', '_')}"
    output_dir.mkdir(parents=True, exist_ok=True)
    run_result = subprocess.run(
        [str(executable), str(tau)],
        cwd=output_dir,
        input=f"{FLAG}\n",
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


def main() -> None:
    executable = build_executable()
    y_over_h = np.linspace(0.0, 1.0, NY + 1)

    figure, profile_axis = plt.subplots(figsize=(6.3, 4.35))

    marker_map = {
        0.56: "o",
        3.0: "x",
    }

    for tau in TAU_VALUES:
        data_path = run_case(executable, tau)
        normalized_profile = read_profile(data_path)
        profile_axis.plot(
            normalized_profile,
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
    profile_axis.set_xticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5])
    profile_axis.set_ylim(0.0, 1.0)
    profile_axis.set_xlabel(r"$u_x/u_c$")
    profile_axis.set_ylabel(r"$y/h$")
    profile_axis.set_title("On-grid bounce back 境界での定常速度分布比較")
    profile_axis.grid(True, linestyle=":", alpha=0.6)
    profile_axis.legend(loc="best")

    figure.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "lbmbound_ongrid_tau_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmbound_ongrid_tau_compare.png"
    figure.savefig(output_path, dpi=200)
    figure.savefig(published_path, dpi=200)
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()