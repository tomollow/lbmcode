from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs"
SECTION_OUTPUT_DIR = OUTPUT_DIR / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

GX = 1.0e-5
LBM_NY = 20
LBM_TAU = 0.56
FDLBM_NY = 22
FDLBM_TAU = 0.06
H = 20.0

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

CASES = [
    {
        "name": "LBM",
        "source": ROOT_DIR / "src" / "sec2" / "lbmpoi.c",
        "exe": BUILD_DIR / "lbmpoi.exe",
        "output_dir": SECTION_OUTPUT_DIR / "lbmpoi",
        "output_file": "data",
        "marker": "o",
        "tau": LBM_TAU,
        "h": float(LBM_NY),
        "y_offset": 0.0,
    },
    {
        "name": "FDLBM",
        "source": ROOT_DIR / "src" / "sec2" / "fdlbm.c",
        "exe": BUILD_DIR / "fdlbm.exe",
        "output_dir": SECTION_OUTPUT_DIR / "fdlbm",
        "output_file": "data",
        "marker": "x",
        "tau": FDLBM_TAU,
        "h": float(FDLBM_NY - 2),
        "y_offset": -1.0,
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


def viscosity(case: dict[str, object]) -> float:
    tau = float(case["tau"])
    if case["name"] == "FDLBM":
        return tau / 3.0
    return (tau - 0.5) / 3.0


def physical_profile(case: dict[str, object], normalized_profile: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    h = float(case["h"])
    nu = viscosity(case)
    u_max = GX * h * h / (8.0 * nu)
    y = np.arange(len(normalized_profile), dtype=float) + float(case["y_offset"])
    u = normalized_profile * u_max
    return y, u, u_max


def analytical_profile() -> tuple[np.ndarray, np.ndarray, float]:
    nu = (LBM_TAU - 0.5) / 3.0
    y = np.linspace(0.0, H, 300)
    u = GX * y * (H - y) / (2.0 * nu)
    u_max = GX * H * H / (8.0 * nu)
    return y, u, u_max


def main() -> None:
    profiles = []
    analytical_y, analytical_u, analytical_u_max = analytical_profile()

    for case in CASES:
        build_case(case["source"])
        run_case(case["exe"], case["output_dir"])
        normalized_profile = read_profile(case["output_dir"] / case["output_file"])
        y, profile, _ = physical_profile(case, normalized_profile)
        profiles.append((case, y, profile))

    figure, profile_axis = plt.subplots(figsize=(7.2, 6.8), constrained_layout=True)

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
    profile_axis.set_xlabel(r"$u_x$")
    profile_axis.set_ylabel(r"$y$")
    profile_axis.set_title("Poiseuille flow velocity profile comparison")
    profile_axis.grid(True, linestyle=":", alpha=0.6)
    profile_axis.legend(loc="best")
    profile_axis.set_xlim(left=0.0)
    profile_axis.set_ylim(0.0, H)
    profile_axis.set_yticks(np.arange(0.0, H + 0.1, 5.0))

    top_axis = profile_axis.secondary_xaxis(
        "top",
        functions=(lambda value: value / analytical_u_max, lambda value: value * analytical_u_max),
    )
    top_axis.set_xlabel(r"$u_x/u_{\max}$")
    right_axis = profile_axis.secondary_yaxis(
        "right",
        functions=(lambda value: value / H, lambda value: value * H),
    )
    right_axis.set_ylabel(r"$y/h$")

    SECTION_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = SECTION_OUTPUT_DIR / "poiseuille_profile_comparison.png"
    published_path = PUBLISHED_ASSET_DIR / "poiseuille_profile_comparison.png"
    figure.savefig(output_path, dpi=200)
    figure.savefig(published_path, dpi=200)
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()