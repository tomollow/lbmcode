from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

TARGET_TIMES = [10.0, 20.0, 30.0]
COURANT_CASES = [
    {"dt": 0.5, "label": "C=0.5", "linestyle": "--"},
    {"dt": 1.0, "label": "C=1.0", "linestyle": "-"},
]
SCHEMES = [
    {"flag": 1, "title": "Upwind"},
    {"flag": 2, "title": "FTCS"},
    {"flag": 3, "title": "Lax-Wendroff"},
    {"flag": 4, "title": "Leap-Frog"},
]
TIME_COLORS = {
    10.0: "#1f77b4",
    20.0: "#ff7f0e",
    30.0: "#2ca02c",
}

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
    source = ROOT_DIR / "src" / "sec2" / "fdmadv.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "fdmadv.exe"


def run_case(executable: Path, scheme_flag: int, dt: float) -> Path:
    output_dir = OUTPUT_DIR / f"fdmadv_flag_{scheme_flag}_cfl_{str(dt).replace('.', '_')}"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "run.log").open("w", encoding="utf-8") as log_file:
        subprocess.run(
            [str(executable), str(scheme_flag), str(dt)],
            cwd=output_dir,
            check=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
    return output_dir / "fdmadv_history.dat"


def read_history(file_path: Path) -> dict[float, np.ndarray]:
    history: dict[float, np.ndarray] = {}
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if not stripped:
                continue
            values = np.fromstring(stripped, sep=" ")
            history[float(values[0])] = values[1:]
    return history


def lookup_profile(history: dict[float, np.ndarray], target_time: float) -> np.ndarray:
    for time_value, profile in history.items():
        if np.isclose(time_value, target_time):
            return profile
    raise ValueError(f"time={target_time} was not found in {sorted(history)}")


def profile_x_coordinates(history: dict[float, np.ndarray]) -> np.ndarray:
    first_profile = next(iter(history.values()))
    return np.arange(len(first_profile), dtype=float)


def main() -> None:
    executable = build_executable()
    x: np.ndarray | None = None

    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.4), sharex=True)

    for axis, scheme in zip(axes.flat, SCHEMES):
        y_min = np.inf
        y_max = -np.inf

        for case in COURANT_CASES:
            history_path = run_case(executable, scheme["flag"], case["dt"])
            history = read_history(history_path)
            if x is None:
                x = profile_x_coordinates(history)

            for target_time in TARGET_TIMES:
                profile = lookup_profile(history, target_time)
                y_min = min(y_min, float(np.min(profile)))
                y_max = max(y_max, float(np.max(profile)))
                axis.plot(
                    x,
                    profile,
                    linestyle=case["linestyle"],
                    linewidth=1.8,
                    color=TIME_COLORS[target_time],
                    label=f"t={int(target_time)}, {case['label']}",
                )

        padding = max(0.05, 0.08 * (y_max - y_min if y_max > y_min else 1.0))
        axis.set_ylim(y_min - padding, y_max + padding)
        axis.set_title(scheme["title"])
        axis.grid(True, linestyle=":", alpha=0.6)
        axis.set_xlabel("x")
        axis.set_ylabel("f")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    unique_labels: list[str] = []
    unique_handles = []
    for handle, label in zip(handles, labels):
        if label not in unique_labels:
            unique_labels.append(label)
            unique_handles.append(handle)

    figure.suptitle("fdmadv: Courant 数 0.5 と 1.0 の比較", fontsize=14)
    figure.legend(unique_handles, unique_labels, loc="upper center", ncol=3, frameon=False,
                  bbox_to_anchor=(0.5, 0.98))
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "fdmadv_courant_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "fdmadv_courant_compare.png"
    figure.savefig(output_path, dpi=220)
    figure.savefig(published_path, dpi=220)
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()