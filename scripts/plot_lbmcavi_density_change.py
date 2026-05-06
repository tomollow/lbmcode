from __future__ import annotations

import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec2" / "generated"

UT_VALUES = [0.025, 0.050, 0.075, 0.100]
NX = 256


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
    source = ROOT_DIR / "src" / "sec2" / "lbmcavi.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "lbmcavi.exe"


def run_case(executable: Path, ut: float) -> tuple[float, float, float, float]:
    case_dir = OUTPUT_DIR / f"lbmcavi_nx{NX}_uwx_study" / f"uwx_{ut:.3f}".replace(".", "_")
    case_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(executable), f"{ut:.6f}", str(NX)],
        cwd=case_dir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )

    with (case_dir / "error").open("r", encoding="utf-8") as file:
        line = file.readline().strip()
    values = [float(value.strip()) for value in line.split(",")]
    return values[0], values[1], values[2], values[3]


def write_csv(rows: list[tuple[float, float, float, float]]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmcavi_density_change_vs_uwx.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["uwx", "rho_avg", "delta_rms", "compressibility_error"])
        for row in rows:
            writer.writerow([
                f"{row[0]:.6f}",
                f"{row[1]:.8f}",
                f"{row[2]:.8f}",
                f"{row[3]:.8f}",
            ])
    return csv_path


def main() -> None:
    executable = build_executable()
    with ThreadPoolExecutor(max_workers=len(UT_VALUES)) as executor:
        rows = list(executor.map(lambda ut: run_case(executable, ut), UT_VALUES))
    csv_path = write_csv(rows)

    uwx = [row[0] for row in rows]
    delta = [row[2] for row in rows]

    figure, axis = plt.subplots(figsize=(4.95, 3.45), constrained_layout=True)
    axis.plot(uwx, delta, color="black", linewidth=1.6, marker="o", markersize=6, markerfacecolor="white")
    axis.set_xlabel(r"$u_x^{w}$")
    axis.set_ylabel(r"$\Delta = \bar{\rho}^{-1}\sqrt{N^{-1}\sum (\rho-\bar{\rho})^2}$")
    axis.set_title(rf"Average density change vs. lid speed ($Re=100$, $n_x={NX}$)")
    axis.grid(True, linestyle=":", alpha=0.6)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "lbmcavi_density_change_vs_uwx.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcavi_density_change_vs_uwx.png"
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")
    print(f"Saved CSV to {csv_path}")


if __name__ == "__main__":
    main()