from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec2" / "generated"

UT_VALUES = [0.025, 0.050, 0.075, 0.100]
GRID_SIZES = [51, 101]
TEXTBOOK_DELTA = {
    0.025: 1.6395e-6,
    0.050: 5.9074e-6,
    0.075: 1.2870e-5,
    0.100: 2.2422e-5,
}


def build_executable() -> Path:
    source = ROOT_DIR / "src" / "sec2" / "lbmcavi.c"
    subprocess.run(
        ["cmd.exe", "/c", "scripts\\build_one.cmd", str(source.relative_to(ROOT_DIR))],
        cwd=ROOT_DIR,
        check=True,
    )
    return BUILD_DIR / "lbmcavi.exe"


def run_case(executable: Path, ut: float, nx: int) -> tuple[float, float, float, float]:
    case_dir = OUTPUT_DIR / "lbmcavi_delta_grid_study" / f"nx_{nx}" / f"uwx_{ut:.3f}".replace(".", "_")
    case_dir.mkdir(parents=True, exist_ok=True)
    with (case_dir / "run.log").open("w", encoding="utf-8") as log_file:
        subprocess.run(
            [str(executable), f"{ut:.6f}", str(nx)],
            cwd=case_dir,
            check=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )

    with (case_dir / "error").open("r", encoding="utf-8") as file:
        line = file.readline().strip()
    values = [float(value.strip()) for value in line.split(",")]
    return values[0], values[1], values[2], values[3]


def write_csv(rows: list[list[str]]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmcavi_delta_grid_comparison.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "nx",
                "uwx",
                "rho_avg",
                "delta_rms",
                "textbook_table21",
                "abs_error_against_textbook",
            ]
        )
        writer.writerows(rows)
    return csv_path


def main() -> None:
    executable = build_executable()
    rows: list[list[str]] = []
    for nx in GRID_SIZES:
        for ut in UT_VALUES:
            result = run_case(executable, ut, nx)
            textbook_value = TEXTBOOK_DELTA[ut]
            abs_error = abs(result[3] - textbook_value)
            rows.append(
                [
                    str(nx),
                    f"{result[0]:.6f}",
                    f"{result[1]:.8f}",
                    f"{result[2]:.8e}",
                    f"{textbook_value:.8e}",
                    f"{abs_error:.8e}",
                ]
            )
            print(
                f"nx={nx:3d}, uwx={ut:.3f}, delta_rms={result[2]:.8e}, "
                f"textbook={textbook_value:.8e}, abs_error={abs_error:.8e}"
            )

    csv_path = write_csv(rows)
    print(f"Saved CSV to {csv_path}")


if __name__ == "__main__":
    main()