from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_EXE = ROOT_DIR / "build" / "bin" / "lbmtv.exe"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec1" / "lbmtv_nx_study"
DOCS_DATA_DIR = ROOT_DIR / "docs" / "sec1" / "generated"
DOCS_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec1"
NX_VALUES = [5, 10, 20, 40, 80]
TAU = 0.8


def build_executable() -> None:
    subprocess.run(["cmd.exe", "/c", "scripts\\build_lbmtv.cmd"], cwd=ROOT_DIR, check=True)


def run_case(nx: int) -> dict[str, float]:
    case_dir = OUTPUT_DIR / f"nx_{nx}"
    case_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(BUILD_EXE), str(TAU), str(case_dir), str(nx)],
        cwd=ROOT_DIR,
        check=True,
        capture_output=True,
        text=True,
    )

    (case_dir / "run.log").write_text(completed.stdout, encoding="utf-8")
    error_text = (case_dir / "error").read_text(encoding="utf-8")
    errors = {}
    for key, value in re.findall(r"([uvp])\s+([0-9.eE+-]+)", error_text):
        errors[key] = float(value)

    return {
        "nx": float(nx),
        "erru": errors["u"],
        "errv": errors["v"],
        "errp": errors["p"],
    }


def write_csv(results: list[dict[str, float]]) -> Path:
    csv_path = OUTPUT_DIR / "lbmtv_nx_errors.csv"
    published_csv_path = DOCS_DATA_DIR / "lbmtv_nx_errors.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["nx", "Erru", "Errv", "Errp"])
        for result in results:
            writer.writerow([result["nx"], result["erru"], result["errv"], result["errp"]])
    published_csv_path.write_text(csv_path.read_text(encoding="utf-8"), encoding="utf-8")
    return csv_path


def write_markdown(results: list[dict[str, float]]) -> Path:
    md_path = OUTPUT_DIR / "lbmtv_nx_errors.md"
    published_md_path = DOCS_DATA_DIR / "lbmtv_nx_errors.md"
    lines = [
        "# lbmtv nx study",
        "",
        f"tau = {TAU}",
        "",
        "| nx | Erru | Errv | Errp |",
        "|---:|---:|---:|---:|",
    ]
    for result in results:
        lines.append(
            f"| {int(result['nx'])} | {result['erru']:.8e} | {result['errv']:.8e} | {result['errp']:.8e} |"
        )
    markdown_text = "\n".join(lines) + "\n"
    md_path.write_text(markdown_text, encoding="utf-8")
    published_md_path.write_text(markdown_text, encoding="utf-8")
    return md_path


def plot_results(results: list[dict[str, float]]) -> Path:
    x_values = [int(result["nx"]) for result in results]
    erru = [result["erru"] for result in results]

    figure, axis = plt.subplots(figsize=(8, 5))
    axis.plot(x_values, erru, marker="o", linewidth=1.8, label="Erru")
    axis.set_xscale("log", base=2)
    axis.set_yscale("log")
    axis.set_xlabel("Grid points nx")
    axis.set_ylabel("Relative error")
    axis.set_title("Relationship between grid resolution and Erru in lbmtv")
    axis.grid(True, which="both", linestyle=":", alpha=0.6)
    axis.legend(loc="best")
    axis.set_xticks(x_values)
    axis.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    figure.tight_layout()
    plot_path = OUTPUT_DIR / "lbmtv_nx_errors.png"
    published_plot_path = DOCS_ASSET_DIR / "lbmtv_nx_errors.png"
    figure.savefig(plot_path, dpi=200)
    figure.savefig(published_plot_path, dpi=200)
    return plot_path


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    build_executable()

    results = [run_case(nx) for nx in NX_VALUES]
    csv_path = write_csv(results)
    md_path = write_markdown(results)
    plot_path = plot_results(results)

    print(f"Saved table to {md_path}")
    print(f"Saved csv to {csv_path}")
    print(f"Saved plot to {plot_path}")


if __name__ == "__main__":
    main()