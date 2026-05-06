from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_EXE = ROOT_DIR / "build" / "bin" / "lbmbound.exe"
OUTPUT_DIR = ROOT_DIR / "outputs" / "sec2" / "lbmbound_nx_study"
DOCS_DATA_DIR = ROOT_DIR / "docs" / "sec2" / "generated"
DOCS_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"

TAU = 0.56
NX_VALUES = [5, 10, 20, 40, 80]
CASES = [
    (1, "Equilibrium", "equilibrium"),
    (2, "On-grid BB", "on_grid_bounce_back"),
    (3, "Inamuro", "inamuro_no_slip"),
    (4, "Zou-He", "zou_non_equilibrium"),
    (5, "Half-way BB", "half_way_bounce_back"),
    (6, "Interpolated BB (L)", "interpolated_linear"),
    (7, "Interpolated BB (Q)", "interpolated_quadratic"),
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


def expected_max_time(nx: int) -> int:
    return max(100000, 100 * nx * nx)


def build_executable() -> None:
    subprocess.run(["cmd.exe", "/c", "scripts\\build_one.cmd", "src/sec2/lbmbound.c"], cwd=ROOT_DIR, check=True)


def parse_summary(run_log_text: str, flag: int, nx: int) -> re.Match[str]:
    summary_match = re.search(
        r"FINAL time=(\d+) norm=([0-9.eE+-]+) err=([0-9.eE+-]+) converged=(\d+) max_time=(\d+)",
        run_log_text,
    )
    if summary_match is None:
        raise RuntimeError(f"Could not determine final norm/status for flag={flag}, nx={nx}")
    return summary_match


def run_case(flag: int, output_name: str, nx: int) -> dict[str, float | int | bool]:
    case_dir = OUTPUT_DIR / output_name / f"nx_{nx}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_path = case_dir / "run.log"
    error_path = case_dir / "error"

    if log_path.exists() and error_path.exists():
        existing_log = log_path.read_text(encoding="utf-8", errors="ignore")
        if "FINAL time=" in existing_log:
            summary_match = parse_summary(existing_log, flag, nx)
            if int(summary_match.group(5)) == expected_max_time(nx):
                error_value = float(error_path.read_text(encoding="utf-8").strip())
                print(f"Reuse: flag={flag}, nx={nx}", flush=True)
                return {
                    "error": error_value,
                    "final_norm": float(summary_match.group(2)),
                    "converged": summary_match.group(4) == "1",
                    "time": int(summary_match.group(1)),
                    "max_time": int(summary_match.group(5)),
                }

    print(f"Run: flag={flag}, nx={nx}", flush=True)
    with log_path.open("w", encoding="utf-8", buffering=1) as log_file:
        subprocess.run(
            [str(BUILD_EXE), str(TAU), str(nx), "1"],
            cwd=case_dir,
            input=f"{flag}\n",
            check=True,
            text=True,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
    run_log_text = log_path.read_text(encoding="utf-8")
    if not error_path.exists():
        raise RuntimeError(f"Could not determine error value for flag={flag}, nx={nx}")

    error_value = float(error_path.read_text(encoding="utf-8").strip())
    summary_match = parse_summary(run_log_text, flag, nx)

    return {
        "error": error_value,
        "final_norm": float(summary_match.group(2)),
        "converged": summary_match.group(4) == "1",
        "time": int(summary_match.group(1)),
        "max_time": int(summary_match.group(5)),
    }


def collect_results() -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for flag, label, output_name in CASES:
        metrics_by_nx: dict[int, dict[str, float | int | bool]] = {}
        for nx in NX_VALUES:
            metrics_by_nx[nx] = run_case(flag, output_name, nx)
        results.append({"flag": flag, "label": label, "output_name": output_name, "metrics": metrics_by_nx})
    return results


def format_error_cell(metric: dict[str, float | int | bool]) -> str:
    suffix = "*" if not metric["converged"] else ""
    return f"{metric['error']:.8e}{suffix}"


def format_norm_cell(metric: dict[str, float | int | bool]) -> str:
    suffix = "*" if not metric["converged"] else ""
    return f"{metric['final_norm']:.8e}{suffix}"


def write_csv(results: list[dict[str, object]]) -> Path:
    csv_path = OUTPUT_DIR / "lbmbound_nx_errors.csv"
    published_csv_path = DOCS_DATA_DIR / "lbmbound_nx_errors.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "method",
                *[f"nx={nx}" for nx in NX_VALUES],
                *[f"final_norm_nx={nx}" for nx in NX_VALUES],
                *[f"status_nx={nx}" for nx in NX_VALUES],
            ]
        )
        for result in results:
            metrics = result["metrics"]
            writer.writerow(
                [
                    result["label"],
                    *[f"{metrics[nx]['error']:.8e}" for nx in NX_VALUES],
                    *[f"{metrics[nx]['final_norm']:.8e}" for nx in NX_VALUES],
                    *[("converged" if metrics[nx]["converged"] else "not_converged") for nx in NX_VALUES],
                ]
            )
    published_csv_path.write_text(csv_path.read_text(encoding="utf-8"), encoding="utf-8")
    return csv_path


def write_markdown(results: list[dict[str, object]]) -> Path:
    md_path = OUTPUT_DIR / "lbmbound_nx_errors.md"
    published_md_path = DOCS_DATA_DIR / "lbmbound_nx_errors.md"
    lines = [
        "# lbmbound nx study",
        "",
        f"tau = {TAU}",
        "",
        "## Error",
        "",
        "| error | " + " | ".join([f"nx={nx}" for nx in NX_VALUES]) + " |",
        "|---|" + "---:|" * len(NX_VALUES),
    ]
    for result in results:
        metrics = result["metrics"]
        row = [result["label"], *[format_error_cell(metrics[nx]) for nx in NX_VALUES]]
        lines.append("| " + " | ".join(row) + " |")

    lines.extend(
        [
            "",
            "## Final norm",
            "",
            "| final norm | " + " | ".join([f"nx={nx}" for nx in NX_VALUES]) + " |",
            "|---|" + "---:|" * len(NX_VALUES),
        ]
    )
    for result in results:
        metrics = result["metrics"]
        row = [result["label"], *[format_norm_cell(metrics[nx]) for nx in NX_VALUES]]
        lines.append("| " + " | ".join(row) + " |")

    lines.extend(
        [
            "",
            "- `*` は反復上限 `max_time` までに `norm < 1e-10` を満たさなかった未収束ケースを表します。",
        ]
    )
    markdown_text = "\n".join(lines) + "\n"
    md_path.write_text(markdown_text, encoding="utf-8")
    published_md_path.write_text(markdown_text, encoding="utf-8")
    return md_path


def plot_results(results: list[dict[str, object]]) -> Path:
    figure, axis = plt.subplots(figsize=(8.4, 5.6))
    marker_map = {
        "Equilibrium": "o",
        "On-grid BB": "s",
        "Inamuro": "^",
        "Zou-He": "D",
        "Half-way BB": "x",
        "Interpolated BB (L)": "+",
        "Interpolated BB (Q)": "v",
    }
    color_map = {
        "Equilibrium": "black",
        "On-grid BB": "#444444",
        "Inamuro": "#111111",
        "Zou-He": "#666666",
        "Half-way BB": "#222222",
        "Interpolated BB (L)": "#555555",
        "Interpolated BB (Q)": "#777777",
    }
    unconverged_x: list[int] = []
    unconverged_y: list[float] = []

    for result in results:
        metrics = result["metrics"]
        x_values = NX_VALUES
        y_values = [float(metrics[nx]["error"]) for nx in NX_VALUES]
        label = result["label"]
        axis.plot(
            x_values,
            y_values,
            marker=marker_map[label],
            linewidth=1.5,
            markersize=6,
            markerfacecolor="none" if marker_map[label] not in {"x", "+"} else None,
            markeredgewidth=1.0,
            color=color_map[label],
            label=label,
        )
        for nx in NX_VALUES:
            if not metrics[nx]["converged"]:
                unconverged_x.append(nx)
                unconverged_y.append(float(metrics[nx]["error"]))

    if unconverged_x:
        axis.scatter(
            unconverged_x,
            unconverged_y,
            s=85,
            facecolors="none",
            edgecolors="red",
            linewidths=1.4,
            label="Not converged",
            zorder=5,
        )

    axis.set_xscale("log", base=2)
    axis.set_yscale("log")
    axis.set_xlabel("Grid points nx = ny")
    axis.set_ylabel("Relative error")
    axis.set_title("Grid convergence of seven boundary-condition methods in lbmbound")
    axis.grid(True, which="both", linestyle=":", alpha=0.6)
    axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)
    axis.set_xticks(NX_VALUES)
    axis.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    figure.tight_layout()
    plot_path = OUTPUT_DIR / "lbmbound_nx_errors.png"
    published_plot_path = DOCS_ASSET_DIR / "lbmbound_nx_errors.png"
    figure.savefig(plot_path, dpi=200, bbox_inches="tight")
    figure.savefig(published_plot_path, dpi=200, bbox_inches="tight")
    return plot_path


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    build_executable()

    results = collect_results()
    csv_path = write_csv(results)
    md_path = write_markdown(results)
    plot_path = plot_results(results)

    print(f"Saved table to {md_path}")
    print(f"Saved csv to {csv_path}")
    print(f"Saved plot to {plot_path}")


if __name__ == "__main__":
    main()