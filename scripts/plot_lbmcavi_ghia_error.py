from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT_DIR = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT_DIR / "docs" / "sec2" / "generated" / "lbmcavi_ghia_re100_comparison.csv"
RUN_DIR = ROOT_DIR / "outputs" / "sec2" / "lbmcavi"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec2"


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


def read_rows() -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    u_rows: list[tuple[float, float]] = []
    v_rows: list[tuple[float, float]] = []
    with CSV_PATH.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)
        for row in reader:
            coordinate = float(row["coordinate"])
            abs_error = float(row["abs_error"])
            if row["quantity"] == "u_vertical_centerline":
                u_rows.append((coordinate, abs_error))
            elif row["quantity"] == "v_horizontal_centerline":
                v_rows.append((coordinate, abs_error))
    return u_rows, v_rows


def main() -> None:
    u_rows, v_rows = read_rows()

    u_coord = [row[0] for row in u_rows]
    u_err = [row[1] for row in u_rows]
    v_coord = [row[0] for row in v_rows]
    v_err = [row[1] for row in v_rows]

    figure, (u_axis, v_axis) = plt.subplots(1, 2, figsize=(11.4, 4.4), constrained_layout=True)

    u_axis.plot(u_coord, u_err, color="black", linewidth=1.5, marker="o", markersize=6, markerfacecolor="white")
    u_axis.set_xlabel(r"$y/L$")
    u_axis.set_ylabel(r"$|u_{LBM} - u_{Ghia}|$")
    u_axis.set_title("鉛直中心線上の絶対誤差")
    u_axis.set_xlim(0.0, 1.0)
    u_axis.grid(True, linestyle=":", alpha=0.6)

    v_axis.plot(v_coord, v_err, color="black", linewidth=1.5, marker="s", markersize=5.8, markerfacecolor="white")
    v_axis.set_xlabel(r"$x/L$")
    v_axis.set_ylabel(r"$|v_{LBM} - v_{Ghia}|$")
    v_axis.set_title("水平中心線上の絶対誤差")
    v_axis.set_xlim(0.0, 1.0)
    v_axis.grid(True, linestyle=":", alpha=0.6)

    figure.suptitle("Absolute error against Ghia (1982), Re = 100", fontsize=13.5)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmcavi_ghia_error.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcavi_ghia_error.png"
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")


if __name__ == "__main__":
    main()