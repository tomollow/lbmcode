"""Plot Couette transient profiles from src/sec4/lbmblock.c output.

Reads `dataT` (combined grid, ny+1 points) and `dataTf` (fine grid, nyf+1 points)
for T in {100, 200, 500, 1000, 3000} from outputs/sec4/lbmblock/, and produces:

- Figure 1 (lbmblock_couette.png): center-line u(y) profile vs analytical
  Couette solution u(y) = u_t * y / ny, overlaying all 5 timesteps. Left
  panel shows the combined coarse+fine grid (ny+1 = 33 points), right panel
  shows the fine grid only (nyf+1 = 33 points, covers y in [0, ny/2]).

- Figure 2 (lbmblock_couette_error.png): absolute error |u_num - u_exact|
  vs y for each snapshot, on a log scale, to show convergence to the
  analytical solution as t -> infinity.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec4" / "lbmblock"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec4"

# Match calculation conditions in lbmblock.c
NY = 32
NYF = 32  # fine grid covers physical y in [0, ny/2] with delta=1/2
UT = 0.01

TIMES = [100, 200, 500, 1000, 3000]


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


def read_column(file_path: Path) -> np.ndarray:
    with file_path.open("r", encoding="utf-8") as file:
        return np.array([float(line) for line in file if line.strip()], dtype=float)


def main() -> None:
    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)

    # y axes in physical units (y in [0, ny])
    y_comb = np.arange(NY + 1, dtype=float)  # 0..32
    y_fine = np.arange(NYF + 1, dtype=float) * 0.5  # 0..16, step 0.5
    y_exact = np.linspace(0.0, NY, 401)
    u_exact = UT * y_exact / NY

    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(TIMES)))

    # ---------- Figure 1: profiles ----------
    fig1, axes1 = plt.subplots(1, 2, figsize=(11.5, 5.4), constrained_layout=True)

    ax_c, ax_f = axes1
    for color, t in zip(colors, TIMES):
        u_comb = read_column(RUN_DIR / f"data{t}")
        u_finec = read_column(RUN_DIR / f"data{t}f")
        ax_c.plot(
            u_comb / UT,
            y_comb / NY,
            color=color,
            linewidth=1.5,
            marker="o",
            markersize=3.0,
            markerfacecolor="white",
            label=f"t={t}",
        )
        ax_f.plot(
            u_finec / UT,
            y_fine / NY,
            color=color,
            linewidth=1.5,
            marker="s",
            markersize=3.0,
            markerfacecolor="white",
            label=f"t={t}",
        )

    for ax, title in (
        (ax_c, "全領域 (coarse + fine 結合, $j = 0..n_y$)"),
        (ax_f, "細格子のみ (下半分, $j_f = 0..n_{y,f}$)"),
    ):
        ax.plot(
            u_exact / UT,
            y_exact / NY,
            color="black",
            linewidth=1.4,
            linestyle="--",
            label="解析解 $u_t y/H$",
        )
        ax.axhline(0.5, color="gray", linewidth=0.6, linestyle=":")
        ax.text(0.02, 0.51, "粗 / 細 境界", fontsize=8, color="gray")
        ax.set_xlabel(r"$u / u_t$")
        ax.set_ylabel(r"$y / H$")
        ax.set_title(title)
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend(loc="lower right", fontsize=8.5)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0.0, 1.0)

    fig1.suptitle(
        "lbmblock.c: Multi-block LBM による Couette 過渡応答 (中心鉛直線 $x = n_x/2$)",
        fontsize=13,
    )

    out1 = RUN_DIR / "lbmblock_couette.png"
    pub1 = PUBLISHED_ASSET_DIR / "lbmblock_couette.png"
    fig1.savefig(out1, dpi=200, bbox_inches="tight")
    fig1.savefig(pub1, dpi=200, bbox_inches="tight")
    print(f"Saved {out1}")
    print(f"Saved {pub1}")

    # ---------- Figure 2: absolute error ----------
    fig2, ax2 = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)

    u_exact_comb = UT * y_comb / NY
    for color, t in zip(colors, TIMES):
        u_comb = read_column(RUN_DIR / f"data{t}")
        err = np.abs(u_comb - u_exact_comb)
        err = np.where(err < 1e-16, 1e-16, err)
        ax2.semilogy(
            y_comb / NY,
            err / UT,
            color=color,
            linewidth=1.5,
            marker="o",
            markersize=3.0,
            markerfacecolor="white",
            label=f"t={t}",
        )

    ax2.axvline(0.5, color="gray", linewidth=0.6, linestyle=":")
    ax2.text(0.51, 1.5e-4, "粗 / 細 境界", fontsize=8, color="gray")
    ax2.set_xlabel(r"$y / H$")
    ax2.set_ylabel(r"$|u_{\rm num} - u_{\rm exact}| / u_t$")
    ax2.set_title("解析解からの絶対誤差 (中心鉛直線)")
    ax2.grid(True, which="both", linestyle=":", alpha=0.5)
    ax2.legend(loc="best", fontsize=9)
    ax2.set_xlim(0.0, 1.0)

    out2 = RUN_DIR / "lbmblock_couette_error.png"
    pub2 = PUBLISHED_ASSET_DIR / "lbmblock_couette_error.png"
    fig2.savefig(out2, dpi=200, bbox_inches="tight")
    fig2.savefig(pub2, dpi=200, bbox_inches="tight")
    print(f"Saved {out2}")
    print(f"Saved {pub2}")

    # ---------- Summary print ----------
    print("\nFinal-step (t=3000) check vs analytical:")
    u_final = read_column(RUN_DIR / "data3000")
    err_final = np.abs(u_final - u_exact_comb) / UT
    print(f"  max relative error : {err_final.max():.3e}")
    print(f"  mean relative error: {err_final.mean():.3e}")
    print(f"  err at interface (j=16): {err_final[NY // 2]:.3e}")


if __name__ == "__main__":
    main()
