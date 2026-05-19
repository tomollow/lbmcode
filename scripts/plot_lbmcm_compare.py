"""Compare SRT / MRT / Central Moment collisions for src/sec4/lbmcm.c.

Reads dataCMu / dataCMv / dataCMs from each variant subdirectory under
outputs/sec4/lbmcm/ and produces a side-by-side comparison figure:

- Row 1: stream-function contours for SRT, MRT, CM at Re=100
- Row 2: u(L/2, y) centerline overlay (Re=100), v(x, L/2) centerline overlay (Re=100),
  and stream-function comparison for MRT vs CM at Re=1000 (where SRT diverges).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_BASE = ROOT_DIR / "outputs" / "sec4" / "lbmcm"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec4"


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


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                rows.append([float(value) for value in stripped.split()])
    return np.array(rows, dtype=float)


def load_run(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    run_dir = RUN_BASE / name
    return (
        read_matrix(run_dir / "dataCMu"),
        read_matrix(run_dir / "dataCMv"),
        read_matrix(run_dir / "dataCMs"),
    )


def draw_psi(axis: plt.Axes, psi: np.ndarray, title: str) -> None:
    ny, nx = psi.shape
    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(x, y)
    psi_min = float(psi.min())
    psi_max = float(psi.max())
    main_levels = np.linspace(psi_min, 0.0, 11)[:-1]
    sec_levels = np.linspace(0.0, max(psi_max, 1e-9), 6)
    levels = np.unique(np.concatenate([main_levels, sec_levels]))
    axis.contour(xx, yy, psi, levels=levels, colors="black", linewidths=0.8)
    axis.set_aspect("equal")
    axis.set_xlabel(r"$x/L$")
    axis.set_ylabel(r"$y/L$")
    axis.set_title(rf"{title}" "\n" rf"$\psi_{{\min}}={psi_min:.4f},\ \psi_{{\max}}={psi_max:.4f}$")
    axis.plot([0.0, 1.0], [1.0, 1.0], color="red", linewidth=1.4)


def main() -> None:
    srt_u, srt_v, srt_psi = load_run("srt_re100")
    mrt_u, mrt_v, mrt_psi = load_run("mrt_re100")
    cm_u, cm_v, cm_psi = load_run("cm_re100")
    mrt1k_u, mrt1k_v, mrt1k_psi = load_run("mrt_re1000")
    cm1k_u, cm1k_v, cm1k_psi = load_run("cm_re1000")

    figure, axes = plt.subplots(2, 3, figsize=(14.5, 9.2), constrained_layout=True)

    draw_psi(axes[0, 0], srt_psi, "SRT, Re=100")
    draw_psi(axes[0, 1], mrt_psi, "MRT, Re=100")
    draw_psi(axes[0, 2], cm_psi, "CM, Re=100")

    ny, nx = cm_psi.shape
    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    mid_i = nx // 2
    mid_j = ny // 2

    u_ax = axes[1, 0]
    u_ax.plot(srt_u[:, mid_i], y, color="C0", linewidth=1.6, marker="o", markersize=3.0, markerfacecolor="white", label="SRT")
    u_ax.plot(mrt_u[:, mid_i], y, color="C1", linewidth=1.6, marker="s", markersize=3.0, markerfacecolor="white", label="MRT")
    u_ax.plot(cm_u[:, mid_i], y, color="C3", linewidth=1.6, marker="^", markersize=3.2, markerfacecolor="white", label="CM")
    u_ax.set_xlabel(r"$u(L/2,\,y)/U_{\rm lid}$")
    u_ax.set_ylabel(r"$y/L$")
    u_ax.set_title("中心鉛直線 $u$ プロファイル, Re=100")
    u_ax.grid(True, linestyle=":", alpha=0.6)
    u_ax.legend(loc="best", fontsize=9)

    v_ax = axes[1, 1]
    v_ax.plot(x, srt_v[mid_j, :], color="C0", linewidth=1.6, marker="o", markersize=3.0, markerfacecolor="white", label="SRT")
    v_ax.plot(x, mrt_v[mid_j, :], color="C1", linewidth=1.6, marker="s", markersize=3.0, markerfacecolor="white", label="MRT")
    v_ax.plot(x, cm_v[mid_j, :], color="C3", linewidth=1.6, marker="^", markersize=3.2, markerfacecolor="white", label="CM")
    v_ax.set_xlabel(r"$x/L$")
    v_ax.set_ylabel(r"$v(x,\,L/2)/U_{\rm lid}$")
    v_ax.set_title("中心水平線 $v$ プロファイル, Re=100")
    v_ax.grid(True, linestyle=":", alpha=0.6)
    v_ax.legend(loc="best", fontsize=9)

    high_ax = axes[1, 2]
    high_ax.plot(mrt1k_u[:, mid_i], y, color="C1", linewidth=1.6, marker="s", markersize=3.0, markerfacecolor="white", label="MRT")
    high_ax.plot(cm1k_u[:, mid_i], y, color="C3", linewidth=1.6, marker="^", markersize=3.2, markerfacecolor="white", label="CM")
    high_ax.set_xlabel(r"$u(L/2,\,y)/U_{\rm lid}$")
    high_ax.set_ylabel(r"$y/L$")
    high_ax.set_title("中心鉛直線 $u$ プロファイル, Re=1000\n(SRT は発散)")
    high_ax.grid(True, linestyle=":", alpha=0.6)
    high_ax.legend(loc="best", fontsize=9)

    figure.suptitle(
        "lbmcm.c: SRT / MRT / 中心モーメント衝突の比較 (51×51 lid-driven cavity)",
        fontsize=14,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_BASE / "lbmcm_compare.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmcm_compare.png"
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    figure.savefig(published_path, dpi=200, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")

    def summarize(name: str, u: np.ndarray, v: np.ndarray, psi: np.ndarray) -> None:
        print(
            f"{name:10s}: psi=[{psi.min():.4f}, {psi.max():.4f}],"
            f" u_center_min={u[:, mid_i].min():.4f},"
            f" v_center_min={v[mid_j, :].min():.4f}"
        )

    summarize("SRT  Re100", srt_u, srt_v, srt_psi)
    summarize("MRT  Re100", mrt_u, mrt_v, mrt_psi)
    summarize("CM   Re100", cm_u, cm_v, cm_psi)
    summarize("MRT  Re1k", mrt1k_u, mrt1k_v, mrt1k_psi)
    summarize("CM   Re1k", cm1k_u, cm1k_v, cm1k_psi)


if __name__ == "__main__":
    main()
