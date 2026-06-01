"""Plot the lbmpy reproduction of lbmnc.c (Ra=1e4, Pr=0.71, MRT) and compare
with de Vahl Davis (1983).

Reads the fields written by scripts/lbmnc_lbmpy.py into
outputs/sec3/lbmnc_lbmpy/{datancu,datancv,datance} (same layout as lbmnc.c),
so the helper functions from plot_lbmnc_results.py apply unchanged.

Produces:
  * docs/assets/sec3/lbmnc_lbmpy_results.png   (4-panel: T, psi, u and v centrelines)
  * docs/assets/sec3/lbmnc_lbmpy_nusselt.png   (local Nu and near-wall T)
  * docs/sec3/generated/lbmnc_lbmpy_dvd_ra1e4_comparison.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_lbmnc_results import (
    DVD_BENCHMARK, read_matrix, compute_stream_function,
    compute_nusselt_at_hot_wall,
)


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec3" / "lbmnc_lbmpy"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec3" / "generated"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"

LABEL = "lbmpy (MRT)"


def main() -> None:
    data_u = read_matrix(RUN_DIR / "datancu")
    data_v = read_matrix(RUN_DIR / "datancv")
    data_t = read_matrix(RUN_DIR / "datance")

    ny, nx = data_t.shape
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    xx, yy = np.meshgrid(x, y)

    psi = compute_stream_function(data_u, data_v, x, y)
    mid_i, mid_j = nx // 2, ny // 2
    u_vert = data_u[:, mid_i]
    v_horiz = data_v[mid_j, :]
    j_umax = int(np.argmax(u_vert))
    i_vmax = int(np.argmax(v_horiz))
    u_max, v_max = float(u_vert[j_umax]), float(v_horiz[i_vmax])
    y_at_umax, x_at_vmax = float(y[j_umax]), float(x[i_vmax])
    nu_local, nu_avg = compute_nusselt_at_hot_wall(data_t, y)
    ref = DVD_BENCHMARK[1.0e4]

    # --- 4-panel results figure --------------------------------------------
    fig = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], wspace=0.18, hspace=0.18)
    ax_t = fig.add_subplot(grid[0, 0])
    ax_s = fig.add_subplot(grid[0, 1])
    ax_u = fig.add_subplot(grid[1, 0])
    ax_v = fig.add_subplot(grid[1, 1])

    levels = np.linspace(0.0, 1.0, 11)
    cf = ax_t.contourf(xx, yy, data_t, levels=levels, cmap="coolwarm")
    ax_t.contour(xx, yy, data_t, levels=levels, colors="black", linewidths=0.6, alpha=0.8)
    ax_t.set_aspect("equal"); ax_t.set_xlabel(r"$x/L$"); ax_t.set_ylabel(r"$y/L$")
    ax_t.set_title("(a) 温度等値線 $T$")
    cb = fig.colorbar(cf, ax=ax_t, shrink=0.85, ticks=np.linspace(0.0, 1.0, 6))
    cb.set_label(r"$T$")

    psi_levels = np.linspace(float(psi.min()), float(psi.max()), 15)
    ax_s.contour(xx, yy, psi, levels=psi_levels, colors="black", linewidths=0.95)
    ax_s.set_aspect("equal"); ax_s.set_xlabel(r"$x/L$"); ax_s.set_ylabel(r"$y/L$")
    ax_s.set_title(r"(b) 流れ関数 $\psi$ の等値線")

    ax_u.plot(u_vert, y, color="black", linewidth=1.5, marker="o", markersize=3.8,
              markerfacecolor="white", label=LABEL)
    ax_u.axhline(ref["y_at_umax"], color="0.55", linestyle=":", linewidth=1.0)
    ax_u.axvline(ref["u_max"], color="0.55", linestyle=":", linewidth=1.0)
    ax_u.scatter([ref["u_max"]], [ref["y_at_umax"]], marker="*", s=120,
                 color="crimson", zorder=5, label="de Vahl Davis (1983)")
    ax_u.set_xlabel(r"$u\,h/\chi$ at $x = L/2$"); ax_u.set_ylabel(r"$y/L$")
    ax_u.set_title("(c) 鉛直中心線の x 方向速度"); ax_u.set_ylim(0.0, 1.0)
    ax_u.grid(True, linestyle=":", alpha=0.6); ax_u.legend(loc="lower right", fontsize=9)

    ax_v.plot(x, v_horiz, color="black", linewidth=1.5, marker="s", markersize=3.6,
              markerfacecolor="white", label=LABEL)
    ax_v.axvline(ref["x_at_vmax"], color="0.55", linestyle=":", linewidth=1.0)
    ax_v.axhline(ref["v_max"], color="0.55", linestyle=":", linewidth=1.0)
    ax_v.scatter([ref["x_at_vmax"]], [ref["v_max"]], marker="*", s=120,
                 color="crimson", zorder=5, label="de Vahl Davis (1983)")
    ax_v.set_xlabel(r"$x/L$"); ax_v.set_ylabel(r"$v\,h/\chi$ at $y = L/2$")
    ax_v.set_title("(d) 水平中心線の y 方向速度"); ax_v.set_xlim(0.0, 1.0)
    ax_v.grid(True, linestyle=":", alpha=0.6); ax_v.legend(loc="lower right", fontsize=9)

    fig.suptitle(
        "lbmpy reproduction of natural convection in a square cavity "
        r"($Ra = 10^{4}$, $Pr = 0.71$, MRT)", fontsize=13)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(RUN_DIR / "lbmnc_lbmpy_results.png", dpi=220, bbox_inches="tight")
    fig.savefig(PUBLISHED_ASSET_DIR / "lbmnc_lbmpy_results.png", dpi=220, bbox_inches="tight")

    # --- Nusselt figure -----------------------------------------------------
    nu_fig, (nu_ax, nt_ax) = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)
    nu_ax.plot(nu_local, y, color="black", linewidth=1.5, marker="o", markersize=3.5,
               markerfacecolor="white", label=f"{LABEL} Nu(y)")
    nu_ax.axvline(nu_avg, color="crimson", linestyle="--", linewidth=1.2,
                  label=fr"$\overline{{Nu}} = {nu_avg:.3f}$")
    nu_ax.axvline(ref["Nu_avg"], color="0.45", linestyle=":", linewidth=1.2,
                  label=fr"DVD $\overline{{Nu}} = {ref['Nu_avg']:.3f}$")
    nu_ax.set_xlabel(r"$Nu(y) = -\partial T/\partial \tilde{x}\,|_{x=0}$")
    nu_ax.set_ylabel(r"$y/L$"); nu_ax.set_title("(a) 高温壁の局所 Nusselt 数")
    nu_ax.set_ylim(0.0, 1.0); nu_ax.grid(True, linestyle=":", alpha=0.6)
    nu_ax.legend(loc="upper right", fontsize=9)

    for y_t, c in zip([0.1, 0.5, 0.9], ["#d6604d", "0.20", "#4393c3"]):
        j_idx = int(np.argmin(np.abs(y - y_t)))
        nt_ax.plot(x, data_t[j_idx, :], color=c, linewidth=1.5, marker="o",
                   markersize=3.0, markerfacecolor="white",
                   label=fr"$y/L = {y[j_idx]:.3f}$")
    nt_ax.set_xlabel(r"$x/L$"); nt_ax.set_ylabel(r"$T$")
    nt_ax.set_title("(b) 水平断面の温度分布")
    nt_ax.set_xlim(0.0, 1.0); nt_ax.set_ylim(0.0, 1.0)
    nt_ax.grid(True, linestyle=":", alpha=0.6); nt_ax.legend(loc="upper right", fontsize=9)

    nu_fig.suptitle(r"lbmpy: 高温壁の局所 Nusselt 数と温度境界層 ($Ra = 10^{4}$, $Pr = 0.71$)",
                    fontsize=12)
    nu_fig.savefig(RUN_DIR / "lbmnc_lbmpy_nusselt.png", dpi=220, bbox_inches="tight")
    nu_fig.savefig(PUBLISHED_ASSET_DIR / "lbmnc_lbmpy_nusselt.png", dpi=220, bbox_inches="tight")

    # --- comparison CSV -----------------------------------------------------
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmnc_lbmpy_dvd_ra1e4_comparison.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "lbmpy", "de_vahl_davis_1983", "abs_error", "rel_error_percent"])
        for label, val, refv in (
            ("max |u|h/chi (vertical centerline)", u_max, ref["u_max"]),
            ("y/H at u_max", y_at_umax, ref["y_at_umax"]),
            ("max |v|h/chi (horizontal centerline)", v_max, ref["v_max"]),
            ("x/H at v_max", x_at_vmax, ref["x_at_vmax"]),
            ("average Nusselt at hot wall", nu_avg, ref["Nu_avg"]),
        ):
            err = abs(val - refv)
            w.writerow([label, f"{val:.6f}", f"{refv:.4f}", f"{err:.6f}",
                        f"{err / abs(refv) * 100:.2f}"])

    print("Saved figures to", PUBLISHED_ASSET_DIR)
    print("Saved CSV to", csv_path)
    print(f"u_max={u_max:.4f} (DVD {ref['u_max']}), v_max={v_max:.4f} (DVD {ref['v_max']}), "
          f"Nu_avg={nu_avg:.4f} (DVD {ref['Nu_avg']}), |psi|={abs(psi.min()):.4f} (DVD {ref.get('psi_mid', 5.071)})")


if __name__ == "__main__":
    main()
