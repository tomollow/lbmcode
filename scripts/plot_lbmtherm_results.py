"""Visualise lbmtherm.c results against the analytical advection-diffusion mode.

The C code solves the steady 2D advection-diffusion equation
    u0 dT/dx = chi (d2T/dx2 + d2T/dy2)
in a channel with Dirichlet (T = cos(k x)) or Neumann (q_y = const) walls
located at y = 0 and y = h with sub-grid offset q. The analytical solution
is a complex-exponential mode of the form

    T_a(x, y) = Re[ exp(i k x) * (sinh(beta y) + sinh(beta (h - y))) /
                 sinh(beta h) ]

with beta = k * sqrt(1 + i u0 / (chi k)).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec3" / "lbmtherm"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"

# Defaults that match src/sec3/lbmtherm.c at the time of writing.
NX = 64
NY = 64
PE = 20.0
TAUG = 0.56
Q = 0.7


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def analytical_solution(nx: int, ny: int, pe: float, taug: float,
                        q: float) -> np.ndarray:
    chi = (taug - 0.5) / 3.0
    h = (ny - 2) + 2.0 * q
    u0 = pe * chi / h
    k = 2.0 * np.pi / nx
    beta = k * np.sqrt(1.0 + 1j * u0 / (chi * k))

    # Interior nodes: i = 1..nx-1 along x, j = 1..ny-1 along y.
    # y-distance from bottom wall: y_local = (j - 1) + q.
    i_idx = np.arange(1, nx)              # 1..nx-1   length nx-1
    j_idx = np.arange(1, ny)              # 1..ny-1   length ny-1
    xx = (i_idx)[None, :]
    y_local = ((j_idx - 1) + q)[:, None]

    phase = np.exp(1j * k * xx)
    num = np.sinh(beta * y_local) + np.sinh(beta * (h - y_local))
    den = np.sinh(beta * h)
    return np.real(phase * num / den)


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in file_path.open("r", encoding="utf-8"):
        s = line.strip()
        if s:
            rows.append([float(x) for x in s.split()])
    return np.array(rows, dtype=float)


def main() -> None:
    e_num = read_matrix(RUN_DIR / "datae")
    e_ana = analytical_solution(NX, NY, PE, TAUG, Q)
    assert e_num.shape == e_ana.shape, (e_num.shape, e_ana.shape)

    err_field = e_num - e_ana
    l2_rel = float(np.sqrt(np.sum(err_field**2) / np.sum(e_ana**2)))
    linf = float(np.max(np.abs(err_field)))
    print(f"Grid:      {e_num.shape}")
    print(f"L2 rel err = {l2_rel:.4e}")
    print(f"L_inf err  = {linf:.4e}")

    ny_int, nx_int = e_num.shape
    x = (np.arange(1, nx_int + 1) + 0.0) / NX  # x/L
    # y / h normalised so bottom wall is 0 and top is 1
    chi = (TAUG - 0.5) / 3.0
    h = (NY - 2) + 2.0 * Q
    y_local = ((np.arange(1, ny_int + 1) - 1) + Q) / h
    xx, yy = np.meshgrid(x, y_local)

    fig = plt.figure(figsize=(12.0, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1])

    levels = np.linspace(-1.0, 1.0, 21)

    ax_num = fig.add_subplot(gs[0, 0])
    cs1 = ax_num.contourf(xx, yy, e_num, levels=levels, cmap="coolwarm")
    ax_num.set_aspect("equal")
    ax_num.set_xlabel(r"$x / L_x$")
    ax_num.set_ylabel(r"$y / h$")
    ax_num.set_title(r"(a) 数値解 $T$ (D2Q5 LBM)")
    fig.colorbar(cs1, ax=ax_num, shrink=0.9)

    ax_ana = fig.add_subplot(gs[0, 1])
    cs2 = ax_ana.contourf(xx, yy, e_ana, levels=levels, cmap="coolwarm")
    ax_ana.set_aspect("equal")
    ax_ana.set_xlabel(r"$x / L_x$")
    ax_ana.set_ylabel(r"$y / h$")
    ax_ana.set_title(r"(b) 解析解 $T_a$")
    fig.colorbar(cs2, ax=ax_ana, shrink=0.9)

    ax_err = fig.add_subplot(gs[0, 2])
    err_max = max(linf, 1e-12)
    err_levels = np.linspace(-err_max, err_max, 21)
    cs3 = ax_err.contourf(xx, yy, err_field, levels=err_levels, cmap="RdBu_r")
    ax_err.set_aspect("equal")
    ax_err.set_xlabel(r"$x / L_x$")
    ax_err.set_ylabel(r"$y / h$")
    ax_err.set_title(
        fr"(c) 誤差 $T - T_a$  ($L_2 = {l2_rel*100:.2f}\%$)"
    )
    fig.colorbar(cs3, ax=ax_err, shrink=0.9)

    # y-profile at x = L/4 (where cos(kx) is small) and at x = 0 (cos peak).
    target_x_fracs = [0.00, 0.25, 0.50]
    ax_prof = fig.add_subplot(gs[1, :2])
    colors = ["#d6604d", "0.20", "#4393c3"]
    for x_frac, c in zip(target_x_fracs, colors):
        i_idx = int(round(x_frac * NX))
        i_idx = max(min(i_idx, nx_int - 1), 0)
        ax_prof.plot(e_num[:, i_idx], y_local, color=c, linewidth=1.6,
                     marker="o", markersize=3.2, markerfacecolor="white",
                     label=fr"数値 $x = {x[i_idx]:.3f}\,L$")
        ax_prof.plot(e_ana[:, i_idx], y_local, color=c, linewidth=0.9,
                     linestyle="--", alpha=0.85,
                     label=fr"解析 $x = {x[i_idx]:.3f}\,L$")
    ax_prof.set_xlabel(r"$T$")
    ax_prof.set_ylabel(r"$y / h$")
    ax_prof.set_xlim(-1.05, 1.05)
    ax_prof.set_ylim(0.0, 1.0)
    ax_prof.grid(True, linestyle=":", alpha=0.6)
    ax_prof.set_title("(d) 鉛直断面の温度プロファイル（数値: 実線, 解析: 破線）")
    ax_prof.legend(fontsize=8, loc="lower right", framealpha=0.9, ncol=2)

    # Pointwise error histogram for quick sanity.
    ax_hist = fig.add_subplot(gs[1, 2])
    ax_hist.hist(err_field.ravel(), bins=40, color="0.30", edgecolor="white")
    ax_hist.set_xlabel(r"$T - T_a$")
    ax_hist.set_ylabel("頻度")
    ax_hist.set_title("(e) 誤差ヒストグラム")
    ax_hist.grid(True, linestyle=":", alpha=0.5)

    chi_val = chi
    u0 = PE * chi_val / h
    fig.suptitle(
        f"lbmtherm.c: D2Q5 thermal LBM の advection-diffusion 検証"
        f"   ($Pe = {PE:.1f}$, $\\tau_g = {TAUG}$, $q = {Q}$, "
        f"$u_0 = {u0:.3e}$, $\\chi = {chi_val:.3e}$)",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    out_local = RUN_DIR / "lbmtherm_results.png"
    out_pub = PUBLISHED_ASSET_DIR / "lbmtherm_results.png"
    fig.savefig(out_local, dpi=220, bbox_inches="tight")
    fig.savefig(out_pub, dpi=220, bbox_inches="tight")
    print(f"Saved {out_local}")
    print(f"Saved {out_pub}")


if __name__ == "__main__":
    main()
