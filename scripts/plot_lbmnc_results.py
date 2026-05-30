from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT_DIR / "outputs" / "sec3" / "lbmnc"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec3" / "generated"


# de Vahl Davis (1983) benchmark for natural convection in a square cavity
# Reference: G. de Vahl Davis, Int. J. Numer. Methods Fluids 3, 249 (1983).
DVD_BENCHMARK = {
    1.0e3: {"u_max": 3.649, "y_at_umax": 0.813,
            "v_max": 3.697, "x_at_vmax": 0.178,
            "Nu_avg": 1.118},
    1.0e4: {"u_max": 16.178, "y_at_umax": 0.823,
            "v_max": 19.617, "x_at_vmax": 0.119,
            "Nu_avg": 2.243},
    1.0e5: {"u_max": 34.73, "y_at_umax": 0.855,
            "v_max": 68.59, "x_at_vmax": 0.066,
            "Nu_avg": 4.519},
    1.0e6: {"u_max": 64.63, "y_at_umax": 0.850,
            "v_max": 219.36, "x_at_vmax": 0.0379,
            "Nu_avg": 8.800},
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
plt.rcParams["mathtext.fontset"] = "dejavusans"


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with file_path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if stripped:
                rows.append([float(value) for value in stripped.split()])
    return np.array(rows, dtype=float)


def compute_stream_function(u: np.ndarray, v: np.ndarray,
                            x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ny, nx = u.shape
    psi = np.zeros_like(u)
    dy = float(y[1] - y[0])
    for j in range(1, ny):
        psi[j, :] = psi[j - 1, :] + 0.5 * (u[j, :] + u[j - 1, :]) * dy
    return psi


def compute_nusselt_at_hot_wall(temp: np.ndarray,
                                y: np.ndarray) -> tuple[np.ndarray, float]:
    """Local Nusselt number Nu(y) on the hot wall (x=0) and its average.

    For Boussinesq natural convection in a square cavity with hot wall at x=0,
    cold wall at x=L, and temperature non-dimensionalised so T_h=1, T_c=0,
    the local Nusselt number at the hot wall is

        Nu(y) = -L/(T_h - T_c) * dT/dx(0, y) = -dT/dx_nd(0, y)

    where x_nd = x/L. The wall sits at x_nd = 0, half a grid spacing to the
    left of the first interior column. We use the value at the wall T = 1
    together with the two nearest interior columns and a quadratic fit to
    extrapolate the gradient to the wall. Specifically, for an equally
    spaced stencil with the wall at x = 0 and interior nodes at
    x_1 = h/2, x_2 = 3h/2 (h = 1/N grid units in normalised coordinates),
    the wall gradient is

        dT/dx|_wall = (-8 T_w + 9 T_1 - T_2) / (3 h)

    which is a 3-point one-sided second-order approximation that respects
    the half-cell offset of the bounce-back wall.
    """
    n_cells_x = temp.shape[1]
    h = 1.0 / n_cells_x  # interior cell width in normalised x.
    t_wall = 1.0
    t1 = temp[:, 0]
    t2 = temp[:, 1]
    dTdx_wall = (-8.0 * t_wall + 9.0 * t1 - t2) / (3.0 * h)
    nu_local = -dTdx_wall
    # Trapezoidal integration over y_nd in [0, 1]. The interior y nodes are
    # offset by h/2 from the walls, so we extend the profile to the walls
    # with the nearest interior value (adiabatic walls => dNu/dy != 0 but
    # zero contribution is a reasonable closure for the mean).
    y_full = np.concatenate(([0.0], y, [1.0]))
    nu_full = np.concatenate(([nu_local[0]], nu_local, [nu_local[-1]]))
    nu_avg = float(np.trapezoid(nu_full, x=y_full))
    return nu_local, nu_avg


def write_benchmark_csv(metrics: dict[str, float]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmnc_dvd_ra1e4_comparison.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["quantity", "lbmnc", "de_vahl_davis_1983", "abs_error",
                         "rel_error_percent"])
        for key, label, ref in (
            ("u_max", "max |u|h/chi (vertical centerline)", 16.178),
            ("y_at_umax", "y/H at u_max", 0.823),
            ("v_max", "max |v|h/chi (horizontal centerline)", 19.617),
            ("x_at_vmax", "x/H at v_max", 0.119),
            ("Nu_avg", "average Nusselt at hot wall", 2.243),
        ):
            val = metrics[key]
            err = abs(val - ref)
            rel = err / abs(ref) * 100.0
            writer.writerow([label, f"{val:.6f}", f"{ref:.4f}",
                             f"{err:.6f}", f"{rel:.2f}"])
    return csv_path


def main() -> None:
    data_u = read_matrix(RUN_DIR / "datancu")
    data_v = read_matrix(RUN_DIR / "datancv")
    data_t = read_matrix(RUN_DIR / "datance")

    # The C code writes one row per j with i varying fastest.
    # File row index 0 -> j=1 (bottom interior), last row -> j=ny-1 (top interior).
    # Walls sit at i=0.5 and i=nx-0.5 (half-way bounce-back), so the physical
    # fluid centres of interior cell i lie at x/L = (i - 0.5) / (ny - 1).
    ny, nx = data_t.shape
    # The number of interior cells equals ny (= ny_grid - 1). Cell-centre
    # normalised coordinate: (k + 0.5) / ny  for k = 0..ny-1.
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    xx, yy = np.meshgrid(x, y)

    psi = compute_stream_function(data_u, data_v, x, y)

    # Centerline profiles.
    mid_i = nx // 2
    mid_j = ny // 2
    u_vert = data_u[:, mid_i]
    v_horiz = data_v[mid_j, :]

    # Benchmark quantities follow DVD's convention: the positive peaks of
    # u on the vertical centerline (clockwise circulation puts rightward
    # motion at the top, so the positive peak sits near y/L ~ 0.82) and of
    # v on the horizontal centerline (positive peak near the hot wall).
    j_umax = int(np.argmax(u_vert))
    i_vmax = int(np.argmax(v_horiz))
    i_vmin = int(np.argmin(v_horiz))
    u_max = float(u_vert[j_umax])
    v_max = float(v_horiz[i_vmax])
    v_min = float(v_horiz[i_vmin])
    y_at_umax = float(y[j_umax])
    x_at_vmax = float(x[i_vmax])
    x_at_vmin = float(x[i_vmin])

    nu_local, nu_avg = compute_nusselt_at_hot_wall(data_t, y)

    metrics = {
        "u_max": u_max,
        "y_at_umax": y_at_umax,
        "v_max": v_max,
        "x_at_vmax": x_at_vmax,
        "Nu_avg": nu_avg,
    }

    figure = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.05, 1.0], wspace=0.18, hspace=0.18)
    temp_axis = figure.add_subplot(grid[0, 0])
    stream_axis = figure.add_subplot(grid[0, 1])
    u_axis = figure.add_subplot(grid[1, 0])
    v_axis = figure.add_subplot(grid[1, 1])

    # (a) Temperature contour.
    temp_levels = np.linspace(0.0, 1.0, 11)
    cs_fill = temp_axis.contourf(xx, yy, data_t, levels=temp_levels, cmap="coolwarm")
    temp_axis.contour(xx, yy, data_t, levels=temp_levels, colors="black", linewidths=0.6, alpha=0.8)
    temp_axis.set_aspect("equal")
    temp_axis.set_xlabel(r"$x/L$")
    temp_axis.set_ylabel(r"$y/L$")
    temp_axis.set_title("(a) 温度等値線 $T$")
    cbar = figure.colorbar(cs_fill, ax=temp_axis, shrink=0.85, ticks=np.linspace(0.0, 1.0, 6))
    cbar.set_label(r"$T$")

    # (b) Stream function.
    psi_levels = np.linspace(float(psi.min()), float(psi.max()), 15)
    stream_axis.contour(xx, yy, psi, levels=psi_levels, colors="black", linewidths=0.95)
    stream_axis.set_aspect("equal")
    stream_axis.set_xlabel(r"$x/L$")
    stream_axis.set_ylabel(r"$y/L$")
    stream_axis.set_title(r"(b) 流れ関数 $\psi$ の等値線")

    # (c) Vertical centerline u profile.
    u_axis.plot(u_vert, y, color="black", linewidth=1.5, marker="o",
                markersize=3.8, markerfacecolor="white", label="本コード")
    ref_u = DVD_BENCHMARK[1.0e4]
    u_axis.axhline(ref_u["y_at_umax"], color="0.55", linestyle=":", linewidth=1.0)
    u_axis.axvline(ref_u["u_max"], color="0.55", linestyle=":", linewidth=1.0)
    u_axis.scatter([ref_u["u_max"]], [ref_u["y_at_umax"]],
                   marker="*", s=120, color="crimson", zorder=5,
                   label="de Vahl Davis (1983)")
    u_axis.set_xlabel(r"$u\,h/\chi$ at $x = L/2$")
    u_axis.set_ylabel(r"$y/L$")
    u_axis.set_title("(c) 鉛直中心線の x 方向速度")
    u_axis.set_ylim(0.0, 1.0)
    u_axis.grid(True, linestyle=":", alpha=0.6)
    u_axis.legend(loc="lower right", fontsize=9, framealpha=0.9)

    # (d) Horizontal centerline v profile.
    v_axis.plot(x, v_horiz, color="black", linewidth=1.5, marker="s",
                markersize=3.6, markerfacecolor="white", label="本コード")
    ref_v = DVD_BENCHMARK[1.0e4]
    v_axis.axvline(ref_v["x_at_vmax"], color="0.55", linestyle=":", linewidth=1.0)
    v_axis.axhline(ref_v["v_max"], color="0.55", linestyle=":", linewidth=1.0)
    v_axis.scatter([ref_v["x_at_vmax"]], [ref_v["v_max"]],
                   marker="*", s=120, color="crimson", zorder=5,
                   label="de Vahl Davis (1983)")
    v_axis.set_xlabel(r"$x/L$")
    v_axis.set_ylabel(r"$v\,h/\chi$ at $y = L/2$")
    v_axis.set_title("(d) 水平中心線の y 方向速度")
    v_axis.set_xlim(0.0, 1.0)
    v_axis.grid(True, linestyle=":", alpha=0.6)
    v_axis.legend(loc="lower right", fontsize=9, framealpha=0.9)

    figure.suptitle(
        "Natural convection in a square cavity by double-population LBM "
        r"($Ra = 10^{4}$, $Pr = 0.71$, MRT)",
        fontsize=13,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RUN_DIR / "lbmnc_results.png"
    published_path = PUBLISHED_ASSET_DIR / "lbmnc_results.png"
    csv_path = write_benchmark_csv(metrics)

    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    figure.savefig(published_path, dpi=220, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    print(f"Saved plot to {published_path}")
    print(f"Saved CSV to {csv_path}")

    # --- Local Nusselt-number figure ------------------------------------
    nu_figure, (nu_ax, nu_t_ax) = plt.subplots(1, 2, figsize=(11.0, 4.6),
                                               constrained_layout=True)
    # (a) Nu(y) on the hot wall.
    nu_ax.plot(nu_local, y, color="black", linewidth=1.5, marker="o",
               markersize=3.5, markerfacecolor="white", label="本コード Nu(y)")
    nu_ax.axvline(nu_avg, color="crimson", linestyle="--", linewidth=1.2,
                  label=fr"$\overline{{Nu}} = {nu_avg:.3f}$")
    nu_ax.axvline(DVD_BENCHMARK[1.0e4]["Nu_avg"], color="0.45", linestyle=":",
                  linewidth=1.2,
                  label=fr"DVD $\overline{{Nu}} = {DVD_BENCHMARK[1.0e4]['Nu_avg']:.3f}$")
    nu_ax.set_xlabel(r"$Nu(y) = -\partial T/\partial \tilde{x}\,|_{x=0}$")
    nu_ax.set_ylabel(r"$y/L$")
    nu_ax.set_title("(a) 高温壁の局所 Nusselt 数")
    nu_ax.set_ylim(0.0, 1.0)
    nu_ax.grid(True, linestyle=":", alpha=0.6)
    nu_ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # (b) Wall temperature gradient illustrated via near-wall T profiles.
    # Plot T(x, y) at three y/L slices to visualise the boundary layer.
    target_y = [0.1, 0.5, 0.9]
    colors_y = ["#d6604d", "0.20", "#4393c3"]
    for y_t, c in zip(target_y, colors_y):
        j_idx = int(np.argmin(np.abs(y - y_t)))
        nu_t_ax.plot(x, data_t[j_idx, :], color=c, linewidth=1.5,
                     marker="o", markersize=3.0, markerfacecolor="white",
                     label=fr"$y/L = {y[j_idx]:.3f}$")
    nu_t_ax.set_xlabel(r"$x/L$")
    nu_t_ax.set_ylabel(r"$T$")
    nu_t_ax.set_title("(b) 水平断面の温度分布")
    nu_t_ax.set_xlim(0.0, 1.0)
    nu_t_ax.set_ylim(0.0, 1.0)
    nu_t_ax.grid(True, linestyle=":", alpha=0.6)
    nu_t_ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    nu_figure.suptitle(
        r"高温壁の局所 Nusselt 数と温度境界層 ($Ra = 10^{4}$, $Pr = 0.71$)",
        fontsize=12,
    )

    nu_output_path = RUN_DIR / "lbmnc_nusselt.png"
    nu_published_path = PUBLISHED_ASSET_DIR / "lbmnc_nusselt.png"
    nu_figure.savefig(nu_output_path, dpi=220, bbox_inches="tight")
    nu_figure.savefig(nu_published_path, dpi=220, bbox_inches="tight")
    print(f"Saved Nusselt plot to {nu_output_path}")
    print(f"Saved Nusselt plot to {nu_published_path}")

    # Also export Nu(y) numeric data.
    nu_csv = GENERATED_DOC_DIR / "lbmnc_nusselt_local.csv"
    with nu_csv.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["y_over_L", "Nu_local"])
        for yi, nui in zip(y, nu_local):
            writer.writerow([f"{yi:.6f}", f"{nui:.6f}"])
    print(f"Saved Nu(y) CSV to {nu_csv}")

    print("\n=== 主要指標 ===")
    print(f"max |u*| = {u_max:8.4f}  at y/L = {y_at_umax:.4f}  "
          f"(DVD: 16.178 at 0.823)")
    print(f"max v*   = {v_max:8.4f}  at x/L = {x_at_vmax:.4f}  "
          f"(DVD: 19.617 at 0.119)")
    print(f"min v*   = {v_min:8.4f}  at x/L = {x_at_vmin:.4f}  (cold-wall side)")
    print(f"Nu_avg   = {nu_avg:8.4f}            "
          f"(DVD: 2.243)")
    print(f"min T   = {float(data_t.min()):.4f},  max T = {float(data_t.max()):.4f}")
    print(f"min psi = {float(psi.min()):.4e},  max psi = {float(psi.max()):.4e}")


if __name__ == "__main__":
    main()
