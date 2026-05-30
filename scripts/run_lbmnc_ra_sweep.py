"""Run lbmnc.c at multiple Rayleigh numbers and aggregate results.

For each Ra in {1e3, 1e4, 1e5, 1e6} this script:
  1. Patches the Ra line in src/sec3/lbmnc.c into a temporary source file.
  2. Compiles it via scripts/build_one.cmd to a Ra-tagged executable.
  3. Runs it from outputs/sec3/lbmnc_ra<exp>/ so the data files land there.
  4. Reads back datancu, datancv, datance and computes benchmark metrics.

The script also produces a comparison plot vs de Vahl Davis (1983) and a
CSV summary in docs/sec3/generated/lbmnc_ra_sweep.csv.

For Ra = 1e5, 1e6 the original grid (nx = 46) is marginal. We optionally
boost the inner-loop count via a separate patch so the simulation has more
time to settle. The grid itself is left at 46 to keep wall-clock manageable.
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT_DIR / "src" / "sec3" / "lbmnc.c"
BUILD_BIN_DIR = ROOT_DIR / "build" / "bin"
OUTPUT_BASE = ROOT_DIR / "outputs" / "sec3"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec3" / "generated"


RA_LIST = [1.0e3, 1.0e4, 1.0e5, 1.0e6]
# Outer-loop count override for high-Ra cases (default in code is 50).
# Higher Ra needs more time to settle; this keeps both convergence and runtime
# in a reasonable range.
OUTER_LOOP_OVERRIDE = {1.0e3: 50, 1.0e4: 50, 1.0e5: 100, 1.0e6: 200}


DVD_BENCHMARK = {
    1.0e3: {"u_max": 3.649, "y_at_umax": 0.813,
            "v_max": 3.697, "x_at_vmax": 0.178,
            "Nu_avg": 1.118, "psi_mid": 1.174},
    1.0e4: {"u_max": 16.178, "y_at_umax": 0.823,
            "v_max": 19.617, "x_at_vmax": 0.119,
            "Nu_avg": 2.243, "psi_mid": 5.071},
    1.0e5: {"u_max": 34.73, "y_at_umax": 0.855,
            "v_max": 68.59, "x_at_vmax": 0.066,
            "Nu_avg": 4.519, "psi_mid": 9.111},
    1.0e6: {"u_max": 64.63, "y_at_umax": 0.850,
            "v_max": 219.36, "x_at_vmax": 0.0379,
            "Nu_avg": 8.800, "psi_mid": 16.32},
}


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def patch_source(src_text: str, ra_value: float, outer_loops: int) -> str:
    """Replace the Ra assignment and outer-loop count in the source text."""
    patched = re.sub(
        r"(?m)^\s*ra\s*=\s*[^;]+;",
        f"  ra =   {ra_value:.6e};",
        src_text,
        count=1,
    )
    patched = re.sub(
        r"for\(loop1\s*=\s*0;\s*loop1\s*<\s*\d+;\s*loop1\+\+\)\{",
        f"for(loop1 = 0; loop1 < {outer_loops}; loop1++)" "{",
        patched,
        count=1,
    )
    return patched


def build_variant(ra_value: float, outer_loops: int) -> Path:
    """Compile a Ra-tagged binary and return its path."""
    exp = int(round(np.log10(ra_value)))
    tag = f"lbmnc_ra{exp}"
    variant_src = SOURCE_PATH.parent / f"{tag}.c"
    try:
        variant_src.write_text(
            patch_source(SOURCE_PATH.read_text(encoding="utf-8"),
                         ra_value, outer_loops),
            encoding="utf-8",
        )
        cmd = ["cmd", "/c", str(ROOT_DIR / "scripts" / "build_one.cmd"),
               str(variant_src.relative_to(ROOT_DIR))]
        result = subprocess.run(cmd, cwd=ROOT_DIR, capture_output=True,
                                text=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"Build failed for Ra={ra_value}")
        exe = BUILD_BIN_DIR / f"{tag}.exe"
        if not exe.exists():
            raise RuntimeError(f"Binary not produced: {exe}")
        return exe
    finally:
        variant_src.unlink(missing_ok=True)


def run_variant(ra_value: float, exe: Path, timeout_s: float = 900.0) -> Path:
    exp = int(round(np.log10(ra_value)))
    run_dir = OUTPUT_BASE / f"lbmnc_ra{exp}"
    run_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([str(exe)], cwd=run_dir, capture_output=True,
                            text=True, timeout=timeout_s)
    log_path = run_dir / "run.log"
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"Run failed for Ra={ra_value}")
    return run_dir


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in file_path.open("r", encoding="utf-8"):
        s = line.strip()
        if s:
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                return np.full((1, 1), np.nan)
    return np.array(rows, dtype=float)


def stream_function(u: np.ndarray, dy: float) -> np.ndarray:
    psi = np.zeros_like(u)
    for j in range(1, u.shape[0]):
        psi[j, :] = psi[j - 1, :] + 0.5 * (u[j, :] + u[j - 1, :]) * dy
    return psi


def nu_avg(temp: np.ndarray) -> float:
    n_cells_x = temp.shape[1]
    h = 1.0 / n_cells_x
    t_wall = 1.0
    t1 = temp[:, 0]
    t2 = temp[:, 1]
    dTdx_wall = (-8.0 * t_wall + 9.0 * t1 - t2) / (3.0 * h)
    nu_local = -dTdx_wall
    y = (np.arange(temp.shape[0]) + 0.5) / temp.shape[0]
    y_full = np.concatenate(([0.0], y, [1.0]))
    nu_full = np.concatenate(([nu_local[0]], nu_local, [nu_local[-1]]))
    return float(np.trapezoid(nu_full, x=y_full))


def analyse(run_dir: Path) -> dict[str, float]:
    u = read_matrix(run_dir / "datancu")
    v = read_matrix(run_dir / "datancv")
    t = read_matrix(run_dir / "datance")
    if np.any(np.isnan(u)) or np.any(np.isnan(v)) or np.any(np.isnan(t)):
        return {"diverged": 1.0}
    ny, nx = u.shape
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    mid_i = nx // 2
    mid_j = ny // 2
    u_vert = u[:, mid_i]
    v_horiz = v[mid_j, :]
    j_u = int(np.argmax(u_vert))
    i_v = int(np.argmax(v_horiz))
    dy = float(y[1] - y[0])
    psi = stream_function(u, dy)
    nu = nu_avg(t)
    return {
        "diverged": 0.0,
        "u_max": float(u_vert[j_u]),
        "y_at_umax": float(y[j_u]),
        "v_max": float(v_horiz[i_v]),
        "x_at_vmax": float(x[i_v]),
        "Nu_avg": nu,
        "psi_min": float(psi.min()),
        "psi_max": float(psi.max()),
        "T_min": float(t.min()),
        "T_max": float(t.max()),
    }


def write_csv(rows: list[dict]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmnc_ra_sweep.csv"
    fields = ["Ra", "u_max", "u_max_DVD", "y_at_umax", "y_at_umax_DVD",
              "v_max", "v_max_DVD", "x_at_vmax", "x_at_vmax_DVD",
              "Nu_avg", "Nu_avg_DVD", "abs_psi_min", "psi_mid_DVD", "diverged"]
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for r in rows:
            ra = r["Ra"]
            ref = DVD_BENCHMARK[ra]
            if r["metrics"].get("diverged", 0.0) == 1.0:
                writer.writerow([f"{ra:.0e}", "NaN", ref["u_max"], "NaN",
                                 ref["y_at_umax"], "NaN", ref["v_max"],
                                 "NaN", ref["x_at_vmax"], "NaN",
                                 ref["Nu_avg"], "NaN", ref["psi_mid"], "1"])
                continue
            m = r["metrics"]
            writer.writerow([
                f"{ra:.0e}", f"{m['u_max']:.4f}", f"{ref['u_max']:.3f}",
                f"{m['y_at_umax']:.4f}", f"{ref['y_at_umax']:.3f}",
                f"{m['v_max']:.4f}", f"{ref['v_max']:.3f}",
                f"{m['x_at_vmax']:.4f}", f"{ref['x_at_vmax']:.3f}",
                f"{m['Nu_avg']:.4f}", f"{ref['Nu_avg']:.3f}",
                f"{abs(m['psi_min']):.4f}", f"{ref['psi_mid']:.3f}", "0",
            ])
    return csv_path


def make_sweep_plot(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 4, figsize=(14.0, 3.6), constrained_layout=True)
    ra_vals = []
    code_u, code_v, code_nu, code_psi = [], [], [], []
    ref_u, ref_v, ref_nu, ref_psi = [], [], [], []
    for r in rows:
        if r["metrics"].get("diverged", 0.0) == 1.0:
            continue
        ra = r["Ra"]
        m = r["metrics"]
        ref = DVD_BENCHMARK[ra]
        ra_vals.append(ra)
        code_u.append(m["u_max"]); ref_u.append(ref["u_max"])
        code_v.append(m["v_max"]); ref_v.append(ref["v_max"])
        code_nu.append(m["Nu_avg"]); ref_nu.append(ref["Nu_avg"])
        code_psi.append(abs(m["psi_min"])); ref_psi.append(ref["psi_mid"])

    panel_data = [
        (axes[0], "$u_{\\max}\\, h/\\chi$", code_u, ref_u),
        (axes[1], "$v_{\\max}\\, h/\\chi$", code_v, ref_v),
        (axes[2], r"$\overline{Nu}$", code_nu, ref_nu),
        (axes[3], r"$|\psi|_{\max}/\chi$", code_psi, ref_psi),
    ]
    for ax, title, code_y, ref_y in panel_data:
        ax.plot(ra_vals, ref_y, color="0.40", marker="s", markersize=8,
                markerfacecolor="white", linewidth=1.0, linestyle="--",
                label="de Vahl Davis (1983)")
        ax.plot(ra_vals, code_y, color="crimson", marker="o", markersize=7,
                linewidth=1.4, label="本コード")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$Ra$")
        ax.set_title(title)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, loc="upper left")

    fig.suptitle(r"Ra スイープ: 本コード vs de Vahl Davis (1983)", fontsize=12)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    out = PUBLISHED_ASSET_DIR / "lbmnc_ra_sweep.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    fig.savefig(OUTPUT_BASE / "lbmnc_ra_sweep.png", dpi=220, bbox_inches="tight")
    return out


def cleanup_variant_exe(ra_value: float) -> None:
    exp = int(round(np.log10(ra_value)))
    exe = BUILD_BIN_DIR / f"lbmnc_ra{exp}.exe"
    obj = ROOT_DIR / "build" / "obj" / f"lbmnc_ra{exp}.obj"
    exe.unlink(missing_ok=True)
    obj.unlink(missing_ok=True)


def main() -> None:
    rows: list[dict] = []
    for ra in RA_LIST:
        outer = OUTER_LOOP_OVERRIDE.get(ra, 50)
        print(f"\n=== Ra = {ra:.0e}, outer loops = {outer} ===")
        try:
            exe = build_variant(ra, outer)
            run_dir = run_variant(ra, exe)
            m = analyse(run_dir)
            if m.get("diverged", 0.0) == 1.0:
                print(f"  diverged (NaN in output)")
            else:
                ref = DVD_BENCHMARK[ra]
                print(f"  u_max  = {m['u_max']:8.3f}  (DVD {ref['u_max']:.3f})")
                print(f"  v_max  = {m['v_max']:8.3f}  (DVD {ref['v_max']:.3f})")
                print(f"  Nu_avg = {m['Nu_avg']:8.3f}  (DVD {ref['Nu_avg']:.3f})")
                print(f"  |psi|  = {abs(m['psi_min']):8.3f}  (DVD {ref['psi_mid']:.3f})")
            rows.append({"Ra": ra, "metrics": m})
        finally:
            cleanup_variant_exe(ra)

    csv_path = write_csv(rows)
    plot_path = make_sweep_plot(rows)
    print(f"\nWrote CSV  : {csv_path}")
    print(f"Wrote plot : {plot_path}")


if __name__ == "__main__":
    main()
