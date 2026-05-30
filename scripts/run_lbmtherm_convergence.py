"""Grid-refinement (spatial) convergence test for lbmtherm.c.

For each grid size nx = ny in NX_LIST, we
  1. Patch the source to set nx, ny and scale up the outer-loop count
     so that the iteration count grows with h^2 (the diffusion timescale).
  2. Compile via scripts/build_one.cmd to an nx-tagged executable.
  3. Run it from outputs/sec3/lbmtherm_n<nx>/ so the data file lands there.
  4. Read datae, evaluate the analytical solution on the same interior grid,
     and record relative L2 and Linf errors.

Result: a CSV (docs/sec3/generated/lbmtherm_convergence.csv) and a log-log
plot (docs/assets/sec3/lbmtherm_convergence.png) showing L2 vs 1/nx, the
fitted convergence order, and the theoretical reference line of slope -2
that the interpolated bounce-back of Bouzidi-Firdaouss-Lallemand
(Phys Fluids 13, 3452, 2001) is expected to produce.

The temporary nx-tagged source and binaries are cleaned up in try/finally
so an exception never leaves an orphan file.
"""

from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT_DIR / "src" / "sec3" / "lbmtherm.c"
BUILD_BIN_DIR = ROOT_DIR / "build" / "bin"
BUILD_OBJ_DIR = ROOT_DIR / "build" / "obj"
OUTPUT_BASE = ROOT_DIR / "outputs" / "sec3"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec3" / "generated"


# Fix all other parameters at code defaults so the only difference is grid.
PE_REF = 20.0
TAUG_REF = 0.56
Q_REF = 0.7
NX_LIST = [32, 48, 64, 80, 96]  # ny = nx
NX_REF = 64  # the value used to scale outer-loop count
OUTER_REF = 100
INNER_LOOPS = 500
BUILD_TIMEOUT_S = 600.0
RUN_TIMEOUT_S = 1200.0


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def outer_loops_for(nx: int) -> int:
    """Choose outer loops so the temporal residual is small at every grid.

    The slowest mode decays as exp(-t / tau) with tau = h^2 / (pi^2 * chi)
    in lattice units. To get a residual ~ exp(-8) we need t >= 8 * tau.
    We set the baseline so that nx=64 uses 4x the original code's count
    (200 outer loops => 100k steps), and scale as h^2 for finer grids.
    """
    base = 200
    factor = (nx / NX_REF) ** 2
    return max(base, int(round(base * factor)))


def analytical_solution(nx: int, ny: int, pe: float, taug: float,
                        q: float) -> np.ndarray:
    chi = (taug - 0.5) / 3.0
    h = (ny - 2) + 2.0 * q
    u0 = pe * chi / h
    k = 2.0 * np.pi / nx
    beta = k * np.sqrt(1.0 + 1j * u0 / (chi * k))

    i_idx = np.arange(1, nx)
    j_idx = np.arange(1, ny)
    xx = (i_idx)[None, :]
    y_local = ((j_idx - 1) + q)[:, None]

    phase = np.exp(1j * k * xx)
    num = np.sinh(beta * y_local) + np.sinh(beta * (h - y_local))
    den = np.sinh(beta * h)
    return np.real(phase * num / den)


def patch_source(src_text: str, nx: int, outer: int) -> str:
    """Patch the grid size and outer-loop count.

    The anchor `(?m)^[ \t]*` (tab/space at line start only, not \\n) keeps the
    replacement on a single line and avoids running into an adjacent statement.
    """
    patched = re.sub(
        r"(?m)^[ \t]*int\s+nx\s*=\s*\d+\s*,\s*ny\s*=\s*\d+",
        f"  int    nx = {nx}, ny = {nx}",
        src_text,
        count=1,
    )
    patched = re.sub(
        r"for\s*\(\s*loop1\s*=\s*0;\s*loop1\s*<\s*\d+\s*;\s*loop1\+\+\s*\)\s*\{",
        f"for(loop1 = 0; loop1 < {outer}; loop1++){{",
        patched,
        count=1,
    )
    return patched


def build_variant(nx: int, outer: int) -> Path:
    tag = f"lbmtherm_n{nx}"
    variant_src = SOURCE_PATH.parent / f"{tag}.c"
    try:
        variant_src.write_text(
            patch_source(SOURCE_PATH.read_text(encoding="utf-8"), nx, outer),
            encoding="utf-8",
        )
        cmd = ["cmd", "/c", str(ROOT_DIR / "scripts" / "build_one.cmd"),
               str(variant_src.relative_to(ROOT_DIR))]
        result = subprocess.run(cmd, cwd=ROOT_DIR, capture_output=True,
                                text=True, timeout=BUILD_TIMEOUT_S)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"Build failed for nx={nx}")
        exe = BUILD_BIN_DIR / f"{tag}.exe"
        if not exe.exists():
            raise RuntimeError(f"Binary not produced: {exe}")
        return exe
    finally:
        variant_src.unlink(missing_ok=True)


def run_variant(nx: int, exe: Path) -> Path:
    run_dir = OUTPUT_BASE / f"lbmtherm_n{nx}"
    run_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([str(exe)], cwd=run_dir, capture_output=True,
                            text=True, timeout=RUN_TIMEOUT_S)
    (run_dir / "run.log").write_text(result.stdout + "\n" + result.stderr,
                                     encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"Run failed for nx={nx}")
    return run_dir


def read_matrix(file_path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in file_path.open("r", encoding="utf-8"):
        s = line.strip()
        if s:
            rows.append([float(x) for x in s.split()])
    return np.array(rows, dtype=float)


def measure_error(nx: int, run_dir: Path) -> dict[str, float]:
    e_num = read_matrix(run_dir / "datae")
    e_ana = analytical_solution(nx, nx, PE_REF, TAUG_REF, Q_REF)
    if e_num.shape != e_ana.shape:
        raise RuntimeError(
            f"shape mismatch for nx={nx}: numeric {e_num.shape} vs analytic {e_ana.shape}"
        )
    err = e_num - e_ana
    l2 = float(np.sqrt(np.sum(err**2) / np.sum(e_ana**2)))
    linf = float(np.max(np.abs(err)))
    return {"l2": l2, "linf": linf}


def cleanup_variant(nx: int) -> None:
    tag = f"lbmtherm_n{nx}"
    (BUILD_BIN_DIR / f"{tag}.exe").unlink(missing_ok=True)
    (BUILD_OBJ_DIR / f"{tag}.obj").unlink(missing_ok=True)


def write_csv(rows: list[dict]) -> Path:
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmtherm_convergence.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["nx", "ny", "h_lattice", "outer_loops",
                         "total_steps", "L2_rel", "Linf"])
        for r in rows:
            writer.writerow([
                r["nx"], r["nx"],
                f"{(r['nx'] - 2) + 2.0 * Q_REF:.4f}",
                r["outer"],
                r["outer"] * INNER_LOOPS,
                f"{r['l2']:.6e}", f"{r['linf']:.6e}",
            ])
    return csv_path


def make_plot(rows: list[dict]) -> Path:
    nxs = np.array([r["nx"] for r in rows], dtype=float)
    l2 = np.array([r["l2"] for r in rows], dtype=float)
    linf = np.array([r["linf"] for r in rows], dtype=float)
    h_inv = 1.0 / nxs  # h_phys = L_x/nx; in lattice units dx = 1, so 1/nx is the
                       # natural normalised grid spacing per wavelength

    # Power-law fit slope from log-log regression.
    log_h = np.log(h_inv)
    log_l2 = np.log(l2)
    slope_l2, intercept_l2 = np.polyfit(log_h, log_l2, 1)
    log_linf = np.log(linf)
    slope_linf, _ = np.polyfit(log_h, log_linf, 1)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)

    # L2 panel.
    ax = axes[0]
    ax.loglog(h_inv, l2, "o-", color="crimson", markersize=8,
              linewidth=1.6, label=fr"$L_2$ 相対誤差 (傾き ${slope_l2:.2f}$)")
    ax.loglog(h_inv, linf, "s--", color="0.40", markersize=6,
              linewidth=1.0, alpha=0.85,
              label=fr"$L_\infty$ 絶対誤差 (傾き ${slope_linf:.2f}$)")
    # Reference lines anchored at nx=NX_REF.
    anchor_idx = list(nxs).index(NX_REF) if NX_REF in nxs else len(nxs) // 2
    ref_x = h_inv
    ref_y_2 = l2[anchor_idx] * (ref_x / h_inv[anchor_idx]) ** 2.0
    ref_y_1 = l2[anchor_idx] * (ref_x / h_inv[anchor_idx]) ** 1.0
    ax.loglog(ref_x, ref_y_1, ":", color="0.20", linewidth=1.0,
              label=r"参考: 1 次精度 (傾き $-1$)")
    ax.loglog(ref_x, ref_y_2, "-.", color="0.20", linewidth=1.0,
              label=r"参考: 2 次精度 (傾き $-2$)")
    ax.set_xlabel(r"$1 / n_x$  (= $\Delta x / L_x$)")
    ax.set_ylabel(r"相対誤差")
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(fontsize=9, loc="lower right")
    ax.invert_xaxis()
    ax.set_title("(a) 格子細分化による誤差収束")

    # Order-per-step panel (consecutive ratios).
    ax2 = axes[1]
    ratios_l2 = np.log(l2[1:] / l2[:-1]) / np.log(h_inv[1:] / h_inv[:-1])
    ratios_linf = np.log(linf[1:] / linf[:-1]) / np.log(h_inv[1:] / h_inv[:-1])
    pair_labels = [f"{int(nxs[i])}→{int(nxs[i+1])}" for i in range(len(nxs) - 1)]
    xs_pair = np.arange(len(pair_labels))
    ax2.plot(xs_pair, ratios_l2, "o-", color="crimson",
             markersize=9, linewidth=1.4, label=r"$L_2$ 観測次数")
    ax2.plot(xs_pair, ratios_linf, "s--", color="0.40",
             markersize=7, linewidth=1.0, label=r"$L_\infty$ 観測次数")
    ax2.axhline(2.0, color="0.20", linestyle="-.", linewidth=1.0,
                label="理論 2 次精度")
    ax2.axhline(1.0, color="0.20", linestyle=":", linewidth=1.0,
                label="参考 1 次精度")
    ax2.set_xticks(xs_pair)
    ax2.set_xticklabels(pair_labels)
    ax2.set_xlabel(r"隣接ペア ($n_x \to n_x'$)")
    ax2.set_ylabel("局所収束次数 $p$")
    ax2.set_ylim(0.0, 2.5)
    ax2.grid(True, which="both", linestyle=":", alpha=0.5)
    ax2.legend(fontsize=9, loc="upper right")
    ax2.set_title("(b) 隣接 2 解像度の収束次数")

    fig.suptitle(
        r"lbmtherm.c 補間 bounce-back の格子収束 ($Pe = 20$, $q = 0.7$, $\tau_g = 0.56$)",
        fontsize=12,
    )

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    out_pub = PUBLISHED_ASSET_DIR / "lbmtherm_convergence.png"
    out_local = OUTPUT_BASE / "lbmtherm_convergence.png"
    fig.savefig(out_pub, dpi=220, bbox_inches="tight")
    fig.savefig(out_local, dpi=220, bbox_inches="tight")
    return out_pub


def main() -> None:
    rows: list[dict] = []
    for nx in NX_LIST:
        outer = outer_loops_for(nx)
        steps = outer * INNER_LOOPS
        print(f"\n=== nx = {nx:3d}, outer loops = {outer:4d}  ({steps} steps) ===")
        try:
            exe = build_variant(nx, outer)
            run_dir = run_variant(nx, exe)
            metrics = measure_error(nx, run_dir)
            print(f"  L2  rel = {metrics['l2']:.4e}")
            print(f"  Linf    = {metrics['linf']:.4e}")
            rows.append({"nx": nx, "outer": outer, **metrics})
        finally:
            cleanup_variant(nx)

    csv_path = write_csv(rows)
    plot_path = make_plot(rows)
    print(f"\nWrote CSV  : {csv_path}")
    print(f"Wrote plot : {plot_path}")


if __name__ == "__main__":
    main()
