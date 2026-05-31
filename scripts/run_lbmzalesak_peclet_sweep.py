"""Peclet-number sweep for src/sec5/lbmzalesak.c (Zalesak's rotating disk).

For each Peclet number Pe, this script patches only the ``pe = 400.0;`` line of a
temporary copy of lbmzalesak.c, builds and runs it (one full revolution = 2500
steps), reads the final order-parameter field ``dataphi``, and measures the
shape-preservation metrics against the analytically reconstructed initial disk:

  * geometric error E1 = sum|H_final - H_init| / sum H_init  (H = indicator phi>0)
  * filled-area change (cells with phi>0)
  * phi over/undershoot (max, min)

The numerical model (collision/forcing/streaming) is left untouched; only the
Peclet number is varied. Non-finite / blown-up fields are flagged as diverged
and skipped. Output:

    docs/sec5/generated/lbmzalesak_peclet.csv
    docs/assets/sec5/lbmzalesak_peclet.png
    outputs/sec5/lbmzalesak/lbmzalesak_peclet.png
"""

from __future__ import annotations

import csv
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC = ROOT_DIR / "src" / "sec5" / "lbmzalesak.c"
BUILD_SCRIPT = ROOT_DIR / "scripts" / "build_one.cmd"
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmzalesak"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"
CSV_PATH = ROOT_DIR / "docs" / "sec5" / "generated" / "lbmzalesak_peclet.csv"

NX = NY = 50
SLOT_I = (93 * NX // 200, 107 * NX // 200)   # = (23, 26) for nx=50 (integer div)
SLOT_J = (2, NY // 2)
PECLETS = [100.0, 200.0, 400.0, 800.0, 1600.0]

BUILD_TIMEOUT = 180
RUN_TIMEOUT = 300

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def initial_field() -> np.ndarray:
    """Reconstruct the exact initial slotted disk (matches the C integer math)."""
    phi = -np.ones((NY + 1, NX + 1))
    for j in range(NY + 1):
        for i in range(NX + 1):
            if np.hypot(i - NX * 0.5, j - NY * 0.5) <= NX * 0.4:
                phi[j, i] = 1.0
    for i in range(SLOT_I[0], SLOT_I[1] + 1):
        for j in range(SLOT_J[0], SLOT_J[1] + 1):
            phi[j, i] = -1.0
    return phi


def load_grid(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        vals = [float(v) for v in line.split()]
        if vals:
            rows.append(vals)
    return np.array(rows)


def patch_peclet(source_text: str, pe: float) -> str:
    """Replace the ``pe = 400.0;`` line only (line-anchored, single match)."""
    pattern = re.compile(r"(?m)^([ \t]*pe = )400\.0;")
    patched, n = pattern.subn(rf"\g<1>{pe:.1f};", source_text, count=1)
    if n != 1:
        raise RuntimeError(f"expected exactly one pe assignment, replaced {n}")
    return patched


def run_case(pe: float, src_text: str, phi0: np.ndarray, h0: np.ndarray) -> dict:
    """Build and run a temporary copy with the given Peclet number."""
    diverged = {"pe": pe, "diverged": 1.0}
    tmp_c = SRC.with_name(f"lbmzalesak_pe{int(round(pe)):05d}.c")
    exe = ROOT_DIR / "build" / "bin" / f"{tmp_c.stem}.exe"
    try:
        tmp_c.write_text(patch_peclet(src_text, pe), encoding="utf-8")
        try:
            build = subprocess.run(
                ["cmd", "/c", str(BUILD_SCRIPT), str(tmp_c)],
                cwd=ROOT_DIR, capture_output=True, text=True,
                timeout=BUILD_TIMEOUT, shell=False,
            )
        except subprocess.TimeoutExpired:
            print(f"  build timed out for Pe={pe}")
            return diverged
        if build.returncode != 0 or not exe.exists():
            print(f"  build failed for Pe={pe}: {build.stderr[-300:]}")
            return diverged

        with tempfile.TemporaryDirectory(prefix="lbmzalesak_sweep_") as run_cwd:
            try:
                subprocess.run(
                    [str(exe)], cwd=run_cwd, capture_output=True, text=True,
                    timeout=RUN_TIMEOUT, shell=False,
                )
            except subprocess.TimeoutExpired:
                print(f"  run timed out for Pe={pe}")
                return diverged

            dataphi = Path(run_cwd) / "dataphi"
            if not dataphi.exists():
                print(f"  no dataphi for Pe={pe}")
                return diverged
            phi_f = load_grid(dataphi)

        if phi_f.shape != phi0.shape or not np.all(np.isfinite(phi_f)):
            print(f"  diverged (non-finite) for Pe={pe}")
            return diverged
        # Reject blow-ups: a sane advected field stays within a few of +/-1.
        if np.abs(phi_f).max() > 5.0:
            print(f"  diverged (|phi|>5) for Pe={pe}")
            return diverged

        hf = (phi_f > 0.0).astype(float)
        e1 = float(np.abs(hf - h0).sum() / h0.sum())
        area0 = float(h0.sum())
        areaf = float(hf.sum())
        return {
            "pe": pe,
            "e1": e1,
            "area_change_pct": (areaf - area0) / area0 * 100.0,
            "phi_max": float(phi_f.max()),
            "phi_min": float(phi_f.min()),
            "mass_change_pct": (float(phi_f.sum()) - float(phi0.sum()))
            / abs(float(phi0.sum())) * 100.0,
            "diverged": 0.0,
        }
    finally:
        if tmp_c.exists():
            tmp_c.unlink()
        if exe.exists():
            exe.unlink()


def main() -> None:
    src_text = SRC.read_text(encoding="utf-8")
    phi0 = initial_field()
    h0 = (phi0 > 0.0).astype(float)

    rows = []
    for pe in PECLETS:
        print(f"Running Pe = {pe} ...")
        rows.append(run_case(pe, src_text, phi0, h0))

    ok = [r for r in rows if r.get("diverged", 1.0) == 0.0]

    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pe", "geometric_error_E1", "area_change_pct",
                         "mass_change_pct", "phi_max", "phi_min"])
        for r in ok:
            writer.writerow([f"{r['pe']:.1f}", f"{r['e1']:.6f}",
                             f"{r['area_change_pct']:.3f}",
                             f"{r['mass_change_pct']:.3f}",
                             f"{r['phi_max']:.4f}", f"{r['phi_min']:.4f}"])
    print(f"Saved CSV to {CSV_PATH}")
    if len(rows) - len(ok) > 0:
        print(f"  ({len(rows) - len(ok)} case(s) diverged and were skipped)")

    if len(ok) >= 2:
        pe = np.array([r["pe"] for r in ok])
        e1 = np.array([r["e1"] for r in ok])
        darea = np.array([r["area_change_pct"] for r in ok])
        pmax = np.array([r["phi_max"] for r in ok])
        pmin = np.array([r["phi_min"] for r in ok])

        fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.0, 4.6),
                                       constrained_layout=True)

        axl.plot(pe, e1 * 100, "o-", color="#1f4e79", label=r"幾何誤差 $E_1$")
        axl.set_xscale("log")
        axl.set_xlabel(r"Peclet 数 $\mathrm{Pe}$")
        axl.set_ylabel(r"$E_1$ [%]", color="#1f4e79")
        axl.tick_params(axis="y", labelcolor="#1f4e79")
        axl.axvline(400.0, color="0.6", ls="dotted", lw=1.0)
        axl.grid(alpha=0.3, which="both")
        axl2 = axl.twinx()
        axl2.plot(pe, darea, "s--", color="#c0392b", label="面積変化")
        axl2.set_ylabel("充填面積変化 [%]", color="#c0392b")
        axl2.tick_params(axis="y", labelcolor="#c0392b")
        axl.set_title(r"(左) 形状誤差と面積保存")

        axr.plot(pe, pmax, "o-", color="#7a1f12", label=r"$\phi_{\max}$")
        axr.plot(pe, pmin, "o-", color="#1f5d8c", label=r"$\phi_{\min}$")
        axr.axhline(1.0, color="0.6", ls="dotted", lw=1.0)
        axr.axhline(-1.0, color="0.6", ls="dotted", lw=1.0)
        axr.set_xscale("log")
        axr.set_xlabel(r"Peclet 数 $\mathrm{Pe}$")
        axr.set_ylabel(r"$\phi$ の極値")
        axr.axvline(400.0, color="0.6", ls="dotted", lw=1.0)
        axr.grid(alpha=0.3, which="both")
        axr.legend(fontsize=9)
        axr.set_title(r"(右) over/undershoot ($\phi$ の極値)")

        fig.suptitle("図 5.2 Zalesak 円板：Peclet 数による形状保持のトレードオフ",
                     fontsize=13)

        PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
        RUN_DIR.mkdir(parents=True, exist_ok=True)
        for path in (RUN_DIR / "lbmzalesak_peclet.png",
                     PUBLISHED_ASSET_DIR / "lbmzalesak_peclet.png"):
            fig.savefig(path, dpi=220, bbox_inches="tight")
            print(f"Saved figure to {path}")
    else:
        print("Not enough converged cases to plot.")


if __name__ == "__main__":
    main()
