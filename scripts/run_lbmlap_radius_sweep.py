"""Laplace's-law radius sweep for src/sec5/lbmlap.c.

For each droplet radius R = nx * frac, this script patches only the radius
factor in the initial tanh profile of a temporary copy of lbmlap.c, builds and
runs it, parses the converged pressure jump Delta p from stdout, and compares it
against the Laplace prediction Delta p = sigma / R.

The numerical model (collision/forcing/streaming) is left untouched; only the
initial droplet radius is varied. Output:

    docs/sec5/generated/lbmlap_laplace_radius.csv
    docs/assets/sec5/lbmlap_laplace_radius.png
    outputs/sec5/lbmlap/lbmlap_laplace_radius.png
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
SRC = ROOT_DIR / "src" / "sec5" / "lbmlap.c"
BUILD_SCRIPT = ROOT_DIR / "scripts" / "build_one.cmd"
RUN_DIR = ROOT_DIR / "outputs" / "sec5" / "lbmlap"
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec5"
CSV_PATH = ROOT_DIR / "docs" / "sec5" / "generated" / "lbmlap_laplace_radius.csv"

# nx is fixed at 50 in the source; sigma and the radius factor drive Laplace.
NX = 50
SIGMA = 0.0001
# Radius fraction of nx used in the initial tanh profile (default in source: 0.25).
# Keep R + a few interface widths well inside nx/2 = 25 to avoid self-overlap.
FRACTIONS = [0.16, 0.20, 0.25, 0.30, 0.34]

BUILD_TIMEOUT = 180
RUN_TIMEOUT = 600

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic",
                                   "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"

LAPLACE_RE = re.compile(r"laplace's law\s*:\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)")


def patch_radius(source_text: str, fraction: float) -> str:
    """Replace the radius factor 0.25 in the initial tanh profile only."""
    pattern = re.compile(r"(tmp - \(double\)nx\*)0\.25")
    patched, n = pattern.subn(rf"\g<1>{fraction:.5f}", source_text, count=1)
    if n != 1:
        raise RuntimeError(f"expected exactly one radius factor, replaced {n}")
    return patched


def run_case(fraction: float, src_text: str) -> dict:
    """Build and run a temporary lbmlap copy with the given radius fraction."""
    radius = NX * fraction
    tmp_c = SRC.with_name(f"lbmlap_r{int(round(radius*100)):04d}.c")
    exe = ROOT_DIR / "build" / "bin" / f"{tmp_c.stem}.exe"
    try:
        tmp_c.write_text(patch_radius(src_text, fraction), encoding="utf-8")
        build = subprocess.run(
            ["cmd", "/c", str(BUILD_SCRIPT), str(tmp_c)],
            cwd=ROOT_DIR, capture_output=True, text=True,
            timeout=BUILD_TIMEOUT, shell=False,
        )
        if build.returncode != 0 or not exe.exists():
            print(f"  build failed for frac={fraction}: {build.stderr[-300:]}")
            return {"fraction": fraction, "radius": radius, "diverged": 1.0}

        run = subprocess.run(
            [str(exe)], cwd=ROOT_DIR, capture_output=True, text=True,
            timeout=RUN_TIMEOUT, shell=False,
        )
        matches = LAPLACE_RE.findall(run.stdout)
        if not matches:
            print(f"  no laplace output for frac={fraction}")
            return {"fraction": fraction, "radius": radius, "diverged": 1.0}

        dp_measured = float(matches[-1][1])
        dp_theory = SIGMA / radius
        if not np.isfinite(dp_measured) or dp_measured <= 0:
            return {"fraction": fraction, "radius": radius, "diverged": 1.0}
        return {
            "fraction": fraction,
            "radius": radius,
            "inv_radius": 1.0 / radius,
            "dp_theory": dp_theory,
            "dp_measured": dp_measured,
            "rel_error_pct": abs(dp_measured - dp_theory) / dp_theory * 100.0,
            "diverged": 0.0,
        }
    finally:
        if tmp_c.exists():
            tmp_c.unlink()
        if exe.exists():
            exe.unlink()


def main() -> None:
    src_text = SRC.read_text(encoding="utf-8")
    rows = []
    for frac in FRACTIONS:
        print(f"Running radius fraction {frac} (R = {NX*frac:.2f}) ...")
        rows.append(run_case(frac, src_text))

    ok = [r for r in rows if r.get("diverged", 1.0) == 0.0]

    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["fraction", "radius", "inv_radius",
                         "dp_theory", "dp_measured", "rel_error_pct"])
        for r in rows:
            if r.get("diverged", 1.0) == 0.0:
                writer.writerow([f"{r['fraction']:.5f}", f"{r['radius']:.4f}",
                                 f"{r['inv_radius']:.6e}", f"{r['dp_theory']:.6e}",
                                 f"{r['dp_measured']:.6e}", f"{r['rel_error_pct']:.4f}"])
    print(f"Saved CSV to {CSV_PATH}")

    if len(ok) >= 2:
        inv_r = np.array([r["inv_radius"] for r in ok])
        dp = np.array([r["dp_measured"] for r in ok])
        # Fit Delta p = slope / R (slope = effective surface tension), no intercept.
        slope = float(np.sum(inv_r * dp) / np.sum(inv_r * inv_r))

        fig, ax = plt.subplots(figsize=(6.4, 5.2), constrained_layout=True)
        xr = np.linspace(0.0, inv_r.max() * 1.08, 100)
        ax.plot(xr, SIGMA * xr, "--", color="0.4",
                label=rf"理論 $\Delta p = \sigma/R$ ($\sigma={SIGMA:g}$)")
        ax.plot(xr, slope * xr, "-", color="#c0392b",
                label=rf"回帰 $\sigma_{{\mathrm{{eff}}}}={slope:.4e}$")
        ax.plot(inv_r, dp, "o", color="#1f4e79", ms=8, label="LBM 測定値")
        ax.set_xlabel(r"$1/R$ (lattice$^{-1}$)")
        ax.set_ylabel(r"$\Delta p = (\rho_{\mathrm{in}}-\rho_{\mathrm{out}})/3$")
        ax.set_title("図 5.2 Laplace の法則：圧力差と曲率の線形関係")
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.3)
        ax.legend()

        PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
        RUN_DIR.mkdir(parents=True, exist_ok=True)
        for path in (RUN_DIR / "lbmlap_laplace_radius.png",
                     PUBLISHED_ASSET_DIR / "lbmlap_laplace_radius.png"):
            fig.savefig(path, dpi=220, bbox_inches="tight")
            print(f"Saved figure to {path}")
        print(f"Fitted sigma_eff = {slope:.6e}  (input sigma = {SIGMA:g}, "
              f"error {abs(slope-SIGMA)/SIGMA*100:.2f}%)")
    else:
        print("Not enough converged cases to plot.")


if __name__ == "__main__":
    main()
