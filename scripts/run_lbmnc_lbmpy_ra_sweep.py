"""Rayleigh-number sweep of the lbmpy reproduction of lbmnc.c.

For each Ra in {1e3, 1e4, 1e5, 1e6} this runs the coupled hydro+thermal lbmpy
scenario (MRT hydrodynamic collision) defined in scripts/lbmnc_lbmpy.py, then
compares the benchmark metrics against de Vahl Davis (1983), writing
  * docs/sec3/generated/lbmnc_lbmpy_ra_sweep.csv
  * docs/assets/sec3/lbmnc_lbmpy_ra_sweep.png

Mirrors scripts/run_lbmnc_ra_sweep.py (which sweeps the C code) so the two
sweeps can be compared side by side. As with lbmnc.c on the 46x46 grid, Ra=1e6
is expected to diverge (lattice Mach number too high).
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import lbmnc_lbmpy as L


ROOT_DIR = Path(__file__).resolve().parents[1]
PUBLISHED_ASSET_DIR = ROOT_DIR / "docs" / "assets" / "sec3"
GENERATED_DOC_DIR = ROOT_DIR / "docs" / "sec3" / "generated"

RA_LIST = [1.0e3, 1.0e4, 1.0e5, 1.0e6]
# Higher Ra needs more iterations to settle (and a tighter Mach margin).
STEPS_OVERRIDE = {1.0e3: 40000, 1.0e4: 60000, 1.0e5: 120000, 1.0e6: 120000}
METHOD = "mrt"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP", "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"


def write_csv(rows):
    GENERATED_DOC_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = GENERATED_DOC_DIR / "lbmnc_lbmpy_ra_sweep.csv"
    fields = ["Ra", "u_max", "u_max_DVD", "y_at_umax", "y_at_umax_DVD",
              "v_max", "v_max_DVD", "x_at_vmax", "x_at_vmax_DVD",
              "Nu_avg", "Nu_avg_DVD", "abs_psi", "psi_mid_DVD", "diverged"]
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for r in rows:
            ra = r["Ra"]; ref = L.DVD[ra]
            if r["diverged"]:
                w.writerow([f"{ra:.0e}", "NaN", ref["u_max"], "NaN", ref["y_um"],
                            "NaN", ref["v_max"], "NaN", ref["x_vm"], "NaN",
                            ref["Nu"], "NaN", ref["psi"], "1"])
                continue
            m = r["m"]
            w.writerow([f"{ra:.0e}", f"{m['u_max']:.4f}", f"{ref['u_max']:.3f}",
                        f"{m['y_um']:.4f}", f"{ref['y_um']:.3f}",
                        f"{m['v_max']:.4f}", f"{ref['v_max']:.3f}",
                        f"{m['x_vm']:.4f}", f"{ref['x_vm']:.3f}",
                        f"{m['Nu']:.4f}", f"{ref['Nu']:.3f}",
                        f"{m['psi']:.4f}", f"{ref['psi']:.3f}", "0"])
    return csv_path


def make_plot(rows):
    ra_vals, cu, cv, cn, cp = [], [], [], [], []
    ru, rv, rn, rp = [], [], [], []
    for r in rows:
        if r["diverged"]:
            continue
        ra = r["Ra"]; m = r["m"]; ref = L.DVD[ra]
        ra_vals.append(ra)
        cu.append(m["u_max"]); ru.append(ref["u_max"])
        cv.append(m["v_max"]); rv.append(ref["v_max"])
        cn.append(m["Nu"]); rn.append(ref["Nu"])
        cp.append(m["psi"]); rp.append(ref["psi"])

    fig, axes = plt.subplots(1, 4, figsize=(14.0, 3.6), constrained_layout=True)
    panels = [
        (axes[0], r"$u_{\max}\, h/\chi$", cu, ru),
        (axes[1], r"$v_{\max}\, h/\chi$", cv, rv),
        (axes[2], r"$\overline{Nu}$", cn, rn),
        (axes[3], r"$|\psi|_{\max}/\chi$", cp, rp),
    ]
    for ax, title, cy, ry in panels:
        ax.plot(ra_vals, ry, color="0.40", marker="s", markersize=8,
                markerfacecolor="white", linewidth=1.0, linestyle="--",
                label="de Vahl Davis (1983)")
        ax.plot(ra_vals, cy, color="crimson", marker="o", markersize=7,
                linewidth=1.4, label="lbmpy (MRT)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$Ra$"); ax.set_title(title)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(r"Ra スイープ: lbmpy (MRT) vs de Vahl Davis (1983)", fontsize=12)

    PUBLISHED_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    out = PUBLISHED_ASSET_DIR / "lbmnc_lbmpy_ra_sweep.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    return out


def main():
    rows = []
    for ra in RA_LIST:
        steps = STEPS_OVERRIDE.get(ra, 60000)
        print(f"\n=== Ra = {ra:.0e}  ({METHOD.upper()}, steps<= {steps}) ===")
        p = L.build_params(N=46, Ra=ra, Pr=0.71, tauf=0.8)
        dh, info = L.run(p, method=METHOD, steps=steps, check_every=2000,
                         tol=1e-9, verbose=False)
        if info["diverged"]:
            print(f"  diverged after {info['steps']} steps")
            rows.append({"Ra": ra, "diverged": True})
            continue
        m = L.analyse(dh, p)
        ref = L.DVD[ra]
        print(f"  steps  = {info['steps']}")
        print(f"  u_max  = {m['u_max']:8.3f}  (DVD {ref['u_max']:.3f})")
        print(f"  v_max  = {m['v_max']:8.3f}  (DVD {ref['v_max']:.3f})")
        print(f"  Nu_avg = {m['Nu']:8.3f}  (DVD {ref['Nu']:.3f})")
        print(f"  |psi|  = {m['psi']:8.3f}  (DVD {ref['psi']:.3f})")
        rows.append({"Ra": ra, "diverged": False, "m": m})

    csv_path = write_csv(rows)
    plot_path = make_plot(rows)
    print(f"\nWrote CSV  : {csv_path}")
    print(f"Wrote plot : {plot_path}")


if __name__ == "__main__":
    main()
