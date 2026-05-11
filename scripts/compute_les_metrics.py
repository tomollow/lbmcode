"""One-off helper: extract LES vs pure LBM vs k-eps metrics from sec4 outputs.

Mirrors the comparison values in docs/sec4/keps_summary.md so les_summary.md
can quote concrete suppression ratios.
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "sec4"


def kelbm():
    # u_max (LES) vs laminar Poiseuille prediction
    # params from src/sec4/kelbm_les.c
    TAU = 0.55; NY = 60; FORCE_X = 5e-6
    nu0 = (TAU - 0.5) / 3.0
    u_max_lam = FORCE_X * NY**2 / (8.0 * nu0)

    df = pd.read_csv(OUT / "kelbm_les" / "kelbm_les_output.csv")
    u_max = float(df["u"].abs().max())
    nut = df["nut"].values
    nut_max = float(nut.max())
    nut_mean = float(nut.mean())
    # k-eps version for reference
    df_keps = pd.read_csv(OUT / "kelbm" / "kelbm_output.csv")
    u_max_keps = float(df_keps["u"].abs().max())
    print(f"[kelbm] u_max_LES = {u_max:.5f}, u_max_keps = {u_max_keps:.5f}, u_max_lam = {u_max_lam:.5f}")
    print(f"[kelbm] LES ratio = {u_max/u_max_lam:.3f}, keps ratio = {u_max_keps/u_max_lam:.3f}")
    print(f"[kelbm] nut/nu0 max = {nut_max/nu0:.3f}, mean = {nut_mean/nu0:.4f}")


def taylor_green():
    TAU = 1.0
    nu0 = (TAU - 0.5) / 3.0
    # Compare u_max at end step to analytical pure decay
    hist_les = pd.read_csv(OUT / "taylor_green_les" / "tg_les_history.csv")
    hist_keps = pd.read_csv(OUT / "taylor_green_keps" / "tg_keps_history.csv")
    end_les = hist_les.iloc[-1]
    end_keps = hist_keps.iloc[-1]
    print(f"[TG] end step u_max_LES = {end_les['u_max']:.5f}, pure_theory = {end_les['u_max_pure_theory']:.5f}, ratio = {end_les['u_max']/end_les['u_max_pure_theory']:.3f}")
    print(f"[TG] end step u_max_keps = {end_keps['u_max']:.5f}, pure_theory = {end_keps['u_max_pure_theory']:.5f}, ratio = {end_keps['u_max']/end_keps['u_max_pure_theory']:.3f}")
    print(f"[TG] LES nut_mean / nu0 (end) = {end_les['nut_mean']/nu0:.5f}")
    # Mean nut over history
    print(f"[TG] LES nut_mean / nu0 (history mean) = {hist_les['nut_mean'].mean()/nu0:.5f}")


def kelvin_helmholtz():
    TAU = 0.52
    nu0 = (TAU - 0.5) / 3.0
    hist_pure = pd.read_csv(OUT / "kelvin_helmholtz" / "kh_pure_history.csv") if (OUT / "kelvin_helmholtz" / "kh_pure_history.csv").exists() else None
    if hist_pure is None:
        # try alternative path
        cands = list((OUT / "kelvin_helmholtz").glob("*history*.csv"))
        hist_pure = pd.read_csv(cands[0]) if cands else None
    hist_keps = pd.read_csv(OUT / "kelvin_helmholtz_keps" / "kh_keps_history.csv")
    hist_les = pd.read_csv(OUT / "kelvin_helmholtz_les" / "kh_les_history.csv")
    vmax_pure = float(hist_pure["v_max"].max()) if hist_pure is not None else None
    vmax_keps = float(hist_keps["v_max"].max())
    vmax_les = float(hist_les["v_max"].max())
    print(f"[KH] v_max peak: pure={vmax_pure}, keps={vmax_keps:.5f}, les={vmax_les:.5f}")
    if vmax_pure:
        print(f"[KH] keps/pure = {vmax_keps/vmax_pure:.3f}, les/pure = {vmax_les/vmax_pure:.3f}")
    print(f"[KH] LES nut_mean / nu0 mean = {hist_les['nut_mean'].mean()/nu0:.5f}, max = {hist_les['nut_mean'].max()/nu0:.5f}")


def cavity():
    TAU = 0.55
    nu0 = (TAU - 0.5) / 3.0
    hist_pure = None
    for cand in (OUT / "cavity").glob("*history*.csv"):
        hist_pure = pd.read_csv(cand)
        break
    hist_keps = pd.read_csv(OUT / "cavity_keps" / "cavity_keps_history.csv")
    hist_les = pd.read_csv(OUT / "cavity_les" / "cavity_les_history.csv")
    psi_pure = float(hist_pure["psi_min"].iloc[-1]) if hist_pure is not None else None
    psi_keps = float(hist_keps["psi_min"].iloc[-1])
    psi_les = float(hist_les["psi_min"].iloc[-1])
    print(f"[Cavity] psi_min final: pure={psi_pure}, keps={psi_keps:.5f}, les={psi_les:.5f}")
    if psi_pure:
        print(f"[Cavity] keps/pure = {psi_keps/psi_pure:.3f}, les/pure = {psi_les/psi_pure:.3f}")
    print(f"[Cavity] LES nut_mean / nu0 mean = {hist_les['nut_mean'].mean()/nu0:.5f}, end = {hist_les['nut_mean'].iloc[-1]/nu0:.5f}")


def backward_step():
    TAU = 0.55; STEP_HEIGHT = 30
    nu0 = (TAU - 0.5) / 3.0
    hist_pure = None
    for cand in (OUT / "backward_step").glob("*history*.csv"):
        hist_pure = pd.read_csv(cand)
        break
    hist_keps = pd.read_csv(OUT / "backward_step_keps" / "step_keps_history.csv")
    hist_les = pd.read_csv(OUT / "backward_step_les" / "step_les_history.csv")
    # use reattach_x at final step in lattice units, divided by STEP_HEIGHT and shifted by STEP_LENGTH=30
    def xr(hist):
        if hist is None:
            return None
        end = hist.iloc[-1]
        ra = end["reattach_x"]
        if ra <= 0:
            return None
        return (float(ra) - 30) / STEP_HEIGHT  # x_R / H downstream of step
    xr_pure = xr(hist_pure)
    xr_keps = xr(hist_keps)
    xr_les = xr(hist_les)
    print(f"[Step] x_R/H final: pure={xr_pure}, keps={xr_keps}, les={xr_les}")
    if xr_pure and xr_keps and xr_les:
        print(f"[Step] keps/pure = {xr_keps/xr_pure:.3f}, les/pure = {xr_les/xr_pure:.3f}")
    print(f"[Step] LES nut_mean / nu0 mean = {hist_les['nut_mean'].mean()/nu0:.5f}, end = {hist_les['nut_mean'].iloc[-1]/nu0:.5f}")


def karman():
    TAU = 0.55
    nu0 = (TAU - 0.5) / 3.0
    probe_pure = None
    for cand in (OUT / "karman").glob("*probe*.csv"):
        probe_pure = pd.read_csv(cand)
        break
    probe_keps = pd.read_csv(OUT / "karman_keps" / "karman_keps_probe.csv")
    probe_les = pd.read_csv(OUT / "karman_les" / "karman_les_probe.csv")
    def amp(p):
        if p is None:
            return None
        # Use second half to skip transient
        v = p["v_probe"].values[len(p)//2:]
        return float(v.max() - v.min()) / 2
    a_pure = amp(probe_pure)
    a_keps = amp(probe_keps)
    a_les = amp(probe_les)
    print(f"[Karman] probe-v amplitude (half-PtoP): pure={a_pure}, keps={a_keps:.6f}, les={a_les:.6f}")
    if a_pure:
        print(f"[Karman] keps/pure = {a_keps/a_pure:.4f}, les/pure = {a_les/a_pure:.4f}")
    print(f"[Karman] LES nut_mean / nu0 mean = {probe_les['nut_mean'].mean()/nu0:.5f}, max = {probe_les['nut_mean'].max()/nu0:.5f}")


if __name__ == "__main__":
    kelbm()
    print()
    taylor_green()
    print()
    kelvin_helmholtz()
    print()
    cavity()
    print()
    backward_step()
    print()
    karman()
