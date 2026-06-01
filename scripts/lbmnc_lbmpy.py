"""Reproduce the side-heated square-cavity natural convection of src/sec3/lbmnc.c
inside the walberla code-generation ecosystem, using lbmpy / pystencils.

This is the "lbmpy standalone" path: the same lbmpy that generates walberla's
LBM kernels drives a pure-Python (numpy + MSVC-compiled kernels) simulation.

Two coupled lattice-Boltzmann fields share one pystencils data handling:
  * hydro  : D2Q9, BGK (SRT) or MRT, compressible, Boussinesq buoyancy force
             fy = rbetag*(T - 0.5)
  * thermal: D2Q9 advection-diffusion (first-order equilibrium), advected by the
             hydro velocity field.  (lbmnc.c uses D2Q5; lbmpy 1.4 ships no D2Q5
             stencil, so we use a D2Q9 AD field, which is macroscopically
             equivalent for this passive-scalar transport.)

Boundary conditions (matching lbmnc.c):
  * velocity : NoSlip on all four walls (half-way bounce-back).
  * temperature: Dirichlet T=1 (left, hot), T=0 (right, cold);
                 adiabatic (Neumann, dT/dy=0) on top and bottom.

Setup (one-time, in this repo's .venv):
  1. Install the code-generation stack:
         .venv\\Scripts\\python.exe -m pip install lbmpy pystencils sympy
  2. pystencils compiles its CPU kernels with MSVC (cl.exe).  Its bundled
     Visual-Studio detection breaks for VS 2017+ install layouts (it looks for
     the legacy VC\\vcvarsall.bat and decodes the env dump as UTF-16).  Two small
     fixes in .venv\\Lib\\site-packages\\pystencils\\cpu\\msvc_detection.py are
     needed (already applied in this environment):
       * get_environment_from_vc_vars_file: wrap the whole command in an extra
         pair of quotes -> f'cmd /c ""{vc_vars_file}" {arch} && set"', and
         decode the output as 'mbcs' instead of 'utf-16le'.
       * get_vc_vars_path_via_environment_variable: fall back to the modern
         VC\\Auxiliary\\Build\\vcvarsall.bat (and a filesystem search) when the
         legacy path is absent.

Run:
    .venv\\Scripts\\python.exe scripts\\lbmnc_lbmpy.py
(If you ever disable the patch above, run from a "Developer Command Prompt for
VS" / a shell where vcvars64.bat has been sourced, so cl.exe is on PATH.)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import sympy as sp

from lbmpy import LBStencil, Stencil, Method, LBMConfig, LatticeBoltzmannStep
from lbmpy.boundaries import NoSlip, DiffusionDirichlet, NeumannByCopy
from pystencils import create_data_handling, make_slice, Target


ROOT_DIR = Path(__file__).resolve().parents[1]

METHOD_MAP = {"srt": Method.SRT, "mrt": Method.MRT}


# -----------------------------------------------------------------------------
# Physical / numerical parameters (identical to src/sec3/lbmnc.c)
# -----------------------------------------------------------------------------
def build_params(N=46, Ra=1.0e4, Pr=0.71, tauf=0.8):
    nu = (tauf - 0.5) / 3.0          # kinematic viscosity (lattice)
    chi = nu / Pr                    # thermal diffusivity (lattice)
    taug = 3.0 * chi + 0.5           # thermal relaxation time (cs^2 = 1/3)
    L = float(N)                     # cavity width in lattice units (half-way BB)
    rbetag = Ra * nu * chi / L**3    # Boussinesq coefficient
    return dict(N=N, Ra=Ra, Pr=Pr, tauf=tauf, nu=nu, chi=chi, taug=taug,
                L=L, rbetag=rbetag, omega_f=1.0 / tauf, omega_g=1.0 / taug)


def build_scenario(p, method="srt"):
    """Construct the coupled hydro + thermal scenario sharing one data handling."""
    N = p["N"]
    dh = create_data_handling((N, N), periodicity=False, default_target=Target.CPU,
                              default_ghost_layers=1)

    # Shared macroscopic fields.
    dh.add_array("u", values_per_cell=2)     # hydro velocity  (advects temperature)
    dh.add_array("rho", values_per_cell=1)   # hydro density
    dh.add_array("T", values_per_cell=1)     # temperature (= AD "density")
    T_field = dh.fields["T"]

    # --- Hydrodynamic step: D2Q9 SRT/MRT with Boussinesq buoyancy -----------
    buoyancy = (0, p["rbetag"] * (T_field.center - sp.Rational(1, 2)))
    hydro_cfg = LBMConfig(stencil=LBStencil(Stencil.D2Q9), method=METHOD_MAP[method],
                          relaxation_rate=p["omega_f"], compressible=True,
                          force=buoyancy)
    hydro = LatticeBoltzmannStep(data_handling=dh, name="hydro", lbm_config=hydro_cfg,
                                 velocity_data_name="u", density_data_name="rho",
                                 compute_velocity_in_every_step=True,
                                 compute_density_in_every_step=True)

    # --- Thermal step: D2Q9 advection-diffusion (linear equilibrium, SRT) ---
    thermal_cfg = LBMConfig(stencil=LBStencil(Stencil.D2Q9), method=Method.SRT,
                            relaxation_rate=p["omega_g"], compressible=True,
                            zero_centered=False, equilibrium_order=1,
                            velocity_input=dh.fields["u"])
    thermal = LatticeBoltzmannStep(data_handling=dh, name="thermal", lbm_config=thermal_cfg,
                                   density_data_name="T",
                                   compute_density_in_every_step=True)

    # --- Boundary conditions ------------------------------------------------
    hbh = hydro.boundary_handling
    for sl in (make_slice[0, :], make_slice[-1, :], make_slice[:, 0], make_slice[:, -1]):
        hbh.set_boundary(NoSlip(), sl)

    tbh = thermal.boundary_handling
    tbh.set_boundary(NeumannByCopy(), make_slice[:, 0])    # bottom adiabatic
    tbh.set_boundary(NeumannByCopy(), make_slice[:, -1])   # top adiabatic
    tbh.set_boundary(DiffusionDirichlet(1.0), make_slice[0, :])    # left hot
    tbh.set_boundary(DiffusionDirichlet(0.0), make_slice[-1, :])   # right cold

    # --- Initial condition: u=0, rho=1, T linear (hot left -> cold right) ---
    dh.fill("u", 0.0, ghost_layers=True, inner_ghost_layers=True)
    dh.fill("rho", 1.0, ghost_layers=True, inner_ghost_layers=True)
    Tprof = np.linspace(1.0, 0.0, N)                       # along x (i index)
    for b in dh.iterate(ghost_layers=False):
        arrT = b["T"]
        for i in range(arrT.shape[0]):
            arrT[i, :] = Tprof[i]
    hydro.set_pdf_fields_from_macroscopic_values()
    thermal.set_pdf_fields_from_macroscopic_values()
    return dh, hydro, thermal


def run(p, method="srt", steps=60000, check_every=2000, tol=1e-9, verbose=True):
    dh, hydro, thermal = build_scenario(p, method=method)
    u_prev = None
    diverged = False
    ran = 0
    for it in range(0, steps, check_every):
        for _ in range(check_every):
            hydro.time_step()
            thermal.time_step()
        ran += check_every
        hydro.post_run(); thermal.post_run()
        u = dh.gather_array("u").copy()
        if not np.all(np.isfinite(u)):
            diverged = True
            if verbose:
                print(f"step {ran:6d}  DIVERGED (non-finite velocity)")
            break
        if u_prev is not None:
            norm = float(np.max(np.abs(u - u_prev)))
            if verbose:
                print(f"step {ran:6d}  max|du| = {norm:.3e}")
            if norm < tol:
                break
        u_prev = u
    hydro.post_run(); thermal.post_run()
    return dh, dict(steps=ran, diverged=diverged)


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------
def fields_ij(dh, p):
    """Return (u_star, v_star, T) as (N, N) arrays indexed [i, j]."""
    u = np.asarray(dh.gather_array("u"))
    T = np.asarray(dh.gather_array("T"))
    if T.ndim == 3:
        T = T[..., 0]
    L, chi = p["L"], p["chi"]
    return u[..., 0] * L / chi, u[..., 1] * L / chi, T


def save_datanc(dh, p, out_dir):
    """Write datancu / datancv / datance in lbmnc.c's layout: one row per j
    (bottom -> top), columns i (left -> right), so the repo's read_matrix and
    plot helpers apply unchanged."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    us, vs, T = fields_ij(dh, p)
    # transpose [i, j] -> [j, i] so rows are j, columns i
    for name, arr in (("datancu", us.T), ("datancv", vs.T), ("datance", T.T)):
        with (out_dir / name).open("w", encoding="utf-8") as f:
            for row in arr:
                f.write(" " + " ".join(f"{v:10.8e}" for v in row) + "\n")
    return out_dir


def write_vtk(dh, name, out_dir):
    """Write a single VTK ImageData (.vti) snapshot with u, T, rho — the same
    format walberla emits, readable in ParaView."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    writer = dh.create_vtk_writer(str(out_dir / name), ["u", "T", "rho"])
    writer(0)
    return out_dir


def analyse(dh, p):
    """Quick in-script metrics (the canonical figure/CSV come from the plot
    script, which reuses the repo's Nusselt convention)."""
    N, L, chi = p["N"], p["L"], p["chi"]
    us, vs, T = fields_ij(dh, p)
    xc = (np.arange(N) + 0.5) / N
    yc = (np.arange(N) + 0.5) / N
    ic = jc = N // 2

    col = us[ic, :]
    j_um = int(np.argmax(col))
    u_max, y_um = col[j_um], yc[j_um]
    row = vs[:, jc]
    i_vm = int(np.argmax(row))
    v_max, x_vm = row[i_vm], xc[i_vm]

    # stream function over chi: cumsum of lattice u over y (dy=1) / chi
    # = cumsum(us)/L  since us = u_lattice * L / chi
    psi = np.cumsum(us, axis=1) / L
    psi_absmax = float(np.max(np.abs(psi)))

    # local Nusselt at hot wall (3-point one-sided, half-cell offset), then
    # trapezoidal average over y in [0,1] with wall extension (repo convention)
    h = 1.0 / N
    nu_local = (8.0 * 1.0 - 9.0 * T[0, :] + T[1, :]) / (3.0 * h)
    y_full = np.concatenate(([0.0], yc, [1.0]))
    nu_full = np.concatenate(([nu_local[0]], nu_local, [nu_local[-1]]))
    Nu_avg = float(np.trapezoid(nu_full, x=y_full))
    return dict(u_max=float(u_max), y_um=float(y_um), v_max=float(v_max),
                x_vm=float(x_vm), psi=psi_absmax, Nu=Nu_avg)


DVD = {
    1.0e3: dict(u_max=3.649, y_um=0.813, v_max=3.697, x_vm=0.178, Nu=1.118, psi=1.174),
    1.0e4: dict(u_max=16.178, y_um=0.823, v_max=19.617, x_vm=0.119, Nu=2.243, psi=5.071),
    1.0e5: dict(u_max=34.730, y_um=0.855, v_max=68.590, x_vm=0.066, Nu=4.519, psi=9.111),
    1.0e6: dict(u_max=64.630, y_um=0.850, v_max=219.36, x_vm=0.0379, Nu=8.800, psi=16.32),
}


def print_table(m, Ra):
    dvd = DVD[Ra]
    print("\n  quantity                 lbmpy        DVD(1983)    rel.err")
    print("  " + "-" * 58)
    rows = [("max u*  (vert. centre)", "u_max"), ("  y/L of u peak", "y_um"),
            ("max v*  (horiz. centre)", "v_max"), ("  x/L of v peak", "x_vm"),
            ("avg Nusselt (hot wall)", "Nu"), ("|psi|max / chi", "psi")]
    for label, key in rows:
        a, b = m[key], dvd[key]
        err = abs(a - b) / abs(b) * 100 if b else 0.0
        print(f"  {label:24s} {a:9.4f}  {b:9.4f}   {err:6.2f}%")


def main():
    method = "mrt"   # lbmnc.c selects MRT by default
    Ra = 1.0e4
    p = build_params(N=46, Ra=Ra, Pr=0.71, tauf=0.8)
    print(f"Method = {method.upper()}")
    for k in ("N", "Ra", "Pr", "nu", "chi", "taug", "L", "rbetag", "omega_f", "omega_g"):
        print(f"  {k:8s} = {p[k]}")
    dh, info = run(p, method=method, steps=60000, check_every=2000, tol=1e-9)
    print(f"\nconverged after {info['steps']} steps (diverged={info['diverged']})")

    out_dir = ROOT_DIR / "outputs" / "sec3" / "lbmnc_lbmpy"
    save_datanc(dh, p, out_dir)
    write_vtk(dh, "lbmnc_lbmpy", out_dir / "vtk")
    print(f"fields + VTK written to {out_dir}")

    print_table(analyse(dh, p), Ra)


if __name__ == "__main__":
    main()
