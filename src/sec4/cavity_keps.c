// cavity_keps.c
// 2D lid-driven cavity flow + standard k-epsilon turbulence model (D2Q9, BGK)
//
// Same geometry and boundary conditions as cavity.c (4-wall halfway BB with
// the top wall moving at (U_LID, 0)). The k-eps transport equations supply
// nu_t = Cmu * k^2 / eps, fed back into a locally varying tau_eff in the
// BGK collision so the eddy viscosity slows the primary vortex spin-up.
//
// Wall treatment for k, eps:
//   - All four walls use Dirichlet wall functions on the first off-wall cell.
//     k_wall = u_tau^2 / sqrt(Cmu),  eps_wall = u_tau^3 / (kappa * Delta y)
//     with u_tau computed locally from |du_tangential/dn| at the wall.
//   - For top wall: tangential velocity gradient uses (U_LID - u_first_cell).
//   - For other walls: gradient uses just the first-cell tangential component.
//
// Caveats:
//   - At Re ~ 400 the cavity is laminar, so the strict physical role of k-eps
//     is questionable. This file documents how the model behaves on a
//     classical wall-bounded recirculating flow.
//   - Wall functions are formally valid for y+ >= 30; here y+_1 << 30, so
//     they act more as a "k injector" than a strict log-law match.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 128
#define NY 128
#define NDIR 9
#define NSTEPS 30000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 100
#define TAU 0.55
#define OMEGA (1.0/TAU)
#define U_LID 0.05
#define KAPPA 0.41
#define WALL_DY 0.5
// k-eps integration time step
#define KEPS_DT 0.05

const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
const int opp[NDIR] = {0,3,4,1,2,7,8,5,6};

#define IDX(x,y) ((x) + NX*(y))
#define nu0 ((TAU - 0.5)/3.0)

#define Cmu 0.09
#define Ce1 1.44
#define Ce2 1.92
#define sig_k 1.0
#define sig_e 1.3

static double f_buf_a[NX*NY*NDIR];
static double f_buf_b[NX*NY*NDIR];
static double *f = f_buf_a;
static double *f2 = f_buf_b;
static double u[NX*NY], v[NX*NY], rho[NX*NY];
static double k[NX*NY], eps[NX*NY];
static double k_new[NX*NY], eps_new[NX*NY];
static double nut_field[NX*NY];
static double vort[NX*NY], psi[NX*NY];

void initialize() {
    // Seed k, eps so that nu_t/nu_0 ~ 0.05 initially.
    double k_seed = 0.005 * U_LID * U_LID;
    double nut_target = 0.05 * nu0;
    double eps_seed = Cmu * k_seed * k_seed / nut_target;

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            k[i] = k_seed;
            eps[i] = eps_seed;
            nut_field[i] = nut_target;
            for (int d = 0; d < NDIR; ++d) {
                f[i*NDIR + d] = w[d] * rho[i];
            }
        }
    }
}

void stream_collide_with_keps() {
    // Same geometry handling as cavity.c, but tau is local: tau_eff = 0.5 + 3*(nu0+nu_t)
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double tau_eff = 0.5 + 3.0 * (nu0 + nut_field[i]);
            double omega_eff = 1.0 / tau_eff;
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double post = f[i*NDIR + d] - omega_eff*(f[i*NDIR + d] - feq);
                int xp = x + cx[d];
                int yp = y + cy[d];
                if (yp >= NY) {
                    f2[i*NDIR + opp[d]] = post - 6.0 * w[d] * rho[i] * cx[d] * U_LID;
                } else if (yp < 0 || xp < 0 || xp >= NX) {
                    f2[i*NDIR + opp[d]] = post;
                } else {
                    f2[IDX(xp, yp)*NDIR + d] = post;
                }
            }
        }
    }
    double *tmp = f; f = f2; f2 = tmp;
}

void macroscopic() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double rr = 0, ru = 0, rv = 0;
            for (int d = 0; d < NDIR; ++d) {
                double ff = f[i*NDIR + d];
                rr += ff;
                ru += ff * cx[d];
                rv += ff * cy[d];
            }
            rho[i] = rr;
            u[i] = ru / rr;
            v[i] = rv / rr;
        }
    }
}

static void apply_wall_function(int i, double u_tangent_jump) {
    // Local wall shear -> u_tau, then standard log-law equilibrium values
    double shear = fabs(u_tangent_jump) / WALL_DY;
    double u_tau = sqrt(nu0 * shear + 1e-12);
    double k_wall = u_tau * u_tau / sqrt(Cmu);
    double eps_wall = u_tau * u_tau * u_tau / (KAPPA * WALL_DY);
    // Take the larger of the existing (already-set) and the new estimate
    // so corner cells (touched by two walls) get the dominant shear.
    if (k_new[i] < k_wall)   k_new[i]   = k_wall;
    if (eps_new[i] < eps_wall) eps_new[i] = eps_wall;
}

void update_kepsilon() {
    double dx = 1.0, dy = 1.0, dt = KEPS_DT;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            // Mirror at walls (Neumann zero-gradient) for the bulk transport.
            int ixp = (x+1 < NX) ? IDX(x+1, y) : i;
            int ixm = (x-1 >= 0) ? IDX(x-1, y) : i;
            int iyp = (y+1 < NY) ? IDX(x, y+1) : i;
            int iym = (y-1 >= 0) ? IDX(x, y-1) : i;
            // Velocity gradients with one-sided differences at walls
            double dudx = (x == 0)      ? (u[IDX(1, y)] - u[i])
                       : (x == NX-1)    ? (u[i] - u[IDX(NX-2, y)])
                       : 0.5 * (u[IDX(x+1, y)] - u[IDX(x-1, y)]);
            double dvdx = (x == 0)      ? (v[IDX(1, y)] - v[i])
                       : (x == NX-1)    ? (v[i] - v[IDX(NX-2, y)])
                       : 0.5 * (v[IDX(x+1, y)] - v[IDX(x-1, y)]);
            double dudy = (y == 0)      ? (u[IDX(x, 1)] - u[i])
                       : (y == NY-1)    ? (u[i] - u[IDX(x, NY-2)])
                       : 0.5 * (u[IDX(x, y+1)] - u[IDX(x, y-1)]);
            double dvdy = (y == 0)      ? (v[IDX(x, 1)] - v[i])
                       : (y == NY-1)    ? (v[i] - v[IDX(x, NY-2)])
                       : 0.5 * (v[IDX(x, y+1)] - v[IDX(x, y-1)]);
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double nut = Cmu * k[i]*k[i] / (eps[i] + 1e-12);
            double Pk = nut * S2;
            double Dk = nu0 + nut/sig_k;
            double De = nu0 + nut/sig_e;
            double lap_k = (k[ixp] + k[ixm] + k[iyp] + k[iym] - 4.0*k[i])/(dx*dx);
            double lap_e = (eps[ixp] + eps[ixm] + eps[iyp] + eps[iym] - 4.0*eps[i])/(dx*dx);
            // 1st-order upwind convection (cavity flow has reverse zones)
            double uk = u[i]*((u[i]>0 ? (k[i] - k[ixm]) : (k[ixp] - k[i]))/dx);
            double vk = v[i]*((v[i]>0 ? (k[i] - k[iym]) : (k[iyp] - k[i]))/dy);
            double ue = u[i]*((u[i]>0 ? (eps[i] - eps[ixm]) : (eps[ixp] - eps[i]))/dx);
            double ve = v[i]*((v[i]>0 ? (eps[i] - eps[iym]) : (eps[iyp] - eps[i]))/dy);
            k_new[i] = k[i] + dt * (Dk*lap_k - uk - vk + Pk - eps[i]);
            eps_new[i] = eps[i] + dt * (De*lap_e - ue - ve
                            + Ce1*Pk*eps[i]/(k[i]+1e-12)
                            - Ce2*eps[i]*eps[i]/(k[i]+1e-12));
        }
    }
    // Wall functions on the four wall layers (overwrite k_new, eps_new).
    // Top wall (y=NY-1): tangential gradient uses U_LID - u(NY-1)
    for (int x = 0; x < NX; ++x) {
        int i = IDX(x, NY-1);
        // Initialize for fresh top-wall calc
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, U_LID - u[i]);
    }
    // Bottom wall (y=0): u_wall = 0
    for (int x = 0; x < NX; ++x) {
        int i = IDX(x, 0);
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, u[i]);
    }
    // Left wall (x=0): v_wall = 0; tangential is v
    for (int y = 0; y < NY; ++y) {
        int i = IDX(0, y);
        // Don't reset corners — top/bottom already wrote them
        if (y == 0 || y == NY-1) continue;
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, v[i]);
    }
    // Right wall (x=NX-1)
    for (int y = 0; y < NY; ++y) {
        int i = IDX(NX-1, y);
        if (y == 0 || y == NY-1) continue;
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, v[i]);
    }
    // Apply non-negativity floor and commit
    for (int i = 0; i < NX*NY; ++i) {
        if (k_new[i] < 1e-10) k_new[i] = 1e-10;
        if (eps_new[i] < 1e-12) eps_new[i] = 1e-12;
        k[i] = k_new[i];
        eps[i] = eps_new[i];
        nut_field[i] = Cmu * k[i]*k[i] / (eps[i] + 1e-12);
    }
}

void compute_vorticity() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            double dvdx, dudy;
            if (x == 0)         dvdx = v[IDX(1, y)] - v[IDX(0, y)];
            else if (x == NX-1) dvdx = v[IDX(NX-1, y)] - v[IDX(NX-2, y)];
            else                dvdx = 0.5 * (v[IDX(x+1, y)] - v[IDX(x-1, y)]);
            if (y == 0)         dudy = u[IDX(x, 1)] - u[IDX(x, 0)];
            else if (y == NY-1) dudy = u[IDX(x, NY-1)] - u[IDX(x, NY-2)];
            else                dudy = 0.5 * (u[IDX(x, y+1)] - u[IDX(x, y-1)]);
            vort[IDX(x, y)] = dvdx - dudy;
        }
    }
}

void compute_streamfunction() {
    for (int x = 0; x < NX; ++x) psi[IDX(x, 0)] = 0.0;
    for (int y = 1; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            psi[IDX(x, y)] = psi[IDX(x, y-1)] + 0.5 * (u[IDX(x, y)] + u[IDX(x, y-1)]);
        }
    }
}

int output_snapshot(int step) {
    char fname[64];
    snprintf(fname, sizeof(fname), "cavity_keps_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    compute_streamfunction();
    fprintf(fp, "x,y,u,v,vorticity,psi,k,eps,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], psi[i], k[i], eps[i], nut_field[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main() {
    initialize();
    int snap_steps[SNAPSHOTS];
    for (int i = 0; i < SNAPSHOTS; ++i) {
        if (i == 0) snap_steps[i] = 0;
        else snap_steps[i] = (int)(NSTEPS * pow((double)i / (SNAPSHOTS - 1), 1.5));
    }

    FILE* hist = fopen("cavity_keps_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open cavity_keps_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,v_max,psi_min,k_mean,eps_mean,nut_mean\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();   // update u,v BEFORE snapshot/history
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, vmax = 0, k_sum = 0, eps_sum = 0, nut_sum = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                if (fabs(v[i]) > vmax) vmax = fabs(v[i]);
                k_sum += k[i];
                eps_sum += eps[i];
                nut_sum += nut_field[i];
            }
            double inv_n = 1.0 / (NX*NY);
            compute_streamfunction();
            double psi_min = 0;
            for (int i = 0; i < NX*NY; ++i) if (psi[i] < psi_min) psi_min = psi[i];
            fprintf(hist, "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    t, umax, vmax, psi_min, k_sum*inv_n, eps_sum*inv_n, nut_sum*inv_n);
        }
        update_kepsilon();
        stream_collide_with_keps();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double Re = U_LID * NX / nu0;
    printf("Done. Snapshots: cavity_keps_snapshot_*.csv, history: cavity_keps_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U_LID=%.3f nu0=%.5f Re=%.0f\n",
           NX, NY, NSTEPS, TAU, U_LID, nu0, Re);
    return 0;
}
