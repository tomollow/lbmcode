// cavity_des.c
// 2D lid-driven cavity + Spalart-Allmaras DES97 (D2Q9, BGK)
//
// Same geometry/BC as cavity_les.c (4-wall halfway BB with top wall at U_LID).
// SA-DES hybrid: a 1-equation transport for the SA working variable nu_tilde
// with the wall distance replaced by a length-scale switch
//
//     d_tilde = min(d_wall, C_DES * Delta),    C_DES = 0.65, Delta = 1 LU
//
// In wall-adjacent cells where d_wall < C_DES*Delta the model behaves like the
// original SA RANS (destruction scales with the true wall distance). In the
// bulk d_tilde saturates at C_DES*Delta and the equilibrium destruction takes
// a Smagorinsky-like form (nu_t ~ (C_DES Delta)^2 |S|) — i.e. the "LES branch".
//
// At Re ~ 384 the cavity is laminar: |S| is small, the SA production is small,
// and nu_tilde decays toward a tiny equilibrium ~ O(1e-3 nu_0) in the bulk
// and ~ 0 next to walls. This file is therefore a methodological consistency
// check: DES sleeps when there is nothing to model, exactly like LES.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sa_closure.h"

#define NX 128
#define NY 128
#define NDIR 9
#define NSTEPS 30000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 100
#define TAU 0.55
#define OMEGA (1.0/TAU)
#define U_LID 0.05

// DES97 length-scale switch
#define C_DES 0.65
#define DELTA_LES 1.0
// SA transport time step (explicit Euler)
#define SA_DT 0.05

const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
const int opp[NDIR] = {0,3,4,1,2,7,8,5,6};

#define IDX(x,y) ((x) + NX*(y))
#define nu0 ((TAU - 0.5)/3.0)

static double f_buf_a[NX*NY*NDIR];
static double f_buf_b[NX*NY*NDIR];
static double *f = f_buf_a;
static double *f2 = f_buf_b;
static double u[NX*NY], v[NX*NY], rho[NX*NY];
static double nu_tilde[NX*NY], nu_tilde_new[NX*NY];
static double nut_field[NX*NY];
static double d_wall[NX*NY], d_tilde[NX*NY];
static double vort[NX*NY], psi[NX*NY];

void initialize() {
    // Standard SA free-stream seed: chi = nu_tilde/nu0 = 3 -> nu_t/nu0 ~ 0.07
    // The SA destruction term will drive nu_tilde toward its natural laminar
    // equilibrium over the run; seeding above the equilibrium lets us watch
    // the model "go to sleep" rather than starting flat-lined.
    double nut_seed = 3.0 * nu0;
    double fv1_seed = sa_fv1(nut_seed / nu0);

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            nu_tilde[i] = nut_seed;
            nut_field[i] = nut_seed * fv1_seed;
            // Wall distance: halfway BB places the wall at the cell-face,
            // so the first off-wall cell is 0.5 LU away.
            double dl = x + 0.5;
            double dr = NX - x - 0.5;
            double db = y + 0.5;
            double dt = NY - y - 0.5;
            double d = dl < dr ? dl : dr;
            if (db < d) d = db;
            if (dt < d) d = dt;
            d_wall[i] = d;
            d_tilde[i] = (d < C_DES*DELTA_LES) ? d : C_DES*DELTA_LES;
            for (int dd = 0; dd < NDIR; ++dd) {
                f[i*NDIR + dd] = w[dd] * rho[i];
            }
        }
    }
}

void stream_collide_with_des() {
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

void update_sa_des() {
    double dx = 1.0, dy = 1.0, dt = SA_DT;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            // Mirror at walls (Neumann zero-gradient) — same idiom as the k-eps
            // version. The destruction term carries the wall info via d_tilde.
            int ixp = (x+1 < NX) ? IDX(x+1, y) : i;
            int ixm = (x-1 >= 0) ? IDX(x-1, y) : i;
            int iyp = (y+1 < NY) ? IDX(x, y+1) : i;
            int iym = (y-1 >= 0) ? IDX(x, y-1) : i;
            // Velocity gradients (one-sided at walls)
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
            double Smag = sqrt(S2);

            // SA closure (Stilde, fw) — see sa_closure.h
            double nt = nu_tilde[i];
            double dl = d_tilde[i];
            sa_terms_t sa = sa_compute_terms(nt, nu0, dl, Smag);

            // Source: production - destruction + diffusion + cross-diffusion
            double prod = C_B1 * sa.Stilde * nt;
            double dest = C_W1 * sa.fw * (nt/dl) * (nt/dl);
            double lap = (nu_tilde[ixp] + nu_tilde[ixm] + nu_tilde[iyp] + nu_tilde[iym]
                          - 4.0*nt) / (dx*dx);
            double dntdx = 0.5 * (nu_tilde[ixp] - nu_tilde[ixm]) / dx;
            double dntdy = 0.5 * (nu_tilde[iyp] - nu_tilde[iym]) / dy;
            double grad2 = dntdx*dntdx + dntdy*dntdy;
            double diff = ((nu0 + nt)*lap + (1.0 + C_B2)*grad2) / SIG_SA;

            // 1st-order upwind convection
            double conv_x = u[i] * ((u[i] > 0) ? (nt - nu_tilde[ixm])/dx
                                              : (nu_tilde[ixp] - nt)/dx);
            double conv_y = v[i] * ((v[i] > 0) ? (nt - nu_tilde[iym])/dy
                                              : (nu_tilde[iyp] - nt)/dy);
            nu_tilde_new[i] = nt + dt * (prod - dest + diff - conv_x - conv_y);
        }
    }
    // Floor and recompute nut = nu_tilde * fv1(chi)
    for (int i = 0; i < NX*NY; ++i) {
        if (nu_tilde_new[i] < 1e-12) nu_tilde_new[i] = 1e-12;
        nu_tilde[i] = nu_tilde_new[i];
        nut_field[i] = nu_tilde[i] * sa_fv1(nu_tilde[i] / nu0);
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
    snprintf(fname, sizeof(fname), "cavity_des_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    compute_streamfunction();
    fprintf(fp, "x,y,u,v,vorticity,psi,nut,nu_tilde,d_wall,d_tilde\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], psi[i],
                    nut_field[i], nu_tilde[i], d_wall[i], d_tilde[i]);
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

    FILE* hist = fopen("cavity_des_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open cavity_des_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,v_max,psi_min,nut_mean,nu_tilde_mean,les_frac\n");

    // LES-branch fraction = cells where d_wall >= C_DES * Delta (precomputed)
    int les_count = 0;
    for (int i = 0; i < NX*NY; ++i) if (d_wall[i] >= C_DES*DELTA_LES) ++les_count;
    double les_frac = (double)les_count / (NX*NY);

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, vmax = 0, nut_sum = 0, nt_sum = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                if (fabs(v[i]) > vmax) vmax = fabs(v[i]);
                nut_sum += nut_field[i];
                nt_sum += nu_tilde[i];
            }
            double inv_n = 1.0 / (NX*NY);
            compute_streamfunction();
            double psi_min = 0;
            for (int i = 0; i < NX*NY; ++i) if (psi[i] < psi_min) psi_min = psi[i];
            fprintf(hist, "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    t, umax, vmax, psi_min, nut_sum*inv_n, nt_sum*inv_n, les_frac);
        }
        update_sa_des();
        stream_collide_with_des();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double Re = U_LID * NX / nu0;
    printf("Done. Snapshots: cavity_des_snapshot_*.csv, history: cavity_des_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U_LID=%.3f C_DES=%.3f Delta=%.1f nu0=%.5f Re=%.0f\n",
           NX, NY, NSTEPS, TAU, U_LID, C_DES, DELTA_LES, nu0, Re);
    printf("DES branch: LES fraction = %.4f (RANS layer = first off-wall cell only)\n", les_frac);
    return 0;
}
