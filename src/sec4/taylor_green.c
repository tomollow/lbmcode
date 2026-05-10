// taylor_green.c
// 2D Taylor-Green vortex (pure LBM, BGK, D2Q9, fully periodic)
// Initial condition:  u = -U0 cos(kx*x) sin(ky*y),  v =  U0 sin(kx*x) cos(ky*y)
//                     rho = rho0 - 0.25 U0^2 [cos(2 kx x) + cos(2 ky y)] / cs^2
// Analytical decay:   u(x,y,t) = u(x,y,0) * exp(-(kx^2 + ky^2) * nu * t)
// Output: snapshots of u, v, vorticity at multiple time steps, plus
//         the centerline u-amplitude history for log-decay validation.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 128
#define NY 128
#define NDIR 9
#define NSTEPS 3000
#define SNAPSHOTS 6                  // number of full-field snapshots
#define HISTORY_INTERVAL 25          // steps between amplitude-history samples
#define TAU 1.0
#define OMEGA (1.0/TAU)
#define U0 0.04                      // peak initial speed (keep |u| << 0.1 for Mach safety)
#define M_PI_F 3.14159265358979323846

// D2Q9 lattice
const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

#define IDX(x,y) ((x) + NX*(y))
#define nu0 ((TAU - 0.5)/3.0)

static double f_buf_a[NX*NY*NDIR];
static double f_buf_b[NX*NY*NDIR];
static double *f = f_buf_a;
static double *f2 = f_buf_b;
static double u[NX*NY], v[NX*NY], rho[NX*NY];

void initialize() {
    double kx = 2.0 * M_PI_F / (double)NX;
    double ky = 2.0 * M_PI_F / (double)NY;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            // analytical Taylor-Green initial state
            u[i] = -U0 * cos(kx * x) * sin(ky * y);
            v[i] =  U0 * sin(kx * x) * cos(ky * y);
            // pressure perturbation translated to density (cs^2 = 1/3)
            rho[i] = 1.0 - 0.75 * U0 * U0 * (cos(2.0 * kx * x) + cos(2.0 * ky * y));
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
            }
        }
    }
}

void stream_collide() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double post = f[i*NDIR + d] - OMEGA*(f[i*NDIR + d] - feq);
                int xp = (x + cx[d] + NX) % NX;
                int yp = (y + cy[d] + NY) % NY;
                f2[IDX(xp, yp)*NDIR + d] = post;
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

// vorticity field (z-component): omega = dv/dx - du/dy, central difference + periodic
static double vort[NX*NY];
void compute_vorticity() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int xp = (x + 1) % NX, xm = (x - 1 + NX) % NX;
            int yp = (y + 1) % NY, ym = (y - 1 + NY) % NY;
            double dvdx = 0.5 * (v[IDX(xp, y)] - v[IDX(xm, y)]);
            double dudy = 0.5 * (u[IDX(x, yp)] - u[IDX(x, ym)]);
            vort[IDX(x, y)] = dvdx - dudy;
        }
    }
}

int output_snapshot(int step) {
    char fname[64];
    snprintf(fname, sizeof(fname), "tg_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,rho\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g\n", x, y, u[i], v[i], vort[i], rho[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main() {
    initialize();
    macroscopic();

    // Schedule snapshot steps logarithmically over [0, NSTEPS]
    int snap_steps[SNAPSHOTS];
    for (int i = 0; i < SNAPSHOTS; ++i) {
        // 0, then geometric progression to NSTEPS
        if (i == 0) snap_steps[i] = 0;
        else snap_steps[i] = (int)(NSTEPS * pow((double)(i) / (SNAPSHOTS - 1), 1.5));
    }

    // amplitude history (peak |u|) sampled every HISTORY_INTERVAL steps
    FILE* hist = fopen("tg_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open tg_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,u_max_theory\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0;
            for (int i = 0; i < NX*NY; ++i) {
                double um = fabs(u[i]);
                if (um > umax) umax = um;
            }
            double k2 = pow(2.0 * M_PI_F / NX, 2) + pow(2.0 * M_PI_F / NY, 2);
            double u_theory = U0 * exp(-k2 * nu0 * t);
            fprintf(hist, "%d,%.9g,%.9g\n", t, umax, u_theory);
        }
        macroscopic();
        stream_collide();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    printf("Done. Snapshots: tg_snapshot_*.csv, history: tg_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U0=%.3f nu0=%.5f\n",
           NX, NY, NSTEPS, TAU, U0, nu0);
    return 0;
}
