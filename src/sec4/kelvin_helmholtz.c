// kelvin_helmholtz.c
// 2D Kelvin-Helmholtz instability (pure LBM, BGK, D2Q9)
//
// Two parallel streams with opposite x-velocity meet across a thin shear layer.
// A small sinusoidal perturbation in v seeds the instability; the linearly
// unstable mode grows exponentially and rolls up into a row of vortices that
// later merge through pairing.
//
// Domain: NX (streamwise) x NY (cross-stream), periodic in BOTH directions.
// Two shear layers are used (top and bottom) so periodicity is consistent —
// each is centered on y = NY/4 and y = 3*NY/4.
//
// Initial mean profile (smooth tanh shear):
//   u(y) = U0 * tanh((y - NY/4) / delta) for y in [0, NY/2]
//   u(y) = -U0 * tanh((y - 3 NY/4) / delta) for y in [NY/2, NY]
// Perturbation in v:
//   v(x, y) = A * sin(2*pi*x/NX) * exp(-((y - y_c) / delta_pert)^2)
//
// Output: snapshots of u, v, vorticity at multiple time steps + a peak |v|
// history that reveals the exponential growth phase before saturation.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 256
#define NY 128
#define NDIR 9
#define NSTEPS 8000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 25
#define TAU 0.52                         // close to 0.5 → low nu so the shear layer doesn't diffuse before roll-up
#define OMEGA (1.0/TAU)
#define U0 0.05                          // half the velocity jump (top has +U0, bottom -U0)
#define SHEAR_DELTA 4.0                  // shear-layer thickness in lattice units
#define PERT_AMP 0.002                   // initial perturbation amplitude (= 0.04 * U0)
#define PERT_SIGMA 8.0                   // Gaussian envelope half-width
#define PERT_MODE 4                      // streamwise mode number (lambda = NX/PERT_MODE = 14*delta is near most-unstable)
#define M_PI_F 3.14159265358979323846

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
static double vort[NX*NY];

static double shear_profile(double y) {
    // Two opposite shear layers in the periodic box, centered at NY/4 and 3*NY/4
    double y1 = NY * 0.25, y2 = NY * 0.75;
    if (y < NY * 0.5) {
        return  U0 * tanh((y - y1) / SHEAR_DELTA);
    } else {
        return -U0 * tanh((y - y2) / SHEAR_DELTA);
    }
}

void initialize() {
    double kx = 2.0 * M_PI_F / (double)NX;
    double y1 = NY * 0.25, y2 = NY * 0.75;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double sx = sin(PERT_MODE * kx * x);
            double dy1 = (y - y1) / PERT_SIGMA;
            double dy2 = (y - y2) / PERT_SIGMA;
            // Perturbation localized to both shear layers, opposite sign on each
            double pert = PERT_AMP * sx * (exp(-dy1*dy1) - exp(-dy2*dy2));
            u[i] = shear_profile((double)y);
            v[i] = pert;
            rho[i] = 1.0;
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
    snprintf(fname, sizeof(fname), "kh_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g\n", x, y, u[i], v[i], vort[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main() {
    initialize();
    macroscopic();

    int snap_steps[SNAPSHOTS];
    for (int i = 0; i < SNAPSHOTS; ++i) {
        if (i == 0) snap_steps[i] = 0;
        else snap_steps[i] = (int)(NSTEPS * pow((double)i / (SNAPSHOTS - 1), 1.2));
    }

    FILE* hist = fopen("kh_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open kh_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,v_max,u_rms\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double vmax = 0;
            double u_sq = 0;
            for (int i = 0; i < NX*NY; ++i) {
                double vm = fabs(v[i]);
                if (vm > vmax) vmax = vm;
                u_sq += u[i]*u[i];
            }
            double u_rms = sqrt(u_sq / (NX*NY));
            fprintf(hist, "%d,%.9g,%.9g\n", t, vmax, u_rms);
        }
        macroscopic();
        stream_collide();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    printf("Done. Snapshots: kh_snapshot_*.csv, history: kh_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U0=%.3f delta=%.1f nu0=%.5f\n",
           NX, NY, NSTEPS, TAU, U0, SHEAR_DELTA, nu0);
    return 0;
}
