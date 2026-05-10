// karman.c
// 2D Karman vortex street (pure LBM, BGK, D2Q9)
//
// A long channel with a circular cylinder obstacle at (CX, CY). Flow is driven
// by a constant body force (Guo) in +x; x is periodic and top/bottom walls are
// halfway bounce-back. The cylinder is centered slightly off the channel axis
// (CY = NY/2 + 1) so its staircase representation is asymmetric, breaking the
// y-symmetry and seeding the vortex-shedding instability that develops above
// the critical Reynolds number Re_D ~ 47.
//
// Domain: NX=360, NY=80, D=20 (R=10), Re_D ~ 100.
//
// Output:
//   * snapshots of u, v, vorticity, solid mask at logarithmically spaced times
//   * a probe time series (u, v at a downstream point) for spectral analysis
//     of the Strouhal frequency.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 360
#define NY 80
#define NDIR 9
#define NSTEPS 50000                 // longer to let the shedding instability fully develop
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 5           // probe sampled often for FFT
#define TAU 0.55
#define OMEGA (1.0/TAU)
// Confined cylinder (blockage D/NY = 0.25) raises the critical Re from the
// free-stream value of 47 to roughly 60-70. F=6e-6 -> Re_D ~ 130, well past
// the confined-cylinder shedding onset. Mach number 0.11/cs ~ 0.19 (safe).
#define FORCE_X 6e-6
#define CX 80
#define CY 41                        // off-center to break y-symmetry
#define R_CYL 10                     // cylinder radius -> D=20
#define PROBE_X 200                  // 6 D downstream
#define PROBE_Y 50                   // 1 R above cylinder centerline (in the upper wake) for stronger v signal

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
static char solid[NX*NY];
static double vort[NX*NY];

void init_geometry() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int dx = x - CX, dy = y - CY;
            solid[IDX(x, y)] = (dx*dx + dy*dy <= R_CYL*R_CYL) ? 1 : 0;
        }
    }
}

void initialize() {
    init_geometry();
    // Anti-symmetric v perturbation localized in the wake region to kick-start
    // the Karman shedding instability. With only the small CY-offset asymmetry
    // and zero initial flow, the cylinder generates a near-symmetric steady
    // wake that takes very long to bifurcate into shedding. A 0.5%-of-U
    // explicit kick gets vortex roll-up going within a few thousand steps.
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            double dx_ = x - (CX + 25);   // wake center, ~1 D downstream
            double dy_ = y - CY;
            double r2 = dx_*dx_ + dy_*dy_;
            double pert = 0.0005 * (y > CY ? 1.0 : -1.0) * exp(-r2 / 200.0);
            v[i] = solid[i] ? 0.0 : pert;
            rho[i] = 1.0;
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu);
            }
        }
    }
}

void stream_collide() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) continue;
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double Fi = (1.0 - 0.5*OMEGA) * w[d] * FORCE_X *
                            (3.0*(cx[d] - u[i]) + 9.0*eu*cx[d]);
                double post = f[i*NDIR + d] - OMEGA*(f[i*NDIR + d] - feq) + Fi;
                int xp = (x + cx[d] + NX) % NX;
                int yp = y + cy[d];
                if (yp < 0 || yp >= NY) {
                    f2[i*NDIR + opp[d]] = post;
                } else if (solid[IDX(xp, yp)]) {
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
            if (solid[i]) {
                u[i] = 0.0; v[i] = 0.0; rho[i] = 1.0;
                continue;
            }
            double rr = 0, ru = 0, rv = 0;
            for (int d = 0; d < NDIR; ++d) {
                double ff = f[i*NDIR + d];
                rr += ff;
                ru += ff * cx[d];
                rv += ff * cy[d];
            }
            rho[i] = rr;
            u[i] = (ru + 0.5*FORCE_X) / rr;
            v[i] = rv / rr;
        }
    }
}

void compute_vorticity() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            if (solid[IDX(x, y)]) { vort[IDX(x, y)] = 0; continue; }
            int xp = (x + 1) % NX, xm = (x - 1 + NX) % NX;
            double dvdx = 0.5 * (v[IDX(xp, y)] - v[IDX(xm, y)]);
            double dudy;
            if (y == 0)         dudy = u[IDX(x, 1)] - u[IDX(x, 0)];
            else if (y == NY-1) dudy = u[IDX(x, NY-1)] - u[IDX(x, NY-2)];
            else                dudy = 0.5 * (u[IDX(x, y+1)] - u[IDX(x, y-1)]);
            vort[IDX(x, y)] = dvdx - dudy;
        }
    }
}

int output_snapshot(int step) {
    char fname[64];
    snprintf(fname, sizeof(fname), "karman_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,solid\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%d\n",
                    x, y, u[i], v[i], vort[i], solid[i]);
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

    FILE* hist = fopen("karman_probe.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open karman_probe.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,u_probe,v_probe\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (!solid[i] && fabs(u[i]) > umax) umax = fabs(u[i]);
            }
            int p = IDX(PROBE_X, PROBE_Y);
            fprintf(hist, "%d,%.9g,%.9g,%.9g\n", t, umax, u[p], v[p]);
        }
        stream_collide();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double umax_actual = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (!solid[i] && fabs(u[i]) > umax_actual) umax_actual = fabs(u[i]);
    }
    double Re_D = umax_actual * 2 * R_CYL / nu0;
    printf("Done. Snapshots: karman_snapshot_*.csv, probe: karman_probe.csv\n");
    printf("Parameters: NX=%d NY=%d D=%d Cyl=(%d,%d) NSTEPS=%d TAU=%.3f F=%g nu0=%.5f u_max=%.4f Re_D=%.0f\n",
           NX, NY, 2*R_CYL, CX, CY, NSTEPS, TAU, FORCE_X, nu0, umax_actual, Re_D);
    return 0;
}
