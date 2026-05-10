// backward_step.c
// 2D backward-facing step (pure LBM, BGK, D2Q9)
//
// A long channel of NX x NY is driven by a constant body force in +x. A solid
// block of size STEP_LENGTH x STEP_HEIGHT sits in the lower-left corner; the
// flow accelerates through the narrow section above the block, then expands
// downstream of the block trailing edge and forms a recirculation zone before
// reattaching to the bottom wall. x is periodic so the flow loops through the
// step on every passage; the long downstream segment lets the recirculation
// fully develop.
//
// Geometry (NX=240, NY=60, STEP=30x30, ER=2):
//   ┌──────────────────────────────────────┐
//   │                                      │   y = NY-1 (top wall)
//   │  flow direction →                    │
//   │█████┐                                │
//   │█████│                                │   y = STEP_HEIGHT
//   │█████│                                │
//   └─────┴──────────────────────────────────┘   y = 0  (bottom wall)
//   x=0   STEP_LENGTH                       NX-1
//
// Output: snapshots of u, v, vorticity, streamfunction, plus a wall-shear
// history along y=STEP_HEIGHT (top of step) and y=0 (downstream of step) so
// the reattachment point can be located as the sign change of du/dy at y=0.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 240
#define NY 60
#define NDIR 9
#define NSTEPS 30000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 100
#define TAU 0.55
#define OMEGA (1.0/TAU)
#define FORCE_X 2e-6                 // body force per unit mass in +x
#define STEP_LENGTH 30
#define STEP_HEIGHT 30               // expansion ratio = NY / (NY - STEP_HEIGHT) = 2

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
static char solid[NX*NY];           // 1 if cell is inside the step block
static double vort[NX*NY], psi[NX*NY];

void init_geometry() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            solid[IDX(x, y)] = (x < STEP_LENGTH && y < STEP_HEIGHT) ? 1 : 0;
        }
    }
}

void initialize() {
    init_geometry();
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            for (int d = 0; d < NDIR; ++d) {
                f[i*NDIR + d] = w[d] * rho[i];
            }
        }
    }
}

void stream_collide() {
    // BGK + Guo forcing, halfway bounce-back at top/bottom walls AND at solid block surfaces.
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) continue;     // skip solid cells entirely
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double Fi = (1.0 - 0.5*OMEGA) * w[d] * FORCE_X *
                            (3.0*(cx[d] - u[i]) + 9.0*eu*cx[d]);
                double post = f[i*NDIR + d] - OMEGA*(f[i*NDIR + d] - feq) + Fi;
                int xp = (x + cx[d] + NX) % NX;     // periodic in x
                int yp = y + cy[d];
                if (yp < 0 || yp >= NY) {
                    // top/bottom domain walls (stationary)
                    f2[i*NDIR + opp[d]] = post;
                } else if (solid[IDX(xp, yp)]) {
                    // streaming into solid block surface
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
            u[i] = (ru + 0.5*FORCE_X) / rr;     // Guo half-step body-force correction
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

void compute_streamfunction() {
    // Integrate u upward from y=0. Solid cells contribute 0.
    for (int x = 0; x < NX; ++x) psi[IDX(x, 0)] = 0.0;
    for (int y = 1; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            double u_avg = 0.5 * (u[IDX(x, y)] + u[IDX(x, y-1)]);
            psi[IDX(x, y)] = psi[IDX(x, y-1)] + u_avg;
        }
    }
}

int output_snapshot(int step) {
    char fname[64];
    snprintf(fname, sizeof(fname), "step_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    compute_streamfunction();
    fprintf(fp, "x,y,u,v,vorticity,psi,solid\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d\n",
                    x, y, u[i], v[i], vort[i], psi[i], solid[i]);
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

    FILE* hist = fopen("step_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open step_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,reattach_x\n");

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
            // Find reattachment: first x downstream of step where u(x, 1) > 0
            int reattach = -1;
            for (int xq = STEP_LENGTH; xq < NX; ++xq) {
                if (u[IDX(xq, 1)] > 0) { reattach = xq; break; }
            }
            fprintf(hist, "%d,%.9g,%d\n", t, umax, reattach);
        }
        stream_collide();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    // Re_H computed from the actually-observed u_max so it matches the
    // headline number quoted in docs and the history CSV.
    double umax_actual = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (!solid[i] && fabs(u[i]) > umax_actual) umax_actual = fabs(u[i]);
    }
    double Re_H = umax_actual * STEP_HEIGHT / nu0;
    printf("Done. Snapshots: step_snapshot_*.csv, history: step_history.csv\n");
    printf("Parameters: NX=%d NY=%d STEP=%dx%d NSTEPS=%d TAU=%.3f F=%g nu0=%.5f u_max=%.4f Re_H=%.0f\n",
           NX, NY, STEP_LENGTH, STEP_HEIGHT, NSTEPS, TAU, FORCE_X, nu0, umax_actual, Re_H);
    return 0;
}
