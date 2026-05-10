// cavity.c
// 2D lid-driven cavity flow (pure LBM, BGK, D2Q9)
//
// Square box [0, NX-1] x [0, NY-1] with halfway bounce-back walls. The top
// wall (lid) moves with velocity (U_LID, 0); the other three walls are
// stationary. Driven from rest, the flow develops a primary recirculating
// vortex below the lid and (at sufficient Re) secondary corner vortices.
//
// Re = U_LID * NX / nu_0. With NX=128, TAU=0.55 -> nu_0=1/60, U_LID=0.05,
// the Reynolds number is about 384, in the range covered by Ghia et al.
// (1982) for benchmark comparison.
//
// Output: snapshots of u, v, vorticity, streamfunction at multiple time
// steps; history of |u|_max, |v|_max, and minimum streamfunction (primary
// vortex strength) sampled every HISTORY_INTERVAL steps.
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
static double vort[NX*NY], psi[NX*NY];

void initialize() {
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
    // BGK collision then push streaming with wall handling.
    // Top wall (y=NY-1, moving): Ladd's halfway bounce-back with momentum
    //   correction f^post[opp(d)] = f^post[d] - 2 w[d] rho (c[d]·u_w) / cs^2
    //   = f^post[d] - 6 w[d] rho cx[d] U_LID  (cs^2 = 1/3, u_w = (U_LID, 0))
    // Other walls (stationary): plain halfway bounce-back.
    // Corner cells where xp and yp would both be out of bounds: the y check
    // wins (treated as top wall), so the moving correction is also applied
    // at the top corners — a small inaccuracy of standard practice.
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double post = f[i*NDIR + d] - OMEGA*(f[i*NDIR + d] - feq);
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

void compute_vorticity() {
    // Central differences in the interior, one-sided at walls
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
    // psi(x, 0) = 0; integrate u upward via trapezoidal rule.
    // psi has units of velocity*length so |psi|_max ~ U_LID * NX / 2.
    for (int x = 0; x < NX; ++x) psi[IDX(x, 0)] = 0.0;
    for (int y = 1; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            psi[IDX(x, y)] = psi[IDX(x, y-1)] + 0.5 * (u[IDX(x, y)] + u[IDX(x, y-1)]);
        }
    }
}

int output_snapshot(int step) {
    char fname[64];
    snprintf(fname, sizeof(fname), "cavity_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    compute_streamfunction();
    fprintf(fp, "x,y,u,v,vorticity,psi\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g\n", x, y, u[i], v[i], vort[i], psi[i]);
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

    FILE* hist = fopen("cavity_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open cavity_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,v_max,psi_min\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();   // update u,v from current f BEFORE snapshot/history
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, vmax = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                if (fabs(v[i]) > vmax) vmax = fabs(v[i]);
            }
            compute_streamfunction();
            double psi_min = 0;
            for (int i = 0; i < NX*NY; ++i) if (psi[i] < psi_min) psi_min = psi[i];
            fprintf(hist, "%d,%.9g,%.9g,%.9g\n", t, umax, vmax, psi_min);
        }
        stream_collide();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double Re = U_LID * NX / nu0;
    printf("Done. Snapshots: cavity_snapshot_*.csv, history: cavity_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U_LID=%.3f nu0=%.5f Re=%.0f\n",
           NX, NY, NSTEPS, TAU, U_LID, nu0, Re);
    return 0;
}
