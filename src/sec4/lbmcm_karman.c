// lbmcm_karman.c
// 2D flow past a rectangular obstacle solved with the central-moment (cascaded)
// D2Q9 LBM of Geier et al. 2006. The geometry/body-force/periodic-x driver
// follows karman.c; the collision operator is the central-moment block of
// lbmcm.c (flag == 3).
//
// Purpose: qualitatively reproduce Fig. 4.5 of Seta's "LBM" book (Sec. 4.6,
// after Geier 2006) — a turbulent wake behind a rectangular obstacle and the
// associated 1D energy spectrum E(k) showing an approximate -5/3 inertial
// range. We cannot reach Re = 1.4e6, but the central-moment operator is stable
// in the tau -> 0.5 regime where SRT/MRT would diverge, so the qualitative
// picture (broadband wake + power-law spectrum) is reachable on a workstation.
//
// Caveat: the body force here is plain Guo forcing (half-force velocity
// correction + post-collision F_i), which is the BGK-era recipe and is NOT
// strictly consistent with central-moment relaxation. The combination is fine
// for the qualitative demonstration in this file but should not be used for
// benchmark-grade work; use a moment-space (or Strang-split) CM forcing
// scheme there.
//
//   6  2  5
//      |
//   3--0--1
//      |
//   7  4  8
//
// Output
//   * lbmcm_karman_snapshot_%05d.csv  — u, v, |U|, solid at log-spaced times
//   * lbmcm_karman_wake_%05d.csv      — wake-window u, v frames for spectrum
//   * lbmcm_karman_probe.csv          — probe time series (FFT validation)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NX 256
#define NY 80
#define NDIR 9
#define NSTEPS 25000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 5

// Wake-window frames for spatial spectrum start once the wake is developed.
#define WAKE_FIRST_STEP   12000
#define WAKE_STRIDE        500

// Central-moment LBM tolerates tau close to 0.5. SRT would diverge here.
#define TAU 0.503
#define OMEGA (1.0/TAU)
#define FORCE_X 4.0e-6

// Rectangular obstacle (vertical bar, similar shape to the bluff body in
// Geier's figure). Width 4 lu, height 16 lu.
#define OBS_X0 48
#define OBS_X1 52
#define OBS_Y0 32
#define OBS_Y1 48

#define PROBE_X 180
#define PROBE_Y 48

// Wake-window for the spatial energy spectrum: downstream of the obstacle,
// leaving room from periodic-x seam.
#define SPEC_X0 96
#define SPEC_X1 (NX - 8)
#define SPEC_Y0 8
#define SPEC_Y1 (NY - 8)

#define IDX(x, y)   ((x) + NX*(y))
#define FIDX(i, d)  ((i)*NDIR + (d))
#define NU0         ((TAU - 0.5)/3.0)

static const int    cx[NDIR]   = { 0, 1, 0,-1, 0, 1,-1,-1, 1};
static const int    cy[NDIR]   = { 0, 0, 1, 0,-1, 1, 1,-1,-1};
static const double w[NDIR]   = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9,
                                  1.0/36, 1.0/36, 1.0/36, 1.0/36};
static const int    opp[NDIR]  = { 0, 3, 4, 1, 2, 7, 8, 5, 6};

static double f_buf_a[NX*NY*NDIR];
static double f_buf_b[NX*NY*NDIR];
static double *f  = f_buf_a;
static double *f2 = f_buf_b;
static double u[NX*NY], v[NX*NY], rho[NX*NY];
static char   solid[NX*NY];

// Raw-moment basis (rows match lbmcm.c flag==3) and its inverse.
static double Mmat[9][9];
static double Minv[9][9];
// Relaxation diagonal: shear modes at 1/tau, others at 1 (Geier-style).
static double Smat[9];

static void init_cm_matrices(void) {
    for (int k = 0; k < 9; ++k) {
        double a = (double)cx[k], b = (double)cy[k];
        Mmat[0][k] = 1.0;
        Mmat[1][k] = a;
        Mmat[2][k] = b;
        Mmat[3][k] = a*a + b*b;
        Mmat[4][k] = a*a - b*b;
        Mmat[5][k] = a*b;
        Mmat[6][k] = a*a*b;
        Mmat[7][k] = a*b*b;
        Mmat[8][k] = a*a*b*b;
    }
    // M^-1 copied verbatim from lbmcm.c flag==3 block.
    double mi[9][9] = {
        { 1.0,   0.0,   0.0,  -1.0,   0.0,   0.0,   0.0,   0.0,   1.0},
        { 0.0,   0.5,   0.0,  0.25,  0.25,   0.0,   0.0,  -0.5,  -0.5},
        { 0.0,   0.0,   0.5,  0.25, -0.25,   0.0,  -0.5,   0.0,  -0.5},
        { 0.0,  -0.5,   0.0,  0.25,  0.25,   0.0,   0.0,   0.5,  -0.5},
        { 0.0,   0.0,  -0.5,  0.25, -0.25,   0.0,   0.5,   0.0,  -0.5},
        { 0.0,   0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25},
        { 0.0,   0.0,   0.0,   0.0,   0.0, -0.25,  0.25, -0.25,  0.25},
        { 0.0,   0.0,   0.0,   0.0,   0.0,  0.25, -0.25, -0.25,  0.25},
        { 0.0,   0.0,   0.0,   0.0,   0.0, -0.25, -0.25,  0.25,  0.25}
    };
    memcpy(Minv, mi, sizeof(mi));

    for (int k = 0; k < 9; ++k) Smat[k] = 1.0;
    Smat[4] = OMEGA;
    Smat[5] = OMEGA;
}

static void init_geometry(void) {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int s = (x >= OBS_X0 && x < OBS_X1 && y >= OBS_Y0 && y < OBS_Y1) ? 1 : 0;
            solid[IDX(x, y)] = (char)s;
        }
    }
}

static inline double feq_dir(int d, double rho_, double ux, double uy) {
    double eu = cx[d]*ux + cy[d]*uy;
    double u2 = ux*ux + uy*uy;
    return w[d] * rho_ * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*u2);
}

static void initialize(void) {
    init_cm_matrices();
    init_geometry();
    // Small anti-symmetric kick localised in the wake to break y-symmetry
    // quickly (same trick as karman.c).
    double cy_wake = 0.5*(OBS_Y0 + OBS_Y1);
    double cx_wake = OBS_X1 + 16.0;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double dx_ = x - cx_wake;
            double dy_ = y - cy_wake;
            double r2 = dx_*dx_ + dy_*dy_;
            double pert = 0.001 * (y > cy_wake ? 1.0 : -1.0) * exp(-r2/200.0);
            u[i]   = 0.0;
            v[i]   = solid[i] ? 0.0 : pert;
            rho[i] = 1.0;
            for (int d = 0; d < NDIR; ++d) {
                f[FIDX(i, d)] = feq_dir(d, 1.0, u[i], v[i]);
            }
        }
    }
}

static void macroscopic(void) {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) { u[i] = 0.0; v[i] = 0.0; rho[i] = 1.0; continue; }
            double rr = 0, ru = 0, rv = 0;
            for (int d = 0; d < 9; ++d) {
                double ff = f[FIDX(i, d)];
                rr += ff; ru += ff*cx[d]; rv += ff*cy[d];
            }
            rho[i] = rr;
            // Guo half-force correction to velocity
            u[i] = (ru + 0.5*FORCE_X) / rr;
            v[i] = rv / rr;
        }
    }
}

static void cm_collide_cell(int i) {
    double fl[9], fe[9];
    double ux = u[i], uy = v[i], rho_ = rho[i];

    for (int d = 0; d < 9; ++d) {
        fl[d] = f[FIDX(i, d)];
        fe[d] = feq_dir(d, rho_, ux, uy);
    }

    double m[9], me[9];
    for (int k = 0; k < 9; ++k) {
        double sm = 0, se = 0;
        for (int d = 0; d < 9; ++d) {
            sm += Mmat[k][d]*fl[d];
            se += Mmat[k][d]*fe[d];
        }
        m[k] = sm; me[k] = se;
    }

    // Shift matrix N(u) and its inverse N^-1 — definitions match lbmcm.c.
    double N[9][9]  = {{0}};
    double Ni[9][9] = {{0}};

    N[0][0] = 1.0;
    N[1][0] = -ux;            N[1][1] = 1.0;
    N[2][0] = -uy;            N[2][2] = 1.0;
    N[3][0] = ux*ux + uy*uy;  N[3][1] = -2*ux; N[3][2] = -2*uy; N[3][3] = 1.0;
    N[4][0] = ux*ux - uy*uy;  N[4][1] = -2*ux; N[4][2] =  2*uy; N[4][4] = 1.0;
    N[5][0] = ux*uy;          N[5][1] = -uy;   N[5][2] = -ux;   N[5][5] = 1.0;

    N[6][0] = -ux*ux*uy;
    N[6][1] =  2*ux*uy;
    N[6][2] =  ux*ux;
    N[6][3] = -uy*0.5;
    N[6][4] = -uy*0.5;
    N[6][5] = -2*ux;
    N[6][6] = 1.0;

    N[7][0] = -uy*uy*ux;
    N[7][1] =  uy*uy;
    N[7][2] =  2*ux*uy;
    N[7][3] = -ux*0.5;
    N[7][4] =  ux*0.5;
    N[7][5] = -2*uy;
    N[7][7] = 1.0;

    N[8][0] =  ux*ux*uy*uy;
    N[8][1] = -2*ux*uy*uy;
    N[8][2] = -2*ux*ux*uy;
    N[8][3] =  uy*uy*0.5 + ux*ux*0.5;
    N[8][4] =  uy*uy*0.5 - ux*ux*0.5;
    N[8][5] =  4*ux*uy;
    N[8][6] = -2*uy;
    N[8][7] = -2*ux;
    N[8][8] = 1.0;

    Ni[0][0] = 1.0;
    Ni[1][0] =  ux;             Ni[1][1] = 1.0;
    Ni[2][0] =  uy;             Ni[2][2] = 1.0;
    Ni[3][0] =  ux*ux + uy*uy;  Ni[3][1] =  2*ux; Ni[3][2] =  2*uy; Ni[3][3] = 1.0;
    Ni[4][0] =  ux*ux - uy*uy;  Ni[4][1] =  2*ux; Ni[4][2] = -2*uy; Ni[4][4] = 1.0;
    Ni[5][0] =  ux*uy;          Ni[5][1] =  uy;   Ni[5][2] =  ux;   Ni[5][5] = 1.0;

    Ni[6][0] =  ux*ux*uy;
    Ni[6][1] =  2*ux*uy;
    Ni[6][2] =  ux*ux;
    Ni[6][3] =  uy*0.5;
    Ni[6][4] =  uy*0.5;
    Ni[6][5] =  2*ux;
    Ni[6][6] = 1.0;

    Ni[7][0] =  uy*uy*ux;
    Ni[7][1] =  uy*uy;
    Ni[7][2] =  2*ux*uy;
    Ni[7][3] =  ux*0.5;
    Ni[7][4] = -ux*0.5;
    Ni[7][5] =  2*uy;
    Ni[7][7] = 1.0;

    Ni[8][0] =  ux*ux*uy*uy;
    Ni[8][1] =  2*ux*uy*uy;
    Ni[8][2] =  2*ux*ux*uy;
    Ni[8][3] =  uy*uy*0.5 + ux*ux*0.5;
    Ni[8][4] =  uy*uy*0.5 - ux*ux*0.5;
    Ni[8][5] =  4*ux*uy;
    Ni[8][6] =  2*uy;
    Ni[8][7] =  2*ux;
    Ni[8][8] = 1.0;

    // tn = N*m, tne = N*me  (raw moments -> central moments)
    double tn[9], tne[9];
    for (int k = 0; k < 9; ++k) {
        double sa = 0, sb = 0;
        for (int j = 0; j < 9; ++j) {
            sa += N[k][j]*m[j];
            sb += N[k][j]*me[j];
        }
        tn[k]  = sa;
        tne[k] = sb;
    }

    for (int k = 0; k < 9; ++k) {
        tn[k] = tn[k] - Smat[k]*(tn[k] - tne[k]);
    }

    // back to raw moments
    double rm[9];
    for (int k = 0; k < 9; ++k) {
        double s = 0;
        for (int j = 0; j < 9; ++j) s += Ni[k][j]*tn[j];
        rm[k] = s;
    }

    // back to distributions
    for (int k = 0; k < 9; ++k) {
        double s = 0;
        for (int j = 0; j < 9; ++j) s += Minv[k][j]*rm[j];
        f[FIDX(i, k)] = s;
    }
}

static inline double guo_force_term(int d, double ux, double uy) {
    // F_i = (1 - 1/(2 tau)) * w_i * [3(c_i - u) + 9(c_i . u) c_i] . g,
    // here g = (FORCE_X, 0).
    double eu = cx[d]*ux + cy[d]*uy;
    double gx = (cx[d] - ux)*3.0 + 9.0*eu*cx[d];
    return (1.0 - 0.5*OMEGA) * w[d] * FORCE_X * gx;
}

static void collide(void) {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) continue;
            cm_collide_cell(i);
            for (int d = 0; d < 9; ++d) {
                f[FIDX(i, d)] += guo_force_term(d, u[i], v[i]);
            }
        }
    }
}

static void stream(void) {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) {
                for (int d = 0; d < 9; ++d) f2[FIDX(i, d)] = f[FIDX(i, d)];
                continue;
            }
            for (int d = 0; d < 9; ++d) {
                int xp = (x + cx[d] + NX) % NX;          // periodic in x
                int yp = y + cy[d];
                double post = f[FIDX(i, d)];
                if (yp < 0 || yp >= NY) {
                    // halfway bounce-back at top/bottom walls
                    f2[FIDX(i, opp[d])] = post;
                } else if (solid[IDX(xp, yp)]) {
                    // halfway bounce-back at solid obstacle
                    f2[FIDX(i, opp[d])] = post;
                } else {
                    f2[FIDX(IDX(xp, yp), d)] = post;
                }
            }
        }
    }
    double *tmp = f; f = f2; f2 = tmp;
}

static int output_snapshot(int step_no) {
    char fname[64];
    snprintf(fname, sizeof(fname), "lbmcm_karman_snapshot_%05d.csv", step_no);
    FILE* fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", fname); return 1; }
    fprintf(fp, "x,y,u,v,speed,solid\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            double sp = sqrt(u[i]*u[i] + v[i]*v[i]);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%d\n",
                    x, y, u[i], v[i], sp, solid[i]);
        }
    }
    fclose(fp);
    return 0;
}

static int output_wake_frame(int step_no) {
    char fname[64];
    snprintf(fname, sizeof(fname), "lbmcm_karman_wake_%05d.csv", step_no);
    FILE* fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", fname); return 1; }
    fprintf(fp, "x,y,u,v\n");
    for (int y = SPEC_Y0; y < SPEC_Y1; ++y) {
        for (int x = SPEC_X0; x < SPEC_X1; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g\n", x, y, u[i], v[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main(void) {
    initialize();

    int snap_steps[SNAPSHOTS];
    for (int s = 0; s < SNAPSHOTS; ++s) {
        if (s == 0) snap_steps[s] = 0;
        else snap_steps[s] = (int)(NSTEPS * pow((double)s / (SNAPSHOTS - 1), 1.5));
    }

    FILE* hist = fopen("lbmcm_karman_probe.csv", "w");
    if (!hist) { fprintf(stderr, "cannot open probe csv\n"); return 1; }
    fprintf(hist, "step,u_max,u_probe,v_probe\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();

        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t >= WAKE_FIRST_STEP && (t - WAKE_FIRST_STEP) % WAKE_STRIDE == 0) {
            if (output_wake_frame(t) != 0) return 1;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (!solid[i] && fabs(u[i]) > umax) umax = fabs(u[i]);
            }
            int p = IDX(PROBE_X, PROBE_Y);
            fprintf(hist, "%d,%.9g,%.9g,%.9g\n", t, umax, u[p], v[p]);
        }

        collide();
        stream();
    }

    macroscopic();
    output_snapshot(NSTEPS);
    fclose(hist);

    double umax = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (!solid[i] && fabs(u[i]) > umax) umax = fabs(u[i]);
    }
    double D = (double)(OBS_Y1 - OBS_Y0);
    double Re_D = umax * D / NU0;
    printf("Done. snapshots/wake-frames/probe written.\n");
    printf("NX=%d NY=%d obstacle=[%d,%d]x[%d,%d] D=%.0f NSTEPS=%d TAU=%.4f F=%g nu=%.5f umax=%.4f Re_D=%.0f\n",
           NX, NY, OBS_X0, OBS_X1, OBS_Y0, OBS_Y1, D, NSTEPS, TAU, FORCE_X, NU0, umax, Re_D);
    return 0;
}
