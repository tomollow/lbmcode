// karman_les.c
// 2D Karman vortex street + standard Smagorinsky LES (D2Q9, BGK)
//
// Same geometry, body force, and asymmetric initial v perturbation as karman.c.
// Smagorinsky SGS: nu_t = (Cs*Delta)^2 * |S|, Cs = 0.16, Delta = 1 (grid
// spacing in LU). nu_t feeds back to BGK via tau_eff = 0.5 + 3*(nu_0 + nu_t).
//
// No wall damping. The cylinder is a staircase obstacle; |S| is computed
// with reflected (zero-gradient) neighbors when adjacent to solid cells.
// At Re_D ~ 130, the wake is laminar shedding and |S| is moderate only in
// the shear layers, so nu_t/nu_0 stays small and shedding survives.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 360
#define NY 80
#define NDIR 9
#define NSTEPS 50000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 5
#define TAU 0.55
#define OMEGA (1.0/TAU)
#define FORCE_X 6e-6
#define CX 80
#define CY 41
#define R_CYL 10
#define PROBE_X 200
#define PROBE_Y 50
#define CS_SMAG 0.16
#define DELTA_LES 1.0

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
static double nut_field[NX*NY];
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
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            double dx_ = x - (CX + 25);
            double dy_ = y - CY;
            double r2 = dx_*dx_ + dy_*dy_;
            double pert = 0.0005 * (y > CY ? 1.0 : -1.0) * exp(-r2 / 200.0);
            v[i] = solid[i] ? 0.0 : pert;
            rho[i] = 1.0;
            nut_field[i] = 0.0;
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu);
            }
        }
    }
}

void stream_collide_with_les() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) continue;
            double tau_eff = 0.5 + 3.0 * (nu0 + nut_field[i]);
            double omega_eff = 1.0 / tau_eff;
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                double Fi = (1.0 - 0.5*omega_eff) * w[d] * FORCE_X *
                            (3.0*(cx[d] - u[i]) + 9.0*eu*cx[d]);
                double post = f[i*NDIR + d] - omega_eff*(f[i*NDIR + d] - feq) + Fi;
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

void update_les() {
    double dx = 1.0, dy = 1.0;
    double Cs2 = (CS_SMAG * DELTA_LES) * (CS_SMAG * DELTA_LES);
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) {
                nut_field[i] = 0.0;
                continue;
            }
            int ixp = IDX((x+1)%NX, y);
            int ixm = IDX((x-1+NX)%NX, y);
            if (solid[ixp]) ixp = i;
            if (solid[ixm]) ixm = i;
            int iyp = (y+1 < NY && !solid[IDX(x, y+1)]) ? IDX(x, y+1) : i;
            int iym = (y-1 >= 0 && !solid[IDX(x, y-1)]) ? IDX(x, y-1) : i;
            double dudx = (u[ixp] - u[ixm]) / (2.0*dx);
            double dvdx = (v[ixp] - v[ixm]) / (2.0*dx);
            double dudy = (u[iyp] - u[iym]) / (2.0*dy);
            double dvdy = (v[iyp] - v[iym]) / (2.0*dy);
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double Smag = sqrt(S2);
            nut_field[i] = Cs2 * Smag;
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
    snprintf(fname, sizeof(fname), "karman_les_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,solid,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%d,%.9g\n",
                    x, y, u[i], v[i], vort[i], solid[i], nut_field[i]);
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

    FILE* hist = fopen("karman_les_probe.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open karman_les_probe.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,u_probe,v_probe,nut_mean\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, nut_sum = 0;
            int n_fluid = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (solid[i]) continue;
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                nut_sum += nut_field[i];
                ++n_fluid;
            }
            int p = IDX(PROBE_X, PROBE_Y);
            double inv_n = 1.0 / n_fluid;
            fprintf(hist, "%d,%.9g,%.9g,%.9g,%.9g\n",
                    t, umax, u[p], v[p], nut_sum*inv_n);
        }
        update_les();
        stream_collide_with_les();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double umax_actual = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (!solid[i] && fabs(u[i]) > umax_actual) umax_actual = fabs(u[i]);
    }
    double Re_D = umax_actual * 2 * R_CYL / nu0;
    printf("Done. Snapshots: karman_les_snapshot_*.csv, probe: karman_les_probe.csv\n");
    printf("Parameters: NX=%d NY=%d D=%d Cyl=(%d,%d) NSTEPS=%d TAU=%.3f F=%g Cs=%.3f nu0=%.5f u_max=%.4f Re_D=%.0f\n",
           NX, NY, 2*R_CYL, CX, CY, NSTEPS, TAU, FORCE_X, CS_SMAG, nu0, umax_actual, Re_D);
    return 0;
}
