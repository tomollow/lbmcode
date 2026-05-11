// kelvin_helmholtz_les.c
// 2D Kelvin-Helmholtz instability + standard Smagorinsky LES (D2Q9, BGK)
//
// Same setup and initial conditions as kelvin_helmholtz.c. Smagorinsky SGS:
// nu_t = (Cs*Delta)^2 * |S|, Cs = 0.16, Delta = 1 (grid spacing in LU).
// nu_t feeds back via tau_eff = 0.5 + 3*(nu_0 + nu_t) in BGK.
//
// Smagorinsky is a local closure with no transport, so unlike k-eps it
// cannot build up nu_t in the shear layer over time — it just responds
// to instantaneous |S|. Expected outcome: much weaker damping of the
// roll-up than k-eps, closer to pure-LBM behavior.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 256
#define NY 128
#define NDIR 9
#define NSTEPS 8000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 25
#define TAU 0.52
#define OMEGA (1.0/TAU)
#define U0 0.05
#define SHEAR_DELTA 4.0
#define PERT_AMP 0.002
#define PERT_SIGMA 8.0
#define PERT_MODE 4
#define M_PI_F 3.14159265358979323846
#define CS_SMAG 0.16
#define DELTA_LES 1.0

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
static double nut_field[NX*NY];
static double vort[NX*NY];

static double shear_profile(double y) {
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
            double pert = PERT_AMP * sx * (exp(-dy1*dy1) - exp(-dy2*dy2));
            u[i] = shear_profile((double)y);
            v[i] = pert;
            rho[i] = 1.0;
            nut_field[i] = 0.0;
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
            }
        }
    }
}

void stream_collide_with_les() {
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

void update_les() {
    double dx = 1.0, dy = 1.0;
    double Cs2 = (CS_SMAG * DELTA_LES) * (CS_SMAG * DELTA_LES);
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            int ixp = IDX((x+1)%NX, y), ixm = IDX((x-1+NX)%NX, y);
            int iyp = IDX(x, (y+1)%NY), iym = IDX(x, (y-1+NY)%NY);
            double dudx = (u[ixp] - u[ixm]) / (2.0*dx);
            double dudy = (u[iyp] - u[iym]) / (2.0*dy);
            double dvdx = (v[ixp] - v[ixm]) / (2.0*dx);
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
    snprintf(fname, sizeof(fname), "kh_les_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], nut_field[i]);
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

    FILE* hist = fopen("kh_les_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open kh_les_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,v_max,u_rms,nut_mean\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double vmax = 0, u_sq = 0, nut_sum = 0;
            for (int i = 0; i < NX*NY; ++i) {
                double vm = fabs(v[i]);
                if (vm > vmax) vmax = vm;
                u_sq += u[i]*u[i];
                nut_sum += nut_field[i];
            }
            double inv_n = 1.0 / (NX*NY);
            double u_rms = sqrt(u_sq * inv_n);
            fprintf(hist, "%d,%.9g,%.9g,%.9g\n",
                    t, vmax, u_rms, nut_sum*inv_n);
        }
        macroscopic();
        update_les();
        stream_collide_with_les();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    printf("Done. Snapshots: kh_les_snapshot_*.csv, history: kh_les_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U0=%.3f delta=%.1f Cs=%.3f nu0=%.5f\n",
           NX, NY, NSTEPS, TAU, U0, SHEAR_DELTA, CS_SMAG, nu0);
    return 0;
}
