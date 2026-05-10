// kelvin_helmholtz_keps.c
// 2D Kelvin-Helmholtz instability + standard k-epsilon model overlay (D2Q9, BGK)
//
// Same setup and initial conditions as kelvin_helmholtz.c. The k-eps transport
// equations supply nu_t = Cmu*k^2/eps, fed back into a locally varying
// tau_eff = 0.5 + 3*(nu_0 + nu_t). Because the shear layer has high |du/dy|,
// the production term Pk = nu_t * S^2 is expected to grow k near each shear
// layer; the resulting eddy viscosity then *suppresses* the linear roll-up
// compared to pure DNS.
//
// Physical caveat: standard k-eps is a RANS closure. It can over-predict
// nu_t in transition regions, smearing out the very vortices we want to see.
// This is exactly the point of the comparison — the side-by-side with the
// pure-LBM run shows when k-eps is too dissipative for instability physics.
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
// k-eps integration time step (smooth periodic fields; central convection)
#define KEPS_DT 0.05

const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

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
    // Seed k, eps so that nu_t/nu_0 ~ 0.05 initially — small enough not to
    // suppress the linear instability immediately, but enough that k-eps source
    // terms are well-posed.
    double k_seed = 0.002 * U0 * U0;
    double nut_target = 0.05 * nu0;
    double eps_seed = Cmu * k_seed * k_seed / nut_target;

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
            k[i] = k_seed;
            eps[i] = eps_seed;
            nut_field[i] = nut_target;
            double usqr = u[i]*u[i] + v[i]*v[i];
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
            }
        }
    }
}

void stream_collide_with_keps() {
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

void update_kepsilon() {
    // Periodic everywhere — no wall functions
    double dx = 1.0, dy = 1.0, dt = KEPS_DT;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            int ixp = IDX((x+1)%NX, y), ixm = IDX((x-1+NX)%NX, y);
            int iyp = IDX(x, (y+1)%NY), iym = IDX(x, (y-1+NY)%NY);
            double dudx = (u[ixp] - u[ixm])/(2.0*dx);
            double dudy = (u[iyp] - u[iym])/(2.0*dy);
            double dvdx = (v[ixp] - v[ixm])/(2.0*dx);
            double dvdy = (v[iyp] - v[iym])/(2.0*dy);
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double nut = Cmu * k[i]*k[i] / (eps[i] + 1e-12);
            double Pk = nut * S2;
            double Dk = nu0 + nut/sig_k;
            double De = nu0 + nut/sig_e;
            double lap_k = (k[ixp] + k[ixm] + k[iyp] + k[iym] - 4.0*k[i])/(dx*dx);
            double lap_e = (eps[ixp] + eps[ixm] + eps[iyp] + eps[iym] - 4.0*eps[i])/(dx*dx);
            double uk = u[i] * (k[ixp] - k[ixm]) / (2.0*dx)
                      + v[i] * (k[iyp] - k[iym]) / (2.0*dy);
            double ue = u[i] * (eps[ixp] - eps[ixm]) / (2.0*dx)
                      + v[i] * (eps[iyp] - eps[iym]) / (2.0*dy);
            k_new[i] = k[i] + dt * (Dk*lap_k - uk + Pk - eps[i]);
            eps_new[i] = eps[i] + dt * (De*lap_e - ue
                            + Ce1*Pk*eps[i]/(k[i]+1e-12)
                            - Ce2*eps[i]*eps[i]/(k[i]+1e-12));
        }
    }
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
    snprintf(fname, sizeof(fname), "kh_keps_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,k,eps,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], k[i], eps[i], nut_field[i]);
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

    FILE* hist = fopen("kh_keps_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open kh_keps_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,v_max,u_rms,k_mean,eps_mean,nut_mean\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double vmax = 0, u_sq = 0, k_sum = 0, eps_sum = 0, nut_sum = 0;
            for (int i = 0; i < NX*NY; ++i) {
                double vm = fabs(v[i]);
                if (vm > vmax) vmax = vm;
                u_sq += u[i]*u[i];
                k_sum += k[i];
                eps_sum += eps[i];
                nut_sum += nut_field[i];
            }
            double inv_n = 1.0 / (NX*NY);
            double u_rms = sqrt(u_sq * inv_n);
            fprintf(hist, "%d,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    t, vmax, u_rms, k_sum*inv_n, eps_sum*inv_n, nut_sum*inv_n);
        }
        macroscopic();
        update_kepsilon();
        stream_collide_with_keps();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    printf("Done. Snapshots: kh_keps_snapshot_*.csv, history: kh_keps_history.csv\n");
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f U0=%.3f delta=%.1f nu0=%.5f\n",
           NX, NY, NSTEPS, TAU, U0, SHEAR_DELTA, nu0);
    return 0;
}
