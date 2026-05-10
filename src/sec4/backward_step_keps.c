// backward_step_keps.c
// 2D backward-facing step + standard k-epsilon model overlay (D2Q9, BGK)
//
// Same channel + step geometry, body force, and periodic-x boundary as
// backward_step.c. The k-eps transport equations supply nu_t = Cmu*k^2/eps,
// fed back into a locally varying tau_eff in the BGK collision.
//
// Wall functions (Dirichlet) are applied to:
//   - top wall (y = NY-1)
//   - bottom wall (y = 0), but only where there is fluid (i.e., x >= STEP_LENGTH)
//   - top of step (y = STEP_HEIGHT, 0 <= x < STEP_LENGTH)
//   - downstream face of step (x = STEP_LENGTH, 0 <= y < STEP_HEIGHT)
//
// The local-shear u_tau formula (sqrt(nu_0 * |du_t/dn|)) is reused; corner
// cells take the larger of two-wall estimates as in cavity_keps.c.
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
#define FORCE_X 2e-6
#define STEP_LENGTH 30
#define STEP_HEIGHT 30
#define KAPPA 0.41
#define WALL_DY 0.5
#define KEPS_DT 0.05

const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
const int opp[NDIR] = {0,3,4,1,2,7,8,5,6};

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
static char solid[NX*NY];
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
    double k_seed = 0.005 * 0.05 * 0.05;     // typical u^2 scale
    double nut_target = 0.05 * nu0;
    double eps_seed = Cmu * k_seed * k_seed / nut_target;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            k[i] = k_seed;
            eps[i] = eps_seed;
            nut_field[i] = nut_target;
            for (int d = 0; d < NDIR; ++d) {
                f[i*NDIR + d] = w[d] * rho[i];
            }
        }
    }
}

void stream_collide_with_keps() {
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

static void apply_wall_function(int i, double tangent_jump) {
    if (solid[i]) return;
    double shear = fabs(tangent_jump) / WALL_DY;
    double u_tau = sqrt(nu0 * shear + 1e-12);
    double k_wall = u_tau * u_tau / sqrt(Cmu);
    double eps_wall = u_tau * u_tau * u_tau / (KAPPA * WALL_DY);
    if (k_new[i] < k_wall) k_new[i] = k_wall;
    if (eps_new[i] < eps_wall) eps_new[i] = eps_wall;
}

void update_kepsilon() {
    double dx = 1.0, dy = 1.0, dt = KEPS_DT;
    // Bulk transport
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) {
                k_new[i] = 1e-10;
                eps_new[i] = 1e-12;
                continue;
            }
            // Mirror boundary handling: if neighbor is a wall (out of domain) or solid,
            // reflect to current cell (Neumann zero-gradient).
            int ixp = IDX((x+1)%NX, y);
            int ixm = IDX((x-1+NX)%NX, y);
            if (solid[ixp]) ixp = i;
            if (solid[ixm]) ixm = i;
            int iyp_idx, iym_idx;
            if (y+1 < NY && !solid[IDX(x, y+1)]) iyp_idx = IDX(x, y+1);
            else                                  iyp_idx = i;
            if (y-1 >= 0 && !solid[IDX(x, y-1)])  iym_idx = IDX(x, y-1);
            else                                  iym_idx = i;
            double dudx = (u[ixp] - u[ixm])/(2.0*dx);
            double dvdx = (v[ixp] - v[ixm])/(2.0*dx);
            double dudy = (u[iyp_idx] - u[iym_idx])/(2.0*dy);
            double dvdy = (v[iyp_idx] - v[iym_idx])/(2.0*dy);
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double nut = Cmu * k[i]*k[i] / (eps[i] + 1e-12);
            double Pk = nut * S2;
            double Dk = nu0 + nut/sig_k;
            double De = nu0 + nut/sig_e;
            double lap_k = (k[ixp] + k[ixm] + k[iyp_idx] + k[iym_idx] - 4.0*k[i])/(dx*dx);
            double lap_e = (eps[ixp] + eps[ixm] + eps[iyp_idx] + eps[iym_idx] - 4.0*eps[i])/(dx*dx);
            double uk = u[i]*((u[i]>0 ? (k[i] - k[ixm]) : (k[ixp] - k[i]))/dx);
            double vk = v[i]*((v[i]>0 ? (k[i] - k[iym_idx]) : (k[iyp_idx] - k[i]))/dy);
            double ue = u[i]*((u[i]>0 ? (eps[i] - eps[ixm]) : (eps[ixp] - eps[i]))/dx);
            double ve = v[i]*((v[i]>0 ? (eps[i] - eps[iym_idx]) : (eps[iyp_idx] - eps[i]))/dy);
            k_new[i] = k[i] + dt * (Dk*lap_k - uk - vk + Pk - eps[i]);
            eps_new[i] = eps[i] + dt * (De*lap_e - ue - ve
                            + Ce1*Pk*eps[i]/(k[i]+1e-12)
                            - Ce2*eps[i]*eps[i]/(k[i]+1e-12));
        }
    }

    // Wall functions
    // Top wall y=NY-1
    for (int x = 0; x < NX; ++x) {
        int i = IDX(x, NY-1);
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, u[i]);
    }
    // Bottom wall y=0 (only where fluid: x >= STEP_LENGTH)
    for (int x = STEP_LENGTH; x < NX; ++x) {
        int i = IDX(x, 0);
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, u[i]);
    }
    // Top of step: y = STEP_HEIGHT, 0 <= x < STEP_LENGTH
    for (int x = 0; x < STEP_LENGTH; ++x) {
        int i = IDX(x, STEP_HEIGHT);
        if (solid[i]) continue;     // shouldn't happen but safe
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, u[i]);
    }
    // Downstream face of step: x = STEP_LENGTH, 0 <= y < STEP_HEIGHT
    for (int y = 0; y < STEP_HEIGHT; ++y) {
        int i = IDX(STEP_LENGTH, y);
        if (solid[i]) continue;
        k_new[i] = 0; eps_new[i] = 0;
        apply_wall_function(i, v[i]);
    }

    // Apply non-negativity floor and commit
    for (int i = 0; i < NX*NY; ++i) {
        if (k_new[i] < 1e-10) k_new[i] = 1e-10;
        if (eps_new[i] < 1e-12) eps_new[i] = 1e-12;
        k[i] = k_new[i];
        eps[i] = eps_new[i];
        nut_field[i] = solid[i] ? 0.0 : Cmu * k[i]*k[i] / (eps[i] + 1e-12);
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
    snprintf(fname, sizeof(fname), "step_keps_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    compute_streamfunction();
    fprintf(fp, "x,y,u,v,vorticity,psi,solid,k,eps,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], psi[i], solid[i],
                    k[i], eps[i], nut_field[i]);
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

    FILE* hist = fopen("step_keps_history.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open step_keps_history.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,reattach_x,k_mean,eps_mean,nut_mean\n");

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, k_sum = 0, eps_sum = 0, nut_sum = 0;
            int n_fluid = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (solid[i]) continue;
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                k_sum += k[i];
                eps_sum += eps[i];
                nut_sum += nut_field[i];
                ++n_fluid;
            }
            int reattach = -1;
            for (int xq = STEP_LENGTH; xq < NX; ++xq) {
                if (u[IDX(xq, 1)] > 0) { reattach = xq; break; }
            }
            double inv_n = 1.0 / n_fluid;
            fprintf(hist, "%d,%.9g,%d,%.9g,%.9g,%.9g\n",
                    t, umax, reattach, k_sum*inv_n, eps_sum*inv_n, nut_sum*inv_n);
        }
        update_kepsilon();
        stream_collide_with_keps();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    printf("Done. Snapshots: step_keps_snapshot_*.csv, history: step_keps_history.csv\n");
    printf("Parameters: NX=%d NY=%d STEP=%dx%d NSTEPS=%d TAU=%.3f F=%g nu0=%.5f\n",
           NX, NY, STEP_LENGTH, STEP_HEIGHT, NSTEPS, TAU, FORCE_X, nu0);
    return 0;
}
