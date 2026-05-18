// karman_des_hires.c
// 2D Karman vortex street + Spalart-Allmaras DES97 (D2Q9, BGK)
//
// Same geometry, body force, and asymmetric initial v perturbation as
// karman.c / karman_les.c. SA-DES hybrid: 1-equation transport for nu_tilde
// with length-scale switch
//
//     d_tilde = min(d_wall, C_DES * Delta),   C_DES = 0.65, Delta = 1 LU
//
// Wall distance combines the top/bottom channel walls (halfway BB) and the
// staircase cylinder, treated as an idealized circle of radius R_CYL:
//
//     d_wall = min( y+0.5, NY-y-0.5, sqrt((x-CX)^2+(y-CY)^2) - R_CYL )
//
// The first off-wall layer (one cell away from any wall or cylinder surface)
// has d_wall < C_DES*Delta and runs in SA-RANS mode; the rest of the fluid
// (~99%) runs in the LES branch. At Re_D ~ 130 the wake shedding is laminar
// and |S| in the wake shear layer is moderate, so SA production cannot
// overcome destruction 窶・the model decays toward chi << c_v1 and nu_t goes
// to essentially zero, exactly as in cavity_des. The educational value is in
// watching the length-scale switch operate and confirming SA-DES "correctly
// sleeps" in a 2D laminar shedding regime.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 360
#define NY 80
#define NDIR 9
#define NSTEPS 80000
#define SNAPSHOTS 6
#define HISTORY_INTERVAL 5
#define TAU 0.51
#define OMEGA (1.0/TAU)
#define FORCE_X 1.2e-6
#define CX 80
#define CY 41
#define R_CYL 10
#define PROBE_X 200
#define PROBE_Y 50

// Spalart-Allmaras constants
#define KAPPA 0.41
#define C_B1  0.1355
#define C_B2  0.622
#define SIG_SA (2.0/3.0)
#define C_V1  7.1
#define C_W1  (C_B1/(KAPPA*KAPPA) + (1.0 + C_B2)/SIG_SA)
#define C_W2  0.3
#define C_W3  2.0
// DES97 length-scale switch
#define C_DES 0.65
#define DELTA_LES 1.0
// SA transport time step (explicit Euler)
#define SA_DT 0.05

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
static double nu_tilde[NX*NY], nu_tilde_new[NX*NY];
static double nut_field[NX*NY];
static double d_wall[NX*NY], d_tilde[NX*NY];
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

void init_wall_distance() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) { d_wall[i] = 0.0; d_tilde[i] = 0.0; continue; }
            double dtop = NY - y - 0.5;
            double dbot = y + 0.5;
            double dx_ = x - CX, dy_ = y - CY;
            double dcyl = sqrt(dx_*dx_ + dy_*dy_) - (double)R_CYL;
            if (dcyl < 0.5) dcyl = 0.5;  // floor at the first off-wall layer
            double d = dtop < dbot ? dtop : dbot;
            if (dcyl < d) d = dcyl;
            d_wall[i] = d;
            d_tilde[i] = (d < C_DES*DELTA_LES) ? d : C_DES*DELTA_LES;
        }
    }
}

void initialize() {
    init_geometry();
    init_wall_distance();
    double nut_seed = 3.0 * nu0;
    double chi_seed = nut_seed / nu0;
    double chi3_seed = chi_seed * chi_seed * chi_seed;
    double fv1_seed = chi3_seed / (chi3_seed + C_V1*C_V1*C_V1);
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
            nu_tilde[i] = solid[i] ? 0.0 : nut_seed;
            nut_field[i] = solid[i] ? 0.0 : nut_seed * fv1_seed;
            for (int d = 0; d < NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                f[i*NDIR + d] = w[d] * rho[i] * (1.0 + 3.0*eu);
            }
        }
    }
}

void stream_collide_with_des() {
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

void update_sa_des() {
    double dx = 1.0, dy = 1.0, dt = SA_DT;
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            if (solid[i]) { nut_field[i] = 0.0; continue; }
            // Periodic x, mirror at top/bottom walls and at solid neighbors
            int ixp = IDX((x+1)%NX, y);
            int ixm = IDX((x-1+NX)%NX, y);
            if (solid[ixp]) ixp = i;
            if (solid[ixm]) ixm = i;
            int iyp = (y+1 < NY && !solid[IDX(x, y+1)]) ? IDX(x, y+1) : i;
            int iym = (y-1 >= 0 && !solid[IDX(x, y-1)]) ? IDX(x, y-1) : i;
            // Velocity gradients
            double dudx = 0.5 * (u[ixp] - u[ixm]) / dx;
            double dvdx = 0.5 * (v[ixp] - v[ixm]) / dx;
            double dudy = 0.5 * (u[iyp] - u[iym]) / dy;
            double dvdy = 0.5 * (v[iyp] - v[iym]) / dy;
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double Smag = sqrt(S2);

            // SA closure
            double nt = nu_tilde[i];
            double chi = nt / nu0;
            double chi3 = chi*chi*chi;
            double fv1 = chi3 / (chi3 + C_V1*C_V1*C_V1);
            double fv2 = 1.0 - chi / (1.0 + chi*fv1);
            double dl = d_tilde[i];
            double inv_kdt2 = 1.0 / (KAPPA*KAPPA * dl*dl);
            double Stilde = Smag + nt * fv2 * inv_kdt2;
            if (Stilde < 1e-12) Stilde = 1e-12;
            double r = nt * inv_kdt2 / Stilde;
            if (r > 10.0) r = 10.0;
            double r6 = r*r*r*r*r*r;
            double g = r + C_W2*(r6 - r);
            double g6 = g*g*g*g*g*g;
            double cw36 = C_W3*C_W3*C_W3*C_W3*C_W3*C_W3;
            double fw = g * pow((1.0 + cw36) / (g6 + cw36), 1.0/6.0);

            double prod = C_B1 * Stilde * nt;
            double dest = C_W1 * fw * (nt/dl) * (nt/dl);
            double lap = (nu_tilde[ixp] + nu_tilde[ixm] + nu_tilde[iyp] + nu_tilde[iym]
                          - 4.0*nt) / (dx*dx);
            double dntdx = 0.5 * (nu_tilde[ixp] - nu_tilde[ixm]) / dx;
            double dntdy = 0.5 * (nu_tilde[iyp] - nu_tilde[iym]) / dy;
            double grad2 = dntdx*dntdx + dntdy*dntdy;
            double diff = ((nu0 + nt)*lap + (1.0 + C_B2)*grad2) / SIG_SA;

            double conv_x = u[i] * ((u[i] > 0) ? (nt - nu_tilde[ixm])/dx
                                              : (nu_tilde[ixp] - nt)/dx);
            double conv_y = v[i] * ((v[i] > 0) ? (nt - nu_tilde[iym])/dy
                                              : (nu_tilde[iyp] - nt)/dy);
            nu_tilde_new[i] = nt + dt * (prod - dest + diff - conv_x - conv_y);
        }
    }
    for (int i = 0; i < NX*NY; ++i) {
        if (solid[i]) continue;
        if (nu_tilde_new[i] < 1e-12) nu_tilde_new[i] = 1e-12;
        nu_tilde[i] = nu_tilde_new[i];
        double chi = nu_tilde[i] / nu0;
        double chi3 = chi*chi*chi;
        double fv1 = chi3 / (chi3 + C_V1*C_V1*C_V1);
        nut_field[i] = nu_tilde[i] * fv1;
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
    snprintf(fname, sizeof(fname), "karman_des_hires_snapshot_%05d.csv", step);
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_snapshot: cannot open %s\n", fname);
        return 1;
    }
    compute_vorticity();
    fprintf(fp, "x,y,u,v,vorticity,solid,nut,nu_tilde,d_wall,d_tilde\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%d,%.9g,%.9g,%.9g,%.9g\n",
                    x, y, u[i], v[i], vort[i], solid[i],
                    nut_field[i], nu_tilde[i], d_wall[i], d_tilde[i]);
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

    FILE* hist = fopen("karman_des_hires_probe.csv", "w");
    if (!hist) {
        fprintf(stderr, "main: cannot open karman_des_hires_probe.csv\n");
        return 1;
    }
    fprintf(hist, "step,u_max,u_probe,v_probe,nut_mean,nu_tilde_mean,les_frac\n");

    // LES-branch fraction over fluid cells
    int n_fluid = 0, n_les = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (solid[i]) continue;
        ++n_fluid;
        if (d_wall[i] >= C_DES*DELTA_LES) ++n_les;
    }
    double les_frac = (double)n_les / (double)n_fluid;

    int snap_idx = 0;
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        if (snap_idx < SNAPSHOTS && t == snap_steps[snap_idx]) {
            if (output_snapshot(t) != 0) return 1;
            ++snap_idx;
        }
        if (t % HISTORY_INTERVAL == 0) {
            double umax = 0, nut_sum = 0, nt_sum = 0;
            int nf = 0;
            for (int i = 0; i < NX*NY; ++i) {
                if (solid[i]) continue;
                if (fabs(u[i]) > umax) umax = fabs(u[i]);
                nut_sum += nut_field[i];
                nt_sum += nu_tilde[i];
                ++nf;
            }
            int p = IDX(PROBE_X, PROBE_Y);
            double inv_n = 1.0 / nf;
            fprintf(hist, "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    t, umax, u[p], v[p], nut_sum*inv_n, nt_sum*inv_n, les_frac);
        }
        update_sa_des();
        stream_collide_with_des();
    }
    macroscopic();
    if (output_snapshot(NSTEPS) != 0) return 1;
    fclose(hist);

    double umax_actual = 0;
    for (int i = 0; i < NX*NY; ++i) {
        if (!solid[i] && fabs(u[i]) > umax_actual) umax_actual = fabs(u[i]);
    }
    double Re_D = umax_actual * 2 * R_CYL / nu0;
    printf("Done. Snapshots: karman_des_hires_snapshot_*.csv, probe: karman_des_hires_probe.csv\n");
    printf("Parameters: NX=%d NY=%d D=%d Cyl=(%d,%d) NSTEPS=%d TAU=%.3f F=%g C_DES=%.3f nu0=%.5f u_max=%.4f Re_D=%.0f\n",
           NX, NY, 2*R_CYL, CX, CY, NSTEPS, TAU, FORCE_X, C_DES, nu0, umax_actual, Re_D);
    printf("DES branch: LES fraction = %.4f over fluid cells (RANS layer = 1 cell at walls and around cylinder)\n", les_frac);
    return 0;
}
