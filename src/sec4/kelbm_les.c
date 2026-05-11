// kelbm_les.c
// 2D LBM + Smagorinsky LES (channel flow, D2Q9, BGK)
//
// Same geometry, body force, and halfway bounce-back walls as kelbm.c. Instead
// of solving k-epsilon transport equations, we use the standard Smagorinsky
// sub-grid-scale model: nu_t = (Cs * Delta)^2 * |S|, with |S| = sqrt(2 S_ij S_ij)
// and Cs = 0.16, Delta = 1 (grid spacing in LU). nu_t is fed back into a
// locally varying tau_eff = 0.5 + 3*(nu_0 + nu_t) in the BGK collision.
// No wall damping (standard Smagorinsky), so nu_t stays finite at the wall.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 200
#define NY 60
#define NDIR 9
#define NSTEPS 30000
#define TAU 0.55
#define OMEGA (1.0/TAU)
#define FORCE_X 5e-6
#define CS_SMAG 0.16                 // Smagorinsky constant
#define DELTA_LES 1.0                // grid spacing in LU

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

void initialize() {
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            nut_field[i] = 0.0;
            for (int d = 0; d < NDIR; ++d) {
                f[i*NDIR + d] = w[d] * rho[i];
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
                double Fi = (1.0 - 0.5*omega_eff) * w[d] * FORCE_X *
                            (3.0*(cx[d] - u[i]) + 9.0*eu*cx[d]);
                double post = f[i*NDIR + d] - omega_eff*(f[i*NDIR + d] - feq) + Fi;
                int xp = (x + cx[d] + NX) % NX;
                int yp = y + cy[d];
                if (yp < 0 || yp >= NY) {
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
            u[i] = (ru + 0.5*FORCE_X) / rr;
            v[i] = rv / rr;
        }
    }
}

void update_les() {
    // Smagorinsky: nu_t = (Cs*Delta)^2 * sqrt(2 S_ij S_ij)
    double dx = 1.0, dy = 1.0;
    double Cs2 = (CS_SMAG * DELTA_LES) * (CS_SMAG * DELTA_LES);
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            int ixp = IDX((x+1) % NX, y);
            int ixm = IDX((x-1+NX) % NX, y);
            double dudx = (u[ixp] - u[ixm]) / (2.0*dx);
            double dvdx = (v[ixp] - v[ixm]) / (2.0*dx);
            double dudy, dvdy;
            if (y == 0) {
                dudy = (u[IDX(x, 1)] - u[i]) / dy;
                dvdy = (v[IDX(x, 1)] - v[i]) / dy;
            } else if (y == NY-1) {
                dudy = (u[i] - u[IDX(x, NY-2)]) / dy;
                dvdy = (v[i] - v[IDX(x, NY-2)]) / dy;
            } else {
                dudy = (u[IDX(x, y+1)] - u[IDX(x, y-1)]) / (2.0*dy);
                dvdy = (v[IDX(x, y+1)] - v[IDX(x, y-1)]) / (2.0*dy);
            }
            double S11 = dudx, S22 = dvdy, S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            double Smag = sqrt(S2);
            nut_field[i] = Cs2 * Smag;
        }
    }
}

int output_csv(const char* fname) {
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_csv: cannot open %s for writing\n", fname);
        return 1;
    }
    fprintf(fp, "x,y,u,v,nut\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int i = IDX(x, y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g\n", x, y, u[i], v[i], nut_field[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main() {
    initialize();
    for (int t = 0; t < NSTEPS; ++t) {
        macroscopic();
        update_les();
        stream_collide_with_les();
        if (t % 200 == 0) printf("step %d\n", t);
    }
    macroscopic();
    if (output_csv("kelbm_les_output.csv") != 0) {
        return 1;
    }
    printf("Done. Output: kelbm_les_output.csv\n");
    double Re = 0;
    double umax = 0;
    for (int i = 0; i < NX*NY; ++i) if (fabs(u[i]) > umax) umax = fabs(u[i]);
    Re = umax * NY / nu0;
    printf("Parameters: NX=%d NY=%d NSTEPS=%d TAU=%.3f Cs=%.3f nu0=%.5f u_max=%.4f Re_max=%.0f\n",
           NX, NY, NSTEPS, TAU, CS_SMAG, nu0, umax, Re);
    return 0;
}
