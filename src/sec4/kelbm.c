// kelbm.c
// 2D LBM + k-epsilon turbulence model (channel flow test case)
// Domain: NX x NY, periodic in x, halfway bounce-back walls at y=0 and y=NY-1
// Driven by constant body force in x direction (Guo forcing 2002)
// Wall functions (Dirichlet) for k, eps at the first off-wall fluid cells
// Output: u,v,k,epsilon as CSV
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 200
#define NY 60
#define NDIR 9
#define NSTEPS 30000
#define TAU 0.55                    // BGK relaxation time (closer to 0.5 -> lower nu)
#define OMEGA (1.0/TAU)             // 1/tau, must be < 2 for stability
#define FORCE_X 5e-6                // body force per unit mass in +x direction
#define KAPPA 0.41                  // von Karman constant
#define WALL_DY 0.5                 // distance from first fluid node to wall (halfway BB)

// D2Q9 Lattice velocities
const int cx[NDIR] = {0,1,0,-1,0,1,-1,-1,1};
const int cy[NDIR] = {0,0,1,0,-1,1,1,-1,-1};
const double w[NDIR] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
// Opposite-direction index for halfway bounce-back
const int opp[NDIR] = {0,3,4,1,2,7,8,5,6};

// Macros
#define IDX(x,y) ((x) + NX*(y))

// Distribution function buffers — swapped via pointers each step (no memcpy)
static double f_buf_a[NX*NY*NDIR];
static double f_buf_b[NX*NY*NDIR];
static double *f = f_buf_a;
static double *f2 = f_buf_b;

// Main variables
static double u[NX*NY], v[NX*NY], rho[NX*NY];
static double k[NX*NY], eps[NX*NY];

// k-epsilon model parameters (standard values)
#define Cmu 0.09
#define Ce1 1.44
#define Ce2 1.92
#define sig_k 1.0
#define sig_e 1.3
// LBM分子動粘性: nu = c_s^2 * (tau - 0.5) = (TAU - 0.5)/3
#define nu0 ((TAU - 0.5)/3.0)

// 一時配列
static double k_new[NX*NY], eps_new[NX*NY];

void initialize() {
    // Start from rest; body force drives flow.
    // Seed k, eps with non-trivial values so that k-eps can self-sustain
    // (otherwise k stays at the floor and nu_t never grows).
    // u_tau from force balance: u_tau^2 = F * delta, delta = NY/2
    double u_tau = sqrt(FORCE_X * (NY/2.0));
    double k_seed = u_tau*u_tau / sqrt(Cmu);          // wall-equilibrium k
    double eps_seed = u_tau*u_tau*u_tau / (KAPPA*WALL_DY); // wall-equilibrium eps
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = 1.0;
            k[i] = k_seed;
            eps[i] = eps_seed;
            for(int d=0; d<NDIR; ++d) {
                f[i*NDIR + d] = w[d]*rho[i];
            }
        }
    }
}

void stream_collide() {
    // BGK collision with Guo forcing, periodic in x, halfway bounce-back at y=0 and y=NY-1
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            double usqr = u[i]*u[i] + v[i]*v[i];
            for(int d=0; d<NDIR; ++d) {
                double eu = cx[d]*u[i] + cy[d]*v[i];
                double feq = w[d]*rho[i]*(1.0 + 3.0*eu + 4.5*eu*eu - 1.5*usqr);
                // Guo forcing source term for body force F = (FORCE_X, 0)
                double Fi = (1.0 - 0.5*OMEGA) * w[d] * FORCE_X *
                            (3.0*(cx[d] - u[i]) + 9.0*eu*cx[d]);
                double post = f[i*NDIR + d] - OMEGA*(f[i*NDIR + d] - feq) + Fi;
                int xp = (x + cx[d] + NX) % NX;
                int yp = y + cy[d];
                if (yp < 0 || yp >= NY) {
                    // halfway bounce-back: reflect to opposite direction at the same node
                    f2[i*NDIR + opp[d]] = post;
                } else {
                    int ip = IDX(xp, yp);
                    f2[ip*NDIR + d] = post;
                }
            }
        }
    }
    // Swap pointers (avoids 2*NX*NY*NDIR memory writes per step)
    double *tmp = f; f = f2; f2 = tmp;
}

void macroscopic() {
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            double ru = 0, rv = 0, rr = 0;
            for(int d=0; d<NDIR; ++d) {
                double ff = f[i*NDIR + d];
                rr += ff;
                ru += ff * cx[d];
                rv += ff * cy[d];
            }
            rho[i] = rr;
            // Half-step body-force correction (Guo)
            u[i] = (ru + 0.5*FORCE_X) / rr;
            v[i] = rv / rr;
        }
    }
}

void update_kepsilon() {
    // 標準k-ε輸送方程式（陽解法）。x方向は周期境界、y方向は壁。
    // バルクの y 境界はミラー（Neumann）で扱い、第1流体セル(y=0, y=NY-1)は
    // 末尾で壁関数値を Dirichlet 上書きする。
    double dx = 1.0, dy = 1.0, dt = 0.01;
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            int ixp = IDX((x+1)%NX, y);
            int ixm = IDX((x-1+NX)%NX, y);
            // y方向は壁境界（ミラー = ゼロ勾配）
            int iyp = IDX(x, (y+1<NY) ? y+1 : y);
            int iym = IDX(x, (y-1>=0) ? y-1 : y);
            // 速度勾配（x: 中心差分、y: 壁では一次片側差分）
            double dudx = (u[ixp] - u[ixm])/(2.0*dx);
            double dvdx = (v[ixp] - v[ixm])/(2.0*dx);
            double dudy, dvdy;
            if (y == 0) {
                dudy = (u[IDX(x,1)] - u[i])/dy;
                dvdy = (v[IDX(x,1)] - v[i])/dy;
            } else if (y == NY-1) {
                dudy = (u[i] - u[IDX(x,NY-2)])/dy;
                dvdy = (v[i] - v[IDX(x,NY-2)])/dy;
            } else {
                dudy = (u[iyp] - u[iym])/(2.0*dy);
                dvdy = (v[iyp] - v[iym])/(2.0*dy);
            }
            // S_ij S_ij
            double S11 = dudx, S22 = dvdy;
            double S12 = 0.5*(dudy + dvdx);
            double S2 = 2.0*(S11*S11 + S22*S22) + 4.0*S12*S12;
            // 乱流粘性
            double nut = Cmu * k[i]*k[i]/(eps[i]+1e-12);
            // 生成項
            double Pk = nut * S2;
            // k, eps の拡散係数
            double Dk = nu0 + nut/sig_k;
            double De = nu0 + nut/sig_e;
            // k, eps の拡散項（2次中心差分、壁ではミラーによりゼロ勾配）
            double lap_k = (k[ixp] + k[ixm] + k[iyp] + k[iym] - 4.0*k[i])/(dx*dx);
            double lap_e = (eps[ixp] + eps[ixm] + eps[iyp] + eps[iym] - 4.0*eps[i])/(dx*dx);
            // k, eps の対流項（1次風上差分）
            double uk = u[i]*(
                (k[i] - k[ixm])/dx * (u[i]>0?1:0) + (k[ixp] - k[i])/dx * (u[i]<=0?1:0)
            );
            double vk = v[i]*(
                (k[i] - k[iym])/dy * (v[i]>0?1:0) + (k[iyp] - k[i])/dy * (v[i]<=0?1:0)
            );
            double ue = u[i]*(
                (eps[i] - eps[ixm])/dx * (u[i]>0?1:0) + (eps[ixp] - eps[i])/dx * (u[i]<=0?1:0)
            );
            double ve = v[i]*(
                (eps[i] - eps[iym])/dy * (v[i]>0?1:0) + (eps[iyp] - eps[i])/dy * (v[i]<=0?1:0)
            );
            // k, eps の時間発展
            k_new[i] = k[i] + dt * (
                Dk*lap_k - uk - vk + Pk - eps[i]
            );
            eps_new[i] = eps[i] + dt * (
                De*lap_e - ue - ve + Ce1*Pk*eps[i]/(k[i]+1e-12) - Ce2*eps[i]*eps[i]/(k[i]+1e-12)
            );
        }
    }
    // 壁関数（halfway bounce-back の第1流体セルに k, eps を Dirichlet 設定）
    // 摩擦速度 u_tau は力学平衡 u_tau^2 = F*delta から評価
    double u_tau = sqrt(FORCE_X * (NY/2.0));
    double k_wall = u_tau*u_tau / sqrt(Cmu);
    double eps_wall = u_tau*u_tau*u_tau / (KAPPA * WALL_DY);
    for(int x=0; x<NX; ++x) {
        int i_bot = IDX(x, 0);
        int i_top = IDX(x, NY-1);
        k_new[i_bot] = k_wall;
        eps_new[i_bot] = eps_wall;
        k_new[i_top] = k_wall;
        eps_new[i_top] = eps_wall;
    }
    // 非負制約と値の更新
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            if(k_new[i]<1e-8) k_new[i]=1e-8;
            if(eps_new[i]<1e-8) eps_new[i]=1e-8;
            k[i] = k_new[i];
            eps[i] = eps_new[i];
        }
    }
}

int output_csv(const char* fname) {
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "output_csv: cannot open %s for writing\n", fname);
        return 1;
    }
    fprintf(fp, "x,y,u,v,k,epsilon\n");
    for(int y=0; y<NY; ++y) {
        for(int x=0; x<NX; ++x) {
            int i = IDX(x,y);
            fprintf(fp, "%d,%d,%.9g,%.9g,%.9g,%.9g\n", x, y, u[i], v[i], k[i], eps[i]);
        }
    }
    fclose(fp);
    return 0;
}

int main() {
    initialize();
    for(int t=0; t<NSTEPS; ++t) {
        macroscopic();
        update_kepsilon();
        stream_collide();
        if(t%200==0) printf("step %d\n", t);
    }
    macroscopic();
    if (output_csv("kelbm_output.csv") != 0) {
        return 1;
    }
    printf("Done. Output: kelbm_output.csv\n");
    return 0;
}
