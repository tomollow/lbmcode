// sa_closure.h
// Spalart-Allmaras (1992) standard closure constants and helper routines
// shared by the SA-DES97 implementations in this directory:
//   - cavity_des.c, karman_des.c, karman_des_hires.c, backward_step_des.c
//
// Each file additionally defines its own DES length-scale switch and
// transport time step (C_DES, DELTA_LES, SA_DT) — those are kept per-file
// so they can be tuned without touching this header.
#ifndef SEC4_SA_CLOSURE_H
#define SEC4_SA_CLOSURE_H

#include <math.h>

#define KAPPA 0.41
#define C_B1  0.1355
#define C_B2  0.622
#define SIG_SA (2.0/3.0)
#define C_V1  7.1
#define C_W1  (C_B1/(KAPPA*KAPPA) + (1.0 + C_B2)/SIG_SA)
#define C_W2  0.3
#define C_W3  2.0

// SA viscosity ratio: nu_t = nu_tilde * f_v1(chi), chi = nu_tilde / nu0.
// For chi << c_v1 = 7.1 the model is "asleep" (f_v1 ~ 1e-8 at chi = 0.1).
static inline double sa_fv1(double chi) {
    double chi3 = chi*chi*chi;
    return chi3 / (chi3 + C_V1*C_V1*C_V1);
}

// SA-DES production/destruction closure: given nu_tilde nt, molecular
// viscosity nu0, length scale dl = d_tilde, and strain-rate magnitude Smag,
// returns the modified strain rate Stilde and wall function fw used by
// the source terms (prod = C_B1*Stilde*nt, dest = C_W1*fw*(nt/dl)^2).
typedef struct { double Stilde; double fw; } sa_terms_t;

static inline sa_terms_t sa_compute_terms(double nt, double nu0, double dl, double Smag) {
    double chi = nt / nu0;
    double fv1 = sa_fv1(chi);
    double fv2 = 1.0 - chi / (1.0 + chi*fv1);
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
    sa_terms_t out = { Stilde, fw };
    return out;
}

#endif // SEC4_SA_CLOSURE_H
