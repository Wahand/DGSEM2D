#include "DGSEM2D/equationFlux.h"

void equationFluxX(const double* const restrict U,
                   #ifdef parabol
                     const double* const restrict gradUx, const double* const restrict gradUy,
                   #endif
                   double* const restrict flux) {

  // Linear Advection
  #if (eqSys == 1)
    // advection-velocity
    flux[0] = U[0] * eqVars_ax;
  #endif

  // Euler
  #if (eqSys == 2)
    const double gamma = eqVars_gamma;
    const double pressure = (gamma - 1.0) * (U[3] - 0.5 * (U[1] * U[1] / U[0] + U[2] * U[2] / U[0]));
    flux[0] = U[1];
    flux[1] = U[1] * U[1] / U[0] + pressure;
    flux[2] = U[1] * U[2] / U[0];
    flux[3] = U[1] / U[0] * (U[3] + pressure);
  #endif

  // Advection-Diffusion
  #if (eqSys == 3)
//    flux[0] = U[0] * eqVars_ax - gradUx[0] * eqVars_b;
    flux[0] = U[0] * eqVars_ax;
  #endif

}

void equationFluxY(const double* const restrict U,
                   #ifdef parabol
                     const double* const restrict gradUx, const double* const restrict gradUy,
                   #endif
                   double* const restrict flux) {

  // Linear Advection
  #if (eqSys == 1)
    flux[0] = U[0] * eqVars_ay;
  #endif

  // Euler
  #if (eqSys == 2)
    const double gamma = eqVars_gamma;
    const double pressure = (gamma - 1.0) * (U[3] - 0.5 * (U[1] * U[1] / U[0] + U[2] * U[2] / U[0]));
    flux[0] = U[2];
    flux[1] = U[1] * U[2] / U[0];
    flux[2] = U[2] * U[2] / U[0] + pressure;
    flux[3] = U[2] / U[0] * (U[3] + pressure);
  #endif

  // Advection-Diffusion
  #if (eqSys == 3)
//    flux[0] = U[0] * eqVars_ay - gradUy[0] * eqVars_b;
    flux[0] = U[0] * eqVars_ay;
  #endif

}
