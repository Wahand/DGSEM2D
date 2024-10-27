#include "DGSEM2D/riemannSolver.h"

void riemannSolver(const int N, const int nCellsX, const int nCellsY,
                   const double* const restrict ULeft, const double* const restrict URight,
                   #ifdef parabol
                     const double* const restrict gradUxLeft, const double* const restrict gradUxRight,
                     const double* const restrict gradUyLeft, const double* const restrict gradUyRight,
                   #endif
                   const double normalX, const double normalY,
                   double* const restrict flux) {

  // Linear advection
  #if (eqSys == 1)
    // ax, ay > 0 --> ULeft nehmen
    if (eqVars_ax > 0.0) {
      flux[0] += normalX * eqVars_ax * ULeft[0];
    } else {
      flux[0] += normalX * eqVars_ax * URight[0];
    }
    if (eqVars_ay > 0.0) {
      flux[0] += normalY * eqVars_ay * ULeft[0];
    } else {
      flux[0] += normalY * eqVars_ay * URight[0];
    }
  #endif

  // Euler-equations
  #if (eqSys == 2)

    // variables in local coordinate-system
    double localLeftState[nVars];
    double localRightState[nVars];
    double localFlux[nVars];

    // transform left and right state into local coordinate system
    localLeftState[0] = ULeft[0];
    localLeftState[1] = normalX * ULeft[1] + normalY * ULeft[2];
    localLeftState[2] = -normalY * ULeft[1] + normalX * ULeft[2];
    localLeftState[3] = ULeft[3];

    localRightState[0] = URight[0];
    localRightState[1] = normalX * URight[1] + normalY * URight[2];
    localRightState[2] = -normalY * URight[1] + normalX * URight[2];
    localRightState[3] = URight[3];


    #if (rSolver == 1)
      eulerLaxFriedrichs(localLeftState, localRightState, localFlux);
    #endif

    // transform flux back into global coordinate system
    flux[0] = localFlux[0];
    flux[1] = normalX * localFlux[1] - normalY * localFlux[2];
    flux[2] = normalY * localFlux[1] + normalX * localFlux[2];
    flux[3] = localFlux[3];
  #endif

  // Advection-Diffusion
  #if (eqSys == 3)
    // ax, ay > 0 --> ULeft nehmen
    if (eqVars_ax > 0.0) {
      flux[0] += normalX * eqVars_ax * ULeft[0];
    } else {
      flux[0] += normalX * eqVars_ax * URight[0];
    }
    if (eqVars_ay > 0.0) {
      flux[0] += normalY * eqVars_ay * ULeft[0];
    } else {
      flux[0] += normalY * eqVars_ay * URight[0];
    }
  #endif

}
