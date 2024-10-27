#include "DGSEM2D/calcTimestep.h"

void calcTimestep(const int N,
                  const int nCellsX, const int nCellsY,
                  const double* const restrict U, const double* const restrict sJ,
                  const double* const restrict Sx, const double* const restrict Sy,
                  const double* const restrict ACell,
                  const double CFL,
                  double* const timestep) {

  // Linear Advection
  #if (eqSys == 1)
    // Calculation as in DGSEM
    const double ax = eqVars_ax;
    const double ay = eqVars_ay;

    double maxLambda = sJ[0] * ax;

    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPY = 0; iGPY < N + 1; iGPY++) {
          for (int iGPX = 0; iGPX < N + 1; iGPX++) {
            // advection
            //  maxLambda = max(sJ * ax, sJ * ay)
            maxLambda = fmax(maxLambda, fmax(sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * ax, sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * ay));
          }
        }
      }
    }

    *timestep = CFL * 2.0 / (maxLambda * (2.0 * N + 1)) ;
  #endif


  // Euler
  #if (eqSys == 2)
    // Calculation as in DGSEM

    const double gamma = eqVars_gamma;

    double maxLambda = 0.0;

    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPY = 0; iGPY < N + 1; iGPY++) {
          for (int iGPX = 0; iGPX < N + 1; iGPX++) {
            int index = iCellY * (N + 1) * nCellsX * (N + 1) * nVars + iGPY * nCellsX * (N + 1) * nVars + iCellX * (N + 1) * nVars + iGPX;
            double v1 = U[index + 1] / U[index];
            double v2 = U[index + 2] / U[index];
            double p = (gamma - 1.0) * (U[index + 3] - 0.5 * (U[index + 1] * v1 + U[index + 2] * v2));
            double c = sqrt(gamma * p / U[index]);
            maxLambda = fmax(maxLambda, sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX* (N + 1) + iCellX * (N + 1) + iGPX] * fmax(fmax(fabs(v1 - c), fabs(v1 + c)), fmax(fabs(v2 - c), fabs(v2 + c))));
          }
        }
      }
    }

    *timestep = CFL * 2.0 / (maxLambda * (2.0 * N + 1)) ;
  #endif


  // Advection-Diffusion
  #if (eqSys == 3)

    // Calculation as in DGSEM
    double ax = eqVars_ax;
    double ay = eqVars_ay;
    double b = eqVars_b;

    double maxLambda = sJ[0] * ax;

    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPY = 0; iGPY < N + 1; iGPY++) {
          for (int iGPX = 0; iGPX < N + 1; iGPX++) {
            // advection
            //  maxLambda = max(sJ * ax, sJ * ay)
            maxLambda = fmax(maxLambda, fmax(sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * ax, sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * ay));
            // diffusion
            //  maxLambda = sJ^2 * b
            maxLambda = fmax(maxLambda, sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
                                        sJ[iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
                                        b);
          }
        }
      }
    }

    *timestep = CFL * 2.0 / (maxLambda * (2.0 * N + 1)) ;
  #endif



}
