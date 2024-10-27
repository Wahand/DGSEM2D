#include "DGSEM2D/eulerLaxFriedrichs.h"


void eulerLaxFriedrichs(const double* const restrict ULeft, const double* const restrict URight, double* const restrict flux) {

//TODO unschön
#if (eqSys == 2)

  const double gamma = eqVars_gamma;
  const double pLeft = (gamma - 1.0) * (ULeft[3] - 0.5 * (ULeft[1] * ULeft[1] / ULeft[0] + ULeft[2] * ULeft[2] / ULeft[0]));
  const double pRight = (gamma - 1.0) * (URight[3] - 0.5 * (URight[1] * URight[1] / URight[0] + URight[2] * URight[2] / URight[0]));

  const double cLeft = sqrt(gamma * pLeft / ULeft[0]);
  const double cRight = sqrt(gamma * pRight / URight[0]);
  const double a = fmax(fmax(fabs(ULeft[1] / ULeft[0]) + cLeft, fabs(ULeft[2] / ULeft[0]) + cLeft),
                        fmax(fabs(URight[1] / URight[0]) + cRight, fabs(URight[2] / URight[0]) + cRight));

  double leftFlux[nVars];
  double rightFlux[nVars];

  equationFluxX(ULeft,
                #ifdef parabol
                  gradUxLeft, gradUyLeft,
                #endif
                leftFlux);
  equationFluxX(URight,
                #ifdef parabol
                  gradUxRight, gradUyRight,
                #endif
                rightFlux);

  for (int iVar = 0; iVar < nVars; iVar++) {
    flux[iVar] = (leftFlux[iVar] + rightFlux[iVar] - (URight[iVar] - ULeft[iVar]) * a) * 0.5;
  }

#endif

}
