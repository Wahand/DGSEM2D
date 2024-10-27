#include "DGSEM2D/interpolation.h"

void applyJacobianToPhysical(double* const restrict U, const double* const restrict sJ, const int size) {
  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int i = 0; i < size; i++) {
      U[iVar * size + i] = U[iVar * size +i] * sJ[i];
    }
  }
}

void applyJacobianToReference(double* const restrict U, double* const restrict JU, const double* const restrict sJ, const int size) {
  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int i = 0; i < size; i++) {
      JU[iVar * size + i] = U[iVar * size + i] / sJ[i];
    }
  }
}
