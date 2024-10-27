#ifndef RIEMANNSOLVER_INCLUDED
#define RIEMANNSOLVER_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include "DGSEM2D/eulerLaxFriedrichs.h"
#include "DGSEM2D/equationVars.h"

void riemannSolver(const int N, const int nCellsX, const int nCellsY,
                   const double* const restrict ULeft, const double* const restrict URight,
                   #ifdef parabol
                     const double* const restrict gradUxLeft, const double* const restrict gradUxRight,
                     const double* const restrict gradUyLeft, const double* const restrict gradUyRight,
                   #endif
                   const double normalX, const double normalY,
                   double* const restrict flux);

#endif
