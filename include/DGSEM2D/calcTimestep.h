#ifndef CALCTIMESTEP_INCLUDED
#define CALCTIMESTEP_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "DGSEM2D/equationVars.h"

void calcTimestep(const int N,
                  const int nCellsX, const int nCellsY,
                  const double* const restrict U, const double* const restrict sJ,
                  const double* const restrict Sx, const double* const restrict Sy,
                  const double* const restrict ACell,
                  const double CFL,
                  double* const restrict timestep);
                  

#endif
