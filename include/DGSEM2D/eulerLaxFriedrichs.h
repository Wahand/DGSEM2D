#ifndef EULERLAXFRIEDRICHS_INCLUDED
#define EULERLAXFRIEDRICHS_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "DGSEM2D/equationFlux.h"
#include "DGSEM2D/equationVars.h"

void eulerLaxFriedrichs(const double* const restrict ULeft, const double* const restrict URight, double* const restrict flux);

#endif
