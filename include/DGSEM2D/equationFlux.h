#ifndef EQUATIONFLUX_INCLUDED
#define EQUATIONFLUX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include "DGSEM2D/equationVars.h"

void equationFluxX(const double* const restrict U,
                   #ifdef parabol
                     const double* const restrict gradUx, const double* const restrict gradUy,
                   #endif
                   double* const restrict flux);

void equationFluxY(const double* const restrict U,
                   #ifdef parabol
                     const double* const restrict gradUx, const double* const restrict gradUy,
                   #endif
                   double* const restrict flux);

#endif
