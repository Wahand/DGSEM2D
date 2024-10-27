#ifndef INTERPOLATION_INCLUDED
#define INTERPOLATION_INCLUDED

#include <stdlib.h>
#include <stdio.h>

void applyJacobianToPhysical(double* const restrict U, const double* const restrict sJ, const int size);

void applyJacobianToReference(double* const restrict U, double* const restrict JU, const double* const restrict sJ, const int size);

#endif
