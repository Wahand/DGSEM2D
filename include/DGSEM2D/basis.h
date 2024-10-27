#ifndef BASIS_INCLUDED
#define BASIS_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

void legendreGaussNodesAndWeights(const int N, double* const restrict xGP, double* const restrict wGP);

void helpLegendrePolynomialAndDerivative(const int N, const double* const restrict x, double* const restrict L, double* const restrict Lder);

void barycentricWeights(const int N, const double* const restrict xGP, double* const restrict wBary);

void lagrangeInterpolationPolys(const int N, const double x, const double* const restrict xGP, const double* const restrict wBary, double* const restrict L);

void polynomialDerivativeMatrix(const int N, const double* const restrict xGP, const double* const restrict wBary, double* const restrict D);

#endif
