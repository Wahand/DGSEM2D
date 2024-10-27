#include "DGSEM2D/basis.h"

void legendreGaussNodesAndWeights(const int N, double* const restrict xGP, double* const restrict wGP) {

  if (xGP == NULL && wGP == NULL) {
    printf("no arguments given\n");
    return;
  }

  double* xGP_ = malloc((N + 1) * sizeof(double));
  double* wGP_ = malloc((N + 1) * sizeof(double));

  if (N == 0) {
    xGP_[0] = 0.0;
    wGP_[0] = 2.0;
  } else if (N == 1) {
    xGP_[0] = -sqrt(1.0 / 3.0);
    xGP_[1] = -xGP_[0];
    wGP_[0] = 1.0;
    wGP_[1] = 1.0;
  } else {
    double tol = 1.0e-15;
    int nIter = 10;
    double cheb_tmp = M_PI / (2.0 * N + 2.0);

    for (int iGP = 0; iGP < (N + 1) / 2; iGP++) {
      xGP_[iGP] = -cos(cheb_tmp * (2.0 * iGP + 1.0));
      for (int iter = 0; iter < nIter; iter++) {
        double L_Np1, Lder_Np1;
        helpLegendrePolynomialAndDerivative(N + 1, &xGP_[iGP], &L_Np1, &Lder_Np1);
        double dx = -L_Np1 / Lder_Np1;
        xGP_[iGP] += dx;
        if (fabs(dx) < tol * fabs(xGP_[iGP])) {
          break;
        }
      }

      double L_Np1;
      double Lder_Np1;
      helpLegendrePolynomialAndDerivative(N + 1, &xGP_[iGP], &L_Np1, &Lder_Np1);
      xGP_[N - iGP] = -xGP_[iGP];
      wGP_[iGP] = (2.0 * N + 3) / ((1.0 - xGP_[iGP] * xGP_[iGP]) * Lder_Np1 * Lder_Np1);
      wGP_[N - iGP] = wGP_[iGP];
    }

    if (N % 2 == 0) {
      xGP_[N / 2] = 0.0;
      double L_Np1;
      double Lder_Np1;
      helpLegendrePolynomialAndDerivative(N + 1, &xGP_[N / 2], &L_Np1, &Lder_Np1);
      wGP_[N / 2] = (2.0 * N + 3) / (Lder_Np1 * Lder_Np1);
    }
  }

  for (int i = 0; i < N + 1; i++) {
    if (xGP != NULL) {
      xGP[i] = xGP_[i];
    }
    if (wGP != NULL) {
      wGP[i] = wGP_[i];
    }
  }
  free(xGP_);
  free(wGP_);
}

void helpLegendrePolynomialAndDerivative(const int N, const double* const restrict x, double* const restrict L, double* const restrict Lder) {

  if (N == 0) {
    *L = 1.0;
    *Lder = 0.0;
  } else if (N == 1) {
    *L = *x;
    *Lder = 1.0;
  } else {
    double L_Nm2 = 1.0;
    double L_Nm1 = *x;
    double Lder_Nm2 = 0.0;
    double Lder_Nm1 = 1.0;

    *L = 0.0;
    *Lder = 0.0;

    for (int iLegendre = 2; iLegendre <= N; iLegendre++) {
      *L = ((2.0 * iLegendre - 1.0) * (*x) * L_Nm1 - (iLegendre - 1.0) * L_Nm2) / iLegendre;
      *Lder = Lder_Nm2 + (2.0 * iLegendre - 1) * L_Nm1;
      L_Nm2 = L_Nm1;
      L_Nm1 = *L;
      Lder_Nm2 = Lder_Nm1;
      Lder_Nm1 = *Lder;
    }

    *L = (*L) * sqrt(N + 0.5);
    *Lder = (*Lder) * sqrt(N + 0.5);
  }
}

void barycentricWeights(const int N, const double* const restrict xGP, double* const restrict wBary) {
  for (int i = 0; i < N + 1; i++) {
    wBary[i] = 1.0;
    for (int j = 0; j < N + 1; j++) {
      if (j != i) {
        wBary[i] = wBary[i] / (xGP[i] - xGP[j]);
      }
    }
  }
}

void lagrangeInterpolationPolys(const int N, const double x, const double* const restrict xGP, const double* const restrict wBary, double* const restrict L) {
  // initialize with zeros
  for (int i = 0; i < N + 1; i++) {
    L[i] = 0.0;
  }
  // check if x coincides with any of xGP
  int found = 0;
  for (int i = 0; i < N + 1; i++) {
    if (x == xGP[i]) {
      found = 1;
      L[i] = 1.0;
      break;
    }
  }
  // otherwise compute L
  if (!found) {
    for (int j = 0; j < N + 1; j++) {
      L[j] = 0.0;
      for (int i = 0; i < N + 1; i++) {
        L[j] = L[j] + (wBary[i] / (x - xGP[i]));
      }
      L[j] = wBary[j] / ((x - xGP[j]) * L[j]);
    }
  }
}

void polynomialDerivativeMatrix(const int N, const double* const restrict xGP, const double* const restrict wBary, double* const restrict D) {
  for (int iGP = 0; iGP < N + 1; iGP++) {
    // D[i, i] = 0.0
    D[iGP * (N + 1) + iGP] = 0.0;
    for (int iLagrange = 0; iLagrange < N + 1; iLagrange++) {
      if (iLagrange != iGP) {
        // D[iGP, iLagrange]
        D[(N + 1) * iGP + iLagrange] = wBary[iLagrange] / (wBary[iGP] * (xGP[iGP] - xGP[iLagrange]));
        // D[iGP, iGP] = D[iGP, iGP] - D[iGP, iLagrange];
        D[iGP * (N + 1) + iGP] = D[iGP * (N + 1) + iGP] - D[iGP * (N + 1) + iLagrange];
//        D[(N + 1) * iGP + iGP] -= D[(N + 1) * iGP + iLagrange];
      }
    }
  }
}
