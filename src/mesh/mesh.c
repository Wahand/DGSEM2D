#include "DGSEM2D/mesh.h"

void createMesh(const int N,
                const double* const restrict nodesX, const double* const restrict nodesY,
                const int nCellsX, const int nCellsY,
                double* const restrict xGP, double* const restrict yGP,
                double* const restrict Sx, double* const restrict Sy, double* const restrict ACell,
                double* const restrict XiBoundX, double* const restrict XiBoundY,
                double* const restrict EtaBoundX, double* const restrict EtaBoundY,
                double* const restrict XiBoundNormalX, double* const restrict XiBoundNormalY,
                double* const restrict EtaBoundNormalX, double* const restrict EtaBoundNormalY,
                double* const restrict XiBoundMetricsPlus, double* const restrict XiBoundMetricsMinus,
                double* const restrict EtaBoundMetricsPlus, double* const restrict EtaBoundMetricsMinus,
                double* const restrict XiBoundScalingPlus, double* const restrict XiBoundScalingMinus,
                double* const restrict EtaBoundScalingPlus, double* const restrict EtaBoundScalingMinus,
                double* const restrict dXdXi, double* const restrict dXdEta,
                double* const restrict dYdXi, double* const restrict dYdEta,
                double* const restrict sJ) {

  // allocate memory for temporary coordinates
  double* xGP_ref = malloc((N + 1) * sizeof(double));

  // get coordinates in reference element
  legendreGaussNodesAndWeights(N, xGP_ref, NULL);

  // local coordinates and metrics in cells
  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    for (int iCellX = 0; iCellX < nCellsX; iCellX++) {

      // node-coordinates
      double x1 = nodesX[iCellY * (nCellsX + 1) + iCellX];
      double x2 = nodesX[iCellY * (nCellsX + 1) + iCellX + 1];
      double x3 = nodesX[(iCellY + 1) * (nCellsX + 1) + iCellX + 1];
      double x4 = nodesX[(iCellY + 1) * (nCellsX + 1) + iCellX];

      double y1 = nodesY[iCellY * (nCellsX + 1) + iCellX];
      double y2 = nodesY[iCellY * (nCellsX + 1) + iCellX + 1];
      double y3 = nodesY[(iCellY + 1) * (nCellsX + 1) + iCellX + 1];
      double y4 = nodesY[(iCellY + 1) * (nCellsX + 1) + iCellX];

      // cell-dimensions in coordinate-directions
      Sx[iCellY * nCellsX + iCellX] = fabs(fmax(fmax(x1, x2), fmax(x3, x4)) - fmin(fmin(x1, x2), fmin(x3, x4)));
      Sy[iCellY * nCellsX + iCellX] = fabs(fmax(fmax(y1, y2), fmax(y3, y4)) - fmin(fmin(y1, y2), fmin(y3, y4)));

      // cell-areas
      ACell[iCellY * nCellsX + iCellX] = 0.5 * fabs((y1 - y3) * (x4 - x2) + (y2 - y4) * (x1 - x3));

      for (int iGPY = 0; iGPY < N + 1; iGPY++) {
        for (int iGPX = 0; iGPX < N + 1; iGPX++) {

          int index = iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX;

          // physical coordinates in cells
          double xRef = xGP_ref[iGPX];
          double yRef = xGP_ref[iGPY];

          xGP[index] = x1 + 0.5 * (yRef + 1) * (x4 - x1) + 0.5 * (xRef + 1) * (x2 + 0.5 * (yRef + 1) * (x3 - x2) - (x1 + 0.5 * (yRef + 1) * (x4 - x1)));
          yGP[index] = y1 + 0.5 * (yRef + 1) * (y4 - y1) + 0.5 * (xRef + 1) * (y2 + 0.5 * (yRef + 1) * (y3 - y2) - (y1 + 0.5 * (yRef + 1) * (y4 - y1)));

          // metrics in cells
          quadMapMetrics(x1, y1,
                         x2, y2,
                         x3, y3,
                         x4, y4,
                         xRef, yRef,
                         &dXdXi[index], &dXdEta[index],
                         &dYdXi[index], &dYdEta[index]);

          // inverse of the jacobian in cells
          sJ[index] = 1.0 / (dXdXi[index] * dYdEta[index] - dXdEta[index] * dYdXi[index]);
        }
      }
    }
  }

  // metrics on boundary
  double xXi, xEta, yXi, yEta;

  for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {

      // node-coordinates
      double x1 = nodesX[iCellY * (nCellsX + 1) + iCellX];
      double x2 = nodesX[iCellY * (nCellsX + 1) + iCellX + 1];
      double x3 = nodesX[(iCellY + 1) * (nCellsX + 1) + iCellX + 1];
      double x4 = nodesX[(iCellY + 1) * (nCellsX + 1) + iCellX];

      double y1 = nodesY[iCellY * (nCellsX + 1) + iCellX];
      double y2 = nodesY[iCellY * (nCellsX + 1) + iCellX + 1];
      double y3 = nodesY[(iCellY + 1) * (nCellsX + 1) + iCellX + 1];
      double y4 = nodesY[(iCellY + 1) * (nCellsX + 1) + iCellX];

      for (int iGPX = 0; iGPX < N + 1; iGPX++) {
        quadMapMetrics(x1, y1,
                       x2, y2,
                       x3, y3,
                       x4, y4,
                       xGP_ref[iGPX], -1.0,
                       &xXi, &xEta,
                       &yXi, &yEta);
        EtaBoundMetricsMinus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 0] = xXi;
        EtaBoundMetricsMinus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 1] = xEta;
        EtaBoundMetricsMinus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 2] = yXi;
        EtaBoundMetricsMinus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 3] = yEta;
        EtaBoundScalingMinus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] = sqrt(yXi * yXi + xXi * xXi);
        quadMapMetrics(x1, y1,
                       x2, y2,
                       x3, y3,
                       x4, y4,
                       xGP_ref[iGPX], 1.0,
                       &xXi, &xEta,
                       &yXi, &yEta);
        EtaBoundMetricsPlus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 0] = xXi;
        EtaBoundMetricsPlus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 1] = xEta;
        EtaBoundMetricsPlus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 2] = yXi;
        EtaBoundMetricsPlus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 3] = yEta;
        EtaBoundScalingPlus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] = sqrt(yXi * yXi + xXi * xXi);
      }
      for (int iGPY = 0; iGPY < N + 1; iGPY++) {
        quadMapMetrics(x1, y1,
                       x2, y2,
                       x3, y3,
                       x4, y4,
                       -1.0, xGP_ref[iGPY],
                       &xXi, &xEta,
                       &yXi, &yEta);
        XiBoundMetricsMinus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 0] = xXi;
        XiBoundMetricsMinus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 1] = xEta;
        XiBoundMetricsMinus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 2] = yXi;
        XiBoundMetricsMinus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 3] = yEta;
        XiBoundScalingMinus[iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] = sqrt(yEta * yEta + xEta * xEta);
        quadMapMetrics(x1, y1,
                       x2, y2,
                       x3, y3,
                       x4, y4,
                       1.0, xGP_ref[iGPY],
                       &xXi, &xEta,
                       &yXi, &yEta);
        XiBoundMetricsPlus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 0] = xXi;
        XiBoundMetricsPlus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 1] = xEta;
        XiBoundMetricsPlus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 2] = yXi;
        XiBoundMetricsPlus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 3] = yEta;
        XiBoundScalingPlus[iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] = sqrt(yEta * yEta + xEta * xEta);
      }
    }
  }


  // save coordinates and normals on boundary

  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    for (int iNodeX = 0; iNodeX < nCellsX + 1; iNodeX++) {

      // node-coordinates on that boundary
      double x1 = nodesX[iCellY * (nCellsX + 1) + iNodeX];
      double x2 = nodesX[(iCellY + 1) * (nCellsX + 1) + iNodeX];
      double y1 = nodesY[iCellY * (nCellsX + 1) + iNodeX];
      double y2 = nodesY[(iCellY + 1) * (nCellsX + 1) + iNodeX];

      for (int iGP = 0; iGP < N + 1; iGP++) {

        int index = iNodeX * nCellsY * (N + 1) + iCellY * (N + 1) + iGP;

        XiBoundX[index] = x1 + 0.5 * (xGP_ref[iGP] + 1.0) * (x2 - x1);
        XiBoundY[index] = y1 + 0.5 * (xGP_ref[iGP] + 1.0) * (y2 - y1);

        XiBoundNormalX[index] = (y2 - y1) / sqrt((y2 - y1) * (y2 - y1) + (x1 - x2) * (x1 - x2));
        XiBoundNormalY[index] = (x1 - x2) / sqrt((y2 - y1) * (y2 - y1) + (x1 - x2) * (x1 - x2));
      }
    }
  }

  for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
    for (int iNodeY = 0; iNodeY < nCellsY + 1; iNodeY++) {

      // node-coordinates on that boundary
      double x1 = nodesX[iNodeY * (nCellsX + 1) + iCellX];
      double x2 = nodesX[iNodeY * (nCellsX + 1) + iCellX + 1];
      double y1 = nodesY[iNodeY * (nCellsX + 1) + iCellX];
      double y2 = nodesY[iNodeY * (nCellsX + 1) + iCellX + 1];

      for (int iGP = 0; iGP < N + 1; iGP++) {

        int index = iNodeY * nCellsX * (N + 1) + iCellX * (N + 1) + iGP;

        EtaBoundX[index] = x1 + 0.5 * (xGP_ref[iGP] + 1.0) * (x2 - x1);
        EtaBoundY[index] = y1 + 0.5 * (xGP_ref[iGP] + 1.0) * (y2 - y1);

        EtaBoundNormalX[index] = (y1 - y2) / sqrt((y1 - y2) * (y1 - y2) + (x2 - x1) * (x2 - x1));
        EtaBoundNormalY[index] = (x2 - x1) / sqrt((y1 - y2) * (y1 - y2) + (x2 - x1) * (x2 - x1));
      }
    }
  }

  free(xGP_ref);

}

void createCartesianNodes(const double xMin, const double xMax,
                          const double yMin, const double yMax,
                          const int nCellsX, const int nCellsY,
                          double** restrict nodesX, double** restrict nodesY) {

  // allocate memory for node-coordinates
  *nodesX = malloc((nCellsX + 1) * (nCellsY + 1) * sizeof(double));
  *nodesY = malloc((nCellsX + 1) * (nCellsY + 1) * sizeof(double));

  double dx = (xMax - xMin) / nCellsX;
  double dy = (yMax - yMin) / nCellsY;

  for (int iNodeY = 0; iNodeY < nCellsY + 1; iNodeY++) {
    for (int iNodeX = 0; iNodeX < nCellsX + 1; iNodeX++) {
      (*nodesX)[iNodeY * (nCellsX + 1) + iNodeX] = xMin + dx * iNodeX;
      (*nodesY)[iNodeY * (nCellsX + 1) + iNodeX] = yMin + dy * iNodeY;
    }
  }

}

void quadMapMetrics(const double x1, const double y1,
                    const double x2, const double y2,
                    const double x3, const double y3,
                    const double x4, const double y4,
                    const double Xi, const double Eta,
                    double* const restrict dXdXi, double* const restrict dXdEta,
                    double* const restrict dYdXi, double* const restrict dYdEta) {

  *dXdXi = 0.25 * ((1.0 - Eta) * (x2 - x1) + (1.0 + Eta) * (x3 - x4));
  *dXdEta = 0.25 * ((1.0 - Xi) * (x4 - x1) + (1.0 + Xi) * (x3 - x2));
  *dYdXi = 0.25 * ((1.0 - Eta) * (y2 - y1) + (1.0 + Eta) * (y3 - y4));
  *dYdEta = 0.25 * ((1.0 - Xi) * (y4 - y1) + (1.0 + Xi) * (y3 - y2));
}

int sign(const double val) {
  return (0.0 < val) - (val < 0.0);
}
