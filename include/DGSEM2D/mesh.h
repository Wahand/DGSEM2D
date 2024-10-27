#ifndef MESH_INCLUDED
#define MESH_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "DGSEM2D/basis.h"

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
                double* const restrict sJ);

void createCartesianNodes(const double xMin, const double xMax,
                          const double yMin, const double yMax,
                          const int nCellsX, const int nCellsY,
                          double** restrict nodesX, double** restrict nodesY);

//! Kopriva algorithm 100
/*!
  Computation of the metric terms on a straight sided quadrilateral
*/
void quadMapMetrics(const double x1, const double y1,
                    const double x2, const double y2,
                    const double x3, const double y3,
                    const double x4, const double y4,
                    const double Xi, const double Eta,
                    double* const restrict dXdXi, double* const restrict dXdEta,
                    double* const restrict dYdXi, double* const restrict dYdEta);

int sign(const double val);

#endif
