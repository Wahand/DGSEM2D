#ifndef LSERK_INCLUDED
#define LSERK_INCLUDED

#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "DGSEM2D/dg.h"

void timeStepByLSERK(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                     const double t, const double dt,
                     const double* const restrict RKa, const double* const restrict RKb, const double* const restrict RKc,
                     const int N, const int nCellsX, const int nCellsY,
                     const int BCtop, const int BCbottom,
                     const int BCleft, const int BCright,
                     const double* const restrict xGP, const double* const restrict yGP,
                     double* const restrict U, double* const restrict Ut, double* const restrict dU, const double* const restrict sJ,
                     double* const restrict UXiMinus, double* const restrict UXiPlus,
                     double* const restrict UEtaMinus, double* const restrict UEtaPlus,
                     #ifdef parabol
                       double* const restrict gradUx, double* const restrict gradUy,
                       double* const restrict gradUxXiMinus, double* const restrict gradUxXiPlus,
                       double* const restrict gradUxEtaMinus, double* const restrict gradUxEtaPlus,
                       double* const restrict gradUyXiMinus, double* const restrict gradUyXiPlus,
                       double* const restrict gradUyEtaMinus, double* const restrict gradUyEtaPlus,
                     #endif
                     const double* const restrict dXdXi, const double* const restrict dXdEta,
                     const double* const restrict dYdXi, const double* const restrict dYdEta,
                     const double* const restrict XiBoundNormalX, const double* const restrict XiBoundNormalY,
                     const double* const restrict EtaBoundNormalX, const double* const restrict EtaBoundNormalY,
                     const double* const restrict XiBoundMetricsMinus, const double* const restrict XiBoundMetricsPlus,
                     const double* const restrict EtaBoundMetricsMinus, const double* const restrict EtaBoundMetricsPlus,
                     const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus,
                     const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus,
                     double* const restrict innerFluxXi, double* const restrict innerFluxEta,
                     double* const restrict sideFluxXi, double* const restrict sideFluxEta,
                     const double* const restrict LMinus, const double* const restrict LPlus,
                     const double* const restrict LHatMinus, const double* const restrict LHatPlus,
                     const double* const restrict DHat);

void initLSERK(double* const restrict RKa, double* const restrict RKb, double* const restrict RKc);

#endif
