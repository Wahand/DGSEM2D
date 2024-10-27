#ifndef BR1_INCLUDED
#define BR1_INCLUDED

#include <stdlib.h>
#include "DGSEM2D/interpolation.h"
#include "DGSEM2D/mpiCommunication.h"
#include "DGSEM2D/output.h"

void lifting(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
             MPI_Request* const restrict reqHandles,
             const int N, const int nCellsX, const int nCellsY,
             const int BCtop, const int BCbottom,
             const int BCleft, const int BCright,
             const double* const restrict U, const double* const restrict sJ,
             double* const restrict UXiMinus, double* const restrict UXiPlus,
             double* const restrict UEtaMinus, double* const restrict UEtaPlus,
             double* const restrict gradUx, double* const restrict gradUy,
             double* const restrict gradUxXiMinus, double* const restrict gradUxXiPlus,
             double* const restrict gradUxEtaMinus, double* const restrict gradUxEtaPlus,
             double* const restrict gradUyXiMinus, double* const restrict gradUyXiPlus,
             double* const restrict gradUyEtaMinus, double* const restrict gradUyEtaPlus,
             const double* const restrict dXdXi, const double* const restrict dXdEta,
             const double* const restrict dYdXi, const double* const restrict dYdEta,
             const double* const restrict XiBoundNormalX, const double* const restrict XiBoundNormalY,
             const double* const restrict EtaBoundNormalX, const double* const restrict EtaBoundNormalY,
             const double* const restrict XiBoundMetricsMinus, const double* const restrict XiBoundMetricsPlus,
             const double* const restrict EtaBoundMetricsMinus, const double* const restrict EtaBoundMetricsPlus,
             const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus,
             const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus,
             double* const restrict sideFluxXi, double* const restrict sideFluxEta,
             const double* const restrict LMinus, const double* const restrict LPlus,
             const double* const restrict LHatMinus, const double* const restrict LHatPlus,
             const double* const restrict DHat);

void liftingSurfaceIntegral(const int N, const int nCellsX, const int nCellsY,
                            double* const restrict gradUx, double* const restrict gradUy,
                            const double* const restrict LHatMinus, const double* const restrict LHatPlus,
                            const double* const restrict sideFluxXi, const double* const restrict sideFluxEta,
                            const double* const restrict XiBoundNormalX, const double* const restrict XiBoundNormalY,
                            const double* const restrict EtaBoundNormalX, const double* const restrict EtaBoundNormalY,
                            const double* const restrict XiBoundMetricsMinus, const double* const restrict XiBoundMetricsPlus,
                            const double* const restrict EtaBoundMetricsMinus, const double* const restrict EtaBoundMetricsPlus,
                            const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus,
                            const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus);

void liftingVolumeIntegral(const int N, const int nCellsX, const int nCellsY,
                           const double* const restrict U, double* const restrict gradUx, double* const restrict gradUy,
                           const double* const restrict dXdXi, const double* const restrict dXdEta,
                           const double* const restrict dYdXi, const double* const restrict dYdEta,
                           const double* const restrict DHat);

#endif
