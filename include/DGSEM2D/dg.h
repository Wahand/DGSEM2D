#ifndef DG_INCLUDED
#define DG_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "DGSEM2D/equationFlux.h"
#include "DGSEM2D/riemannSolver.h"
#include "DGSEM2D/interpolation.h"
#include "DGSEM2D/mpiCommunication.h"
#include "DGSEM2D/output.h"

#ifdef parabol
  #include "src/dg/br1/br1.h"
#endif

void fillU(const char* iniState, const double* const restrict xGP, const int sizeX, const double* const restrict yGP, const int sizeY, double* const restrict U);

void DGTimeDerivativeWeakForm(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                              const int N, const int nCellsX, const int nCellsY,
                              const int BCtop, const int BCbottom,
                              const int BCleft, const int BCright,
                              const double* const restrict xGP, const double* const restrict yGP,
                              const double* const restrict U, double* const restrict Ut, const double* const restrict sJ,
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

void innerFlux(const int N, const int nCellsX, const int nCellsY,
               const double* const restrict U,
               #ifdef parabol
                 const double* const restrict gradUx, const double* const restrict gradUy,
               #endif
               const double* const restrict dXdXi, const double* const restrict dXdEta,
               const double* const restrict dYdXi, const double* const restrict dYdEta,
               double* restrict const innerFluxXi, double* const restrict innerFluxEta);

void volumeIntegral(const int N, const int nCellsX, const int nCellsY, double* const restrict Ut, const double* const restrict innerFluxXi, const double* const restrict innerFluxEta, const double* const restrict DHat);

void surfaceIntegral(int N, int nCellsX, int nCellsY, double* const restrict Ut, const double* const restrict LHatPlus, const double* const restrict LHatMinus, const double* const restrict sideFluxXi, const double* const restrict sideFluxEta, const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus, const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus);

void sideFlux(const int N, const int nCellsX, const int nCellsY,
              const double* const restrict UXiMinus, const double* const restrict UXiPlus,
              const double* const restrict UEtaMinus, const double* const restrict UEtaPlus,
              #ifdef parabol
                const double* const restrict gradUxXiMinus, const double* const restrict gradUxXiPlus,
                const double* const restrict gradUxEtaMinus, const double* const restrict gradUxEtaPlus,
                const double* const restrict gradUyXiMinus, const double* const restrict gradUyXiPlus,
                const double* const restrict gradUyEtaMinus, const double* const restrict gradUyEtaPlus,
              #endif
              const double* const restrict XiBoundNormalX, const double* const restrict XiBoundNormalY,
              const double* const restrict EtaBoundNormalX, const double* const restrict EtaBoundNormalY,
              double* const restrict sideFluxXi, double* const restrict sideFluxEta);

#endif
