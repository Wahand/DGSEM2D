#include "DGSEM2D/br1.h"

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
             const double* const restrict DHat) {

  // clear variables
  for (int i = 0; i < nVars * (N + 1) * nCellsX * (N + 1) * nCellsY; i++) {
    gradUx[i] = 0.0;
    gradUy[i] = 0.0;
  }
  for (int i = 0; i < nVars * (nCellsY + 1) * (N + 1) * nCellsX; i++) {
    gradUxEtaMinus[i] = 0.0;
    gradUxEtaPlus[i] = 0.0;
    gradUyEtaMinus[i] = 0.0;
    gradUyEtaPlus[i] = 0.0;
  }
  for (int i = 0; i < nVars * (N + 1) * nCellsY * (nCellsX + 1); i++) {
    gradUxXiMinus[i] = 0.0;
    gradUxXiPlus[i] = 0.0;
    gradUyXiMinus[i] = 0.0;
    gradUyXiPlus[i] = 0.0;
  }

writeMatrix("ULifting.out", nCellsX * (N + 1), nCellsY * (N + 1), U);
  // volume integral
  liftingVolumeIntegral(N, nCellsX, nCellsY,
                        U, gradUx, gradUy,
                        dXdXi, dXdEta,
                        dYdXi, dYdEta,
                        DHat);

writeMatrix("gradUxAfterVolInt.out", nCellsX * (N + 1), nCellsY * (N + 1), gradUx);
  // Receive U (-> tagOffset = 0)
  receiveMPIData(cartComm, neighbourRanks,
                 N, nCellsX, nCellsY,
                 BCtop, BCbottom,
                 BCleft, BCright,
                 UXiMinus, UXiPlus,
                 UEtaMinus, UEtaPlus,
                 0);
  MPI_Waitall(8, reqHandles, MPI_STATUSES_IGNORE);

  // UPlus / Minus are already set
  // calculate fluxes on sides
  for (int iBoundX = 0; iBoundX < nCellsX + 1; iBoundX++) {
    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iGPY = 0; iGPY < N + 1; iGPY++) {
        sideFluxXi[iBoundX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] = 0.5 *
          (UXiPlus[iBoundX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] +
          UXiMinus[iBoundX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY]);
      }
    }
  }
  for (int iBoundY = 0; iBoundY < nCellsY + 1; iBoundY++) {
    for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
      for (int iGPX = 0; iGPX < N + 1; iGPX++) {
        sideFluxEta[iBoundY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] = 0.5 *
          (UEtaPlus[iBoundY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +
          UEtaMinus[iBoundY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX]);
      }
    }
  }
printf("UXiPlus:\n");
for (int i = 0; i < (nCellsX + 1) * nCellsY * (N + 1); i++) {
  printf("%4.2f  ", UXiPlus[i]);
}
printf("\n");
printf("UXiMinus:\n");
for (int i = 0; i < (nCellsX + 1) * nCellsY * (N + 1); i++) {
  printf("%4.2f  ", UXiMinus[i]);
}
printf("\n");
printf("SideFluxXi:\n");
for (int i = 0; i < (nCellsX + 1) * nCellsY * (N + 1); i++) {
  printf("%4.2f  ", sideFluxXi[i]);
}
printf("\n\n");
  // surface integral
  liftingSurfaceIntegral(N, nCellsX, nCellsY,
                         gradUx, gradUy,
                         LHatMinus, LHatPlus,
                         sideFluxXi, sideFluxEta,
                         XiBoundNormalX, XiBoundNormalY,
                         EtaBoundNormalX, EtaBoundNormalY,
                         XiBoundMetricsMinus, XiBoundMetricsPlus,
                         EtaBoundMetricsMinus, EtaBoundMetricsPlus,
                         XiBoundScalingMinus, XiBoundScalingPlus,
                         EtaBoundScalingMinus, EtaBoundScalingPlus);
writeMatrix("gradUxAfterSurfInt.out", nCellsX * (N + 1), nCellsY * (N + 1), gradUx);

  // transform gradients to physical space
  applyJacobianToPhysical(gradUx, sJ, nCellsX * nCellsY * (N + 1) * (N + 1));
  applyJacobianToPhysical(gradUy, sJ, nCellsX * nCellsY * (N + 1) * (N + 1));

  // prolong gradients to sides
//TODO
/*
  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
      for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
        for (int iGPX = 0; iGPX < N + 1; iGPX++) {
          // prolong to EtaPlus
          prolong(&gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPX],
                  &gradUxEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX],
                  LPlus, N + 1, nCellsX * (N + 1));
          prolong(&gradUy[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPX],
                  &gradUyEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX],
                  LPlus, N + 1, nCellsX * (N + 1));
          // prolong to EtaMinus
          prolong(&gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPX],
                  &gradUxEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX],
                  LMinus, N + 1, nCellsX * (N + 1));
          prolong(&gradUy[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPX],
                  &gradUyEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX],
                  LMinus, N + 1, nCellsX * (N + 1));
        }
        for (int iGPY = 0; iGPY < N + 1; iGPY++) {
          // prolong to XiPlus
          prolong(&gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPY * nCellsX * (N + 1)],
                  &gradUxXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + (iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY],
                  LPlus, N + 1, 1);
          prolong(&gradUy[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPY * nCellsX * (N + 1)],
                  &gradUyXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + (iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY],
                  LPlus, N + 1, 1);
          // prolong to XiMinus
          prolong(&gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPY * nCellsX * (N + 1)],
                  &gradUxXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY],
                  LMinus, N + 1, 1);
          prolong(&gradUy[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + (iCellX * (N + 1)) + (iCellY * nCellsX * (N + 1) * (N + 1)) + iGPY * nCellsX * (N + 1)],
                  &gradUyXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY],
                  LMinus, N + 1, 1);
        }
      }
    }
  }
*/

  // update boundary conditions for gradients
  if (BCtop == 1) {
    // von-Neumann boundaries
    double refGradX[nVars];
    double refGradY[nVars];
    memset(refGradX, 0.0, sizeof(refGradX));
    memset(refGradY, 0.0, sizeof(refGradY));
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iXi = 0; iXi < nCellsX * (N + 1); iXi++) {
        gradUxEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi] = refGradX[iVar];
        gradUyEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi] = refGradY[iVar];
      }
    }
  } else {
    // periodic boundaries
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iXi = 0; iXi < nCellsX * (N + 1); iXi++) {
        gradUxEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi] = gradUxEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi];
        gradUyEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi] = gradUyEtaMinus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi];
      }
    }
  }

  if (BCbottom == 1) {
    // von-Neumann boundaries
    double refGradX[nVars];
    double refGradY[nVars];
    memset(refGradX, 0.0, sizeof(refGradX));
    memset(refGradY, 0.0, sizeof(refGradY));
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iXi = 0; iXi < nCellsX * (N + 1); iXi++) {
        gradUxEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi] = refGradX[iVar];
        gradUyEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi] = refGradY[iVar];
      }
    }
  } else {
    // periodic boundaries
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iXi = 0; iXi < nCellsX * (N + 1); iXi++) {
        gradUxEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi] = gradUxEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi];
        gradUyEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iXi] = gradUyEtaPlus[iVar * (nCellsY + 1) * nCellsX * (N + 1) + nCellsY * nCellsX * (N + 1) + iXi];
      }
    }
  }

  if (BCleft == 1) {
    // von-Neumann boundaries
    double refGradX[nVars];
    double refGradY[nVars];
    memset(refGradX, 0.0, sizeof(refGradX));
    memset(refGradY, 0.0, sizeof(refGradY));
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iEta = 0; iEta < nCellsY * (N + 1); iEta++) {
        gradUxXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta] = refGradX[iVar];
        gradUyXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta] = refGradY[iVar];
      }
    }
  } else {
    // periodic boundaries
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iEta = 0; iEta < nCellsY * (N + 1); iEta++) {
        gradUxXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta] = gradUxXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta];
        gradUyXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta] = gradUyXiPlus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta];
      }
    }
  }

  if (BCright == 1) {
    // von-Neumann boundaries
    double refGradX[nVars];
    double refGradY[nVars];
    memset(refGradX, 0.0, sizeof(refGradX));
    memset(refGradY, 0.0, sizeof(refGradY));
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iEta = 0; iEta < nCellsY * (N + 1); iEta++) {
        gradUxXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta] = refGradX[iVar];
        gradUyXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta] = refGradY[iVar];
      }
    }
  } else {
    // periodic boundaries
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iEta = 0; iEta < nCellsY * (N + 1); iEta++) {
        gradUxXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta] = gradUxXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta];
        gradUyXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + nCellsX * nCellsY * (N + 1) + iEta] = gradUyXiMinus[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iEta];
      }
    }
  }

  // Send gradU
  sendMPIData(cartComm, neighbourRanks,
              N, nCellsX, nCellsY,
              BCtop, BCbottom,
              BCleft, BCright,
              gradUxXiMinus, gradUxXiPlus,
              gradUxEtaMinus, gradUxEtaPlus,
              reqHandles, 1);
  sendMPIData(cartComm, neighbourRanks,
              N, nCellsX, nCellsY,
              BCtop, BCbottom,
              BCleft, BCright,
              gradUyXiMinus, gradUyXiPlus,
              gradUyEtaMinus, gradUyEtaPlus,
              &reqHandles[3], 2);
writeMatrix("gradUxAfterEnd.out", nCellsX * (N + 1), nCellsY * (N + 1), gradUx);
}

void liftingSurfaceIntegral(const int N, const int nCellsX, const int nCellsY,
                            double* const restrict gradUx, double* const restrict gradUy,
                            const double* const restrict LHatMinus, const double* const restrict LHatPlus,
                            const double* const restrict sideFluxXi, const double* const restrict sideFluxEta,
                            const double* const restrict XiBoundNormalX, const double* const restrict XiBoundNormalY,
                            const double* const restrict EtaBoundNormalX, const double* const restrict EtaBoundNormalY,
                            const double* const restrict XiBoundMetricsMinus, const double* const restrict XiBoundMetricsPlus,
                            const double* const restrict EtaBoundMetricsMinus, const double* const restrict EtaBoundMetricsPlus,
                            const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus,
                            const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus) {

  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
      for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
        for (int iGPX = 0; iGPX < N + 1; iGPX++) {
          for (int iGPY = 0; iGPY < N + 1; iGPY++) {
/*
            gradUx[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
             // flux on xi plus
             LHatPlus[iGPX] * XiBoundScalingPlus[iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + (iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] -
             // flux on xi minus
             LHatMinus[iGPX] * XiBoundScalingMinus[iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] +
             LHatPlus[iGPY] * EtaBoundScalingPlus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] - 
             LHatMinus[iGPY] * EtaBoundScalingMinus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX];
*/
/*
            gradUx[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
             // flux on xi plus
             LHatPlus[iGPX] * XiBoundNormalX[(iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + (iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] -
             // flux on xi minus
             LHatMinus[iGPX] * XiBoundNormalX[iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] +
             // flux on eta plus
             LHatPlus[iGPY] * EtaBoundNormalX[(iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] - 
             // flux on eta minus
             LHatMinus[iGPY] * EtaBoundNormalX[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX];
*/
            gradUx[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
             LHatPlus[iGPX] * XiBoundMetricsPlus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 3] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + (iCellX + 1) * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] -
             LHatMinus[iGPX] * XiBoundMetricsMinus[iCellX * nCellsY * (N + 1) * 4 + iCellY * (N + 1) * 4 + iGPY * 4 + 3] * sideFluxXi[iVar * (nCellsX + 1) * nCellsY * (N + 1) + iCellX * nCellsY * (N + 1) + iCellY * (N + 1) + iGPY] -
             LHatPlus[iGPY] * EtaBoundMetricsPlus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 2] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +
             LHatMinus[iGPY] * EtaBoundMetricsMinus[iCellY * nCellsX * (N + 1) * 4 + iCellX * (N + 1) * 4 + iGPX * 4 + 2] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX];

            gradUy[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
             // flux on eta plus
             LHatPlus[iGPY] * EtaBoundScalingPlus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + (iCellY + 1) * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] -
             // flux on eta minus
             LHatMinus[iGPY] * EtaBoundScalingMinus[iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] * sideFluxEta[iVar * (nCellsY + 1) * nCellsX * (N + 1) + iCellY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX];
          }
        }
      }
    }
  }
}

void liftingVolumeIntegral(const int N, const int nCellsX, const int nCellsY,
                           const double* const restrict U, double* const restrict gradUx, double* const restrict gradUy,
                           const double* const restrict dXdXi, const double* const restrict dXdEta,
                           const double* const restrict dYdXi, const double* const restrict dYdEta,
                           const double* const restrict DHat) {

  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iGPY = 0; iGPY < N + 1; iGPY++) {
        for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
          for (int iGPX = 0; iGPX < N + 1; iGPX++) {
            for (int j = 0; j < N + 1; j++) {
/*
Ut +=
innerFluxXi[j, iGPY] * DHat[iGPX, j] +
innerFluxEta[iGPX, j] * DHat[iGPY, j]

innerFluxXi = dYdETa * Fx - dXdEta * Fy
innerFluxEta = dXdXI * Fy - dYdXi * Fx

Ut +=
dYdEta[j, iGPY] * U[j, iGPY] * DHat[iGPX, j] -
dYdXi[iGPX, j] * U[iGPX, j] * DHat[iGPY, j]
printf("index: %d\n", iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX);
gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
dYdEta[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
DHat[iGPX * (N + 1) + j] -
dYdXi[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
DHat[iGPY * (N + 1) + j];
*/
/*
              gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
              (dYdEta[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] -
               dYdXi[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX]) *
               DHat[iGPX * (N + 1) + j];
*/
              gradUx[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
               dYdEta[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
               DHat[iGPX * (N + 1) + j] -
               dYdXi[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
               DHat[iGPY * (N + 1) + j];

              gradUy[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] +=
              (dXdXi[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + j * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX] -
               dYdXi[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j] *
               U[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + iCellY * (N + 1) * nCellsX * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + j]) *
               DHat[iGPY * (N + 1) + j];
            }
          }
        }
      }
    }
  }
}
