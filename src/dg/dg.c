#include "DGSEM2D/dg.h"

void fillU(const char* iniState, const double* restrict xGP, const int sizeX, const double* restrict yGP, const int sizeY, double* restrict U) {

  double* URef = malloc(nVars * sizeof(double));
  #if (eqSys == 1)
    URef[0] = 1.0;
  #endif
  #if (eqSys == 2)
    const double gamma = 1.4;

    double UPrim[nVars];
    UPrim[0] = 1.0;
    UPrim[1] = 0.0;
    UPrim[2] = 0.0;
    UPrim[3] = 1.0;

    URef[0] = UPrim[0];
    URef[1] = UPrim[1] * UPrim[0];
    URef[2] = UPrim[2] * UPrim[0];
    URef[3] = UPrim[3] / (gamma - 1.0) + 0.5 * (UPrim[0] * UPrim[1] * UPrim[1] + UPrim[0] * UPrim[2] * UPrim[2]);
  #endif
  #if (eqSys == 3)
    URef[0] = 1.0;
  #endif


  // linear advection
  #if (eqSys == 1)
    if (strcmp(iniState, "1") == 0) {
      for (int iY = 0; iY < sizeY; iY++) {
        for (int iX = 0; iX < sizeX; iX++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            U[iY * sizeX * nVars + iX * nVars + iVar] = URef[iVar];
          }
        }
      }
    } else if (strcmp(iniState, "2") == 0) {
      double xp = 0.0;
      double yp = 0.0;
      double r = 0.8;
      double l = 1.0;
      double u = 2.0;
      for (int iY = 0; iY < sizeY; iY++) {
        for (int iX = 0; iX < sizeX; iX++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            U[iY*sizeX + iX] = l + (u-l) * exp(-pow(sqrt(pow(xGP[iY*sizeX+iX]-xp,2) + pow(yGP[iY*sizeX+iX]-yp,2)), 2) / (2*pow(r,2)));
          }
        }
      }
    }
  #endif

  // euler
  #if (eqSys == 2)
//TODO parse iniState
    if (!strcmp(iniState, "0")) {
      for (int iVar = 0; iVar < nVars; iVar++) {
        for (int iX = 0; iX < sizeX * sizeY; iX++) {
          U[iVar * sizeX * sizeY + iX] = iVar + 1;
        }
      }
    } else if (!strcmp(iniState, "3")) {
      for (int iX = 0; iX < sizeX; iX++) {
        for (int iY = 0; iY < sizeY; iY++) {
          // rho
          U[iY * sizeX * nVars + iX * nVars] = 1.0 + 0.5 * exp(-(0.5 * (xGP[iY * sizeX + iX] * xGP[iY * sizeX + iX] + yGP[iY * sizeX + iX] * yGP[iY * sizeX + iX])) / (2 * 0.2 * 0.2));
          // v1
          U[iY * sizeX * nVars + iX * nVars + 1] = URef[1];
          // v2
          U[iY * sizeX * nVars + iX * nVars + 2] = URef[2];
          // p
          U[iY * sizeX * nVars + iX * nVars + 3] = URef[3];
        }
      }
    } else if (!strcmp(iniState, "1")) {
      for (int iX = 0; iX < sizeX; iX++) {
        for (int iY = 0; iY < sizeY; iY++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            U[iY * sizeX * nVars + iX * nVars + iVar] = URef[iVar];
          }
        }
      }
    } else if (!strcmp(iniState, "2")) {
      for (int iX = 0; iX < sizeX; iX++) {
        for (int iY = 0; iY < sizeY; iY++) {
          U[iY * sizeX * nVars + iX * nVars] = 1.0;
          U[iY * sizeX * nVars + iX * nVars + 1] = 0.0;
          U[iY * sizeX * nVars + iX * nVars + 2] = 0.0;
          U[iY * sizeX * nVars + iX * nVars + 3] = 2.5 + 0.5 * exp(-(0.5 * (xGP[iY * sizeX + iX] * xGP[iY * sizeX + iX] + yGP[iY * sizeX + iX] * yGP[iY * sizeX + iX])) / (2 * 0.2 * 0.2));
        }
      }
    }
  #endif

  // advection-diffusion
  #if (eqSys == 3)
//TODO parse iniState
    if (!strcmp(iniState, "1")) {
      for (int iVar = 0; iVar < nVars; iVar++) {
        for (int iX = 0; iX < sizeX; iX++) {
          for (int iY = 0; iY < sizeY; iY++) {
            U[iY * sizeX * nVars + iX * nVars + iVar] = URef[iVar];
          }
        }
      }
    } else if (!strcmp(iniState, "2")) {
      for (int iVar = 0; iVar < nVars; iVar++) {
        for (int iX = 0; iX < sizeX; iX++) {
          for (int iY = 0; iY < sizeY; iY++) {
            U[iY * sizeX * nVars + iX * nVars + iVar] = 1.0 + 0.5 * exp(-(0.5 * (xGP[iY * sizeX + iX] * xGP[iY * sizeX + iX] + yGP[iY * sizeX + iX] * yGP[iY * sizeX + iX])) / (2 * 0.2 * 0.2));
          }
        }
      }
    } else if (!strcmp(iniState, "3")) {
      for (int iVar = 0; iVar < nVars; iVar++) {
        for (int iX = 0; iX < sizeX; iX++) {
          for (int iY = 0; iY < sizeY; iY++) {
            U[iY * sizeX * nVars + iX * nVars + iVar] = 1.0 + 0.2 * xGP[iY * sizeX + iX];
          }
        }
      }
    }
  #endif

  free(URef);
}

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
                              const double* const restrict DHat) {

  int Np1 = N + 1;

  //TODO myRank, cartCommDims, cartCommCoords etc. überall übergeben?
  // get information from communicator
  int myRank, numProcs;
  int cartCommDims[2];
  int cartCommPeriodicity[2];
  int myCoords[2];
  MPI_Comm_rank(*cartComm, &myRank);
  MPI_Comm_size(*cartComm, &numProcs);
  MPI_Cart_get(*cartComm, 2, cartCommDims, cartCommPeriodicity, myCoords);


  #if (outputLevel >= 3)
    printf("Calling DGTimeDerivativeWeakForm from processor %d...\n", myRank);
  #endif

  // clear variables
  for (int i = 0; i < nVars * Np1 * nCellsX * Np1 * nCellsY; i++) {
    Ut[i] = 0.0;
  }
  for (int i = 0; i < nVars * (nCellsY + 1) * Np1 * nCellsX; i++) {
    sideFluxEta[i] = 0.0;
    UEtaMinus[i] = 0.0;
    UEtaPlus[i] = 0.0;
  }
  for (int i = 0; i < nVars * Np1 * nCellsY * (nCellsX + 1); i++) {
    sideFluxXi[i] = 0.0;
    UXiMinus[i] = 0.0;
    UXiPlus[i] = 0.0;
  }


  // prolong variables
  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    for (int iGPY = 0; iGPY < Np1; iGPY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPX = 0; iGPX < Np1; iGPX++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            double currentU = U[iCellY * Np1 * nCellsX * Np1 * nVars + iGPY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar];
            UEtaPlus[(iCellY + 1) * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar] +=
             currentU * LPlus[iGPY];
            UEtaMinus[iCellY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar] +=
             currentU * LMinus[iGPY];
            UXiPlus[(iCellX + 1) * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars + iVar] +=
             currentU * LPlus[iGPX];
            UXiMinus[iCellX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars + iVar] +=
             currentU * LMinus[iGPX];
          }
        }
      }
    }
  }

  //TODO fill BC-data (non-MPI)
  if (myCoords[0] == cartCommDims[0] - 1) {
    if (BCtop == 1) {
      //TODO dirichlet-BC
    } else if (BCtop == 2) {
      // periodic BC
      if (cartCommDims[0] == 1) {
        // Copy values for periodic boundaries if only one processor in Y-direction
        for (int iXi = 0; iXi < nCellsX * Np1; iXi++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            UEtaMinus[nCellsY * nCellsX * Np1 * nVars + iXi * nVars + iVar] = UEtaMinus[iXi * nVars + iVar];
          }
        }
      }
    } else if (BCtop == 3) {
      // wall-BC
      #if (eqSys == 2)
        // Euler
        for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
          for (int iGPX = 0; iGPX < Np1; iGPX++) {
            int UIndex = nCellsX * nCellsY * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars;
            int normalIndex = nCellsX * nCellsY * Np1 + iCellX * Np1 + iGPX;
            UEtaMinus[UIndex] = UEtaPlus[UIndex];
            UEtaMinus[UIndex + 1] = -fabs(EtaBoundNormalX[normalIndex]) * UEtaPlus[UIndex + 1];
            UEtaMinus[UIndex + 2] = -fabs(EtaBoundNormalY[normalIndex]) * UEtaPlus[UIndex + 2];
            UEtaMinus[UIndex + 3] = UEtaPlus[UIndex + 3];
          }
        }
      #endif
    }
  }

  // bottom
  if (myCoords[0] == 0) {
    if (BCbottom == 1) {
      //TODO dirichlet-BC
    } else if (BCbottom == 2) {
      if (cartCommDims[0] == 1) {
        for (int iXi = 0; iXi < nCellsX * Np1; iXi++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            UEtaPlus[iXi * nVars + iVar] = UEtaPlus[nCellsY * nCellsX * Np1 * nVars + iXi * nVars + iVar];
          }
        }
      }
    } else if (BCbottom == 3) {
      #if (eqSys == 2)
        for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
          for (int iGPX = 0; iGPX < Np1; iGPX++) {
            int UIndex = iCellX * Np1 * nVars + iGPX * nVars;
            int normalIndex = iCellX * Np1 + iGPX;
            UEtaPlus[UIndex] = UEtaMinus[UIndex];
            UEtaPlus[UIndex + 1] = -fabs(EtaBoundNormalX[normalIndex]) * UEtaMinus[UIndex + 1];
            UEtaPlus[UIndex + 2] = -fabs(EtaBoundNormalY[normalIndex]) * UEtaMinus[UIndex + 2];
            UEtaPlus[UIndex + 3] = UEtaMinus[UIndex + 3];
          }
        }
      #endif
    }
  }
 
  // left
  if (myCoords[1] == 0) {
    if (BCleft == 1) {
      //TODO dirichlet-BC
    } else if (BCleft == 2) {
      if (cartCommDims[1] == 1) {
        for (int iEta = 0; iEta < nCellsY * Np1; iEta++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            UXiPlus[iEta * nVars + iVar] = UXiPlus[nCellsX * nCellsY * Np1 * nVars + iEta * nVars + iVar];
          }
        }
      }
    } else if (BCleft == 3) {
      #if (eqSys == 2)
        for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
          for (int iGPY = 0; iGPY < Np1; iGPY++) {
            int UIndex = iCellY * Np1 * nVars + iGPY * nVars;
            int normalIndex = iCellY * Np1 + iGPY;
            UXiPlus[UIndex] = UXiMinus[UIndex];
            UXiPlus[UIndex + 1] = -fabs(XiBoundNormalX[normalIndex]) * UXiMinus[UIndex + 1];
            UXiPlus[UIndex + 2] = -fabs(XiBoundNormalY[normalIndex]) * UXiMinus[UIndex + 2];
            UXiPlus[UIndex + 3] = UXiMinus[UIndex + 3];
          }
        }
      #endif
    }
  }

  // right
  if (myCoords[1] == cartCommDims[1] - 1) {
    if (BCright == 1) {
      //TODO dirichlet-BC
    } else if (BCright == 2) {
      if (cartCommDims[1] == 1) {
        for (int iEta = 0; iEta < nCellsY * Np1; iEta++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            UXiMinus[nCellsX * nCellsY * Np1 * nVars + iEta * nVars + iVar] = UXiMinus[iEta * nVars + iVar];
          }
        }
      }
    } else if (BCright == 3) {
      #if (eqSys == 2)
        for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
          for (int iGPY = 0; iGPY < Np1; iGPY++) {
            int UIndex = nCellsX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars;
            int normalIndex = nCellsX * nCellsY * Np1 + iCellY * Np1 + iGPY;
            UXiMinus[UIndex] = UXiPlus[UIndex];
            UXiMinus[UIndex + 1] = -fabs(XiBoundNormalX[normalIndex]) * UXiPlus[UIndex + 1];
            UXiMinus[UIndex + 2] = -fabs(XiBoundNormalY[normalIndex]) * UXiPlus[UIndex + 2];
            UXiMinus[UIndex + 3] = UXiPlus[UIndex + 3];
          }
        }
      #endif
    }
  }

  // Request handles
  MPI_Request reqHandles[8];
  for (int i = 0; i < 8; i++) {
    reqHandles[i] = MPI_REQUEST_NULL;
  }

  // Send U (-> tagOffset = 0)
  sendMPIData(cartComm, neighbourRanks,
              N, nCellsX, nCellsY,
              BCtop, BCbottom,
              BCleft, BCright,
              UXiMinus, UXiPlus,
              UEtaMinus, UEtaPlus,
              reqHandles, 0);

  // Lifting
  #ifdef parabol
    lifting(cartComm, neighbourRanks,
            reqHandles,
            N, nCellsX, nCellsY,
            BCtop, BCbottom,
            BCleft, BCright,
            U, sJ,
            UXiMinus, UXiPlus,
            UEtaMinus, UEtaPlus,
            gradUx, gradUy,
            gradUxXiMinus, gradUxXiPlus,
            gradUxEtaMinus, gradUxEtaPlus,
            gradUyXiMinus, gradUyXiPlus,
            gradUyEtaMinus, gradUyEtaPlus,
            dXdXi, dXdEta,
            dYdXi, dYdEta,
            XiBoundNormalX, XiBoundNormalY,
            EtaBoundNormalX, EtaBoundNormalY,
            XiBoundMetricsMinus, XiBoundMetricsPlus,
            EtaBoundMetricsMinus, EtaBoundMetricsPlus,
            XiBoundScalingMinus, XiBoundScalingPlus,
            EtaBoundScalingMinus, EtaBoundScalingPlus,
            sideFluxXi, sideFluxEta,
            LMinus, LPlus,
            LHatMinus, LHatPlus,
            DHat);

    // clear side-flux variables again
    for (int i = 0; i < nVars * (nCellsY + 1) * Np1 * nCellsX; i++) {
      sideFluxEta[i] = 0.0;
    }
    for (int i = 0; i < nVars * Np1 * nCellsY * (nCellsX + 1); i++) {
      sideFluxXi[i] = 0.0;
    }
  #endif

  // fill inner flux on gauss-points
  #ifdef parabol
    innerFlux(N, nCellsX, nCellsY, U, gradUx, gradUy, dXdXi, dXdEta, dYdXi, dYdEta, innerFluxXi, innerFluxEta);
  #else
    innerFlux(N, nCellsX, nCellsY, U, dXdXi, dXdEta, dYdXi, dYdEta, innerFluxXi, innerFluxEta);
  #endif

  // volume integral
  volumeIntegral(N, nCellsX, nCellsY, Ut, innerFluxXi, innerFluxEta, DHat);
//writeMatrix("UtAfterVolInt.out", nCellsX * (N + 1), nCellsY * (N + 1), Ut);

  // Receive data (gradU if parabol, U otherwise)
  #ifdef parabol
    // Receive gradU (-> tagOffset = 1/2)
    receiveMPIData(cartComm, neighbourRanks,
                   N, nCellsX, nCellsY,
                   BCtop, BCbottom,
                   BCleft, BCright,
                   gradUxXiMinus, gradUxXiPlus,
                   gradUxEtaMinus, gradUxEtaPlus,
                   1);
    receiveMPIData(cartComm, neighbourRanks,
                   N, nCellsX, nCellsY,
                   BCtop, BCbottom,
                   BCleft, BCright,
                   gradUyXiMinus, gradUyXiPlus,
                   gradUyEtaMinus, gradUyEtaPlus,
                   2);
  #else
    // Receive U (-> tagOffset = 0)
    receiveMPIData(cartComm, neighbourRanks,
                   N, nCellsX, nCellsY,
                   BCtop, BCbottom,
                   BCleft, BCright,
                   UXiMinus, UXiPlus,
                   UEtaMinus, UEtaPlus,
                   0);
  #endif

  // Make sure sending is complete
  MPI_Waitall(8, reqHandles, MPI_STATUSES_IGNORE);

  // flux on sides
  sideFlux(N, nCellsX, nCellsY,
           UXiMinus, UXiPlus,
           UEtaMinus, UEtaPlus,
           #ifdef parabol
             gradUxXiMinus, gradUxXiPlus,
             gradUxEtaMinus, gradUxEtaPlus,
             gradUyXiMinus, gradUyXiPlus,
             gradUyEtaMinus, gradUyEtaPlus,
           #endif
           XiBoundNormalX, XiBoundNormalY,
           EtaBoundNormalX, EtaBoundNormalY,
           sideFluxXi, sideFluxEta);

  // surface integral
  surfaceIntegral(N, nCellsX, nCellsY, Ut, LHatPlus, LHatMinus, sideFluxXi, sideFluxEta, XiBoundScalingMinus, XiBoundScalingPlus, EtaBoundScalingMinus, EtaBoundScalingPlus);
//writeMatrix("UtAfterSurfInt.out", nCellsX * (N + 1), nCellsY * (N + 1), Ut);

  // transform Ut to physical space
  applyJacobianToPhysical(Ut, sJ, nCellsX * nCellsY * (N + 1) * (N + 1));
//writeMatrix("UtAfterEnd.out", nCellsX * (N + 1), nCellsY * (N + 1), Ut);
}

void surfaceIntegral(const int N, const int nCellsX, const int nCellsY, double* const restrict Ut, const double* const restrict LHatPlus, const double* const restrict LHatMinus, const double* const restrict sideFluxXi, const double* const restrict sideFluxEta, const double* const restrict XiBoundScalingMinus, const double* const restrict XiBoundScalingPlus, const double* const restrict EtaBoundScalingMinus, const double* const restrict EtaBoundScalingPlus) {
  int Np1 = N + 1;
  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    for (int iGPY = 0; iGPY < Np1; iGPY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPX = 0; iGPX < Np1; iGPX++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            Ut[iCellY * Np1 * nCellsX * Np1 * nVars + iGPY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar] +=
             LHatPlus[iGPX] * XiBoundScalingPlus[iCellX * nCellsY * Np1 + iCellY * Np1 + iGPX] * sideFluxXi[(iCellX + 1) * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars + iVar] -
             LHatMinus[iGPX] * XiBoundScalingMinus[iCellX * nCellsY * Np1 + iCellY * Np1 + iGPX] * sideFluxXi[iCellX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars + iVar] +
             LHatPlus[iGPY] * EtaBoundScalingPlus[iCellY * nCellsX * Np1 + iCellX * Np1 + iGPX] * sideFluxEta[(iCellY + 1) * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar] -
             LHatMinus[iGPY] * EtaBoundScalingMinus[iCellY * nCellsX * Np1 + iCellX * Np1 + iGPX] * sideFluxEta[iCellY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars + iVar];
          }
        }
      }
    }
  }
}

void innerFlux(const int N, const int nCellsX, const int nCellsY,
               const double* const restrict U,
               #ifdef parabol
                 const double* const restrict gradUx, const double* const restrict gradUy,
               #endif
               const double* const restrict dXdXi, const double* const restrict dXdEta,
               const double* const restrict dYdXi, const double* const restrict dYdEta,
               double* const restrict innerFluxXi, double* const restrict innerFluxEta) {

  int Np1 = N + 1;
  double tmpFluxX[nVars];
  double tmpFluxY[nVars];
  for (int j = 0; j < nCellsY * Np1; j++) {
    for (int i = 0; i < nCellsX * Np1; i++) {
      int UIndex = j * nCellsX * Np1 * nVars + i * nVars;
      int index = j * nCellsX * Np1 + i;
      #ifdef parabol
        equationFluxX(&U[UIndex], &gradUx[UIndex], &gradUy[UIndex], tmpFluxX);
        equationFluxY(&U[UIndex], &gradUx[UIndex], &gradUy[UIndex], tmpFluxY);
      #else
        equationFluxX(&U[UIndex], tmpFluxX);
        equationFluxY(&U[UIndex], tmpFluxY);
      #endif
      for (int iVar = 0; iVar < nVars; iVar++) {
        innerFluxXi[UIndex + iVar] = dYdEta[index] * tmpFluxX[iVar] - dXdEta[index] * tmpFluxY[iVar];
        innerFluxEta[UIndex + iVar] = dXdXi[index] * tmpFluxY[iVar] - dYdXi[index] * tmpFluxX[iVar];
      }
    }
  }
}

void volumeIntegral(const int N, const int nCellsX, const int nCellsY, double* const restrict Ut, const double* const restrict innerFluxXi, const double* const restrict innerFluxEta, const double* const restrict DHat) {
  int Np1 = N + 1;
  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    int iCellYOffset = iCellY * Np1 * nCellsX * Np1 * nVars;
    for (int iGPY = 0; iGPY < N + 1; iGPY++) {
      int iGPYOffset = nCellsX * Np1 * nVars;
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        int iCellXOffset = iCellX * Np1 * nVars;
        for (int iGPX = 0; iGPX < N + 1; iGPX++) {
          for (int iVar = 0; iVar < nVars; iVar++) {
            double* currentUt = &Ut[iCellYOffset + iGPY * iGPYOffset + iCellXOffset + iGPX * nVars + iVar];
            for (int j = 0; j < N + 1; j++) {
              *currentUt +=
               innerFluxXi[iCellYOffset + iGPY * iGPYOffset + iCellXOffset + j * nVars + iVar] *
               DHat[iGPX * Np1 + j] +
               innerFluxEta[iCellYOffset + j * iGPYOffset + iCellXOffset + iGPX * nVars + iVar] *
               DHat[iGPY * Np1 + j];
            }
          }
        }
      }
    }
  }
}

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
              double* const restrict sideFluxXi, double* const restrict sideFluxEta) {

  int Np1 = N + 1;

  // flux in xi
  for (int iBoundX = 0; iBoundX < nCellsX + 1; iBoundX++) {
    for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
      for (int iGPY = 0; iGPY < Np1; iGPY++) {
        riemannSolver(N, nCellsX, nCellsY,
                      &UXiPlus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                      &UXiMinus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                      #ifdef parabol
                        &gradUxXiPlus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                        &gradUxXiMinus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                        &gradUyXiPlus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                        &gradUyXiMinus[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars],
                      #endif
                      XiBoundNormalX[iBoundX * nCellsY * Np1 + iCellY * Np1 + iGPY],
                      XiBoundNormalY[iBoundX * nCellsY * Np1 + iCellY * Np1 + iGPY],
                      &sideFluxXi[iBoundX * nCellsY * Np1 * nVars + iCellY * Np1 * nVars + iGPY * nVars]);
      }
    }
  }

  // flux in eta
  for (int iBoundY = 0; iBoundY < nCellsY + 1; iBoundY++) {
    for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
      for (int iGPX = 0; iGPX < N + 1; iGPX++) {
        riemannSolver(N, nCellsX, nCellsY,
                      &UEtaPlus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                      &UEtaMinus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                      #ifdef parabol
                        &gradUxEtaPlus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                        &gradUxEtaMinus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                        &gradUyEtaPlus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                        &gradUyEtaMinus[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars],
                      #endif
                      EtaBoundNormalX[iBoundY * nCellsX * Np1 + iCellX * Np1 + iGPX],
                      EtaBoundNormalY[iBoundY * nCellsX * Np1 + iCellX * Np1 + iGPX],
                      &sideFluxEta[iBoundY * nCellsX * Np1 * nVars + iCellX * Np1 * nVars + iGPX * nVars]);
      }
    }
  }
}
