#include "DGSEM2D/LSERK.h"

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
                     const double* const restrict DHat) {

//  double tStage = t;

  // initialize dU with zeros
  for (int i = 0; i < nVars * nCellsX * nCellsY * (N + 1) * (N + 1); i++) {
    dU[i] = 0.0;
  }

  // RK-loop
  for (int iStage = 0; iStage < 5; iStage++) {

    #if (outputLevel >= 3)
      int myRank;
      MPI_Comm_rank(*cartComm, &myRank);
      if (myRank == 0) {
        printf("\nLSERK-stage no %d...\n", iStage + 1);
      }
    #endif

    MPI_Barrier(*cartComm);

    // execute spatial operator
    DGTimeDerivativeWeakForm(cartComm, neighbourRanks,
                             N, nCellsX, nCellsY,
                             BCtop, BCbottom,
                             BCleft, BCright,
                             xGP, yGP,
                             U, Ut, sJ,
                             UXiMinus, UXiPlus,
                             UEtaMinus, UEtaPlus,
                             #ifdef parabol
                               gradUx, gradUy,
                               gradUxXiMinus, gradUxXiPlus,
                               gradUxEtaMinus, gradUxEtaPlus,
                               gradUyXiMinus, gradUyXiPlus,
                               gradUyEtaMinus, gradUyEtaPlus,
                             #endif
                             dXdXi, dXdEta,
                             dYdXi, dYdEta,
                             XiBoundNormalX, XiBoundNormalY,
                             EtaBoundNormalX, EtaBoundNormalY,
                             XiBoundMetricsMinus, XiBoundMetricsPlus,
                             EtaBoundMetricsMinus, EtaBoundMetricsPlus,
                             XiBoundScalingMinus, XiBoundScalingPlus,
                             EtaBoundScalingMinus, EtaBoundScalingPlus,
                             innerFluxXi, innerFluxEta,
                             sideFluxXi, sideFluxEta,
                             LMinus, LPlus,
                             LHatMinus, LHatPlus,
                             DHat);

    // dU = dU * RKa - Ut * dt;
    for (int i = 0; i < nVars * nCellsX * nCellsY * (N + 1) * (N + 1); i++) {
      dU[i] = dU[i] * RKa[iStage] - Ut[i] * dt;
    }

    // U = U + dU * RKb;
    for (int i = 0; i < nVars * nCellsX * nCellsY * (N + 1) * (N + 1); i++) {
      U[i] += dU[i] * RKb[iStage];
    }

    // update tStage
//    tStage = t + RKc[iStage] * dt;

    #ifdef stepIteration
      MPI_Barrier(*cartComm);
      int myRank;
      MPI_Comm_rank(*cartComm, &myRank);
      if (myRank == 0) {
        printf("\n\n---- Finished RK-stage %d of 5! Press Enter to continue ----\n\n", iStage + 1);
        getchar();
      }
      MPI_Barrier(*cartComm);
    #endif
  }

}

void initLSERK(double* const restrict RKa, double* const restrict RKb, double* const restrict RKc) {

  RKa[0] = 0.0;
  RKa[1] = -567301805773.0 / 1357537059087.0;
  RKa[2] = -2404267990393.0 / 2016746695238.0;
  RKa[3] = -3550918686646.0 / 2091501179385.0;
  RKa[4] = -1275806237668.0 / 842570457699.0;

  RKb[0] = 1432997174477.0 / 9575080441755.0;
  RKb[1] = 5161836677717.0 / 13612068292358.0;
  RKb[2] = 1720146321549.0 / 2090206949498.0;
  RKb[3] = 3134564353537.0 / 4481467310338.0;
  RKb[4] = 2277821191437.0 / 14882151754819.0;

  RKc[0] = 0.0;
  RKc[1] = 1432997174477.0 / 9575080441755.0;
  RKc[2] = 2526269341429.0 / 6820363962896.0;
  RKc[3] = 2006345519317.0 / 3224310063776.0;
  RKc[4] = 2802321613138.0 / 2924317926251.0;
}
