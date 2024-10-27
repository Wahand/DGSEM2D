#include "DGSEM2D/mpiCommunication.h"

void sendMPIData(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                 const int N, const int nCellsX, const int nCellsY,
                 const int BCtop, const int BCbottom,
                 const int BCleft, const int BCright,
                 const double* const restrict dataXiMinus, const double* const restrict dataXiPlus,
                 const double* const restrict dataEtaMinus, const double* const restrict dataEtaPlus,
                 MPI_Request* const restrict reqHandles, const int tagOffset) {

  int Np1 = N + 1;

  // get information from communicator
  int myRank;
  int cartCommDims[2];
  int cartCommPeriodicity[2];
  int myCoords[2];
  MPI_Comm_rank(*cartComm, &myRank);
  MPI_Cart_get(*cartComm, 2, cartCommDims, cartCommPeriodicity, myCoords);

  // top
  // check if the process is connected to the physical boundary
  if (myCoords[0] == cartCommDims[0] - 1) {
    // send data only when the BC is periodic
    if (BCtop == 2) {
      MPI_Isend(&dataEtaPlus[nCellsY * nCellsX * (N + 1) * nVars], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[0], 0 + tagOffset, *cartComm, &reqHandles[0]);
    }
  } else {
    // if the process is an inner one, send data
    MPI_Isend(&dataEtaPlus[nCellsY * nCellsX * (N + 1) * nVars], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[0], 0 + tagOffset, *cartComm, &reqHandles[0]);
  }

  // bottom
  if (myCoords[0] == 0) {
    if (BCbottom == 2) {
      MPI_Isend(&dataEtaMinus[0], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[1], 1 + tagOffset, *cartComm, &reqHandles[1]);
    }
  } else {
    MPI_Isend(&dataEtaMinus[0], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[1], 1 + tagOffset, *cartComm, &reqHandles[1]);
  }

  // left
  if (myCoords[1] == 0) {
    if (BCleft == 2) {
      MPI_Isend(&dataXiMinus[0], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[2], 2 + tagOffset, *cartComm, &reqHandles[2]);
    }
  } else {
    MPI_Isend(&dataXiMinus[0], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[2], 2 + tagOffset, *cartComm, &reqHandles[2]);
  }

  // right
  if (myCoords[1] == cartCommDims[1] - 1) {
    if (BCright == 2) {
      MPI_Isend(&dataXiPlus[nCellsX * nCellsY * (N + 1) * nVars], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[3], 3 + tagOffset, *cartComm, &reqHandles[3]);
    }
  } else {
      MPI_Isend(&dataXiPlus[nCellsX * nCellsY * (N + 1) * nVars], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[3], 3 + tagOffset, *cartComm, &reqHandles[3]);
  }
}

void receiveMPIData(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                    const int N, const int nCellsX, const int nCellsY,
                    const int BCtop, const int BCbottom,
                    const int BCleft, const int BCright,
                    double* const restrict dataXiMinus, double* const restrict dataXiPlus,
                    double* const restrict dataEtaMinus, double* const restrict dataEtaPlus,
                    const int tagOffset) {

  int Np1 = N + 1;

  // get information from communicator
  int myRank;
  int cartCommDims[2];
  int cartCommPeriodicity[2];
  int myCoords[2];
  MPI_Comm_rank(*cartComm, &myRank);
  MPI_Cart_get(*cartComm, 2, cartCommDims, cartCommPeriodicity, myCoords);

  // top
  // Receive data if needed
  if (myCoords[0] == cartCommDims[0] - 1) {
    if (BCtop == 2) {
      // periodic BC
      // check if there is more than one processor in Y-direction
      if (cartCommDims[0] > 1) {
        MPI_Recv(&dataEtaMinus[nCellsY * nCellsX * (N + 1) * nVars], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[0], 1 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
      }
    }
  } else {
    // if the top side is an inner one, receive data
    if (cartCommDims[0] > 1) {
      MPI_Recv(&dataEtaMinus[nCellsY * nCellsX * (N + 1) * nVars], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[0], 1 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
    }
  }

  // bottom
  if (myCoords[0] == 0) {
    if (BCbottom == 2) {
      if (cartCommDims[0] > 1) {
        MPI_Recv(&dataEtaPlus[0], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[1], 0 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
      }
    }
  } else {
    if (cartCommDims[0] > 1) {
      MPI_Recv(&dataEtaPlus[0], nCellsX * Np1 * nVars, MPI_DOUBLE, neighbourRanks[1], 0 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
    }
  }

  // left
  if (myCoords[1] == 0) {
    if (BCleft == 2) {
      if (cartCommDims[1] > 1) {
        MPI_Recv(&dataXiPlus[0], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[2], 3 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
      }
    }
  } else {
    if (cartCommDims[0] > 1) {
      MPI_Recv(&dataXiPlus[0], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[2], 3 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
    }
  }

  // right
  if (myCoords[1] == cartCommDims[1] - 1) {
    if (BCright == 2) {
      if (cartCommDims[1] > 1) {
        MPI_Recv(&dataXiMinus[nCellsX * nCellsY * (N + 1) * nVars], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[3], 2 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
      }
    }
  } else {
    if (cartCommDims[1] > 1) {
      MPI_Recv(&dataXiMinus[nCellsX * nCellsY * (N + 1) * nVars], nCellsY * Np1 * nVars, MPI_DOUBLE, neighbourRanks[3], 2 + tagOffset, *cartComm, MPI_STATUS_IGNORE);
    }
  }
}
