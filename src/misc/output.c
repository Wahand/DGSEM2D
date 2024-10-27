#include "DGSEM2D/output.h"

void writeOutput(const MPI_Comm comm, const int myRank, const int numProcs, const char* const restrict filename, const int N, const int nCellsX, const int nCellsY, const double* const restrict data) {

  FILE* f;

  // open file
  if (myRank == 0) {
    f = fopen(filename, "w");
    if (f == NULL) {
      printf("Error opening file \"%s\"!\n", filename);
      return;
    }

    // write data of master process
    for (int iVar = 0; iVar < nVars; iVar++) {
      for (int iY = 0; iY < nCellsY * (N + 1); iY++) {
        for (int iX = 0; iX < nCellsX * (N + 1); iX++) {
          fprintf(f, "%f, ", data[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iY * nCellsX * (N + 1) + iX]);
        }
        fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }
  }

  if (myRank == 0) {
    for (int iProc = 1; iProc < numProcs; iProc++) {
      // get data size
      int bufSize;
      MPI_Recv(&bufSize, 1, MPI_INT, iProc, 0, comm, MPI_STATUS_IGNORE);

      // receive  data
      double otherData[nVars * bufSize];
      MPI_Recv(otherData, nVars*bufSize, MPI_DOUBLE, iProc, 3, comm, MPI_STATUS_IGNORE);
      
      // print data
      for (int i = 0; i < nVars*bufSize; i++) {
        fprintf(f, "%f, ", otherData[i]);
      }
      fprintf(f, "\n");
    }
  } else {
    // send data size
    int bufSize = nCellsX * (N + 1) * nCellsY * (N + 1);
    MPI_Send(&bufSize, 1, MPI_INT, 0, 0, comm);
    // send data
    MPI_Send(data, bufSize, MPI_DOUBLE, 0, 3, comm);
  }

  // close file
  if (myRank == 0) {
    fclose(f);
  }
}

void printOutput(const int N, const int nCellsX, const int nCellsY, const double* const data) {

  for (int iVar = 0; iVar < nVars; iVar++) {
    for (int iY = 0; iY < nCellsY * (N + 1); iY++) {
      for (int iX = 0; iX < nCellsX * (N + 1); iX++) {
        printf("%f, ", data[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iY * nCellsX * (N + 1) + iX]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n\n");
}

void writeOutputWithCoords(const MPI_Comm comm, const int myRank, const int numProcs, const char* const restrict filename, const int N, const int nCellsX, const int nCellsY, const double* const restrict xGP, const double* const restrict yGP, const double* const restrict data) {

  FILE* f;

  if (myRank == 0) {
    f = fopen(filename, "w");
    if (f == NULL) {
      printf("Error opening file \"%s\"!\n", filename);
      return;
    }

    // print header
    fprintf(f, "X, Y");
    for (int iVar = 0; iVar < nVars; iVar++) {
      fprintf(f, ", U%d", iVar);
    }
    fprintf(f, "\n");

    // write data of master process
    for (int i = 0; i < nCellsX * (N + 1) * nCellsY * (N + 1); i++) {
      // print coordinates
      fprintf(f, "%f, %f", xGP[i], yGP[i]);
      // print data
      for (int iVar = 0; iVar < nVars; iVar++) {
        fprintf(f, ", %f", data[i * nVars + iVar]);
/*
        fprintf(f, ", %f", data[iVar * nCellsX * (N + 1) * nCellsY * (N + 1) + i]);
*/
      }
      fprintf(f, "\n");
    }
  }

  if (myRank == 0) {
    for (int iProc = 1; iProc < numProcs; iProc++) {
      // receive data size
      int bufSize;
      MPI_Recv(&bufSize, 1, MPI_INT, iProc, 0, comm, MPI_STATUS_IGNORE);
      // receive coordinates
      double otherxGP[bufSize];
      double otheryGP[bufSize];
      MPI_Recv(otherxGP, bufSize, MPI_DOUBLE, iProc, 1, comm, MPI_STATUS_IGNORE);
      MPI_Recv(otheryGP, bufSize, MPI_DOUBLE, iProc, 2, comm, MPI_STATUS_IGNORE);
      // receive data
      double otherData[nVars*bufSize];
      MPI_Recv(otherData, nVars*bufSize, MPI_DOUBLE, iProc, 3, comm, MPI_STATUS_IGNORE);

      for (int i = 0; i < bufSize; i++) {
        // print coordinates
        fprintf(f, "%f, %f", otherxGP[i], otheryGP[i]);
        // print data
        for (int iVar = 0; iVar < nVars; iVar++) {
          fprintf(f, ", %f", otherData[i * nVars + iVar]);
/*
          fprintf(f, ", %f", otherData[iVar * bufSize + i]);
*/
        }
        fprintf(f, "\n");
      }
    }
      
  } else {
    // send data size
    int bufSize = nCellsX * (N + 1) * nCellsY * (N + 1);
    MPI_Send(&bufSize, 1, MPI_INT, 0, 0, comm);
    // send coordinates
    MPI_Send(xGP, bufSize, MPI_DOUBLE, 0, 1, comm);
    MPI_Send(yGP, bufSize, MPI_DOUBLE, 0, 2, comm);
    // send data
    MPI_Send(data, nVars*bufSize, MPI_DOUBLE, 0, 3, comm);
  }

  if (myRank == 0) {
    fclose(f);
  }

/*
  FILE* f = fopen(filename, "w");
  if (f == NULL) {
    printf("Error opening file \"%s\"!\n", filename);
    return;
  }

  // print header
  fprintf(f, "X, Y");
  for (int iVar = 0; iVar < nVars; iVar++) {
    fprintf(f, ", U%d", iVar);
  }
  fprintf(f, "\n");

  for (int iCellY = 0; iCellY < nCellsY; iCellY++) {
    for (int iGPY = 0; iGPY < N + 1; iGPY++) {
      for (int iCellX = 0; iCellX < nCellsX; iCellX++) {
        for (int iGPX = 0; iGPX < N + 1; iGPX++) {
          // print coordinates
          fprintf(f, "%f, %f", xGP[iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX],
                               yGP[iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX]);

          // print variables
          for (int iVar = 0; iVar < nVars; iVar++) {
            fprintf(f, ", %f", data[iVar * nCellsX * nCellsY * (N + 1) * (N + 1) + iCellY * nCellsX * (N + 1) * (N + 1) + iGPY * nCellsX * (N + 1) + iCellX * (N + 1) + iGPX]);
          }
          fprintf(f, "\n");
        }
      }
    }
  }

  fclose(f);
*/
}

void writeMatrix(const char* const filename, const int sizeX, const int sizeY, const double* const data) {

  FILE* f = fopen(filename, "w");

  for (int j = 0; j < sizeY; j++) {
    for (int i = 0; i < sizeX; i++) {
      fprintf(f, "%f  ", data[j * sizeX + i]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}

void printMatrix(const int sizeX, const int sizeY, const double* const data) {
  for (int j = 0; j < sizeY; j++) {
    for (int i = 0; i < sizeX; i++) {
      printf("%f  ", data[j * sizeX + i]);
    }
    printf("\n");
  }
  printf("\n\n");
}
