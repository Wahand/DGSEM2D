#ifndef OUTPUT_INCLUDED
#define OUTPUT_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

void writeOutput(const MPI_Comm comm, const int myRank, const int numProcs, const char* const restrict filename, const int N, const int nCellsX, const int nCellsY, const double* const restrict data);

void printOutput(const int N, const int nCellsX, const int nCellsY, const double* const restrict data);

void writeOutputWithCoords(const MPI_Comm comm, const int myRank, const int numProcs, const char* const restrict filename, const int N, const int nCellsX, const int nCellsY, const double* const restrict xGP, const double* const restrict yGP, const double* const restrict data);

void writeMatrix(const char* const restrict filename, const int sizeX, const int sizeY, const double* const restrict data);

void printMatrix(const int sizeX, const int sizeY, const double* const restrict data);

#endif
