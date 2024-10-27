#ifndef MPIH_INCLUDED
#define MPIH_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void sendMPIData(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                 const int N, const int nCellsX, const int nCellsY,
                 const int BCtop, const int BCbottom,
                 const int BCleft, const int BCright,
                 const double* const restrict dataXiMinus, const double* const restrict dataXiPlus,
                 const double* const restrict dataEtaMinus, const double* const restrict dataEtaPlus,
                 MPI_Request* const restrict reqHandles, const int tagOffset);

void receiveMPIData(const MPI_Comm* const restrict cartComm, const int* const restrict neighbourRanks,
                    const int N, const int nCellsX, const int nCellsY,
                    const int BCtop, const int BCbottom,
                    const int BCleft, const int BCright,
                    double* const restrict dataXiMinus, double* const restrict dataXiPlus,
                    double* const restrict dataEtaMinus, double* const restrict dataEtaPlus,
                    const int tagOffset);

#endif
