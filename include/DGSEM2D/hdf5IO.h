#ifndef HDF5IO_INCLUDED
#define HDF5IO_INCLUDED

#include <hdf5.h>
#include <mpi.h>

void initOutputFile(const char* const filename, const MPI_Comm comm);

void addBlock(const char* const filename, const MPI_Comm comm, const int IDBlock, const int N, const int nCellsXGlobal, const int nCellsYGlobal, const int nCellsX, const int nCellsY, const double* const restrict nodesX, const double* const restrict nodesY, const double* const restrict xGP, const double* const restrict yGP);

void writeSolution(const char* const filename, const MPI_Comm comm, const int IDBlock, const int N, const int nCellsXGlobal, const int nCellsYGlobal, const int nCellsX, const int nCellsY, const int iOutput, const double time, const double* const restrict U);

#endif
