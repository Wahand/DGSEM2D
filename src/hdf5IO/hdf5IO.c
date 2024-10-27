#include "hdf5IO.h"

void initOutputFile(const char* const filename, const MPI_Comm comm) {

  // Open file
  hid_t IDPList = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(IDPList, comm, MPI_INFO_NULL);
  hid_t IDFile = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, IDPList);
  H5Pclose(IDPList);

  // Add attributes
  hid_t IDDataspaceIntAttribute = H5Screate(H5S_SCALAR);
  hid_t IDAttributeNBlocks = H5Acreate2(IDFile, "nBlocks", H5T_NATIVE_INT, IDDataspaceIntAttribute, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(IDDataspaceIntAttribute);
  int nBlocks = 0;
  H5Awrite(IDAttributeNBlocks, H5T_NATIVE_INT, &nBlocks);
  H5Aclose(IDAttributeNBlocks);

  hsize_t dims[1] = {0};
  hsize_t maxdims[1] = {H5S_UNLIMITED};
  hsize_t chunkdims[1] = {1};
  hid_t IDDataspaceArray = H5Screate_simple(1, dims, maxdims);
  hid_t IDchunkProperty = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(IDchunkProperty, 1, chunkdims);
  hid_t IDDatasetSolutionTimes = H5Dcreate2(IDFile, "solutionTimes", H5T_NATIVE_DOUBLE, IDDataspaceArray, H5P_DEFAULT, IDchunkProperty, H5P_DEFAULT);
  H5Sclose(IDDataspaceArray);
  H5Pclose(IDchunkProperty);
  H5Dclose(IDDatasetSolutionTimes);
  H5Fclose(IDFile);

}

void addBlock(const char* const filename, const MPI_Comm comm, const int IDBlock, const int N, const int nCellsXGlobal, const int nCellsYGlobal, const int nCellsX, const int nCellsY, const double* const restrict nodesX, const double* const restrict nodesY, const double* const restrict xGP, const double* const restrict yGP) {

  char tmpString[255];

  // Open file
  hid_t IDPListCreate = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(IDPListCreate, comm, MPI_INFO_NULL);
  hid_t IDFile = H5Fopen(filename, H5F_ACC_RDWR, IDPListCreate);
  H5Pclose(IDPListCreate);

  // Get data from communicator
  int myRank;
  int commDims[2];
  int commPeriodicity[2];
  int commCoords[2];
  MPI_Comm_rank(comm, &myRank);
  MPI_Cart_get(comm, 2, commDims, commPeriodicity, commCoords);

  // Get offset of local nodes in global nodes
  int nodestartX;
  int nodestartY;
  // nodestartX
  if (nCellsXGlobal % commDims[1] == 0) {
    nodestartX = nCellsX * commCoords[1];
  } else {
    if (commCoords[1] < nCellsXGlobal % commDims[1]) {
      nodestartX = commCoords[1] * (nCellsXGlobal / commDims[1] + 1);
    } else {
      nodestartX = (nCellsXGlobal % commDims[1]) * (nCellsXGlobal / commDims[1] + 1) + (commCoords[1] - (nCellsXGlobal % commDims[1])) * (nCellsXGlobal / commDims[1]);
    }
  }
  // nodestartY
  if (nCellsYGlobal % commDims[0] == 0) {
    nodestartY = nCellsY * commCoords[0];
  } else {
    if (commCoords[0] < nCellsYGlobal % commDims[0]) {
      nodestartY = commCoords[0] * (nCellsYGlobal / commDims[0] + 1);
    } else {
      nodestartY = (nCellsYGlobal % commDims[0]) * (nCellsYGlobal / commDims[0] + 1) + (commCoords[0] - (nCellsYGlobal % commDims[0])) * (nCellsYGlobal / commDims[0]);
    }
  }

  // Add block to nBlocks
  int nBlocks = -1;
  hid_t IDAttributeNBlocks = H5Aopen(IDFile, "nBlocks", H5P_DEFAULT);
  H5Aread(IDAttributeNBlocks, H5T_NATIVE_INT, &nBlocks);
  if (IDBlock >= nBlocks) {
    nBlocks++;
    H5Awrite(IDAttributeNBlocks, H5T_NATIVE_INT, &nBlocks);
  }
  H5Aclose(IDAttributeNBlocks);

  // Create block group
  snprintf(tmpString, 255, "/block%d", IDBlock);
  hid_t IDBlockGroup = H5Gcreate2(IDFile, tmpString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create solution group
  snprintf(tmpString, 255, "/block%d/solution", IDBlock);
  hid_t IDSolutionGroup = H5Gcreate2(IDBlockGroup, tmpString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(IDSolutionGroup);

  // Create attributes
  hid_t IDDataspaceIntAttribute = H5Screate(H5S_SCALAR);
  hid_t IDAttributeN = H5Acreate2(IDBlockGroup, "N", H5T_NATIVE_INT, IDDataspaceIntAttribute, H5P_DEFAULT, H5P_DEFAULT);
  hid_t IDAttributeNCellsX = H5Acreate2(IDBlockGroup, "nCellsX", H5T_NATIVE_INT, IDDataspaceIntAttribute, H5P_DEFAULT, H5P_DEFAULT);
  hid_t IDAttributeNCellsY = H5Acreate2(IDBlockGroup, "nCellsY", H5T_NATIVE_INT, IDDataspaceIntAttribute, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(IDAttributeN, H5T_NATIVE_INT, &N);
  H5Awrite(IDAttributeNCellsX, H5T_NATIVE_INT, &nCellsXGlobal);
  H5Awrite(IDAttributeNCellsY, H5T_NATIVE_INT, &nCellsYGlobal);
  H5Aclose(IDAttributeN);
  H5Aclose(IDAttributeNCellsX);
  H5Aclose(IDAttributeNCellsY);
  H5Sclose(IDDataspaceIntAttribute);

  // Create mesh group
  snprintf(tmpString, 255, "/block%d/mesh", IDBlock);
  hid_t IDMeshGroup = H5Gcreate2(IDBlockGroup, tmpString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write nodes
  // (X- and Y-coordinates are "swapped", so that the X-coordinate is the fastest index)
  snprintf(tmpString, 255, "/block%d/mesh/nodes", IDBlock);
  hsize_t dimsNodes[3];
  dimsNodes[0] = nCellsYGlobal + 1;
  dimsNodes[1] = nCellsXGlobal + 1;
  dimsNodes[2] = 2;
  hid_t IDDataspaceNodes = H5Screate_simple(3, dimsNodes, NULL);
  hid_t IDDatasetNodes = H5Dcreate2(IDMeshGroup, tmpString, H5T_NATIVE_DOUBLE, IDDataspaceNodes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t lengthNodes[2] = {nCellsY + 1, nCellsX + 1};
  hid_t IDMemspace = H5Screate_simple(2, lengthNodes, NULL);
  hid_t IDHyperslab = H5Dget_space(IDDatasetNodes);

  hsize_t offset[3] = {nodestartY, nodestartX, 0};
  hsize_t count[3] = {1, 1, 1};
  hsize_t block[3] = {nCellsY + 1, nCellsX + 1, 1};

  H5Sselect_hyperslab(IDHyperslab, H5S_SELECT_SET, offset, NULL, count, block);
  H5Dwrite(IDDatasetNodes, H5T_NATIVE_DOUBLE, IDMemspace, IDHyperslab, H5P_DEFAULT, nodesX);
  offset[2] = 1;
  H5Sselect_hyperslab(IDHyperslab, H5S_SELECT_SET, offset, NULL, count, block);
  H5Dwrite(IDDatasetNodes, H5T_NATIVE_DOUBLE, IDMemspace, IDHyperslab, H5P_DEFAULT, nodesY);

  // Close nodes
  H5Dclose(IDDatasetNodes);
  H5Sclose(IDDataspaceNodes);

  // Get starting indices of mesh-matrix
  int meshStartX = nCellsXGlobal / commDims[1] * (N + 1) * commCoords[1];
  int meshStartY = nCellsYGlobal / commDims[0] * (N + 1) * commCoords[0];
  if (commCoords[1] >= nCellsXGlobal % commDims[1]) {
    meshStartX += (nCellsXGlobal % commDims[1]) * (N + 1);
  } else {
    meshStartX += commCoords[1] * (N + 1);
  }
  if (commCoords[0] >= nCellsYGlobal % commDims[0]) {
    meshStartY += (nCellsYGlobal % commDims[0]) * (N + 1);
  } else {
    meshStartY += commCoords[0] * (N + 1);
  }

  // Write mesh
  snprintf(tmpString, 255, "/block%d/mesh/meshCoordinates", IDBlock);
  hsize_t dimsCoords[3];
  dimsCoords[0] = nCellsYGlobal * (N + 1);
  dimsCoords[1] = nCellsXGlobal * (N + 1);
  dimsCoords[2] = 2;
  hid_t IDDataspaceCoords = H5Screate_simple(3, dimsCoords, NULL);
  hid_t IDDatasetCoords = H5Dcreate2(IDMeshGroup, tmpString, H5T_NATIVE_DOUBLE, IDDataspaceCoords, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(IDDataspaceCoords);

  hsize_t lengthCoords[2] = {nCellsY * (N + 1), nCellsX * (N + 1)};
  IDMemspace = H5Screate_simple(2, lengthCoords, NULL);
  IDHyperslab = H5Dget_space(IDDatasetCoords);

  offset[0] = meshStartY;
  offset[1] = meshStartX;
  offset[2] = 0;
  block[0] = nCellsY * (N + 1);
  block[1] = nCellsX * (N + 1);

  H5Sselect_hyperslab(IDHyperslab, H5S_SELECT_SET, offset, NULL, count, block);
  H5Dwrite(IDDatasetCoords, H5T_NATIVE_DOUBLE, IDMemspace, IDHyperslab, H5P_DEFAULT, xGP);
  offset[2] = 1;
  H5Sselect_hyperslab(IDHyperslab, H5S_SELECT_SET, offset, NULL, count, block);
  H5Dwrite(IDDatasetCoords, H5T_NATIVE_DOUBLE, IDMemspace, IDHyperslab, H5P_DEFAULT, yGP);

  // Close coords
  H5Dclose(IDDatasetCoords);

  // Close everything
  H5Sclose(IDHyperslab);
  H5Sclose(IDMemspace);
  H5Gclose(IDMeshGroup);
  H5Gclose(IDBlockGroup);
  H5Fclose(IDFile);

}

void writeSolution(const char* const filename, const MPI_Comm comm, const int IDBlock, const int N, const int nCellsXGlobal, const int nCellsYGlobal, const int nCellsX, const int nCellsY, const int iOutput, const double time, const double* const restrict U) {

  char tmpString[255];

  // Open file
  hid_t IDPListCreate = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(IDPListCreate, comm, MPI_INFO_NULL);
  hid_t IDFile = H5Fopen(filename, H5F_ACC_RDWR, IDPListCreate);
  H5Pclose(IDPListCreate);

  // Get data from communicator
  int myRank;
  int commDims[2];
  int commPeriodicity[2];
  int commCoords[2];
  MPI_Comm_rank(comm, &myRank);
  MPI_Cart_get(comm, 2, commDims, commPeriodicity, commCoords);

  // Add solution time if needed
  hid_t IDDatasetSolutionTimes = H5Dopen2(IDFile, "/solutionTimes", H5P_DEFAULT);
  hid_t IDDataspaceSolutionTimes = H5Dget_space(IDDatasetSolutionTimes);
  hsize_t solutionTimesLength;
  H5Sget_simple_extent_dims(IDDataspaceSolutionTimes, &solutionTimesLength, NULL);
  if ((unsigned) iOutput >= solutionTimesLength) {
    // Read solution times
    double solutionTimes[solutionTimesLength + 1];
    H5Dread(IDDatasetSolutionTimes, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, solutionTimes);
    solutionTimes[solutionTimesLength] = time;
    solutionTimesLength++;
    H5Dset_extent(IDDatasetSolutionTimes, &solutionTimesLength);
    H5Dwrite(IDDatasetSolutionTimes, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, solutionTimes);
  }
  H5Sclose(IDDataspaceSolutionTimes);
  H5Dclose(IDDatasetSolutionTimes);

  // Get solution group
  snprintf(tmpString, 255, "/block%d/solution", IDBlock);
  hid_t IDSolutionGroup = H5Gopen(IDFile, tmpString, H5P_DEFAULT);

  // Create solution dataset
  snprintf(tmpString, 255, "/block%d/solution/solution%05d", IDBlock, iOutput);
  hsize_t dimsSolution[3];
  dimsSolution[0] = nCellsXGlobal * (N + 1);
  dimsSolution[1] = nCellsYGlobal * (N + 1);
  dimsSolution[2] = nVars;
  hid_t IDDataspaceSolution = H5Screate_simple(3, dimsSolution, NULL);
  hid_t IDDatasetSolution = H5Dcreate2(IDSolutionGroup, tmpString, H5T_NATIVE_DOUBLE, IDDataspaceSolution, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(IDDataspaceSolution);

  // Get starting indices of solution-matrix
  int solutionStartX = nCellsXGlobal / commDims[1] * (N + 1) * commCoords[1];
  int solutionStartY = nCellsYGlobal / commDims[0] * (N + 1) * commCoords[0];
  if (commCoords[1] >= nCellsXGlobal % commDims[1]) {
    solutionStartX += (nCellsXGlobal % commDims[1]) * (N + 1);
  } else {
    solutionStartX += commCoords[1] * (N + 1);
  }
  if (commCoords[0] >= nCellsYGlobal % commDims[0]) {
    solutionStartY += (nCellsYGlobal % commDims[0]) * (N + 1);
  } else {
    solutionStartY += commCoords[0] * (N + 1);
  }

  hsize_t lengthSolution[3] = {nCellsY * (N + 1), nCellsX * (N + 1), nVars};
  hid_t IDMemspace = H5Screate_simple(3, lengthSolution, NULL);
  hid_t IDHyperslab = H5Dget_space(IDDatasetSolution);

  hsize_t offset[3] = {solutionStartY, solutionStartX, 0};
  hsize_t count[3] = {1, 1, 1};
  hsize_t block[3] = {nCellsY * (N + 1), nCellsX * (N + 1), nVars};

  H5Sselect_hyperslab(IDHyperslab, H5S_SELECT_SET, offset, NULL, count, block);
  H5Dwrite(IDDatasetSolution, H5T_NATIVE_DOUBLE, IDMemspace, IDHyperslab, H5P_DEFAULT, U);

  // Close everything
  H5Dclose(IDDatasetSolution);
  H5Sclose(IDHyperslab);
  H5Sclose(IDMemspace);
  H5Gclose(IDSolutionGroup);
  H5Fclose(IDFile);

}
