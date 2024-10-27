#ifndef READINPUT_INCLUDED
#define READINPUT_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void readIniFile(const char* const restrict filename,
                 char* const restrict meshFile,
                 double* const restrict xMin, double* const restrict xMax,
                 double* const restrict yMin, double* const restrict yMax,
                 int* const restrict nCellsX, int* const nCellsY,
                 int* const restrict polynomeDegree, int* const restrict nodeType,
                 char* const restrict iniState,
                 int* const restrict BCtop, int* const restrict BCbottom,
                 int* const restrict BCleft, int* const restrict BCright,
                 double* const restrict CFL,
                 int* const restrict maxIter, double* const restrict endTime,
                 int* const restrict iOutput, double* const restrict tOutput,
                 char* const restrict outputPrefix);

void readMeshFile(const char* restrict filename,
                  int* const restrict nCellsX, int* const restrict nCellsY,
                  double** restrict nodesX, double** restrict nodesY);

void clearChar(char line[], char rem);

#endif
