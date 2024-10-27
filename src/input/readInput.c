#include "DGSEM2D/readInput.h"

void readIniFile(const char* const restrict filename,
                 char* const restrict meshFile,
                 double* const restrict xMin, double* const restrict xMax,
                 double* const restrict yMin, double* const restrict yMax,
                 int* const restrict nCellsX, int* const restrict nCellsY,
                 int* const restrict polynomeDegree, int* const restrict nodeType,
                 char* const restrict iniState,
                 int* const restrict BCtop, int* const restrict BCbottom,
                 int* const restrict BCleft, int* const restrict BCright,
                 double* const restrict CFL,
                 int* const restrict maxIter, double* const restrict endTime,
                 int* const restrict iOutput, double* const restrict tOutput,
                 char* const restrict outputPrefix) {

  #if (outputLevel >= 2)
    printf("Start reading input from file \"%s\"...\n", filename);
  #endif

  FILE* f = fopen(filename, "r");

  char delim = '=';
  char* delimPos = NULL;

  char line[255] = "";
  char paramName[255] = "";
  char paramVal[255] = "";

  while (fgets(line, sizeof(line), f) != NULL) {
    if (strlen(line) != 1 && line[0] != '#') {
      delimPos = strchr(line, delim);
      *delimPos = '\0';

      strcpy(paramName, line);
      strcpy(paramVal, delimPos + 1);

      clearChar(paramName, ' ');
      clearChar(paramName, '\n');
      clearChar(paramName, '\t');
      clearChar(paramVal, ' ');
      clearChar(paramVal, '\n');
      clearChar(paramVal, '\t');

      if (strcmp(paramName, "meshfile") == 0) {
        memcpy(meshFile, paramVal, strlen(paramVal));
        #if (outputLevel >= 2)
          printf("Set parameter meshFile to \"%s\"\n", meshFile);
        #endif
      }
      else if (strcmp(paramName, "xMin") == 0) {
        *xMin = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter xMin to %4.2f\n", *xMin);
        #endif
      }
      else if (strcmp(paramName, "xMax") == 0) {
        *xMax = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter xMax to %4.2f\n", *xMax);
        #endif
      }
      else if (strcmp(paramName, "yMin") == 0) {
        *yMin = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter yMin to %4.2f\n", *yMin);
        #endif
      }
      else if (strcmp(paramName, "yMax") == 0) {
        *yMax = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter yMax to %4.2f\n", *yMax);
        #endif
      }
      else if (strcmp(paramName, "nCellsX") == 0) {
        *nCellsX = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter nCellsX to %d\n", *nCellsX);
        #endif
      }
      else if (strcmp(paramName, "nCellsY") == 0) {
        *nCellsY = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter nCellsY to %d\n", *nCellsY);
        #endif
      }
      else if (strcmp(paramName, "polynomeDegree") == 0) {
        *polynomeDegree = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter polynomeDegree to %d\n", *polynomeDegree);
        #endif
      }
      else if (strcmp(paramName, "nodeType") == 0) {
        *nodeType = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter nodeType to %d\n", *nodeType);
        #endif
      }
      else if (strcmp(paramName, "iniState") == 0) {
        memcpy(iniState, paramVal, strlen(paramVal));
        #if (outputLevel >= 2)
          printf("Set parameter iniState to \"%s\"\n", *iniState);
        #endif
      }
      else if (strcmp(paramName, "BCtop") == 0) {
        *BCtop = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter BCtop to %d\n", *BCtop);
        #endif
      }
      else if (strcmp(paramName, "BCbottom") == 0) {
        *BCbottom = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter BCbottom to %d\n", *BCbottom);
        #endif
      }
      else if (strcmp(paramName, "BCleft") == 0) {
        *BCleft = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter BCleft to %d\n", *BCleft);
        #endif
      }
      else if (strcmp(paramName, "BCright") == 0) {
        *BCright = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter BCright to %d\n", *BCright);
        #endif
      }
      else if (strcmp(paramName, "CFL") == 0) {
        *CFL = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter CFL to %4.2f\n", *CFL);
        #endif
      }
      else if (strcmp(paramName, "maxIter") == 0) {
        *maxIter = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter maxIter to %d\n", *maxIter);
        #endif
      }
      else if (strcmp(paramName, "endTime") == 0) {
        *endTime = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter endTime to %4.2f\n", *endTime);
        #endif
      }
      else if (strcmp(paramName, "iOutput") == 0) {
        *iOutput = atoi(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter iOutput to %d\n", *iOutput);
        #endif
      }
      else if (strcmp(paramName, "tOutput") == 0) {
        *tOutput = atof(paramVal);
        #if (outputLevel >= 2)
          printf("Set parameter tOutput to %4.2f\n", *tOutput);
        #endif
      }
      else if (strcmp(paramName, "outputPrefix") == 0) {
        memcpy(outputPrefix, paramVal, strlen(paramVal));
        #if (outputLevel >= 2)
          printf("Set parameter outputPrefix to \"%s\"\n", outputPrefix);
        #endif
      }
      else {
        #if (outputLevel >= 1)
          printf("Unknown parameter-pair \"%s\", \"%s\"!\n", paramName, paramVal);
        #endif
      }
    }

    memset(line, '\0', sizeof(line));
  }

  fclose(f);

  // check parameters
  //TODO check if needed input is given
  if (*BCtop == 3 && *BCbottom != 3) {
    printf("\n  +----\n  | Warning: Top BC is periodic, but bottom is not!\n  +----\n\n");
  }
  if (*BCbottom == 3 && *BCtop != 3) {
    printf("\n  +----\n  | Warning: Bottom BC is periodic, but top is not!\n  +----\n\n");
  }
  if (*BCleft == 3 && *BCright != 3) {
    printf("\n  +----\n  | Warning: Left BC is periodic, but right is not!\n  +----\n\n");
  }
  if (*BCright == 3 && *BCleft != 3) {
    printf("\n  +----\n  |Warning: Right BC is periodic, but left is not!\n  +----\n\n");
  }

  #if (outputLevel >= 2)
    printf("Done reading input!\n");
  #endif

}

void readMeshFile(const char* const restrict filename,
                  int* const restrict nCellsX, int* const restrict nCellsY,
                  double** restrict nodesX, double** restrict nodesY) {

  #if (outputLevel >= 2)
    printf("Start reading meshfile \"%s\"...\n", filename);
  #endif

  FILE* f = fopen(filename, "r");

  char delim = ',';
  char* delimPos = NULL;
  char line[255] = "";
  char x[255] = "";
  char y[255] = "";
  int i = 0;

  // read number of nodes in each direction (and save to nCells)
  if (fgets(line, sizeof(line), f) != NULL) {
    delimPos = strchr(line, delim);
    *delimPos = '\0';

    strcpy(x, line);
    strcpy(y, delimPos + 1);

    *nCellsX = atoi(x) - 1;
    #if (outputLevel >= 2)
      printf("Set parameter nCellsX to %d\n", *nCellsX);
    #endif

    *nCellsY = atoi(y) - 1;
    #if (outputLevel >= 2)
      printf("Set parameter nCellsY to %d\n", *nCellsY);
    #endif
  }

  memset(line, '\0', sizeof(line));

  // allocate memory for nodes
  *nodesX = malloc((*nCellsX + 1) * (*nCellsY + 1) * sizeof(double));
  *nodesY = malloc((*nCellsX + 1) * (*nCellsY + 1) * sizeof(double));

  // read node-coordinates
  while (fgets(line, sizeof(line), f) != NULL) {
    if (strlen(line) != 1) {
      delimPos = strchr(line, delim);
      *delimPos = '\0';

      strcpy(x, line);
      strcpy(y, delimPos + 1);

      (*nodesX)[i] = atof(x);
      (*nodesY)[i] = atof(y);
      i++;
    }

    memset(line, '\0', sizeof(line));
  }

  fclose(f);

  #if (outputLevel >= 2)
    printf("Done reading mesh-file!\n");
  #endif

}

void clearChar(char str[], char rem) {

  char buf[strlen(str)];
  memset(buf, '\0', sizeof(buf));
  int i = 0;

  for (unsigned long int j = 0; j < strlen(str); j++) {
    if (str[j] != rem) {
      buf[i++] = str[j];
    }
  }

  memcpy(str, buf, sizeof(buf));
}
