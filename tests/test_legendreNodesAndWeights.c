#include <assert.h>
#include "DGSEM2D/basis.h"


int main() {

  double xGP[5];
  double wGP[5];

  legendreGaussNodesAndWeights(4, xGP, wGP);

  assert(xGP[2] == 0.0);

  return 0;
}
