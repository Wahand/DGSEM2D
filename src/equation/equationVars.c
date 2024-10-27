#include "DGSEM2D/equationVars.h"

// linear advection
#if (eqSys == 1)
  const double eqVars_ax = 1.0;
  const double eqVars_ay = 0.2;
#endif

// Euler
#if (eqSys == 2)
  const double eqVars_gamma = 1.4;
#endif

// advection-diffusion
#if (eqSys == 3)
  const double eqVars_ax = 1.0;
  const double eqVars_ay = 0.5;
  const double eqVars_b = 0.02;
#endif
