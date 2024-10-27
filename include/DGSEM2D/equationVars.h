#ifndef EQUATIONVARS_INCLUDED
#define EQUATIONVARS_INCLUDED

#include <stdlib.h>
#include <stdio.h>

// linear advection
#if (eqSys == 1)
  extern const double eqVars_ax;
  extern const double eqVars_ay;
#endif

// Euler
#if (eqSys == 2)
  extern const double eqVars_gamma;
#endif

// advection-diffusion
#if (eqSys == 3)
  extern const double eqVars_ax;
  extern const double eqVars_ay;
  extern const double eqVars_b;
#endif

#endif
