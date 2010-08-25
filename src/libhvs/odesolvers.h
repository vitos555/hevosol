#define ODESOLVERS_H 1

#ifndef FLOAT_TYPE
#include "gentypes.h"
#endif

FLOAT_TYPE rk4_solve(void (*func)(void,FLOAT_TYPE), FLOAT_TYPE timestep);