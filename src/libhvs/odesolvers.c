#include "odesolvers.h"
#include <math.h>

FLOAT_TYPE rk4_solve(void (*func)(void,FLOAT_TYPE), FLOAT_TYPE tn, FLOAT_TYPE timestep) {
	void k1,k2,k3,k4;
	
}