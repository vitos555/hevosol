#include "hermiteutil.h"
#include <math.h>

#define PI_INV 1/M_PI

#if NMOMENTS > 2
#error "Current implementation doesn't support moments of order higher then 2."
#endif

FLOAT_TYPE he_00(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return PI_INV/lambda_sq*exp(-(x1*x1+x2*x2)/lambda_sq);
}

FLOAT_TYPE he_10(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return he_00(x1,x2,lambda_sq)*(-2*x1)/lambda_sq;
}

FLOAT_TYPE he_01(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return he_00(x1,x2,lambda_sq)*(-2*x2)/lambda_sq;
}

FLOAT_TYPE he_20(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return he_00(x1,x2,lambda_sq)*(-2*lambda_sq+4*x1*x1)/lambda_sq/lambda_sq;
}

FLOAT_TYPE he_02(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return he_00(x1,x2,lambda_sq)*(-2*lambda_sq+4*x2*x2)/lambda_sq/lambda_sq;
}

FLOAT_TYPE he_11(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return he_00(x1,x2,lambda_sq)*(4*x1*x2)/lambda_sq/lambda_sq;
}

FLOAT_TYPE he(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq, unsigned short k1, unsigned short k2) {
	if ( (k1==0) && (k2==0) )
		return he_00(x1,x2,lambda_sq);
	else if ( (k1==1) && (k2==0) )
		return he_10(x1,x2,lambda_sq);
	else if ( (k1==0) && (k2==1) )
		return he_01(x1,x2,lambda_sq);
	else if ( (k1==2) && (k2==0) )
		return he_20(x1,x2,lambda_sq);
	else if ( (k1==1) && (k2==1) )
		return he_11(x1,x2,lambda_sq);
	else if ( (k1==0) && (k2==2) )
		return he_02(x1,x2,lambda_sq);
	else
		return 0.0;
}

