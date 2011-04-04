#include "hermiteutil.h"
#include "factorialutil.h"
#include <math.h>

#if HVS_DEBUG
#include<stdio.h>
#include "errorutil.h"
#endif

#define PI_INV 1/M_PI

#if NMOMENTS > 4
#error "Current implementation doesn't support moments of order higher then 4."
#endif

#if NMOMENTS<=2
#include "he05.c"
#elif NMOMENTS==3
#include "he07.c"
#elif NMOMENTS==4
#include "he10.c"
#endif

FLOAT_TYPE Power(FLOAT_TYPE x, int power) {
	switch(power) {
		case 0:
			return 1.0;
		case 1:
			return x;
		case 2:
			return x*x;
		case 3:
			return x*x*x;
		case 4:
			return x*x*x*x;
		default:
			return M_POW(x,power);
	}
}

FLOAT_TYPE h1(UINT alpha1, UINT alpha2, FLOAT_TYPE lambda_sq) {
	if ((alpha1%2==0) && (alpha2%2==1)) {
		int alpha_temp = (alpha1+alpha2-1)/2;
		int negone = (alpha_temp%2==0?1:-1);
		return
			(-0.5)*PI_INV/Power(2*lambda_sq,alpha_temp+1)*negone*
			factorial(alpha1)/factorial(alpha_temp+1)*
			factorial(alpha2)*binomial(alpha_temp,(int)alpha1/2);
	} else {
#if HVS_DEBUG
		hvsdie("Error. Function h1: Parameters out of range\n");
#endif
		return 0.0;
	}
}

FLOAT_TYPE h2(UINT alpha1, UINT alpha2,FLOAT_TYPE lambda_sq) {
	if ((alpha1%2==1) && (alpha2%2==0)) {
		int alpha_temp = (alpha1+alpha2-1)/2;
		int negone = (alpha_temp%2==0?1:-1);
		return
			(0.5)*PI_INV/Power(2*lambda_sq,alpha_temp+1)*negone*
			factorial(alpha1)/factorial(alpha_temp+1)*
			factorial(alpha2)*binomial(alpha_temp,(UINT)(alpha1-1)/2);
	} else {
#if HVS_DEBUG
		hvsdie("Error. Function h2: Parameters out of range\n");
#endif
		return 0.0;
	}
}

