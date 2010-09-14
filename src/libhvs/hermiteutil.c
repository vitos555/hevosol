#include "hermiteutil.h"
#include "factorialutial.h"
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

FLOAT_TYPE h1(alpha1,alpha2,lambda_sq) {
	if ((alpha_1%2==1) && (alpha_2%2==0)) {
		int alpha_temp = (alpha1+alpha2-1)/2;
		int negone = (alpha_temp%2==0?1:-1);
		return
			(-0.25)*PI_INV/lambda_sq*negone*
			factorial(alpha1)/factorial(alpha_temp)*
			factorial(alpha2)*binomial(alpha_temp,(int)alpha1/2);
	} else {
		return 0.0;
	}
}

FLOAT_TYPE h2(alpha1,alpha2,lambda_sq) {
	if ((alpha_1%2==1) && (alpha_2%2==0)) {
		int alpha_temp = (alpha1+alpha2-1)/2;
		int negone = (alpha_temp%2==0?1:-1);
		return
			(0.25)*PI_INV/lambda_sq*negone*
			factorial(alpha1)/factorial(alpha_temp)*
			factorial(alpha2)*binomial(alpha_temp,(UINT)(alpha1-1)/2);
	} else {
		return 0.0;
	}
}

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
			return log(power*exp(x));
	}
}

FLOAT_TYPE hb1(FLOAT_TYPE x1,FLOAT_TYPE x2,FLOAT_TYPE lambdasq,UINT k1,UINT k2) {
	FLOAT_TYPE he0 = exp(-(x1*x1+x2*x2)/lambda_sq);
	FLOAT_TYPE Pi = M_PI;
	if ( (k1==0) && (k2==0) )
		return ((-1 + he0)*x2)/(4.*Pi*Power(x1,2));
	else if ( (k1==1) && (k2==0) )
		return -((2*(-1 + he0)*lambdasq + he0*Power(x1,2))*x2)/
   (4.*lambdasq*Pi*Power(x1,3));
	else if ( (k1==0) && (k2==1) )
		return -(lambdasq - he0*lambdasq + he0*Power(x2,2))/(4.*lambdasq*Pi*Power(x1,2));
	else if ( (k1==2) && (k2==0) )
		return ((6*(-1 + he0)*Power(lambdasq,2) + 3*he0*lambdasq*Power(x1,2) + 
       he0*Power(x1,4))*x2)/(4.*Power(lambdasq,2)*Pi*Power(x1,4));
	else if ( (k1==1) && (k2==1) )
		return (-(lambdasq*(2*(-1 + he0)*lambdasq + he0*Power(x1,2))) + 
     he0*(2*lambdasq + Power(x1,2))*Power(x2,2))/
   (4.*Power(lambdasq,2)*Pi*Power(x1,3));
	else if ( (k1==0) && (k2==2) )
		return (he0*x2*(-3*lambdasq + Power(x2,2)))/
   (4.*Power(lambdasq,2)*Pi*Power(x1,2));
	else
		return 0.0;
}

FLOAT_TYPE hb2(FLOAT_TYPE x1,FLOAT_TYPE x2,FLOAT_TYPE lambdasq,UINT k1,UINT k2) {
	FLOAT_TYPE he0 = exp(-(x1*x1+x2*x2)/lambda_sq);
	FLOAT_TYPE Pi = M_PI;
	if ( (k1==0) && (k2==0) )
		return (x1 - he0*x1)/(2*Pi*Power(x1,2) + 2*Pi*Power(x2,2));
	else if ( (k1==1) && (k2==0) )
		return ((-1 + he0)*lambdasq*(Power(x1,2) - Power(x2,2)) + 
     he0*Power(x1,2)*(Power(x1,2) + Power(x2,2)))/
   (2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==0) && (k2==1) )
		return (x1*x2*(2*(-1 + he0)*lambdasq + he0*(Power(x1,2) + Power(x2,2))))/
   (2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==2) && (k2==0) )
		return -(x1*(2*(-1 + he0)*Power(lambdasq,2)*(Power(x1,2) - 3*Power(x2,2)) + 
        he0*lambdasq*(Power(x1,2) - 3*Power(x2,2))*
         (Power(x1,2) + Power(x2,2)) + 
        he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),2)))/
   (2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + Power(x2,2),3));
	else if ( (k1==1) && (k2==1) )
		return (x2*(-2*(-1 + he0)*Power(lambdasq,2)*(3*Power(x1,2) - Power(x2,2)) - 
       he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),2) + 
       he0*lambdasq*(-3*Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + 
          Power(x2,4))))/
   (2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + Power(x2,2),3));
	else if ( (k1==0) && (k2==2) )
		return (x1*(2*(-1 + he0)*Power(lambdasq,2)*(Power(x1,2) - 3*Power(x2,2)) + 
       he0*lambdasq*(Power(x1,2) - 3*Power(x2,2))*
        (Power(x1,2) + Power(x2,2)) - 
       he0*Power(x2,2)*Power(Power(x1,2) + Power(x2,2),2)))/
   (2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + Power(x2,2),3));
	else
		return 0.0;
}

