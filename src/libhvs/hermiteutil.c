#include "hermiteutil.h"
#include "factorialutil.h"
#include <math.h>

#define PI_INV 1/M_PI

#if NMOMENTS > 2
#error "Current implementation doesn't support moments of order higher then 2."
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

FLOAT_TYPE he_00(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda_sq) {
	return PI_INV/lambda_sq*M_EXP(-(x1*x1+x2*x2)/lambda_sq);
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

FLOAT_TYPE h1(UINT alpha1, UINT alpha2, FLOAT_TYPE lambda_sq) {
	if ((alpha1%2==0) && (alpha2%2==1)) {
		int alpha_temp = (alpha1+alpha2-1)/2;
		int negone = (alpha_temp%2==0?1:-1);
		return
			(-0.5)*PI_INV/Power(2*lambda_sq,alpha_temp+1)*negone*
			factorial(alpha1)/factorial(alpha_temp+1)*
			factorial(alpha2)*binomial(alpha_temp,(int)alpha1/2);
	} else {
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
		return 0.0;
	}
}

FLOAT_TYPE hb1(FLOAT_TYPE x1,FLOAT_TYPE x2,FLOAT_TYPE lambdasq,UINT k1,UINT k2) {
	FLOAT_TYPE he0 = M_EXP(-(x1*x1+x2*x2)/lambdasq*0.5);
	FLOAT_TYPE Pi = M_PI;
	if ( (k1==0) && (k2==0) )
		return ((-1 + he0)*x2)/(2.*Pi*(Power(x1,2) + Power(x2,2)));
	else if ( (k1==0) && (k2==1) )
		return ((-1 + he0)*lambdasq*(Power(x1,2) - Power(x2,2)) - \
he0*Power(x2,2)*(Power(x1,2) + \
Power(x2,2)))/(2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==0) && (k2==2) )
		return (x2*(-2*(-1 + he0)*Power(lambdasq,2)*(3*Power(x1,2) - \
Power(x2,2)) + he0*Power(x2,2)*Power(Power(x1,2) + Power(x2,2),2) + \
he0*lambdasq*(-3*Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4))))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==0) && (k2==3) )
		return -(-6*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),2) + he0*Power(x2,4)*Power(Power(x1,2) + Power(x2,2),3) \
+ 6*(-1 + he0)*Power(lambdasq,3)*(Power(x1,4) - \
6*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
3*he0*Power(lambdasq,2)*(Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + \
Power(x2,6)))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==1) && (k2==0) )
		return -(x1*x2*(2*(-1 + he0)*lambdasq + he0*(Power(x1,2) + \
Power(x2,2))))/(2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==1) && (k2==1) )
		return (x1*(-2*(-1 + he0)*Power(lambdasq,2)*(Power(x1,2) - \
3*Power(x2,2)) + he0*Power(x2,2)*Power(Power(x1,2) + Power(x2,2),2) + \
he0*lambdasq*(-Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4))))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==1) && (k2==2) )
		return -(x1*x2*(-24*(-1 + he0)*Power(lambdasq,3)*(Power(x1,2) - \
Power(x2,2)) - 3*he0*lambdasq*(Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2) + \
he0*Power(x2,2)*Power(Power(x1,2) + Power(x2,2),3) + \
12*he0*Power(lambdasq,2)*(-Power(x1,4) + \
Power(x2,4))))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==1) && (k2==3) )
		return (x1*(2*he0*lambdasq*Power(x2,2)*(-3*Power(x1,2) + \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3) + \
he0*Power(x2,4)*Power(Power(x1,2) + Power(x2,2),4) + 24*(-1 + \
he0)*Power(lambdasq,4)*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 12*he0*Power(lambdasq,3)*(Power(x1,6) - \
9*Power(x1,4)*Power(x2,2) - 5*Power(x1,2)*Power(x2,4) + \
5*Power(x2,6)) + 3*he0*Power(lambdasq,2)*(Power(x1,8) - \
8*Power(x1,6)*Power(x2,2) - 14*Power(x1,4)*Power(x2,4) + \
5*Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==2) && (k2==0) )
		return (x2*(2*(-1 + he0)*Power(lambdasq,2)*(3*Power(x1,2) - \
Power(x2,2)) + he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),2) + \
he0*lambdasq*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) - \
Power(x2,4))))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==2) && (k2==1) )
		return (-(he0*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),3)) + 6*(-1 + he0)*Power(lambdasq,3)*(Power(x1,4) - \
6*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
he0*lambdasq*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,4) - \
4*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
3*he0*Power(lambdasq,2)*(Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + \
Power(x2,6)))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==2) && (k2==2) )
		return (he0*x2*(Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4) + 24*(-1 + 1/he0)*Power(lambdasq,4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4)) - \
lambdasq*Power(Power(x1,2) + Power(x2,2),3)*(3*Power(x1,4) - \
4*Power(x1,2)*Power(x2,2) + Power(x2,4)) - \
12*Power(lambdasq,3)*(5*Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
9*Power(x1,2)*Power(x2,4) + Power(x2,6)) - \
3*Power(lambdasq,2)*(5*Power(x1,8) - 14*Power(x1,4)*Power(x2,4) - \
8*Power(x1,2)*Power(x2,6) + \
Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==2) && (k2==3) )
		return -(he0*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5) - he0*lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(6*Power(x1,4) - 3*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,6) - 36*Power(x1,4)*Power(x2,2) + \
39*Power(x1,2)*Power(x2,4) - 2*Power(x2,6)) + 120*(-1 + \
he0)*Power(lambdasq,5)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
15*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
60*he0*Power(lambdasq,4)*(Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
14*Power(x1,2)*Power(x2,6) - \
Power(x2,8)))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==3) && (k2==0) )
		return -(x1*x2*(24*(-1 + he0)*Power(lambdasq,3)*(Power(x1,2) - \
Power(x2,2)) + 3*he0*lambdasq*(Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2) + \
he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),3) + \
12*he0*Power(lambdasq,2)*(Power(x1,4) - \
Power(x2,4))))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==3) && (k2==1) )
		return (he0*x1*(Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4) - lambdasq*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 24*(-1 + 1/he0)*Power(lambdasq,4)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) - \
12*Power(lambdasq,3)*(Power(x1,6) - 9*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) - \
3*Power(lambdasq,2)*(Power(x1,8) - 8*Power(x1,6)*Power(x2,2) - \
14*Power(x1,4)*Power(x2,4) + \
5*Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==3) && (k2==2) )
		return -(he0*x1*x2*(Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + 1/he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) - \
30*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 5*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - lambdasq*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 120*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==3) && (k2==3) )
		return (he0*x1*(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6) - 3*lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(2*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + (720*(-1 + he0)*Power(lambdasq,6)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)))/he0 + 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) + \
90*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,6) \
- 21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) + \
3*Power(lambdasq,2)*Power(Power(x1,2) + Power(x2,2),4)*(Power(x1,6) - \
12*Power(x1,4)*Power(x2,2) + 23*Power(x1,2)*Power(x2,4) - \
4*Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else {
		return 0.0;
	}
}

FLOAT_TYPE hb2(FLOAT_TYPE x1,FLOAT_TYPE x2,FLOAT_TYPE lambdasq,UINT k1,UINT k2) {
	FLOAT_TYPE he0 = M_EXP(-(x1*x1+x2*x2)/lambdasq*0.5);
	FLOAT_TYPE Pi = M_PI;
	if ( (k1==0) && (k2==0) )
		return (x1 - he0*x1)/(2*Pi*Power(x1,2) + 2*Pi*Power(x2,2));
	else if ( (k1==0) && (k2==1) )
		return (x1*x2*(2*(-1 + he0)*lambdasq + he0*(Power(x1,2) + \
Power(x2,2))))/(2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==0) && (k2==2) )
		return (x1*(2*(-1 + he0)*Power(lambdasq,2)*(Power(x1,2) - \
3*Power(x2,2)) + he0*lambdasq*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2)) - \
he0*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),2)))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==0) && (k2==3) )
		return (x1*x2*(-24*(-1 + he0)*Power(lambdasq,3)*(Power(x1,2) - \
Power(x2,2)) - 3*he0*lambdasq*(Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2) + \
he0*Power(x2,2)*Power(Power(x1,2) + Power(x2,2),3) + \
12*he0*Power(lambdasq,2)*(-Power(x1,4) + \
Power(x2,4))))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==1) && (k2==0) )
		return ((-1 + he0)*lambdasq*(Power(x1,2) - Power(x2,2)) + \
he0*Power(x1,2)*(Power(x1,2) + \
Power(x2,2)))/(2.*lambdasq*Pi*Power(Power(x1,2) + Power(x2,2),2));
	else if ( (k1==1) && (k2==1) )
		return (x2*(-2*(-1 + he0)*Power(lambdasq,2)*(3*Power(x1,2) - \
Power(x2,2)) - he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),2) + \
he0*lambdasq*(-3*Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4))))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==1) && (k2==2) )
		return -(-(he0*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),3)) + 6*(-1 + he0)*Power(lambdasq,3)*(Power(x1,4) - \
6*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
he0*lambdasq*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,4) - \
4*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
3*he0*Power(lambdasq,2)*(Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + \
Power(x2,6)))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==1) && (k2==3) )
		return (x2*(-(he0*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)) + 24*(-1 + he0)*Power(lambdasq,4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
he0*lambdasq*Power(Power(x1,2) + Power(x2,2),3)*(3*Power(x1,4) - \
4*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
12*he0*Power(lambdasq,3)*(5*Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
9*Power(x1,2)*Power(x2,4) + Power(x2,6)) + \
3*he0*Power(lambdasq,2)*(5*Power(x1,8) - 14*Power(x1,4)*Power(x2,4) - \
8*Power(x1,2)*Power(x2,6) + \
Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==2) && (k2==0) )
		return -(x1*(2*(-1 + he0)*Power(lambdasq,2)*(Power(x1,2) - \
3*Power(x2,2)) + he0*lambdasq*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2)) + \
he0*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),2)))/(2.*Power(lambdasq,2)*Pi*Power(Power(x1,2) + \
Power(x2,2),3));
	else if ( (k1==2) && (k2==1) )
		return (x1*x2*(24*(-1 + he0)*Power(lambdasq,3)*(Power(x1,2) - \
Power(x2,2)) + 3*he0*lambdasq*(Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2) + \
he0*Power(x1,2)*Power(Power(x1,2) + Power(x2,2),3) + \
12*he0*Power(lambdasq,2)*(Power(x1,4) - \
Power(x2,4))))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==2) && (k2==2) )
		return (x1*(-(he0*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)) + he0*lambdasq*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 24*(-1 + he0)*Power(lambdasq,4)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) + \
12*he0*Power(lambdasq,3)*(Power(x1,6) - 9*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) + \
3*he0*Power(lambdasq,2)*(Power(x1,8) - 8*Power(x1,6)*Power(x2,2) - \
14*Power(x1,4)*Power(x2,4) + \
5*Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==2) && (k2==3) )
		return (he0*x1*x2*(Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + 1/he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) - \
30*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 5*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - lambdasq*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 120*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==3) && (k2==0) )
		return (-6*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),2) + he0*Power(x1,4)*Power(Power(x1,2) + Power(x2,2),3) + \
6*(-1 + he0)*Power(lambdasq,3)*(Power(x1,4) - \
6*Power(x1,2)*Power(x2,2) + Power(x2,4)) + \
3*he0*Power(lambdasq,2)*(Power(x1,6) - 5*Power(x1,4)*Power(x2,2) - \
5*Power(x1,2)*Power(x2,4) + \
Power(x2,6)))/(2.*Power(lambdasq,3)*Pi*Power(Power(x1,2) + \
Power(x2,2),4));
	else if ( (k1==3) && (k2==1) )
		return -(x2*(2*he0*lambdasq*Power(x1,2)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3) + \
he0*Power(x1,4)*Power(Power(x1,2) + Power(x2,2),4) + 24*(-1 + \
he0)*Power(lambdasq,4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 12*he0*Power(lambdasq,3)*(5*Power(x1,6) - \
5*Power(x1,4)*Power(x2,2) - 9*Power(x1,2)*Power(x2,4) + Power(x2,6)) \
+ 3*he0*Power(lambdasq,2)*(5*Power(x1,8) - 14*Power(x1,4)*Power(x2,4) \
- 8*Power(x1,2)*Power(x2,6) + \
Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==3) && (k2==2) )
		return -(-(he0*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)) + he0*lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,4) - 3*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(2*Power(x1,6) - 39*Power(x1,4)*Power(x2,2) + \
36*Power(x1,2)*Power(x2,4) - 3*Power(x2,6)) + 120*(-1 + \
he0)*Power(lambdasq,5)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
15*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
60*he0*Power(lambdasq,4)*(Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
14*Power(x1,2)*Power(x2,6) - \
Power(x2,8)))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==3) && (k2==3) )
		return -(he0*x2*(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6) - 3*lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) - Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) - 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(4*Power(x1,6) - 23*Power(x1,4)*Power(x2,2) + \
12*Power(x1,2)*Power(x2,4) - Power(x2,6)) + 720*(-1 + \
1/he0)*Power(lambdasq,6)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) \
+ 21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
360*Power(lambdasq,5)*(Power(x1,2) + Power(x2,2))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 90*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),2)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else {
		return 0.0;
	}
}

