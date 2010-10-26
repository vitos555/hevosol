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
	else if ( (k1==0) && (k2==4) )
		return (x2*(he0*Power(x2,4)*Power(Power(x1,2) + Power(x2,2),4) - \
2*he0*lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(5*Power(x1,2) + Power(x2,2)) + 24*(-1 + \
he0)*Power(lambdasq,4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 12*he0*Power(lambdasq,3)*(5*Power(x1,6) - \
5*Power(x1,4)*Power(x2,2) - 9*Power(x1,2)*Power(x2,4) + Power(x2,6)) \
+ 3*he0*Power(lambdasq,2)*(5*Power(x1,8) - 14*Power(x1,4)*Power(x2,4) \
- 8*Power(x1,2)*Power(x2,6) + \
Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==0) && (k2==5) )
		return (-(he0*Power(x2,6)*Power(Power(x1,2) + Power(x2,2),5)) + \
5*he0*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,2) + Power(x2,2)) - \
5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*Power(-3*Power(x1,2)*x2 + Power(x2,3),2) + 120*(-1 + \
he0)*Power(lambdasq,5)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
15*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
60*he0*Power(lambdasq,4)*(Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
14*Power(x1,2)*Power(x2,6) - \
Power(x2,8)))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
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
	else if ( (k1==1) && (k2==4) )
		return -(x1*x2*(-10*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(\
Power(x1,2) + Power(x2,2),4) + he0*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) + \
30*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 120*he0*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==1) && (k2==5) )
		return (he0*x1*(Power(x2,6)*Power(Power(x1,2) + Power(x2,2),6) - \
3*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,2) + Power(x2,2)) + \
15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 720*(-1 + 1/he0)*Power(lambdasq,6)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) - 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) - \
90*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,6) \
- 21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) - 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
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
	else if ( (k1==2) && (k2==4) )
		return (he0*x2*(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6) - lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(10*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + (720*(-1 + he0)*Power(lambdasq,6)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
90*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),2)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else if ( (k1==2) && (k2==5) )
		return (he0*(-(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),7)) + lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(15*Power(x1,6) - 25*Power(x1,4)*Power(x2,2) + \
17*Power(x1,2)*Power(x2,4) + Power(x2,6)) + (5040*(-1 + \
he0)*Power(lambdasq,7)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)))/he0 + 630*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 24*Power(x1,6)*Power(x2,2) + \
62*Power(x1,4)*Power(x2,4) - 24*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 2520*Power(lambdasq,6)*(Power(x1,10) - \
27*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
42*Power(x1,4)*Power(x2,6) - 27*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
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
	else if ( (k1==3) && (k2==4) )
		return -(he0*x1*x2*(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(10*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,6) - 20*Power(x1,4)*Power(x2,2) + \
29*Power(x1,2)*Power(x2,4) - 2*Power(x2,6)) + (40320*(-1 + \
he0)*Power(lambdasq,7)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)))/he0 + \
5040*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
840*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) \
- 7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 105*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
20160*Power(lambdasq,6)*(Power(x1,8) - 6*Power(x1,6)*Power(x2,2) + \
6*Power(x1,2)*Power(x2,6) - \
Power(x2,8))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==3) && (k2==5) )
		return (he0*x1*(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8) - lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(15*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) + Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(45*Power(x1,6) - 75*Power(x1,4)*Power(x2,2) + \
107*Power(x1,2)*Power(x2,4) + 3*Power(x2,6)) + 40320*(-1 + \
1/he0)*Power(lambdasq,8)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) - 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) - 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) - 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,8) - 120*Power(x1,6)*Power(x2,2) + \
450*Power(x1,4)*Power(x2,4) - 288*Power(x1,2)*Power(x2,6) + \
33*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==4) && (k2==0) )
		return (x2*(2*he0*lambdasq*Power(x1,2)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3) + \
he0*Power(x1,4)*Power(Power(x1,2) + Power(x2,2),4) + 24*(-1 + \
he0)*Power(lambdasq,4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 12*he0*Power(lambdasq,3)*(5*Power(x1,6) - \
5*Power(x1,4)*Power(x2,2) - 9*Power(x1,2)*Power(x2,4) + Power(x2,6)) \
+ 3*he0*Power(lambdasq,2)*(5*Power(x1,8) - 14*Power(x1,4)*Power(x2,4) \
- 8*Power(x1,2)*Power(x2,6) + \
Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==4) && (k2==1) )
		return (-(he0*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
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
	else if ( (k1==4) && (k2==2) )
		return (he0*x2*(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
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
	else if ( (k1==4) && (k2==3) )
		return -(he0*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7) - \
2*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 12*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 31*Power(x1,6)*Power(x2,2) + \
76*Power(x1,4)*Power(x2,4) - 31*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 5040*(-1 + he0)*Power(lambdasq,7)*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
630*he0*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*he0*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 3*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,8) - 10*Power(x1,6)*Power(x2,2) + \
34*Power(x1,4)*Power(x2,4) - 10*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 2520*he0*Power(lambdasq,6)*(Power(x1,10) - \
27*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
42*Power(x1,4)*Power(x2,6) - 27*Power(x1,2)*Power(x2,8) + \
Power(x2,10)))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==4) && (k2==4) )
		return (he0*x2*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(5*Power(x1,4) + 3*Power(x2,4)) + \
6*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(15*Power(x1,8) - 150*Power(x1,6)*Power(x2,2) + \
216*Power(x1,4)*Power(x2,4) - 66*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + (40320*(-1 + he0)*Power(lambdasq,8)*(9*Power(x1,8) - \
84*Power(x1,6)*Power(x2,2) + 126*Power(x1,4)*Power(x2,4) - \
36*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 50*Power(x1,6)*Power(x2,2) + \
138*Power(x1,4)*Power(x2,4) - 18*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==4) && (k2==5) )
		return (he0*(-(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)) + \
3*lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(5*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(15*Power(x1,8) - 20*Power(x1,6)*Power(x2,2) + \
60*Power(x1,4)*Power(x2,4) + Power(x2,8)) + \
9*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(10*Power(x1,10) - 475*Power(x1,8)*Power(x2,2) + \
2200*Power(x1,6)*Power(x2,4) - 2210*Power(x1,4)*Power(x2,6) + \
470*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + (362880*(-1 + \
he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)))/he0 + \
45360*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
7560*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,10) - 105*Power(x1,8)*Power(x2,2) + \
580*Power(x1,6)*Power(x2,4) - 520*Power(x1,4)*Power(x2,6) + \
135*Power(x1,2)*Power(x2,8) + Power(x2,10)) + \
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==5) && (k2==0) )
		return -(x1*x2*(-10*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(\
Power(x1,2) + Power(x2,2),4) + he0*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) + \
30*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 120*he0*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==5) && (k2==1) )
		return (he0*x1*(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6) + 15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) - Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) + 720*(-1 + 1/he0)*Power(lambdasq,6)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) - 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) - \
90*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,6) \
- 21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) - 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else if ( (k1==5) && (k2==2) )
		return -(he0*x1*x2*(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) - 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(2*Power(x1,6) - 29*Power(x1,4)*Power(x2,2) + \
20*Power(x1,2)*Power(x2,4) - 5*Power(x2,6)) + 40320*(-1 + \
1/he0)*Power(lambdasq,7)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
5040*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
840*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) \
- 7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 105*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
20160*Power(lambdasq,6)*(Power(x1,8) - 6*Power(x1,6)*Power(x2,2) + \
6*Power(x1,2)*Power(x2,6) - \
Power(x2,8))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==5) && (k2==3) )
		return (he0*x1*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(3*Power(x1,4) + 5*Power(x2,4)) + \
20160*Power(lambdasq,7)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6)) + (40320*(-1 + \
he0)*Power(lambdasq,8)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)))/he0 + 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 6*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,8) - 66*Power(x1,6)*Power(x2,2) + \
216*Power(x1,4)*Power(x2,4) - 150*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 18*Power(x1,6)*Power(x2,2) + \
138*Power(x1,4)*Power(x2,4) - 50*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==5) && (k2==4) )
		return -(he0*x1*x2*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9) + (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)))/he0 + \
362880*Power(lambdasq,8)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) + \
90720*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 15120*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 1890*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(5*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 12*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 85*Power(x1,6)*Power(x2,2) + \
156*Power(x1,4)*Power(x2,4) - 85*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) + 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,8) - 10*Power(x1,6)*Power(x2,2) + \
66*Power(x1,4)*Power(x2,4) - 10*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==5) && (k2==5) )
		return (he0*x1*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10) - \
5*lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) + 15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,8) - 2*Power(x1,6)*Power(x2,2) + \
18*Power(x1,4)*Power(x2,4) + Power(x2,8)) - \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(4*Power(x1,10) - 295*Power(x1,8)*Power(x2,2) + \
1720*Power(x1,6)*Power(x2,4) - 2438*Power(x1,4)*Power(x2,6) + \
860*Power(x1,2)*Power(x2,8) - 59*Power(x2,10)) + 3628800*(-1 + \
1/he0)*Power(lambdasq,10)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) \
+ 330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) - \
1814400*Power(lambdasq,9)*(Power(x1,2) + Power(x2,2))*(Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 330*Power(x1,6)*Power(x2,4) - \
462*Power(x1,4)*Power(x2,6) + 165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10)) - 453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) - \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) - \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
- 15*Power(x1,8)*Power(x2,2) + 140*Power(x1,6)*Power(x2,4) - \
152*Power(x1,4)*Power(x2,6) + 75*Power(x1,2)*Power(x2,8) - \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else {
		printf("Err:%d,%d\n",k1,k2);
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
	else if ( (k1==0) && (k2==4) )
		return -(x1*(2*he0*lambdasq*Power(x2,2)*(-3*Power(x1,2) + \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3) + \
he0*Power(x2,4)*Power(Power(x1,2) + Power(x2,2),4) + 24*(-1 + \
he0)*Power(lambdasq,4)*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 12*he0*Power(lambdasq,3)*(Power(x1,6) - \
9*Power(x1,4)*Power(x2,2) - 5*Power(x1,2)*Power(x2,4) + \
5*Power(x2,6)) + 3*he0*Power(lambdasq,2)*(Power(x1,8) - \
8*Power(x1,6)*Power(x2,2) - 14*Power(x1,4)*Power(x2,4) + \
5*Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==0) && (k2==5) )
		return (x1*x2*(-10*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(\
x1,2) + Power(x2,2),4) + he0*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) + \
30*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 120*he0*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
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
	else if ( (k1==1) && (k2==4) )
		return (he0*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
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
	else if ( (k1==1) && (k2==5) )
		return (he0*x2*(-(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)) + lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(10*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 720*(-1 + 1/he0)*Power(lambdasq,6)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
90*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),2)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
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
	else if ( (k1==2) && (k2==4) )
		return -(he0*x1*(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
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
	else if ( (k1==2) && (k2==5) )
		return (he0*x1*x2*(Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(10*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,6) - 20*Power(x1,4)*Power(x2,2) + \
29*Power(x1,2)*Power(x2,4) - 2*Power(x2,6)) + (40320*(-1 + \
he0)*Power(lambdasq,7)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)))/he0 + \
5040*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
840*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) \
- 7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 105*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
20160*Power(lambdasq,6)*(Power(x1,8) - 6*Power(x1,6)*Power(x2,2) + \
6*Power(x1,2)*Power(x2,6) - \
Power(x2,8))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
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
	else if ( (k1==3) && (k2==4) )
		return (he0*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),6)*(3*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 12*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 31*Power(x1,6)*Power(x2,2) + \
76*Power(x1,4)*Power(x2,4) - 31*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + (5040*(-1 + he0)*Power(lambdasq,7)*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
630*Power(lambdasq,5)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,8) \
- 28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
105*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,8) \
- 28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
3*Power(lambdasq,2)*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,8) - \
10*Power(x1,6)*Power(x2,2) + 34*Power(x1,4)*Power(x2,4) - \
10*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
2520*Power(lambdasq,6)*(Power(x1,10) - 27*Power(x1,8)*Power(x2,2) + \
42*Power(x1,6)*Power(x2,4) + 42*Power(x1,4)*Power(x2,6) - \
27*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==3) && (k2==5) )
		return -(he0*x2*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(5*Power(x1,4) + 3*Power(x2,4)) + \
6*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(15*Power(x1,8) - 150*Power(x1,6)*Power(x2,2) + \
216*Power(x1,4)*Power(x2,4) - 66*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + (40320*(-1 + he0)*Power(lambdasq,8)*(9*Power(x1,8) - \
84*Power(x1,6)*Power(x2,2) + 126*Power(x1,4)*Power(x2,4) - \
36*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 50*Power(x1,6)*Power(x2,2) + \
138*Power(x1,4)*Power(x2,4) - 18*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==4) && (k2==0) )
		return -(x1*(he0*Power(x1,4)*Power(Power(x1,2) + Power(x2,2),4) - \
2*he0*lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,2) + 5*Power(x2,2)) + 24*(-1 + \
he0)*Power(lambdasq,4)*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 12*he0*Power(lambdasq,3)*(Power(x1,6) - \
9*Power(x1,4)*Power(x2,2) - 5*Power(x1,2)*Power(x2,4) + \
5*Power(x2,6)) + 3*he0*Power(lambdasq,2)*(Power(x1,8) - \
8*Power(x1,6)*Power(x2,2) - 14*Power(x1,4)*Power(x2,4) + \
5*Power(x2,8))))/(2.*Power(lambdasq,4)*Pi*Power(Power(x1,2) + \
Power(x2,2),5));
	else if ( (k1==4) && (k2==1) )
		return (x1*x2*(-10*he0*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(\
x1,2) + Power(x2,2),4) + he0*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5) + 240*(-1 + he0)*Power(lambdasq,5)*(3*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 3*Power(x2,4)) + \
30*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 120*he0*Power(lambdasq,4)*(3*Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) - 7*Power(x1,2)*Power(x2,4) + \
3*Power(x2,6))))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==4) && (k2==2) )
		return (he0*x1*(-(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)) - 15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) - Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) + (720*(-1 + he0)*Power(lambdasq,6)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)))/he0 + 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) + \
90*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,6) \
- 21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else if ( (k1==4) && (k2==3) )
		return (he0*x1*x2*(Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,4) - Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) - 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(2*Power(x1,6) - 29*Power(x1,4)*Power(x2,2) + \
20*Power(x1,2)*Power(x2,4) - 5*Power(x2,6)) + 40320*(-1 + \
1/he0)*Power(lambdasq,7)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
5040*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
840*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) \
- 7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 105*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
7*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
20160*Power(lambdasq,6)*(Power(x1,8) - 6*Power(x1,6)*Power(x2,2) + \
6*Power(x1,2)*Power(x2,6) - \
Power(x2,8))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==4) && (k2==4) )
		return -(he0*x1*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(3*Power(x1,4) + 5*Power(x2,4)) + \
20160*Power(lambdasq,7)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6)) + (40320*(-1 + \
he0)*Power(lambdasq,8)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)))/he0 + 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 36*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 6*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,8) - 66*Power(x1,6)*Power(x2,2) + \
216*Power(x1,4)*Power(x2,4) - 150*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 18*Power(x1,6)*Power(x2,2) + \
138*Power(x1,4)*Power(x2,4) - 50*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==4) && (k2==5) )
		return (he0*x1*x2*(Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9) + (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)))/he0 + \
362880*Power(lambdasq,8)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) + \
90720*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 15120*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 1890*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 2*lambdasq*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(5*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 12*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 85*Power(x1,6)*Power(x2,2) + \
156*Power(x1,4)*Power(x2,4) - 85*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) + 3*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,8) - 10*Power(x1,6)*Power(x2,2) + \
66*Power(x1,4)*Power(x2,4) - 10*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==5) && (k2==0) )
		return (he0*Power(x1,6)*Power(Power(x1,2) + Power(x2,2),5) - \
5*he0*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,2) + 3*Power(x2,2)) + \
5*he0*Power(lambdasq,2)*Power(Power(x1,2) + \
Power(x2,2),3)*Power(Power(x1,3) - 3*x1*Power(x2,2),2) + 120*(-1 + \
he0)*Power(lambdasq,5)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
15*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 15*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - Power(x2,6)) + \
60*he0*Power(lambdasq,4)*(Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
14*Power(x1,2)*Power(x2,6) - \
Power(x2,8)))/(2.*Power(lambdasq,5)*Pi*Power(Power(x1,2) + \
Power(x2,2),6));
	else if ( (k1==5) && (k2==1) )
		return (he0*x2*(-(Power(x1,6)*Power(Power(x1,2) + Power(x2,2),6)) + \
3*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,2) + 5*Power(x2,2)) - \
15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 720*(-1 + 1/he0)*Power(lambdasq,6)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 360*Power(lambdasq,5)*(Power(x1,2) + \
Power(x2,2))*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
90*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),2)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - Power(x2,6)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),3)*(7*Power(x1,6) - 35*Power(x1,4)*Power(x2,2) + \
21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else if ( (k1==5) && (k2==2) )
		return (he0*(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) + 17*Power(x1,4)*Power(x2,2) - \
25*Power(x1,2)*Power(x2,4) + 15*Power(x2,6)) + 5040*(-1 + \
1/he0)*Power(lambdasq,7)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 630*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,8) - 24*Power(x1,6)*Power(x2,2) + \
62*Power(x1,4)*Power(x2,4) - 24*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 2520*Power(lambdasq,6)*(Power(x1,10) - \
27*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
42*Power(x1,4)*Power(x2,6) - 27*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==5) && (k2==3) )
		return (he0*x2*(-(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)) + lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,6) + 107*Power(x1,4)*Power(x2,2) - \
75*Power(x1,2)*Power(x2,4) + 45*Power(x2,6)) + (40320*(-1 + \
he0)*Power(lambdasq,8)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)))/he0 + 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(33*Power(x1,8) - 288*Power(x1,6)*Power(x2,2) + \
450*Power(x1,4)*Power(x2,4) - 120*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==5) && (k2==4) )
		return (he0*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9) - 3*lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(2*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,8) + 60*Power(x1,4)*Power(x2,4) - \
20*Power(x1,2)*Power(x2,6) + 15*Power(x2,8)) + \
9*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(11*Power(x1,10) - 470*Power(x1,8)*Power(x2,2) + \
2210*Power(x1,6)*Power(x2,4) - 2200*Power(x1,4)*Power(x2,6) + \
475*Power(x1,2)*Power(x2,8) - 10*Power(x2,10)) + (362880*(-1 + \
he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)))/he0 + \
45360*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
7560*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
3*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),6)*(Power(x1,10) \
+ 135*Power(x1,8)*Power(x2,2) - 520*Power(x1,6)*Power(x2,4) + \
580*Power(x1,4)*Power(x2,6) - 105*Power(x1,2)*Power(x2,8) + \
5*Power(x2,10)) + 181440*Power(lambdasq,8)*(Power(x1,12) - \
44*Power(x1,10)*Power(x2,2) + 165*Power(x1,8)*Power(x2,4) - \
165*Power(x1,4)*Power(x2,8) + 44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==5) && (k2==5) )
		return (he0*x2*(-(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)) + \
5*lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(2*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,8) + 18*Power(x1,4)*Power(x2,4) - \
2*Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) - \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(59*Power(x1,10) - 860*Power(x1,8)*Power(x2,2) + \
2438*Power(x1,6)*Power(x2,4) - 1720*Power(x1,4)*Power(x2,6) + \
295*Power(x1,2)*Power(x2,8) - 4*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
- 75*Power(x1,8)*Power(x2,2) + 152*Power(x1,6)*Power(x2,4) - \
140*Power(x1,4)*Power(x2,6) + 15*Power(x1,2)*Power(x2,8) - \
Power(x2,10)) + 3628800*(-1 + \
1/he0)*Power(lambdasq,10)*(11*Power(x1,10) - \
165*Power(x1,8)*Power(x2,2) + 462*Power(x1,6)*Power(x2,4) - \
330*Power(x1,4)*Power(x2,6) + 55*Power(x1,2)*Power(x2,8) - \
Power(x2,10)) - 1814400*Power(lambdasq,9)*(Power(x1,2) + \
Power(x2,2))*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else {
		printf("Err:%d,%d\n",k1,k2);
		return 0.0;
	}
}

