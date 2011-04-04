// Hermite functions up until 7th order

FLOAT_TYPE Power(FLOAT_TYPE x, int power);

FLOAT_TYPE he(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambdasq, unsigned short k1, unsigned short k2) {
	FLOAT_TYPE he0 = M_EXP(-(x1*x1+x2*x2)/lambdasq*0.5);
	FLOAT_TYPE Pi = M_PI;
	if ( (k1==0) && (k2==0) )
		return 1/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*lambdasq*Pi);
	else if ( (k1==0) && (k2==1) )
		return (-2*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,2)*Pi);
	else if ( (k1==0) && (k2==2) )
		return (-2*lambdasq + 4*Power(x2,2))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,3)*Pi);
	else if ( (k1==0) && (k2==3) )
		return (12*lambdasq*x2 - 8*Power(x2,3))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,4)*Pi);
	else if ( (k1==0) && (k2==4) )
		return (4*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,5)*Pi);
	else if ( (k1==0) && (k2==5) )
		return (-8*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==0) && (k2==6) )
		return (-8*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x2,2) \
+ 60*lambdasq*Power(x2,4) - 8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==0) && (k2==7) )
		return (16*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==1) && (k2==0) )
		return (-2*x1)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,2)*Pi);
	else if ( (k1==1) && (k2==1) )
		return (4*x1*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,3)*Pi);
	else if ( (k1==1) && (k2==2) )
		return (4*x1*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,4)*Pi);
	else if ( (k1==1) && (k2==3) )
		return (8*x1*x2*(-3*lambdasq + \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,5)*Pi);
	else if ( (k1==1) && (k2==4) )
		return (-8*x1*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==1) && (k2==5) )
		return (16*x1*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==1) && (k2==6) )
		return (16*x1*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==1) && (k2==7) )
		return (32*x1*x2*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x2,2) - 84*lambdasq*Power(x2,4) + \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==2) && (k2==0) )
		return (-2*lambdasq + 4*Power(x1,2))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,3)*Pi);
	else if ( (k1==2) && (k2==1) )
		return (4*(lambdasq - 2*Power(x1,2))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,4)*Pi);
	else if ( (k1==2) && (k2==2) )
		return (4*(lambdasq - 2*Power(x1,2))*(lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,5)*Pi);
	else if ( (k1==2) && (k2==3) )
		return (-8*(lambdasq - 2*Power(x1,2))*x2*(3*lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==2) && (k2==4) )
		return (-8*(lambdasq - 2*Power(x1,2))*(3*Power(lambdasq,2) - \
12*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==2) && (k2==5) )
		return (16*(lambdasq - 2*Power(x1,2))*x2*(15*Power(lambdasq,2) - \
20*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==2) && (k2==6) )
		return (16*(lambdasq - 2*Power(x1,2))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==2) && (k2==7) )
		return (-32*(lambdasq - 2*Power(x1,2))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==3) && (k2==0) )
		return (12*lambdasq*x1 - 8*Power(x1,3))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,4)*Pi);
	else if ( (k1==3) && (k2==1) )
		return (8*x1*(-3*lambdasq + \
2*Power(x1,2))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,5)*Pi);
	else if ( (k1==3) && (k2==2) )
		return (8*x1*(-3*lambdasq + 2*Power(x1,2))*(lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==3) && (k2==3) )
		return (16*x1*(3*lambdasq - 2*Power(x1,2))*x2*(3*lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==3) && (k2==4) )
		return (16*x1*(3*lambdasq - 2*Power(x1,2))*(3*Power(lambdasq,2) - \
12*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==3) && (k2==5) )
		return (32*x1*(-3*lambdasq + \
2*Power(x1,2))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==3) && (k2==6) )
		return (32*x1*(-3*lambdasq + 2*Power(x1,2))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==3) && (k2==7) )
		return (64*x1*(3*lambdasq - \
2*Power(x1,2))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==4) && (k2==0) )
		return (4*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,5)*Pi);
	else if ( (k1==4) && (k2==1) )
		return (-8*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==4) && (k2==2) )
		return (-8*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==4) && (k2==3) )
		return (16*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(3*lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) \
+ Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==4) && (k2==4) )
		return (16*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==4) && (k2==5) )
		return (-32*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==4) && (k2==6) )
		return (-32*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==4) && (k2==7) )
		return (64*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==5) && (k2==0) )
		return (-8*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,6)*Pi);
	else if ( (k1==5) && (k2==1) )
		return (16*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==5) && (k2==2) )
		return (16*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==5) && (k2==3) )
		return (-32*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(3*lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) \
+ Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==5) && (k2==4) )
		return (-32*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==5) && (k2==5) )
		return (64*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==5) && (k2==6) )
		return (64*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==5) && (k2==7) )
		return (-128*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==6) && (k2==0) )
		return (-8*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x1,2) \
+ 60*lambdasq*Power(x1,4) - 8*Power(x1,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,7)*Pi);
	else if ( (k1==6) && (k2==1) )
		return (16*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x1,2) \
+ 60*lambdasq*Power(x1,4) - 8*Power(x1,6))*x2)/(M_EXP(1*(Power(x1,2) \
+ Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==6) && (k2==2) )
		return (16*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x1,2) \
+ 60*lambdasq*Power(x1,4) - 8*Power(x1,6))*(lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==6) && (k2==3) )
		return (-32*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*x2*(3*lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) \
+ Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==6) && (k2==4) )
		return (-32*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==6) && (k2==5) )
		return (64*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x1,2) \
+ 60*lambdasq*Power(x1,4) - 8*Power(x1,6))*x2*(15*Power(lambdasq,2) - \
20*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==6) && (k2==6) )
		return (64*(15*Power(lambdasq,3) - 90*Power(lambdasq,2)*Power(x1,2) \
+ 60*lambdasq*Power(x1,4) - 8*Power(x1,6))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==6) && (k2==7) )
		return (-128*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==7) && (k2==0) )
		return (16*x1*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x1,2) + 84*lambdasq*Power(x1,4) - \
8*Power(x1,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,8)*Pi);
	else if ( (k1==7) && (k2==1) )
		return (32*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==7) && (k2==2) )
		return (32*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==7) && (k2==3) )
		return (64*x1*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x1,2) + 84*lambdasq*Power(x1,4) - \
8*Power(x1,6))*x2*(3*lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) \
+ Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==7) && (k2==4) )
		return (64*x1*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x1,2) + 84*lambdasq*Power(x1,4) - \
8*Power(x1,6))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==7) && (k2==5) )
		return (128*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==7) && (k2==6) )
		return (128*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==7) && (k2==7) )
		return (256*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*x2*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x2,2) - 84*lambdasq*Power(x2,4) + \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else {
#if HVS_DEBUG
		printf("Err:%d,%d\n",k1,k2);
		hvsdie("Error. Function hb1.\n");
#endif
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
	else if ( (k1==0) && (k2==6) )
		return (he0*x2*(Power(x2,6)*Power(Power(x1,2) + Power(x2,2),6) - \
3*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(7*Power(x1,2) + 3*Power(x2,2)) + \
15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(7*Power(x1,4) + Power(x2,4)) + 720*(-1 + \
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
	else if ( (k1==0) && (k2==7) )
		return -(-420*he0*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(\
Power(x1,2) - Power(x2,2),2)*Power(Power(x1,2) + Power(x2,2),4) + \
he0*Power(x2,8)*Power(Power(x1,2) + Power(x2,2),7) - \
14*he0*lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(2*Power(x1,2) + Power(x2,2)) + \
42*he0*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 5040*(-1 + he0)*Power(lambdasq,7)*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
630*he0*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 105*he0*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,8) - 28*Power(x1,6)*Power(x2,2) + \
70*Power(x1,4)*Power(x2,4) - 28*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 2520*he0*Power(lambdasq,6)*(Power(x1,10) - \
27*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
42*Power(x1,4)*Power(x2,6) - 27*Power(x1,2)*Power(x2,8) + \
Power(x2,10)))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
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
	else if ( (k1==1) && (k2==6) )
		return -(he0*x1*x2*(Power(x2,6)*Power(Power(x1,2) + Power(x2,2),7) \
- 7*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,2) + Power(x2,2)) + \
21*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 40320*(-1 + 1/he0)*Power(lambdasq,7)*(Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - Power(x2,6)) \
- 5040*Power(lambdasq,5)*Power(Power(x1,2) + \
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
	else if ( (k1==1) && (k2==7) )
		return (he0*x1*(Power(x2,8)*Power(Power(x1,2) + Power(x2,2),8) - \
4*lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,2) + 3*Power(x2,2)) + \
14*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) + 84*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(-5*Power(x1,6) + 15*Power(x1,4)*Power(x2,2) - \
11*Power(x1,2)*Power(x2,4) + Power(x2,6)) + (40320*(-1 + \
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
9*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
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
	else if ( (k1==2) && (k2==6) )
		return (he0*x2*(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8) - lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 7*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,6) - 9*Power(x1,4)*Power(x2,2) + \
9*Power(x1,2)*Power(x2,4) + Power(x2,6)) + 40320*(-1 + \
1/he0)*Power(lambdasq,8)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) \
+ 126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
66*Power(x1,4)*Power(x2,4) - 16*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==2) && (k2==7) )
		return (he0*(-(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),9)) + lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 6*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) + 15*Power(x1,2)*Power(x2,4) + \
2*Power(x2,6)) + 42*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(10*Power(x1,8) - 35*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 5*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
- 21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,10) - 200*Power(x1,8)*Power(x2,2) + \
950*Power(x1,6)*Power(x2,4) - 940*Power(x1,4)*Power(x2,6) + \
205*Power(x1,2)*Power(x2,8) - 4*Power(x2,10)) + 362880*(-1 + \
1/he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
45360*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
7560*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
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
	else if ( (k1==3) && (k2==6) )
		return (he0*x1*x2*(-(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)) + 3*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
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
5*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
45*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) + \
21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
102*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==3) && (k2==7) )
		return (he0*x1*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(28*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 30*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,6) + 6*Power(x1,2)*Power(x2,4) + \
Power(x2,6)) - 30*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(14*Power(x1,8) - 49*Power(x1,6)*Power(x2,2) + \
105*Power(x1,4)*Power(x2,4) - 19*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) + (3628800*(-1 + he0)*Power(lambdasq,10)*(Power(x1,10) \
- 55*Power(x1,8)*Power(x2,2) + 330*Power(x1,6)*Power(x2,4) - \
462*Power(x1,4)*Power(x2,6) + 165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10)))/he0 + 1814400*Power(lambdasq,9)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,10) - 40*Power(x1,8)*Power(x2,2) + \
250*Power(x1,6)*Power(x2,4) - 344*Power(x1,4)*Power(x2,6) + \
125*Power(x1,2)*Power(x2,8) - \
8*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
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
	else if ( (k1==4) && (k2==6) )
		return (he0*x2*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),9)*(21*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(35*Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
78*Power(x1,4)*Power(x2,4) + 8*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
- 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(8*Power(x1,10) - 125*Power(x1,8)*Power(x2,2) + \
344*Power(x1,6)*Power(x2,4) - 250*Power(x1,4)*Power(x2,6) + \
40*Power(x1,2)*Power(x2,8) - Power(x2,10)) + 3628800*(-1 + \
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
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,10) - 49*Power(x1,8)*Power(x2,2) + \
196*Power(x1,6)*Power(x2,4) - 104*Power(x1,4)*Power(x2,6) + \
29*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==4) && (k2==7) )
		return -(he0*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11) - \
2*he0*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(14*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + he0*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),9)*(210*Power(x1,8) + 28*Power(x1,6)*Power(x2,2) + \
309*Power(x1,4)*Power(x2,4) + 54*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) - \
30*he0*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(14*Power(x1,10) - 42*Power(x1,8)*Power(x2,2) + \
154*Power(x1,6)*Power(x2,4) - 37*Power(x1,4)*Power(x2,6) + \
18*Power(x1,2)*Power(x2,8) + Power(x2,10)) + 39916800*(-1 + \
he0)*Power(lambdasq,11)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 831600*he0*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 103950*he0*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*he0*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 210*he0*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(4*Power(x1,12) - 273*Power(x1,10)*Power(x2,2) + \
2040*Power(x1,8)*Power(x2,4) - 3814*Power(x1,6)*Power(x2,6) + \
2040*Power(x1,4)*Power(x2,8) - 273*Power(x1,2)*Power(x2,10) + \
4*Power(x2,12)) + 15*he0*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,12) - 252*Power(x1,10)*Power(x2,2) + \
2100*Power(x1,8)*Power(x2,4) - 3724*Power(x1,6)*Power(x2,6) + \
2115*Power(x1,4)*Power(x2,8) - 240*Power(x1,2)*Power(x2,10) + \
10*Power(x2,12)) + 19958400*he0*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*he0*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16)))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
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
	else if ( (k1==5) && (k2==6) )
		return -(he0*x1*x2*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(21*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) + 5*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(21*Power(x1,8) + 72*Power(x1,4)*Power(x2,4) + \
8*Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) - \
45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(14*Power(x1,10) - 315*Power(x1,8)*Power(x2,2) + \
1064*Power(x1,6)*Power(x2,4) - 1114*Power(x1,4)*Power(x2,6) + \
290*Power(x1,2)*Power(x2,8) - 19*Power(x2,10)) + 159667200*(-1 + \
1/he0)*Power(lambdasq,11)*(3*Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 198*Power(x1,6)*Power(x2,4) - \
198*Power(x1,4)*Power(x2,6) + 55*Power(x1,2)*Power(x2,8) - \
3*Power(x2,10)) - 19958400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
3326400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
415800*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
41580*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,10) - 35*Power(x1,8)*Power(x2,2) + \
252*Power(x1,6)*Power(x2,4) - 152*Power(x1,4)*Power(x2,6) + \
85*Power(x1,2)*Power(x2,8) + 3*Power(x2,10)) - \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==5) && (k2==7) )
		return (he0*x1*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),12) - \
2*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(14*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(70*Power(x1,8) + 28*Power(x1,6)*Power(x2,2) + \
159*Power(x1,4)*Power(x2,4) + 30*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 60*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,10) - 14*Power(x1,8)*Power(x2,2) + \
105*Power(x1,6)*Power(x2,4) - 26*Power(x1,4)*Power(x2,6) + \
26*Power(x1,2)*Power(x2,8) + 2*Power(x2,10)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)))/he0 + 239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 196*Power(x1,10)*Power(x2,2) + \
2380*Power(x1,8)*Power(x2,4) - 5012*Power(x1,6)*Power(x2,6) + \
4259*Power(x1,4)*Power(x2,8) - 752*Power(x1,2)*Power(x2,10) + \
66*Power(x2,12)) + 90*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,12) - 651*Power(x1,10)*Power(x2,2) + \
5880*Power(x1,8)*Power(x2,4) - 14182*Power(x1,6)*Power(x2,6) + \
10599*Power(x1,4)*Power(x2,8) - 2367*Power(x1,2)*Power(x2,10) + \
106*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==6) && (k2==0) )
		return (he0*x2*(Power(x1,6)*Power(Power(x1,2) + Power(x2,2),6) - \
3*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,2) + 5*Power(x2,2)) + \
15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
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
	else if ( (k1==6) && (k2==1) )
		return (he0*(-(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)) + lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) + 17*Power(x1,4)*Power(x2,2) - \
25*Power(x1,2)*Power(x2,4) + 15*Power(x2,6)) + (5040*(-1 + \
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
	else if ( (k1==6) && (k2==2) )
		return (he0*x2*(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8) - lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,6) + 107*Power(x1,4)*Power(x2,2) - \
75*Power(x1,2)*Power(x2,4) + 45*Power(x2,6)) + 40320*(-1 + \
1/he0)*Power(lambdasq,8)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) \
+ 126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 5040*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),2)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 840*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),3)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),4)*(9*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
126*Power(x1,4)*Power(x2,4) - 36*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(33*Power(x1,8) - 288*Power(x1,6)*Power(x2,2) + \
450*Power(x1,4)*Power(x2,4) - 120*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==6) && (k2==3) )
		return -(he0*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9) - \
3*he0*lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + \
3*he0*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,8) + 60*Power(x1,4)*Power(x2,4) - \
20*Power(x1,2)*Power(x2,6) + 15*Power(x2,8)) + \
9*he0*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(11*Power(x1,10) - 470*Power(x1,8)*Power(x2,2) + \
2210*Power(x1,6)*Power(x2,4) - 2200*Power(x1,4)*Power(x2,6) + \
475*Power(x1,2)*Power(x2,8) - 10*Power(x2,10)) + 362880*(-1 + \
he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
45360*he0*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
7560*he0*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
945*he0*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
3*he0*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,10) + 135*Power(x1,8)*Power(x2,2) - \
520*Power(x1,6)*Power(x2,4) + 580*Power(x1,4)*Power(x2,6) - \
105*Power(x1,2)*Power(x2,8) + 5*Power(x2,10)) + \
181440*he0*Power(lambdasq,8)*(Power(x1,12) - \
44*Power(x1,10)*Power(x2,2) + 165*Power(x1,8)*Power(x2,4) - \
165*Power(x1,4)*Power(x2,8) + 44*Power(x1,2)*Power(x2,10) - \
Power(x2,12)))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==6) && (k2==4) )
		return (he0*x2*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10) - \
5*lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(2*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,8) + 18*Power(x1,4)*Power(x2,4) - \
2*Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) + \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(59*Power(x1,10) - 860*Power(x1,8)*Power(x2,2) + \
2438*Power(x1,6)*Power(x2,4) - 1720*Power(x1,4)*Power(x2,6) + \
295*Power(x1,2)*Power(x2,8) - 4*Power(x2,10)) + \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
- 75*Power(x1,8)*Power(x2,2) + 152*Power(x1,6)*Power(x2,4) - \
140*Power(x1,4)*Power(x2,6) + 15*Power(x1,2)*Power(x2,8) - \
Power(x2,10)) + (3628800*(-1 + \
he0)*Power(lambdasq,10)*(11*Power(x1,10) - \
165*Power(x1,8)*Power(x2,2) + 462*Power(x1,6)*Power(x2,4) - \
330*Power(x1,4)*Power(x2,6) + 55*Power(x1,2)*Power(x2,8) - \
Power(x2,10)))/he0 + 1814400*Power(lambdasq,9)*(Power(x1,2) + \
Power(x2,2))*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(11*Power(x1,10) - 165*Power(x1,8)*Power(x2,2) + \
462*Power(x1,6)*Power(x2,4) - 330*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==6) && (k2==5) )
		return (he0*(-(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)) + lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(15*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - \
5*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,8) + 3*Power(x1,6)*Power(x2,2) + \
76*Power(x1,4)*Power(x2,4) + 3*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,12) - 291*Power(x1,10)*Power(x2,2) + \
1995*Power(x1,8)*Power(x2,4) - 3874*Power(x1,6)*Power(x2,6) + \
1995*Power(x1,4)*Power(x2,8) - 291*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + (39916800*(-1 + he0)*Power(lambdasq,11)*(Power(x1,12) \
- 66*Power(x1,10)*Power(x2,2) + 495*Power(x1,8)*Power(x2,4) - \
924*Power(x1,6)*Power(x2,6) + 495*Power(x1,4)*Power(x2,8) - \
66*Power(x1,2)*Power(x2,10) + Power(x2,12)))/he0 + \
831600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 103950*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,12) - 6*Power(x1,10)*Power(x2,2) + \
165*Power(x1,8)*Power(x2,4) - 184*Power(x1,6)*Power(x2,6) + \
165*Power(x1,4)*Power(x2,8) - 6*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 15*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(59*Power(x1,12) - 3804*Power(x1,10)*Power(x2,2) + \
28605*Power(x1,8)*Power(x2,4) - 53336*Power(x1,6)*Power(x2,6) + \
28605*Power(x1,4)*Power(x2,8) - 3804*Power(x1,2)*Power(x2,10) + \
59*Power(x2,12)) + 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==6) && (k2==6) )
		return -(he0*x2*(-(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)) + \
3*lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(7*Power(x1,4) + 4*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(35*Power(x1,8) + 21*Power(x1,6)*Power(x2,2) + \
172*Power(x1,4)*Power(x2,4) + 25*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 329*Power(x1,10)*Power(x2,2) + \
1253*Power(x1,8)*Power(x2,4) - 1854*Power(x1,6)*Power(x2,6) + \
685*Power(x1,4)*Power(x2,8) - 97*Power(x1,2)*Power(x2,10) - \
Power(x2,12)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(13*Power(x1,12) - \
286*Power(x1,10)*Power(x2,2) + 1287*Power(x1,8)*Power(x2,4) - \
1716*Power(x1,6)*Power(x2,6) + 715*Power(x1,4)*Power(x2,8) - \
78*Power(x1,2)*Power(x2,10) + Power(x2,12)))/he0 + \
239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,12) - 14*Power(x1,10)*Power(x2,2) + \
315*Power(x1,8)*Power(x2,4) - 176*Power(x1,6)*Power(x2,6) + \
197*Power(x1,4)*Power(x2,8) + 6*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(217*Power(x1,12) - 4704*Power(x1,10)*Power(x2,2) + \
21273*Power(x1,8)*Power(x2,4) - 28264*Power(x1,6)*Power(x2,6) + \
11835*Power(x1,4)*Power(x2,8) - 1272*Power(x1,2)*Power(x2,10) + \
19*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==6) && (k2==7) )
		return (he0*(-(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),13)) + lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),12)*(28*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(70*Power(x1,8) + 56*Power(x1,6)*Power(x2,2) + \
229*Power(x1,4)*Power(x2,4) + 50*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 3*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(140*Power(x1,12) - 70*Power(x1,10)*Power(x2,2) + \
2758*Power(x1,8)*Power(x2,4) - 483*Power(x1,6)*Power(x2,6) + \
1225*Power(x1,4)*Power(x2,8) + 105*Power(x1,2)*Power(x2,10) + \
5*Power(x2,12)) - 45*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(217*Power(x1,14) - 19502*Power(x1,12)*Power(x2,2) + \
214767*Power(x1,10)*Power(x2,4) - 644056*Power(x1,8)*Power(x2,6) + \
644231*Power(x1,6)*Power(x2,8) - 214662*Power(x1,4)*Power(x2,10) + \
19537*Power(x1,2)*Power(x2,12) - 212*Power(x2,14)) - \
45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,14) - 1267*Power(x1,12)*Power(x2,2) + \
13272*Power(x1,10)*Power(x2,4) - 40516*Power(x1,8)*Power(x2,6) + \
39991*Power(x1,6)*Power(x2,8) - 13587*Power(x1,4)*Power(x2,10) + \
1162*Power(x1,2)*Power(x2,12) - 22*Power(x2,14)) + 6227020800*(-1 + \
1/he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) \
+ 1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,14) - 112*Power(x1,12)*Power(x2,2) + \
2632*Power(x1,10)*Power(x2,4) - 5936*Power(x1,8)*Power(x2,6) + \
7511*Power(x1,6)*Power(x2,8) - 1792*Power(x1,4)*Power(x2,10) + \
322*Power(x1,2)*Power(x2,12) + 8*Power(x2,14)) - \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==7) && (k2==0) )
		return -(he0*x1*x2*(Power(x1,6)*Power(Power(x1,2) + Power(x2,2),7) \
- 7*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,2) + 3*Power(x2,2)) + \
21*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + (40320*(-1 + he0)*Power(lambdasq,7)*(Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 5040*Power(lambdasq,5)*Power(Power(x1,2) + \
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
	else if ( (k1==7) && (k2==1) )
		return (he0*x1*(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8) - lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) - 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) + 7*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) + 9*Power(x1,4)*Power(x2,2) - \
9*Power(x1,2)*Power(x2,4) + 15*Power(x2,6)) - \
21*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,8) \
- 16*Power(x1,6)*Power(x2,2) + 66*Power(x1,4)*Power(x2,4) - \
40*Power(x1,2)*Power(x2,6) + 5*Power(x2,8)) + 40320*(-1 + \
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
9*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==7) && (k2==2) )
		return (he0*x1*x2*(-(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)) + (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) \
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
5*Power(x2,4)) + 3*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,6) + 45*Power(x1,4)*Power(x2,2) - \
21*Power(x1,2)*Power(x2,4) + 35*Power(x2,6)) + \
21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
102*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==7) && (k2==3) )
		return (he0*x1*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(6*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
78*Power(x1,4)*Power(x2,4) - 14*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + (3628800*(-1 + \
he0)*Power(lambdasq,10)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)))/he0 + \
1814400*Power(lambdasq,9)*(Power(x1,2) + Power(x2,2))*(Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 330*Power(x1,6)*Power(x2,4) - \
462*Power(x1,4)*Power(x2,6) + 165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10)) + 453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,10) - 40*Power(x1,8)*Power(x2,2) + \
250*Power(x1,6)*Power(x2,4) - 344*Power(x1,4)*Power(x2,6) + \
125*Power(x1,2)*Power(x2,8) - 8*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
+ 29*Power(x1,8)*Power(x2,2) - 104*Power(x1,6)*Power(x2,4) + \
196*Power(x1,4)*Power(x2,6) - 49*Power(x1,2)*Power(x2,8) + \
7*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==7) && (k2==4) )
		return -(he0*x1*x2*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),10)*(10*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + 5*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
72*Power(x1,4)*Power(x2,4) + 21*Power(x2,8)) + \
45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,10) - 290*Power(x1,8)*Power(x2,2) + \
1114*Power(x1,6)*Power(x2,4) - 1064*Power(x1,4)*Power(x2,6) + \
315*Power(x1,2)*Power(x2,8) - 14*Power(x2,10)) + (159667200*(-1 + \
he0)*Power(lambdasq,11)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) \
+ 198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)))/he0 + \
19958400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
3326400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
415800*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
41580*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,10) + 85*Power(x1,8)*Power(x2,2) - \
152*Power(x1,6)*Power(x2,4) + 252*Power(x1,4)*Power(x2,6) - \
35*Power(x1,2)*Power(x2,8) + 7*Power(x2,10)) + \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==7) && (k2==5) )
		return -(he0*x1*(-(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)) + \
3*lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(5*Power(x1,4) + 4*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(15*Power(x1,8) + 25*Power(x1,6)*Power(x2,2) + \
172*Power(x1,4)*Power(x2,4) + 21*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,12) + 97*Power(x1,10)*Power(x2,2) - \
685*Power(x1,8)*Power(x2,4) + 1854*Power(x1,6)*Power(x2,6) - \
1253*Power(x1,4)*Power(x2,8) + 329*Power(x1,2)*Power(x2,10) - \
7*Power(x2,12)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,12) + 6*Power(x1,10)*Power(x2,2) + \
197*Power(x1,8)*Power(x2,4) - 176*Power(x1,6)*Power(x2,6) + \
315*Power(x1,4)*Power(x2,8) - 14*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)))/he0 + 239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,12) - 1272*Power(x1,10)*Power(x2,2) + \
11835*Power(x1,8)*Power(x2,4) - 28264*Power(x1,6)*Power(x2,6) + \
21273*Power(x1,4)*Power(x2,8) - 4704*Power(x1,2)*Power(x2,10) + \
217*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==7) && (k2==6) )
		return -(he0*x1*x2*(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(21*Power(x1,4) + 16*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) - (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 - 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(35*Power(x1,8) + 49*Power(x1,6)*Power(x2,2) + \
236*Power(x1,4)*Power(x2,4) + 49*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,12) + 1043*Power(x1,10)*Power(x2,2) - \
3899*Power(x1,8)*Power(x2,4) + 8434*Power(x1,6)*Power(x2,6) - \
3899*Power(x1,4)*Power(x2,8) + 1043*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) - 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),10)*(35*Power(x1,12) + 70*Power(x1,10)*Power(x2,2) + \
1967*Power(x1,8)*Power(x2,4) - 712*Power(x1,6)*Power(x2,6) + \
1967*Power(x1,4)*Power(x2,8) + 70*Power(x1,2)*Power(x2,10) + \
35*Power(x2,12)) - 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(203*Power(x1,12) - 4788*Power(x1,10)*Power(x2,2) + \
27069*Power(x1,8)*Power(x2,4) - 45704*Power(x1,6)*Power(x2,6) + \
27069*Power(x1,4)*Power(x2,8) - 4788*Power(x1,2)*Power(x2,10) + \
203*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==7) && (k2==7) )
		return (he0*x1*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14) - \
7*lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(4*Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
43589145600*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
10897286400*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1816214400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
227026800*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
22702680*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1891890*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
135135*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
7*Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(30*Power(x1,8) + 40*Power(x1,6)*Power(x2,2) + \
135*Power(x1,4)*Power(x2,4) + 36*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) - 21*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),11)*(20*Power(x1,12) + 30*Power(x1,10)*Power(x2,2) + \
510*Power(x1,8)*Power(x2,4) - 5*Power(x1,6)*Power(x2,6) + \
357*Power(x1,4)*Power(x2,8) + 35*Power(x1,2)*Power(x2,10) + \
5*Power(x2,12)) + 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(29*Power(x1,14) - 2800*Power(x1,12)*Power(x2,2) + \
36645*Power(x1,10)*Power(x2,4) - 134120*Power(x1,8)*Power(x2,6) + \
172615*Power(x1,6)*Power(x2,8) - 80472*Power(x1,4)*Power(x2,10) + \
12215*Power(x1,2)*Power(x2,12) - 400*Power(x2,14)) + \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(5*Power(x1,14) + 2100*Power(x1,10)*Power(x2,4) - \
4340*Power(x1,8)*Power(x2,6) + 8655*Power(x1,6)*Power(x2,8) - \
2394*Power(x1,4)*Power(x2,10) + 840*Power(x1,2)*Power(x2,12) + \
30*Power(x2,14)) - 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,14) + 525*Power(x1,12)*Power(x2,2) - \
5880*Power(x1,10)*Power(x2,4) + 22820*Power(x1,8)*Power(x2,6) - \
28215*Power(x1,6)*Power(x2,8) + 13797*Power(x1,4)*Power(x2,10) - \
1890*Power(x1,2)*Power(x2,12) + \
90*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else {
#if HVS_DEBUG
		printf("Err:%d,%d\n",k1,k2);
		hvsdie("Error. Function hb1.\n");
#endif
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
	else if ( (k1==0) && (k2==6) )
		return (he0*x1*(-(Power(x2,6)*Power(Power(x1,2) + Power(x2,2),6)) + \
3*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,2) + Power(x2,2)) - \
15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,4) - 4*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + (720*(-1 + he0)*Power(lambdasq,6)*(Power(x1,6) - \
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
	else if ( (k1==0) && (k2==7) )
		return (he0*x1*x2*(Power(x2,6)*Power(Power(x1,2) + Power(x2,2),7) - \
7*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,2) + Power(x2,2)) + \
21*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 40320*(-1 + 1/he0)*Power(lambdasq,7)*(Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - Power(x2,6)) \
- 5040*Power(lambdasq,5)*Power(Power(x1,2) + \
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
	else if ( (k1==1) && (k2==6) )
		return (he0*(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),7) - lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(15*Power(x1,6) - 25*Power(x1,4)*Power(x2,2) + \
17*Power(x1,2)*Power(x2,4) + Power(x2,6)) + 5040*(-1 + \
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
	else if ( (k1==1) && (k2==7) )
		return (he0*x2*(-(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)) + lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 7*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,6) - 9*Power(x1,4)*Power(x2,2) + \
9*Power(x1,2)*Power(x2,4) + Power(x2,6)) + (40320*(-1 + \
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
Power(x2,8)) + 21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
66*Power(x1,4)*Power(x2,4) - 16*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
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
	else if ( (k1==2) && (k2==6) )
		return (he0*x1*(-(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)) + lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(15*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) - Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(45*Power(x1,6) - 75*Power(x1,4)*Power(x2,2) + \
107*Power(x1,2)*Power(x2,4) + 3*Power(x2,6)) + (40320*(-1 + \
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
9*Power(x2,8)) + 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,8) - 120*Power(x1,6)*Power(x2,2) + \
450*Power(x1,4)*Power(x2,4) - 288*Power(x1,2)*Power(x2,6) + \
33*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==2) && (k2==7) )
		return (he0*x1*x2*(Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9) - 3*lambdasq*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)))/he0 - \
362880*Power(lambdasq,8)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) - \
90720*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 15120*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 1890*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
45*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) - \
21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
102*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
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
	else if ( (k1==3) && (k2==6) )
		return (he0*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9) - 3*lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(5*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(15*Power(x1,8) - 20*Power(x1,6)*Power(x2,2) + \
60*Power(x1,4)*Power(x2,4) + Power(x2,8)) - \
9*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(10*Power(x1,10) - 475*Power(x1,8)*Power(x2,2) + \
2200*Power(x1,6)*Power(x2,4) - 2210*Power(x1,4)*Power(x2,6) + \
470*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + 362880*(-1 + \
1/he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
45360*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
7560*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,10) - 105*Power(x1,8)*Power(x2,2) + \
580*Power(x1,6)*Power(x2,4) - 520*Power(x1,4)*Power(x2,6) + \
135*Power(x1,2)*Power(x2,8) + Power(x2,10)) - \
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==3) && (k2==7) )
		return -(he0*x2*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),9)*(21*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(35*Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
78*Power(x1,4)*Power(x2,4) + 8*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
- 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(8*Power(x1,10) - 125*Power(x1,8)*Power(x2,2) + \
344*Power(x1,6)*Power(x2,4) - 250*Power(x1,4)*Power(x2,6) + \
40*Power(x1,2)*Power(x2,8) - Power(x2,10)) + 3628800*(-1 + \
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
55*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,10) - 49*Power(x1,8)*Power(x2,2) + \
196*Power(x1,6)*Power(x2,4) - 104*Power(x1,4)*Power(x2,6) + \
29*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
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
	else if ( (k1==4) && (k2==6) )
		return (he0*x1*(-(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)) + \
5*lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,4) + Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) - 15*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,8) - 2*Power(x1,6)*Power(x2,2) + \
18*Power(x1,4)*Power(x2,4) + Power(x2,8)) + \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(4*Power(x1,10) - 295*Power(x1,8)*Power(x2,2) + \
1720*Power(x1,6)*Power(x2,4) - 2438*Power(x1,4)*Power(x2,6) + \
860*Power(x1,2)*Power(x2,8) - 59*Power(x2,10)) + (3628800*(-1 + \
he0)*Power(lambdasq,10)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)))/he0 + \
1814400*Power(lambdasq,9)*(Power(x1,2) + Power(x2,2))*(Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 330*Power(x1,6)*Power(x2,4) - \
462*Power(x1,4)*Power(x2,6) + 165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10)) + 453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
- 15*Power(x1,8)*Power(x2,2) + 140*Power(x1,6)*Power(x2,4) - \
152*Power(x1,4)*Power(x2,6) + 75*Power(x1,2)*Power(x2,8) - \
Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==4) && (k2==7) )
		return (he0*x1*x2*(Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(21*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) + 5*Power(lambdasq,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(21*Power(x1,8) + 72*Power(x1,4)*Power(x2,4) + \
8*Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) - \
45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(14*Power(x1,10) - 315*Power(x1,8)*Power(x2,2) + \
1064*Power(x1,6)*Power(x2,4) - 1114*Power(x1,4)*Power(x2,6) + \
290*Power(x1,2)*Power(x2,8) - 19*Power(x2,10)) + 159667200*(-1 + \
1/he0)*Power(lambdasq,11)*(3*Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 198*Power(x1,6)*Power(x2,4) - \
198*Power(x1,4)*Power(x2,6) + 55*Power(x1,2)*Power(x2,8) - \
3*Power(x2,10)) - 19958400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
3326400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
415800*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
41580*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,10) - 35*Power(x1,8)*Power(x2,2) + \
252*Power(x1,6)*Power(x2,4) - 152*Power(x1,4)*Power(x2,6) + \
85*Power(x1,2)*Power(x2,8) + 3*Power(x2,10)) - \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
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
	else if ( (k1==5) && (k2==6) )
		return (he0*(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(15*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + \
5*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,8) + 3*Power(x1,6)*Power(x2,2) + \
76*Power(x1,4)*Power(x2,4) + 3*Power(x1,2)*Power(x2,6) + \
9*Power(x2,8)) - 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,12) - 291*Power(x1,10)*Power(x2,2) + \
1995*Power(x1,8)*Power(x2,4) - 3874*Power(x1,6)*Power(x2,6) + \
1995*Power(x1,4)*Power(x2,8) - 291*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 39916800*(-1 + \
1/he0)*Power(lambdasq,11)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) \
+ 495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 831600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 103950*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,12) - 6*Power(x1,10)*Power(x2,2) + \
165*Power(x1,8)*Power(x2,4) - 184*Power(x1,6)*Power(x2,6) + \
165*Power(x1,4)*Power(x2,8) - 6*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 15*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(59*Power(x1,12) - 3804*Power(x1,10)*Power(x2,2) + \
28605*Power(x1,8)*Power(x2,4) - 53336*Power(x1,6)*Power(x2,6) + \
28605*Power(x1,4)*Power(x2,8) - 3804*Power(x1,2)*Power(x2,10) + \
59*Power(x2,12)) - 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) - 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==5) && (k2==7) )
		return (he0*x2*(-(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)) + \
3*lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(7*Power(x1,4) + 4*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(35*Power(x1,8) + 21*Power(x1,6)*Power(x2,2) + \
172*Power(x1,4)*Power(x2,4) + 25*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 329*Power(x1,10)*Power(x2,2) + \
1253*Power(x1,8)*Power(x2,4) - 1854*Power(x1,6)*Power(x2,6) + \
685*Power(x1,4)*Power(x2,8) - 97*Power(x1,2)*Power(x2,10) - \
Power(x2,12)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(13*Power(x1,12) - \
286*Power(x1,10)*Power(x2,2) + 1287*Power(x1,8)*Power(x2,4) - \
1716*Power(x1,6)*Power(x2,6) + 715*Power(x1,4)*Power(x2,8) - \
78*Power(x1,2)*Power(x2,10) + Power(x2,12)))/he0 + \
239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,12) - 14*Power(x1,10)*Power(x2,2) + \
315*Power(x1,8)*Power(x2,4) - 176*Power(x1,6)*Power(x2,6) + \
197*Power(x1,4)*Power(x2,8) + 6*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(217*Power(x1,12) - 4704*Power(x1,10)*Power(x2,2) + \
21273*Power(x1,8)*Power(x2,4) - 28264*Power(x1,6)*Power(x2,6) + \
11835*Power(x1,4)*Power(x2,8) - 1272*Power(x1,2)*Power(x2,10) + \
19*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==6) && (k2==0) )
		return -(he0*x1*(Power(x1,6)*Power(Power(x1,2) + Power(x2,2),6) - \
3*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(3*Power(x1,2) + 7*Power(x2,2)) + \
15*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,4) + 7*Power(x2,4)) + (720*(-1 + \
he0)*Power(lambdasq,6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)))/he0 + \
360*Power(lambdasq,5)*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6)) + 90*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) + \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) \
- 21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))))/(2.*Power(lambdasq,6)*Pi*Power(Power(x1,2) + \
Power(x2,2),7));
	else if ( (k1==6) && (k2==1) )
		return (he0*x1*x2*(Power(x1,6)*Power(Power(x1,2) + Power(x2,2),7) - \
7*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,2) + 3*Power(x2,2)) + \
21*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) - 2*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + (40320*(-1 + he0)*Power(lambdasq,7)*(Power(x1,6) - \
7*Power(x1,4)*Power(x2,2) + 7*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 5040*Power(lambdasq,5)*Power(Power(x1,2) + \
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
	else if ( (k1==6) && (k2==2) )
		return (he0*x1*(-(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)) + lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) - 7*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) + 9*Power(x1,4)*Power(x2,2) - \
9*Power(x1,2)*Power(x2,4) + 15*Power(x2,6)) + \
21*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,8) \
- 16*Power(x1,6)*Power(x2,2) + 66*Power(x1,4)*Power(x2,4) - \
40*Power(x1,2)*Power(x2,6) + 5*Power(x2,8)) + (40320*(-1 + \
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
9*Power(x2,8))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==6) && (k2==3) )
		return (he0*x1*x2*(Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9) - (725760*(-1 + he0)*Power(lambdasq,9)*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)))/he0 - \
362880*Power(lambdasq,8)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4)) - \
90720*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 15120*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 1890*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 3*lambdasq*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,6) + 45*Power(x1,4)*Power(x2,2) - \
21*Power(x1,2)*Power(x2,4) + 35*Power(x2,6)) - \
21*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
102*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==6) && (k2==4) )
		return -(he0*x1*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(6*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
78*Power(x1,4)*Power(x2,4) - 14*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + (3628800*(-1 + \
he0)*Power(lambdasq,10)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)))/he0 + \
1814400*Power(lambdasq,9)*(Power(x1,2) + Power(x2,2))*(Power(x1,10) - \
55*Power(x1,8)*Power(x2,2) + 330*Power(x1,6)*Power(x2,4) - \
462*Power(x1,4)*Power(x2,6) + 165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10)) + 453600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
75600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
9450*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
330*Power(x1,6)*Power(x2,4) - 462*Power(x1,4)*Power(x2,6) + \
165*Power(x1,2)*Power(x2,8) - 11*Power(x2,10)) + \
105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,10) - 40*Power(x1,8)*Power(x2,2) + \
250*Power(x1,6)*Power(x2,4) - 344*Power(x1,4)*Power(x2,6) + \
125*Power(x1,2)*Power(x2,8) - 8*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,10) \
+ 29*Power(x1,8)*Power(x2,2) - 104*Power(x1,6)*Power(x2,4) + \
196*Power(x1,4)*Power(x2,6) - 49*Power(x1,2)*Power(x2,8) + \
7*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==6) && (k2==5) )
		return (he0*x1*x2*(Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),10)*(10*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + 5*Power(lambdasq,2)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
72*Power(x1,4)*Power(x2,4) + 21*Power(x2,8)) + \
45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,10) - 290*Power(x1,8)*Power(x2,2) + \
1114*Power(x1,6)*Power(x2,4) - 1064*Power(x1,4)*Power(x2,6) + \
315*Power(x1,2)*Power(x2,8) - 14*Power(x2,10)) + (159667200*(-1 + \
he0)*Power(lambdasq,11)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) \
+ 198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)))/he0 + \
19958400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),2)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
3326400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
415800*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
41580*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,10) - 55*Power(x1,8)*Power(x2,2) + \
198*Power(x1,6)*Power(x2,4) - 198*Power(x1,4)*Power(x2,6) + \
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) - \
15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,10) + 85*Power(x1,8)*Power(x2,2) - \
152*Power(x1,6)*Power(x2,4) + 252*Power(x1,4)*Power(x2,6) - \
35*Power(x1,2)*Power(x2,8) + 7*Power(x2,10)) + \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==6) && (k2==6) )
		return (he0*x1*(-(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)) + \
3*lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(5*Power(x1,4) + 4*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(15*Power(x1,8) + 25*Power(x1,6)*Power(x2,2) + \
172*Power(x1,4)*Power(x2,4) + 21*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,12) + 97*Power(x1,10)*Power(x2,2) - \
685*Power(x1,8)*Power(x2,4) + 1854*Power(x1,6)*Power(x2,6) - \
1253*Power(x1,4)*Power(x2,8) + 329*Power(x1,2)*Power(x2,10) - \
7*Power(x2,12)) + 15*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,12) + 6*Power(x1,10)*Power(x2,2) + \
197*Power(x1,8)*Power(x2,4) - 176*Power(x1,6)*Power(x2,6) + \
315*Power(x1,4)*Power(x2,8) - 14*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)))/he0 + 239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,12) - 1272*Power(x1,10)*Power(x2,2) + \
11835*Power(x1,8)*Power(x2,4) - 28264*Power(x1,6)*Power(x2,6) + \
21273*Power(x1,4)*Power(x2,8) - 4704*Power(x1,2)*Power(x2,10) + \
217*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==6) && (k2==7) )
		return (he0*x1*x2*(Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(21*Power(x1,4) + 16*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) - (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 - 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(35*Power(x1,8) + 49*Power(x1,6)*Power(x2,2) + \
236*Power(x1,4)*Power(x2,4) + 49*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,12) + 1043*Power(x1,10)*Power(x2,2) - \
3899*Power(x1,8)*Power(x2,4) + 8434*Power(x1,6)*Power(x2,6) - \
3899*Power(x1,4)*Power(x2,8) + 1043*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) - 3*Power(lambdasq,3)*Power(Power(x1,2) + \
Power(x2,2),10)*(35*Power(x1,12) + 70*Power(x1,10)*Power(x2,2) + \
1967*Power(x1,8)*Power(x2,4) - 712*Power(x1,6)*Power(x2,6) + \
1967*Power(x1,4)*Power(x2,8) + 70*Power(x1,2)*Power(x2,10) + \
35*Power(x2,12)) - 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(203*Power(x1,12) - 4788*Power(x1,10)*Power(x2,2) + \
27069*Power(x1,8)*Power(x2,4) - 45704*Power(x1,6)*Power(x2,6) + \
27069*Power(x1,4)*Power(x2,8) - 4788*Power(x1,2)*Power(x2,10) + \
203*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==7) && (k2==0) )
		return (he0*(-420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(\
Power(x1,2) - Power(x2,2),2)*Power(Power(x1,2) + Power(x2,2),4) + \
Power(x1,8)*Power(Power(x1,2) + Power(x2,2),7) - \
14*lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,2) + 2*Power(x2,2)) + \
42*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + (5040*(-1 + he0)*Power(lambdasq,7)*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
630*Power(lambdasq,5)*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,8) \
- 28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
105*Power(lambdasq,4)*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,8) \
- 28*Power(x1,6)*Power(x2,2) + 70*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
2520*Power(lambdasq,6)*(Power(x1,10) - 27*Power(x1,8)*Power(x2,2) + \
42*Power(x1,6)*Power(x2,4) + 42*Power(x1,4)*Power(x2,6) - \
27*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,7)*Pi*Power(Power(x1,2) + \
Power(x2,2),8));
	else if ( (k1==7) && (k2==1) )
		return -(he0*x2*(Power(x1,8)*Power(Power(x1,2) + Power(x2,2),8) - \
4*lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,2) + 7*Power(x2,2)) + \
14*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + 84*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),5)*(Power(x1,6) - 11*Power(x1,4)*Power(x2,2) + \
15*Power(x1,2)*Power(x2,4) - 5*Power(x2,6)) + (40320*(-1 + \
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
Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==7) && (k2==2) )
		return (he0*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9) - lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) + 6*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(2*Power(x1,6) + 15*Power(x1,4)*Power(x2,2) + \
35*Power(x2,6)) - 42*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),6)*(Power(x1,8) - 5*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 35*Power(x1,2)*Power(x2,6) + \
10*Power(x2,8)) - 21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(4*Power(x1,10) - 205*Power(x1,8)*Power(x2,2) + \
940*Power(x1,6)*Power(x2,4) - 950*Power(x1,4)*Power(x2,6) + \
200*Power(x1,2)*Power(x2,8) - 5*Power(x2,10)) + 362880*(-1 + \
1/he0)*Power(lambdasq,9)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
45360*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
7560*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,10) - 45*Power(x1,8)*Power(x2,2) + \
210*Power(x1,6)*Power(x2,4) - 210*Power(x1,4)*Power(x2,6) + \
45*Power(x1,2)*Power(x2,8) - Power(x2,10)) - \
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==7) && (k2==3) )
		return -(he0*x2*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) + 30*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(Power(x1,6) + 6*Power(x1,4)*Power(x2,2) + \
7*Power(x2,6)) - 30*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,8) - 19*Power(x1,6)*Power(x2,2) + \
105*Power(x1,4)*Power(x2,4) - 49*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(8*Power(x1,10) - 125*Power(x1,8)*Power(x2,2) + \
344*Power(x1,6)*Power(x2,4) - 250*Power(x1,4)*Power(x2,6) + \
40*Power(x1,2)*Power(x2,8) - Power(x2,10)) + 3628800*(-1 + \
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
	else if ( (k1==7) && (k2==4) )
		return (he0*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11) - \
2*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
14*Power(x2,4)) + Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,8) + 54*Power(x1,6)*Power(x2,2) + \
309*Power(x1,4)*Power(x2,4) + 28*Power(x1,2)*Power(x2,6) + \
210*Power(x2,8)) - 30*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(Power(x1,10) + 18*Power(x1,8)*Power(x2,2) - \
37*Power(x1,6)*Power(x2,4) + 154*Power(x1,4)*Power(x2,6) - \
42*Power(x1,2)*Power(x2,8) + 14*Power(x2,10)) + (39916800*(-1 + \
he0)*Power(lambdasq,11)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)))/he0 + 831600*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 103950*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 66*Power(x1,10)*Power(x2,2) + \
495*Power(x1,8)*Power(x2,4) - 924*Power(x1,6)*Power(x2,6) + \
495*Power(x1,4)*Power(x2,8) - 66*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 210*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(4*Power(x1,12) - 273*Power(x1,10)*Power(x2,2) + \
2040*Power(x1,8)*Power(x2,4) - 3814*Power(x1,6)*Power(x2,6) + \
2040*Power(x1,4)*Power(x2,8) - 273*Power(x1,2)*Power(x2,10) + \
4*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(10*Power(x1,12) - 240*Power(x1,10)*Power(x2,2) + \
2115*Power(x1,8)*Power(x2,4) - 3724*Power(x1,6)*Power(x2,6) + \
2100*Power(x1,4)*Power(x2,8) - 252*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) + 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==7) && (k2==5) )
		return -(he0*x2*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12) - \
2*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(5*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
14*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(5*Power(x1,8) + 30*Power(x1,6)*Power(x2,2) + \
159*Power(x1,4)*Power(x2,4) + 28*Power(x1,2)*Power(x2,6) + \
70*Power(x2,8)) - 60*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(2*Power(x1,10) + 26*Power(x1,8)*Power(x2,2) - \
26*Power(x1,6)*Power(x2,4) + 105*Power(x1,4)*Power(x2,6) - \
14*Power(x1,2)*Power(x2,8) + 7*Power(x2,10)) + (479001600*(-1 + \
he0)*Power(lambdasq,12)*(13*Power(x1,12) - \
286*Power(x1,10)*Power(x2,2) + 1287*Power(x1,8)*Power(x2,4) - \
1716*Power(x1,6)*Power(x2,6) + 715*Power(x1,4)*Power(x2,8) - \
78*Power(x1,2)*Power(x2,10) + Power(x2,12)))/he0 + \
239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 90*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(106*Power(x1,12) - 2367*Power(x1,10)*Power(x2,2) + \
10599*Power(x1,8)*Power(x2,4) - 14182*Power(x1,6)*Power(x2,6) + \
5880*Power(x1,4)*Power(x2,8) - 651*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(66*Power(x1,12) - 752*Power(x1,10)*Power(x2,2) + \
4259*Power(x1,8)*Power(x2,4) - 5012*Power(x1,6)*Power(x2,6) + \
2380*Power(x1,4)*Power(x2,8) - 196*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==7) && (k2==6) )
		return (he0*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(15*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) + \
3*Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,8) + 50*Power(x1,6)*Power(x2,2) + \
229*Power(x1,4)*Power(x2,4) + 56*Power(x1,2)*Power(x2,6) + \
70*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(5*Power(x1,12) + 105*Power(x1,10)*Power(x2,2) + \
1225*Power(x1,8)*Power(x2,4) - 483*Power(x1,6)*Power(x2,6) + \
2758*Power(x1,4)*Power(x2,8) - 70*Power(x1,2)*Power(x2,10) + \
140*Power(x2,12)) - 45*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(212*Power(x1,14) - 19537*Power(x1,12)*Power(x2,2) + \
214662*Power(x1,10)*Power(x2,4) - 644231*Power(x1,8)*Power(x2,6) + \
644056*Power(x1,6)*Power(x2,8) - 214767*Power(x1,4)*Power(x2,10) + \
19502*Power(x1,2)*Power(x2,12) - 217*Power(x2,14)) - \
45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(22*Power(x1,14) - 1162*Power(x1,12)*Power(x2,2) + \
13587*Power(x1,10)*Power(x2,4) - 39991*Power(x1,8)*Power(x2,6) + \
40516*Power(x1,6)*Power(x2,8) - 13272*Power(x1,4)*Power(x2,10) + \
1267*Power(x1,2)*Power(x2,12) - 7*Power(x2,14)) + 6227020800*(-1 + \
1/he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) \
+ 1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(8*Power(x1,14) + 322*Power(x1,12)*Power(x2,2) - \
1792*Power(x1,10)*Power(x2,4) + 7511*Power(x1,8)*Power(x2,6) - \
5936*Power(x1,6)*Power(x2,8) + 2632*Power(x1,4)*Power(x2,10) - \
112*Power(x1,2)*Power(x2,12) + 7*Power(x2,14)) - \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==7) && (k2==7) )
		return -(he0*x2*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14) - \
7*lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
4*Power(x2,4)) - (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(3*Power(x1,2) - Power(x2,2))*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 - \
43589145600*Power(lambdasq,13)*(3*Power(x1,2) - \
Power(x2,2))*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
10897286400*Power(lambdasq,12)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1816214400*Power(lambdasq,11)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
227026800*Power(lambdasq,10)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
22702680*Power(lambdasq,9)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1891890*Power(lambdasq,8)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
135135*Power(lambdasq,7)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
7*Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(15*Power(x1,8) + 36*Power(x1,6)*Power(x2,2) + \
135*Power(x1,4)*Power(x2,4) + 40*Power(x1,2)*Power(x2,6) + \
30*Power(x2,8)) - 21*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),11)*(5*Power(x1,12) + 35*Power(x1,10)*Power(x2,2) + \
357*Power(x1,8)*Power(x2,4) - 5*Power(x1,6)*Power(x2,6) + \
510*Power(x1,4)*Power(x2,8) + 30*Power(x1,2)*Power(x2,10) + \
20*Power(x2,12)) - 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(400*Power(x1,14) - 12215*Power(x1,12)*Power(x2,2) + \
80472*Power(x1,10)*Power(x2,4) - 172615*Power(x1,8)*Power(x2,6) + \
134120*Power(x1,6)*Power(x2,8) - 36645*Power(x1,4)*Power(x2,10) + \
2800*Power(x1,2)*Power(x2,12) - 29*Power(x2,14)) - \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(90*Power(x1,14) - 1890*Power(x1,12)*Power(x2,2) + \
13797*Power(x1,10)*Power(x2,4) - 28215*Power(x1,8)*Power(x2,6) + \
22820*Power(x1,6)*Power(x2,8) - 5880*Power(x1,4)*Power(x2,10) + \
525*Power(x1,2)*Power(x2,12) + Power(x2,14)) + \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(30*Power(x1,14) + 840*Power(x1,12)*Power(x2,2) - \
2394*Power(x1,10)*Power(x2,4) + 8655*Power(x1,8)*Power(x2,6) - \
4340*Power(x1,6)*Power(x2,8) + 2100*Power(x1,4)*Power(x2,10) + \
5*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else {
#if HVS_DEBUG
		printf("Err:%d,%d\n",k1,k2);
		hvsdie("Error. Function hb2.\n");
#endif
		return 0.0;
	}
}
