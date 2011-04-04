// Hermite functions up until 10th order

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
	else if ( (k1==0) && (k2==8) )
		return (16*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==0) && (k2==9) )
		return (-32*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==0) && (k2==10) )
		return (-32*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
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
	else if ( (k1==1) && (k2==8) )
		return (-32*x1*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==1) && (k2==9) )
		return (64*x1*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==1) && (k2==10) )
		return (64*x1*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
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
	else if ( (k1==2) && (k2==8) )
		return (-32*(lambdasq - 2*Power(x1,2))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==2) && (k2==9) )
		return (64*(lambdasq - 2*Power(x1,2))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==2) && (k2==10) )
		return (64*(lambdasq - 2*Power(x1,2))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
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
	else if ( (k1==3) && (k2==8) )
		return (64*x1*(3*lambdasq - 2*Power(x1,2))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==3) && (k2==9) )
		return (128*x1*(-3*lambdasq + \
2*Power(x1,2))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==3) && (k2==10) )
		return (128*x1*(-3*lambdasq + 2*Power(x1,2))*(945*Power(lambdasq,5) \
- 9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
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
	else if ( (k1==4) && (k2==8) )
		return (64*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==4) && (k2==9) )
		return (-128*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==4) && (k2==10) )
		return (-128*(3*Power(lambdasq,2) - 12*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
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
	else if ( (k1==5) && (k2==8) )
		return (-128*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==5) && (k2==9) )
		return (256*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else if ( (k1==5) && (k2==10) )
		return (256*x1*(15*Power(lambdasq,2) - 20*lambdasq*Power(x1,2) + \
4*Power(x1,4))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
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
	else if ( (k1==6) && (k2==8) )
		return (-128*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else if ( (k1==6) && (k2==9) )
		return (256*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
	else if ( (k1==6) && (k2==10) )
		return (256*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x1,2) + 60*lambdasq*Power(x1,4) - \
8*Power(x1,6))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,17)*Pi);
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
	else if ( (k1==7) && (k2==8) )
		return (256*x1*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x1,2) + 84*lambdasq*Power(x1,4) - \
8*Power(x1,6))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
	else if ( (k1==7) && (k2==9) )
		return (512*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,17)*Pi);
	else if ( (k1==7) && (k2==10) )
		return (512*x1*(-105*Power(lambdasq,3) + \
210*Power(lambdasq,2)*Power(x1,2) - 84*lambdasq*Power(x1,4) + \
8*Power(x1,6))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,18)*Pi);
	else if ( (k1==8) && (k2==0) )
		return (16*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,9)*Pi);
	else if ( (k1==8) && (k2==1) )
		return (-32*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==8) && (k2==2) )
		return (-32*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*(lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==8) && (k2==3) )
		return (64*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*x2*(3*lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==8) && (k2==4) )
		return (64*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*(3*Power(lambdasq,2) - \
12*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==8) && (k2==5) )
		return (-128*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*x2*(15*Power(lambdasq,2) \
- 20*lambdasq*Power(x2,2) + 4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==8) && (k2==6) )
		return (-128*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else if ( (k1==8) && (k2==7) )
		return (256*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
	else if ( (k1==8) && (k2==8) )
		return (256*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,17)*Pi);
	else if ( (k1==8) && (k2==9) )
		return (-512*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,18)*Pi);
	else if ( (k1==8) && (k2==10) )
		return (-512*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x1,2) + 840*Power(lambdasq,2)*Power(x1,4) \
- 224*lambdasq*Power(x1,6) + 16*Power(x1,8))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,19)*Pi);
	else if ( (k1==9) && (k2==0) )
		return (-32*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,10)*Pi);
	else if ( (k1==9) && (k2==1) )
		return (64*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==9) && (k2==2) )
		return (64*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==9) && (k2==3) )
		return (-128*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(3*lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==9) && (k2==4) )
		return (-128*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==9) && (k2==5) )
		return (256*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else if ( (k1==9) && (k2==6) )
		return (256*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
	else if ( (k1==9) && (k2==7) )
		return (-512*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,17)*Pi);
	else if ( (k1==9) && (k2==8) )
		return (-512*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,18)*Pi);
	else if ( (k1==9) && (k2==9) )
		return (1024*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,19)*Pi);
	else if ( (k1==9) && (k2==10) )
		return (1024*x1*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x1,2) + \
1512*Power(lambdasq,2)*Power(x1,4) - 288*lambdasq*Power(x1,6) + \
16*Power(x1,8))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,20)*Pi);
	else if ( (k1==10) && (k2==0) )
		return (-32*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,11)*Pi);
	else if ( (k1==10) && (k2==1) )
		return (64*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*x2)/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,12)*Pi);
	else if ( (k1==10) && (k2==2) )
		return (64*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*(lambdasq - 2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,13)*Pi);
	else if ( (k1==10) && (k2==3) )
		return (-128*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*x2*(3*lambdasq - \
2*Power(x2,2)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,14)*Pi);
	else if ( (k1==10) && (k2==4) )
		return (-128*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*(3*Power(lambdasq,2) - 12*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,15)*Pi);
	else if ( (k1==10) && (k2==5) )
		return (256*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*x2*(15*Power(lambdasq,2) - 20*lambdasq*Power(x2,2) + \
4*Power(x2,4)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,16)*Pi);
	else if ( (k1==10) && (k2==6) )
		return (256*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*(15*Power(lambdasq,3) - \
90*Power(lambdasq,2)*Power(x2,2) + 60*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,17)*Pi);
	else if ( (k1==10) && (k2==7) )
		return (-512*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*x2*(105*Power(lambdasq,3) - \
210*Power(lambdasq,2)*Power(x2,2) + 84*lambdasq*Power(x2,4) - \
8*Power(x2,6)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,18)*Pi);
	else if ( (k1==10) && (k2==8) )
		return (-512*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*(105*Power(lambdasq,4) - \
840*Power(lambdasq,3)*Power(x2,2) + 840*Power(lambdasq,2)*Power(x2,4) \
- 224*lambdasq*Power(x2,6) + 16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,19)*Pi);
	else if ( (k1==10) && (k2==9) )
		return (1024*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*x2*(945*Power(lambdasq,4) - \
2520*Power(lambdasq,3)*Power(x2,2) + \
1512*Power(lambdasq,2)*Power(x2,4) - 288*lambdasq*Power(x2,6) + \
16*Power(x2,8)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,20)*Pi);
	else if ( (k1==10) && (k2==10) )
		return (1024*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x1,2) + \
12600*Power(lambdasq,3)*Power(x1,4) - \
5040*Power(lambdasq,2)*Power(x1,6) + 720*lambdasq*Power(x1,8) - \
32*Power(x1,10))*(945*Power(lambdasq,5) - \
9450*Power(lambdasq,4)*Power(x2,2) + \
12600*Power(lambdasq,3)*Power(x2,4) - \
5040*Power(lambdasq,2)*Power(x2,6) + 720*lambdasq*Power(x2,8) - \
32*Power(x2,10)))/(M_EXP(1*(Power(x1,2) + \
Power(x2,2))/lambdasq)*Power(lambdasq,21)*Pi);
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
	else if ( (k1==0) && (k2==8) )
		return (he0*x2*(Power(x2,8)*Power(Power(x1,2) + Power(x2,2),8) - \
4*lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(9*Power(x1,2) + 5*Power(x2,2)) + \
14*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(27*Power(x1,4) + 18*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) - 84*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(15*Power(x1,6) - 9*Power(x1,4)*Power(x2,2) + \
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
Power(x2,8)) + 20160*Power(lambdasq,7)*(9*Power(x1,10) - \
75*Power(x1,8)*Power(x2,2) + 42*Power(x1,6)*Power(x2,4) + \
90*Power(x1,4)*Power(x2,6) - 35*Power(x1,2)*Power(x2,8) + \
Power(x2,10))))/(2.*Power(lambdasq,8)*Pi*Power(Power(x1,2) + \
Power(x2,2),9));
	else if ( (k1==0) && (k2==9) )
		return (he0*(-(Power(x2,10)*Power(Power(x1,2) + Power(x2,2),9)) + \
9*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(5*Power(x1,2) + 3*Power(x2,2)) + \
126*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,2) + Power(x2,2))*(5*Power(x1,4) + \
3*Power(x2,4)) - 18*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,4) + 30*Power(x1,2)*Power(x2,2) + \
11*Power(x2,4)) - 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*Power(5*Power(x1,4)*x2 - 10*Power(x1,2)*Power(x2,3) + \
Power(x2,5),2) + (362880*(-1 + he0)*Power(lambdasq,9)*(Power(x1,10) - \
45*Power(x1,8)*Power(x2,2) + 210*Power(x1,6)*Power(x2,4) - \
210*Power(x1,4)*Power(x2,6) + 45*Power(x1,2)*Power(x2,8) - \
Power(x2,10)))/he0 + 45360*Power(lambdasq,7)*Power(Power(x1,2) + \
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
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==0) && (k2==10) )
		return (he0*x2*(Power(x2,10)*Power(Power(x1,2) + Power(x2,2),10) - \
5*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(11*Power(x1,2) + 7*Power(x2,2)) + \
90*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(11*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
4*Power(x2,4)) - 90*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(77*Power(x1,6) + 55*Power(x1,4)*Power(x2,2) + \
55*Power(x1,2)*Power(x2,4) + 13*Power(x2,6)) + \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(55*Power(x1,8) - 88*Power(x1,6)*Power(x2,2) + \
110*Power(x1,4)*Power(x2,4) + 3*Power(x2,8)) + 3628800*(-1 + \
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
	else if ( (k1==1) && (k2==8) )
		return -(he0*x1*x2*(Power(x2,8)*Power(Power(x1,2) + Power(x2,2),9) \
- 18*lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,2) + Power(x2,2)) + (725760*(-1 + \
he0)*Power(lambdasq,9)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 + 362880*Power(lambdasq,8)*(Power(x1,2) + \
Power(x2,2))*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 90720*Power(lambdasq,7)*Power(Power(x1,2) + \
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
5*Power(x2,4)) - \
252*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,4) - 6*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==1) && (k2==9) )
		return (he0*x1*(Power(x2,10)*Power(Power(x1,2) + Power(x2,2),10) - \
5*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,2) + 5*Power(x2,2)) + \
90*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) - 90*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
25*Power(x1,2)*Power(x2,4) + 3*Power(x2,6)) + \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 80*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 3628800*(-1 + \
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
165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==1) && (k2==10) )
		return -(he0*x1*x2*(Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),11) - 11*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(5*Power(x1,2) + 3*Power(x2,2)) + \
110*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - 990*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(7*Power(x1,6) + 3*Power(x1,4)*Power(x2,2) + \
5*Power(x1,2)*Power(x2,4) + Power(x2,6)) + \
495*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,8) - 84*Power(x1,6)*Power(x2,2) + \
114*Power(x1,4)*Power(x2,4) - 20*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 159667200*(-1 + \
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
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
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
	else if ( (k1==2) && (k2==8) )
		return (he0*x2*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(36*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 18*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(21*Power(x1,6) + 8*Power(x1,4)*Power(x2,2) + \
8*Power(x1,2)*Power(x2,4) + Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(14*Power(x1,8) - 21*Power(x1,6)*Power(x2,2) + \
29*Power(x1,4)*Power(x2,4) + Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 80*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + (3628800*(-1 + \
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
	else if ( (k1==2) && (k2==9) )
		return (he0*(-(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),11)) + lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(45*Power(x1,4) + 24*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 5*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(126*Power(x1,6) + 81*Power(x1,4)*Power(x2,2) + \
48*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) + \
90*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(35*Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) + 8*Power(x1,2)*Power(x2,6) + \
2*Power(x2,8)) - 45*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(105*Power(x1,10) - 630*Power(x1,8)*Power(x2,2) + \
1358*Power(x1,6)*Power(x2,4) - 600*Power(x1,4)*Power(x2,6) + \
129*Power(x1,2)*Power(x2,8) + 6*Power(x2,10)) + (39916800*(-1 + \
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
Power(x2,12)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,12) - 180*Power(x1,10)*Power(x2,2) + \
1365*Power(x1,8)*Power(x2,4) - 2536*Power(x1,6)*Power(x2,6) + \
1365*Power(x1,4)*Power(x2,8) - 180*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) + 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==2) && (k2==10) )
		return (he0*x2*(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),12) - 330*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),9)*(Power(x1,4) + Power(x2,4))*(21*Power(x1,4) + \
6*Power(x1,2)*Power(x2,2) + Power(x2,4)) - \
lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(55*Power(x1,4) + 32*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 1485*Power(lambdasq,5)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,2) - 4*x1*x2 + \
Power(x2,2))*(Power(x1,2) + 4*x1*x2 + Power(x2,2))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 33*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(30*Power(x1,6) + 25*Power(x1,4)*Power(x2,2) + \
12*Power(x1,2)*Power(x2,4) + Power(x2,6)) + \
495*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(35*Power(x1,10) - 98*Power(x1,8)*Power(x2,2) + \
194*Power(x1,6)*Power(x2,4) - 40*Power(x1,4)*Power(x2,6) + \
19*Power(x1,2)*Power(x2,8) + 2*Power(x2,10)) + 479001600*(-1 + \
1/he0)*Power(lambdasq,12)*(13*Power(x1,12) - \
286*Power(x1,10)*Power(x2,2) + 1287*Power(x1,8)*Power(x2,4) - \
1716*Power(x1,6)*Power(x2,6) + 715*Power(x1,4)*Power(x2,8) - \
78*Power(x1,2)*Power(x2,10) + Power(x2,12)) - \
239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
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
	else if ( (k1==3) && (k2==8) )
		return -(he0*x1*x2*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(36*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 2*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(189*Power(x1,6) + 72*Power(x1,4)*Power(x2,2) + \
127*Power(x1,2)*Power(x2,4) + 24*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(14*Power(x1,8) - 21*Power(x1,6)*Power(x2,2) + \
51*Power(x1,4)*Power(x2,4) + Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) \
+ 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,10) - 280*Power(x1,8)*Power(x2,2) + \
1134*Power(x1,6)*Power(x2,4) - 1044*Power(x1,4)*Power(x2,6) + \
325*Power(x1,2)*Power(x2,8) - 12*Power(x2,10)) + (159667200*(-1 + \
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
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==3) && (k2==9) )
		return (he0*x1*(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),12) - 3*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(210*Power(x1,6) + 135*Power(x1,4)*Power(x2,2) + \
124*Power(x1,2)*Power(x2,4) + 23*Power(x2,6)) - \
30*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(105*Power(x1,8) - 42*Power(x1,6)*Power(x2,2) + \
234*Power(x1,4)*Power(x2,4) + 46*Power(x1,2)*Power(x2,6) + \
17*Power(x2,8)) + 135*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(35*Power(x1,10) - 210*Power(x1,8)*Power(x2,2) + \
658*Power(x1,6)*Power(x2,4) - 376*Power(x1,4)*Power(x2,6) + \
131*Power(x1,2)*Power(x2,8) + 2*Power(x2,10)) + 479001600*(-1 + \
1/he0)*Power(lambdasq,12)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) \
+ 715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,12) - 78*Power(x1,10)*Power(x2,2) + \
715*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
1287*Power(x1,4)*Power(x2,8) - 286*Power(x1,2)*Power(x2,10) + \
13*Power(x2,12)) - 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,12) - 420*Power(x1,10)*Power(x2,2) + \
3955*Power(x1,8)*Power(x2,4) - 9408*Power(x1,6)*Power(x2,6) + \
7101*Power(x1,4)*Power(x2,8) - 1564*Power(x1,2)*Power(x2,10) + \
73*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==3) && (k2==10) )
		return -(he0*x1*x2*(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),12)*(55*Power(x1,4) + 32*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - (12454041600*(-1 + \
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
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(330*Power(x1,6) + 275*Power(x1,4)*Power(x2,2) + \
184*Power(x1,2)*Power(x2,4) + 31*Power(x2,6)) - \
66*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(105*Power(x1,8) + 30*Power(x1,6)*Power(x2,2) + \
175*Power(x1,4)*Power(x2,4) + 56*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) + 165*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(105*Power(x1,10) - 294*Power(x1,8)*Power(x2,2) + \
894*Power(x1,6)*Power(x2,4) - 224*Power(x1,4)*Power(x2,6) + \
161*Power(x1,2)*Power(x2,8) + 14*Power(x2,10)) - \
1485*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 140*Power(x1,10)*Power(x2,2) + \
833*Power(x1,8)*Power(x2,4) - 1368*Power(x1,6)*Power(x2,6) + \
833*Power(x1,4)*Power(x2,8) - 140*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
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
	else if ( (k1==4) && (k2==8) )
		return (he0*x2*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),12) - \
6*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(6*Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(126*Power(x1,8) + 60*Power(x1,6)*Power(x2,2) + \
139*Power(x1,4)*Power(x2,4) + 30*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 12*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(105*Power(x1,10) - 126*Power(x1,8)*Power(x2,2) + \
591*Power(x1,6)*Power(x2,4) + 4*Power(x1,4)*Power(x2,6) + \
66*Power(x1,2)*Power(x2,8) + 4*Power(x2,10)) + (479001600*(-1 + \
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
Power(x2,12)) + 270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,12) - 791*Power(x1,10)*Power(x2,2) + \
3528*Power(x1,8)*Power(x2,4) - 4734*Power(x1,6)*Power(x2,6) + \
1955*Power(x1,4)*Power(x2,8) - 219*Power(x1,2)*Power(x2,10) + \
2*Power(x2,12)) + 135*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 84*Power(x1,10)*Power(x2,2) + \
476*Power(x1,8)*Power(x2,4) - 548*Power(x1,6)*Power(x2,6) + \
275*Power(x1,4)*Power(x2,8) - 16*Power(x1,2)*Power(x2,10) + \
2*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==4) && (k2==9) )
		return (he0*(-(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),13)) + lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),12)*(45*Power(x1,4) + 25*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(210*Power(x1,8) + 150*Power(x1,6)*Power(x2,2) + \
191*Power(x1,4)*Power(x2,4) + 44*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) + 3*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(1050*Power(x1,10) - 210*Power(x1,8)*Power(x2,2) + \
3675*Power(x1,6)*Power(x2,4) + 763*Power(x1,4)*Power(x2,6) + \
427*Power(x1,2)*Power(x2,8) + 23*Power(x2,10)) - \
15*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(315*Power(x1,12) - 1680*Power(x1,10)*Power(x2,2) + \
7896*Power(x1,8)*Power(x2,4) - 5376*Power(x1,6)*Power(x2,6) + \
3059*Power(x1,4)*Power(x2,8) + 56*Power(x1,2)*Power(x2,10) + \
34*Power(x2,12)) + 135*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(70*Power(x1,14) - 6517*Power(x1,12)*Power(x2,2) + \
71540*Power(x1,10)*Power(x2,4) - 214767*Power(x1,8)*Power(x2,6) + \
214662*Power(x1,6)*Power(x2,8) - 71603*Power(x1,4)*Power(x2,10) + \
6496*Power(x1,2)*Power(x2,12) - 73*Power(x2,14)) + (6227020800*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,14) - 385*Power(x1,12)*Power(x2,2) + \
4550*Power(x1,10)*Power(x2,4) - 13272*Power(x1,8)*Power(x2,6) + \
13587*Power(x1,6)*Power(x2,8) - 4361*Power(x1,4)*Power(x2,10) + \
448*Power(x1,2)*Power(x2,12) + 2*Power(x2,14)) + \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==4) && (k2==10) )
		return (he0*x2*(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),13)*(55*Power(x1,4) + 33*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) - (87178291200*(-1 + \
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
Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(990*Power(x1,8) + 880*Power(x1,6)*Power(x2,2) + \
795*Power(x1,4)*Power(x2,4) + 180*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(2310*Power(x1,10) + 990*Power(x1,8)*Power(x2,2) + \
5885*Power(x1,6)*Power(x2,4) + 2055*Power(x1,4)*Power(x2,6) + \
705*Power(x1,2)*Power(x2,8) + 31*Power(x2,10)) + \
231*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(75*Power(x1,12) - 180*Power(x1,10)*Power(x2,2) + \
900*Power(x1,8)*Power(x2,4) - 250*Power(x1,6)*Power(x2,6) + \
285*Power(x1,4)*Power(x2,8) + 30*Power(x1,2)*Power(x2,10) + \
4*Power(x2,12)) - 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(12*Power(x1,14) - 371*Power(x1,12)*Power(x2,2) + \
2436*Power(x1,10)*Power(x2,4) - 5235*Power(x1,8)*Power(x2,6) + \
4060*Power(x1,6)*Power(x2,8) - 1113*Power(x1,4)*Power(x2,10) + \
84*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
1155*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,14) - 165*Power(x1,12)*Power(x2,2) + \
1278*Power(x1,10)*Power(x2,4) - 2520*Power(x1,8)*Power(x2,6) + \
2125*Power(x1,6)*Power(x2,8) - 501*Power(x1,4)*Power(x2,10) + \
60*Power(x1,2)*Power(x2,12) + \
2*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
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
	else if ( (k1==5) && (k2==8) )
		return -(he0*x1*x2*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),13) - \
2*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(18*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(126*Power(x1,8) + 84*Power(x1,6)*Power(x2,2) + \
211*Power(x1,4)*Power(x2,4) + 50*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 6*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(210*Power(x1,10) - 126*Power(x1,8)*Power(x2,2) + \
1686*Power(x1,6)*Power(x2,4) + 49*Power(x1,4)*Power(x2,6) + \
350*Power(x1,2)*Power(x2,8) + 35*Power(x2,10)) + \
270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,12) - 833*Power(x1,10)*Power(x2,2) + \
4424*Power(x1,8)*Power(x2,4) - 7734*Power(x1,6)*Power(x2,6) + \
4424*Power(x1,4)*Power(x2,8) - 833*Power(x1,2)*Power(x2,10) + \
28*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(63*Power(x1,12) - 588*Power(x1,10)*Power(x2,2) + \
5124*Power(x1,8)*Power(x2,4) - 6684*Power(x1,6)*Power(x2,6) + \
5299*Power(x1,4)*Power(x2,8) - 448*Power(x1,2)*Power(x2,10) + \
98*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==5) && (k2==9) )
		return (he0*x1*(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),13)*(45*Power(x1,4) + 27*Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) - (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 - \
43589145600*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
10897286400*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1816214400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
227026800*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
22702680*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1891890*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
135135*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(630*Power(x1,8) + 540*Power(x1,6)*Power(x2,2) + \
843*Power(x1,4)*Power(x2,4) + 220*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) - 21*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),11)*(150*Power(x1,10) + 30*Power(x1,8)*Power(x2,2) + \
765*Power(x1,6)*Power(x2,4) + 183*Power(x1,4)*Power(x2,6) + \
145*Power(x1,2)*Power(x2,8) + 15*Power(x2,10)) + \
21*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(225*Power(x1,12) - 900*Power(x1,10)*Power(x2,2) + \
7140*Power(x1,8)*Power(x2,4) - 5130*Power(x1,6)*Power(x2,6) + \
4899*Power(x1,4)*Power(x2,8) + 110*Power(x1,2)*Power(x2,10) + \
120*Power(x2,12)) - 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(8*Power(x1,14) - 945*Power(x1,12)*Power(x2,2) + \
12180*Power(x1,10)*Power(x2,4) - 44765*Power(x1,8)*Power(x2,6) + \
57480*Power(x1,6)*Power(x2,8) - 26859*Power(x1,4)*Power(x2,10) + \
4060*Power(x1,2)*Power(x2,12) - 135*Power(x2,14)) - \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,14) - 405*Power(x1,12)*Power(x2,2) + \
6390*Power(x1,10)*Power(x2,4) - 21720*Power(x1,8)*Power(x2,6) + \
29565*Power(x1,6)*Power(x2,8) - 12837*Power(x1,4)*Power(x2,10) + \
2260*Power(x1,2)*Power(x2,12) - \
30*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==5) && (k2==10) )
		return -(he0*x1*x2*(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),15) - \
5*lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14)*(11*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) + 15*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(66*Power(x1,8) + 66*Power(x1,6)*Power(x2,2) + \
75*Power(x1,4)*Power(x2,4) + 20*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 15*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(462*Power(x1,10) + 330*Power(x1,8)*Power(x2,2) + \
1705*Power(x1,6)*Power(x2,4) + 657*Power(x1,4)*Power(x2,6) + \
305*Power(x1,2)*Power(x2,8) + 29*Power(x2,10)) + \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(55*Power(x1,12) - 88*Power(x1,10)*Power(x2,2) + \
880*Power(x1,8)*Power(x2,4) - 220*Power(x1,6)*Power(x2,6) + \
467*Power(x1,4)*Power(x2,8) + 60*Power(x1,2)*Power(x2,10) + \
14*Power(x2,12)) - 17325*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(6*Power(x1,14) - 231*Power(x1,12)*Power(x2,2) + \
1764*Power(x1,10)*Power(x2,4) - 4665*Power(x1,8)*Power(x2,6) + \
4630*Power(x1,6)*Power(x2,8) - 1785*Power(x1,4)*Power(x2,10) + \
224*Power(x1,2)*Power(x2,12) - 7*Power(x2,14)) + 20922789888000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,14) - 45*Power(x1,12)*Power(x2,2) + \
486*Power(x1,10)*Power(x2,4) - 1080*Power(x1,8)*Power(x2,6) + \
1255*Power(x1,6)*Power(x2,8) - 381*Power(x1,4)*Power(x2,10) + \
80*Power(x1,2)*Power(x2,12) + 2*Power(x2,14)) - \
10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) - 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
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
	else if ( (k1==6) && (k2==8) )
		return (he0*x2*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),13)*(36*Power(x1,4) + 23*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(3*Power(x1,2) - Power(x2,2))*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
43589145600*Power(lambdasq,13)*(3*Power(x1,2) - \
Power(x2,2))*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
10897286400*Power(lambdasq,12)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1816214400*Power(lambdasq,11)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
227026800*Power(lambdasq,10)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
22702680*Power(lambdasq,9)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1891890*Power(lambdasq,8)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
135135*Power(lambdasq,7)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(378*Power(x1,8) + 360*Power(x1,6)*Power(x2,2) + \
905*Power(x1,4)*Power(x2,4) + 240*Power(x1,2)*Power(x2,6) + \
45*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(420*Power(x1,12) + 126*Power(x1,10)*Power(x2,2) + \
4590*Power(x1,8)*Power(x2,4) + 475*Power(x1,6)*Power(x2,6) + \
1605*Power(x1,4)*Power(x2,8) + 195*Power(x1,2)*Power(x2,10) + \
5*Power(x2,12)) + 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(45*Power(x1,14) - 2175*Power(x1,12)*Power(x2,2) + \
13032*Power(x1,10)*Power(x2,4) - 29340*Power(x1,8)*Power(x2,6) + \
21845*Power(x1,6)*Power(x2,8) - 6375*Power(x1,4)*Power(x2,10) + \
390*Power(x1,2)*Power(x2,12) - 14*Power(x2,14)) + \
945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(135*Power(x1,14) - 4060*Power(x1,12)*Power(x2,2) + \
26859*Power(x1,10)*Power(x2,4) - 57480*Power(x1,8)*Power(x2,6) + \
44765*Power(x1,6)*Power(x2,8) - 12180*Power(x1,4)*Power(x2,10) + \
945*Power(x1,2)*Power(x2,12) - 8*Power(x2,14)) + \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(45*Power(x1,14) - 240*Power(x1,12)*Power(x2,2) + \
4284*Power(x1,10)*Power(x2,4) - 5580*Power(x1,8)*Power(x2,6) + \
7115*Power(x1,6)*Power(x2,8) - 750*Power(x1,4)*Power(x2,10) + \
300*Power(x1,2)*Power(x2,12) + \
10*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==6) && (k2==9) )
		return (he0*(-(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),15)) + \
15*lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14)*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(42*Power(x1,8) + 45*Power(x1,6)*Power(x2,2) + \
79*Power(x1,4)*Power(x2,4) + 23*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 15*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(210*Power(x1,12) + 168*Power(x1,10)*Power(x2,2) + \
1485*Power(x1,8)*Power(x2,4) + 452*Power(x1,6)*Power(x2,6) + \
440*Power(x1,4)*Power(x2,8) + 60*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,14) - 30*Power(x1,12)*Power(x2,2) + \
588*Power(x1,10)*Power(x2,4) - 375*Power(x1,8)*Power(x2,6) + \
647*Power(x1,6)*Power(x2,8) + 20*Power(x1,4)*Power(x2,10) + \
30*Power(x1,2)*Power(x2,12) + Power(x2,14)) + (1307674368000*(-1 + \
he0)*Power(lambdasq,15)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) \
+ 1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)))/he0 + \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,16) - 549*Power(x1,14)*Power(x2,2) + \
8085*Power(x1,12)*Power(x2,4) - 35889*Power(x1,10)*Power(x2,6) + \
57375*Power(x1,8)*Power(x2,8) - 35903*Power(x1,6)*Power(x2,10) + \
8071*Power(x1,4)*Power(x2,12) - 555*Power(x1,2)*Power(x2,14) + \
2*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,16) - 90*Power(x1,14)*Power(x2,2) + \
2265*Power(x1,12)*Power(x2,4) - 8376*Power(x1,10)*Power(x2,6) + \
15165*Power(x1,8)*Power(x2,8) - 8306*Power(x1,6)*Power(x2,10) + \
2335*Power(x1,4)*Power(x2,12) - 60*Power(x1,2)*Power(x2,14) + \
8*Power(x2,16)) + 14175*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(9*Power(x1,16) - 1072*Power(x1,14)*Power(x2,2) + \
16268*Power(x1,12)*Power(x2,4) - 71568*Power(x1,10)*Power(x2,6) + \
115030*Power(x1,8)*Power(x2,8) - 71568*Power(x1,6)*Power(x2,10) + \
16268*Power(x1,4)*Power(x2,12) - 1072*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) + 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==6) && (k2==10) )
		return -(he0*x2*(-(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),15)*(55*Power(x1,4) + 38*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(66*Power(x1,8) + 77*Power(x1,6)*Power(x2,2) + \
103*Power(x1,4)*Power(x2,4) + 31*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 15*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(462*Power(x1,12) + 528*Power(x1,10)*Power(x2,2) + \
2365*Power(x1,8)*Power(x2,4) + 1044*Power(x1,6)*Power(x2,6) + \
620*Power(x1,4)*Power(x2,8) + 84*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 15*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(1155*Power(x1,14) - 462*Power(x1,12)*Power(x2,2) + \
23892*Power(x1,10)*Power(x2,4) - 3355*Power(x1,8)*Power(x2,6) + \
19771*Power(x1,6)*Power(x2,8) + 3020*Power(x1,4)*Power(x2,10) + \
1006*Power(x1,2)*Power(x2,12) + 29*Power(x2,14)) + \
3465*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(21*Power(x1,16) - 1155*Power(x1,14)*Power(x2,2) + \
9891*Power(x1,12)*Power(x2,4) - 31911*Power(x1,10)*Power(x2,6) + \
39145*Power(x1,8)*Power(x2,8) - 20377*Power(x1,6)*Power(x2,10) + \
3745*Power(x1,4)*Power(x2,12) - 253*Power(x1,2)*Power(x2,14) - \
2*Power(x2,16)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(17*Power(x1,16) - \
680*Power(x1,14)*Power(x2,2) + 6188*Power(x1,12)*Power(x2,4) - \
19448*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
12376*Power(x1,6)*Power(x2,10) + 2380*Power(x1,4)*Power(x2,12) - \
136*Power(x1,2)*Power(x2,14) + Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),2)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),3)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 17325*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(111*Power(x1,16) - 4416*Power(x1,14)*Power(x2,2) + \
40236*Power(x1,12)*Power(x2,4) - 126384*Power(x1,10)*Power(x2,6) + \
158050*Power(x1,8)*Power(x2,8) - 80416*Power(x1,6)*Power(x2,10) + \
15484*Power(x1,4)*Power(x2,12) - 880*Power(x1,2)*Power(x2,14) + \
7*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(33*Power(x1,16) - 330*Power(x1,14)*Power(x2,2) + \
5973*Power(x1,12)*Power(x2,4) - 13728*Power(x1,10)*Power(x2,6) + \
22385*Power(x1,8)*Power(x2,8) - 8066*Power(x1,6)*Power(x2,10) + \
2795*Power(x1,4)*Power(x2,12) + 76*Power(x1,2)*Power(x2,14) + \
14*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
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
	else if ( (k1==7) && (k2==8) )
		return -(he0*x1*x2*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),15) - \
3*lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(12*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + \
21*Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(18*Power(x1,8) + 24*Power(x1,6)*Power(x2,2) + \
59*Power(x1,4)*Power(x2,4) + 18*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 105*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(12*Power(x1,12) + 18*Power(x1,10)*Power(x2,2) + \
174*Power(x1,8)*Power(x2,4) + 41*Power(x1,6)*Power(x2,6) + \
93*Power(x1,4)*Power(x2,8) + 13*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(75*Power(x1,14) - 2478*Power(x1,12)*Power(x2,2) + \
19593*Power(x1,10)*Power(x2,4) - 51000*Power(x1,8)*Power(x2,6) + \
51245*Power(x1,6)*Power(x2,8) - 19446*Power(x1,4)*Power(x2,10) + \
2527*Power(x1,2)*Power(x2,12) - 68*Power(x2,14)) + \
315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,14) - 735*Power(x1,12)*Power(x2,2) + \
4536*Power(x1,10)*Power(x2,4) - 13380*Power(x1,8)*Power(x2,6) + \
12155*Power(x1,6)*Power(x2,8) - 5271*Power(x1,4)*Power(x2,10) + \
490*Power(x1,2)*Power(x2,12) - 38*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
315*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,14) + 336*Power(x1,10)*Power(x2,4) - \
360*Power(x1,8)*Power(x2,6) + 815*Power(x1,6)*Power(x2,8) - \
84*Power(x1,4)*Power(x2,10) + 70*Power(x1,2)*Power(x2,12) + \
4*Power(x2,14)) + 10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) + 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==7) && (k2==9) )
		return -(he0*x1*(-(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),15)*(45*Power(x1,4) + 34*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(42*Power(x1,8) + 57*Power(x1,6)*Power(x2,2) + \
107*Power(x1,4)*Power(x2,4) + 35*Power(x1,2)*Power(x2,6) + \
7*Power(x2,8)) + 105*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(30*Power(x1,12) + 48*Power(x1,10)*Power(x2,2) + \
285*Power(x1,8)*Power(x2,4) + 116*Power(x1,6)*Power(x2,6) + \
124*Power(x1,4)*Power(x2,8) + 20*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 105*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,14) + 30*Power(x1,12)*Power(x2,2) + \
2172*Power(x1,10)*Power(x2,4) - 885*Power(x1,8)*Power(x2,6) + \
3509*Power(x1,6)*Power(x2,8) + 196*Power(x1,4)*Power(x2,10) + \
290*Power(x1,2)*Power(x2,12) + 19*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,16)*(Power(x1,16) - \
136*Power(x1,14)*Power(x2,2) + 2380*Power(x1,12)*Power(x2,4) - \
12376*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
19448*Power(x1,6)*Power(x2,10) + 6188*Power(x1,4)*Power(x2,12) - \
680*Power(x1,2)*Power(x2,14) + 17*Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) \
+ Power(x2,2),2)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,16) - 30*Power(x1,14)*Power(x2,2) + \
2415*Power(x1,12)*Power(x2,4) - 8736*Power(x1,10)*Power(x2,6) + \
21795*Power(x1,8)*Power(x2,8) - 13894*Power(x1,6)*Power(x2,10) + \
6097*Power(x1,4)*Power(x2,12) - 220*Power(x1,2)*Power(x2,14) + \
58*Power(x2,16)) + 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,16) - 2613*Power(x1,14)*Power(x2,2) + \
41685*Power(x1,12)*Power(x2,4) - 223377*Power(x1,10)*Power(x2,6) + \
431295*Power(x1,8)*Power(x2,8) - 350671*Power(x1,6)*Power(x2,10) + \
108871*Power(x1,4)*Power(x2,12) - 12715*Power(x1,2)*Power(x2,14) + \
226*Power(x2,16)) + 1575*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(75*Power(x1,16) - 9696*Power(x1,14)*Power(x2,2) + \
170268*Power(x1,12)*Power(x2,4) - 884688*Power(x1,10)*Power(x2,6) + \
1738410*Power(x1,8)*Power(x2,8) - 1390336*Power(x1,6)*Power(x2,10) + \
442540*Power(x1,4)*Power(x2,12) - 48592*Power(x1,2)*Power(x2,14) + \
1219*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==7) && (k2==10) )
		return (he0*x1*x2*(-(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),17)) + lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),16)*(55*Power(x1,4) + 42*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) + (711374856192000*(-1 + \
he0)*Power(lambdasq,17)*(Power(x1,2) - 3*Power(x2,2))*(3*Power(x1,2) \
- Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 355687428096000*Power(lambdasq,16)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 88921857024000*Power(lambdasq,15)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 14820309504000*Power(lambdasq,14)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1852538688000*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 185253868800*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 15437822400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1102701600*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 68918850*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3828825*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - \
Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(990*Power(x1,8) + 1375*Power(x1,6)*Power(x2,2) + \
2061*Power(x1,4)*Power(x2,4) + 693*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) + 15*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),14)*(462*Power(x1,12) + 792*Power(x1,10)*Power(x2,2) + \
3179*Power(x1,8)*Power(x2,4) + 1660*Power(x1,6)*Power(x2,6) + \
1176*Power(x1,4)*Power(x2,8) + 196*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) - 105*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),13)*(165*Power(x1,14) + 198*Power(x1,12)*Power(x2,2) + \
4356*Power(x1,10)*Power(x2,4) + 275*Power(x1,8)*Power(x2,6) + \
5229*Power(x1,6)*Power(x2,8) + 972*Power(x1,4)*Power(x2,10) + \
426*Power(x1,2)*Power(x2,12) + 27*Power(x2,14)) + \
315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(99*Power(x1,16) - 12573*Power(x1,14)*Power(x2,2) + \
113421*Power(x1,12)*Power(x2,4) - 450153*Power(x1,10)*Power(x2,6) + \
657415*Power(x1,8)*Power(x2,8) - 451623*Power(x1,6)*Power(x2,10) + \
111951*Power(x1,4)*Power(x2,12) - 13203*Power(x1,2)*Power(x2,14) - \
6*Power(x2,16)) + 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(99*Power(x1,16) - 330*Power(x1,14)*Power(x2,2) + \
20097*Power(x1,12)*Power(x2,4) - 41976*Power(x1,10)*Power(x2,6) + \
102245*Power(x1,8)*Power(x2,8) - 39330*Power(x1,6)*Power(x2,10) + \
22743*Power(x1,4)*Power(x2,12) + 804*Power(x1,2)*Power(x2,14) + \
288*Power(x2,16)) + 3465*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(513*Power(x1,16) - 22416*Power(x1,14)*Power(x2,2) + \
237132*Power(x1,12)*Power(x2,4) - 878256*Power(x1,10)*Power(x2,6) + \
1344230*Power(x1,8)*Power(x2,8) - 878256*Power(x1,6)*Power(x2,10) + \
237132*Power(x1,4)*Power(x2,12) - 22416*Power(x1,2)*Power(x2,14) + \
513*Power(x2,16))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==8) && (k2==0) )
		return (he0*x2*(Power(x1,8)*Power(Power(x1,2) + Power(x2,2),8) - \
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
	else if ( (k1==8) && (k2==1) )
		return (he0*(-(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)) + lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) - 6*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(2*Power(x1,6) + 15*Power(x1,4)*Power(x2,2) + \
35*Power(x2,6)) + 42*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),6)*(Power(x1,8) - 5*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 35*Power(x1,2)*Power(x2,6) + \
10*Power(x2,8)) + 21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(4*Power(x1,10) - 205*Power(x1,8)*Power(x2,2) + \
940*Power(x1,6)*Power(x2,4) - 950*Power(x1,4)*Power(x2,6) + \
200*Power(x1,2)*Power(x2,8) - 5*Power(x2,10)) + (362880*(-1 + \
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
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==8) && (k2==2) )
		return (he0*x2*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
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
	else if ( (k1==8) && (k2==3) )
		return -(he0*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11) - \
2*he0*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
14*Power(x2,4)) + he0*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),9)*(3*Power(x1,8) + 54*Power(x1,6)*Power(x2,2) + \
309*Power(x1,4)*Power(x2,4) + 28*Power(x1,2)*Power(x2,6) + \
210*Power(x2,8)) - \
30*he0*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,10) + 18*Power(x1,8)*Power(x2,2) - \
37*Power(x1,6)*Power(x2,4) + 154*Power(x1,4)*Power(x2,6) - \
42*Power(x1,2)*Power(x2,8) + 14*Power(x2,10)) + 39916800*(-1 + \
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
Power(x2,2),7)*(10*Power(x1,12) - 240*Power(x1,10)*Power(x2,2) + \
2115*Power(x1,8)*Power(x2,4) - 3724*Power(x1,6)*Power(x2,6) + \
2100*Power(x1,4)*Power(x2,8) - 252*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12)) + 19958400*he0*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*he0*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16)))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==8) && (k2==4) )
		return (he0*x2*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
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
	else if ( (k1==8) && (k2==5) )
		return (he0*(-(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13)) + lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(15*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) - \
3*Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,8) + 50*Power(x1,6)*Power(x2,2) + \
229*Power(x1,4)*Power(x2,4) + 56*Power(x1,2)*Power(x2,6) + \
70*Power(x2,8)) + 3*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(5*Power(x1,12) + 105*Power(x1,10)*Power(x2,2) + \
1225*Power(x1,8)*Power(x2,4) - 483*Power(x1,6)*Power(x2,6) + \
2758*Power(x1,4)*Power(x2,8) - 70*Power(x1,2)*Power(x2,10) + \
140*Power(x2,12)) + 45*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(212*Power(x1,14) - 19537*Power(x1,12)*Power(x2,2) + \
214662*Power(x1,10)*Power(x2,4) - 644231*Power(x1,8)*Power(x2,6) + \
644056*Power(x1,6)*Power(x2,8) - 214767*Power(x1,4)*Power(x2,10) + \
19502*Power(x1,2)*Power(x2,12) - 217*Power(x2,14)) + \
45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(22*Power(x1,14) - 1162*Power(x1,12)*Power(x2,2) + \
13587*Power(x1,10)*Power(x2,4) - 39991*Power(x1,8)*Power(x2,6) + \
40516*Power(x1,6)*Power(x2,8) - 13272*Power(x1,4)*Power(x2,10) + \
1267*Power(x1,2)*Power(x2,12) - 7*Power(x2,14)) + (6227020800*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(8*Power(x1,14) + 322*Power(x1,12)*Power(x2,2) - \
1792*Power(x1,10)*Power(x2,4) + 7511*Power(x1,8)*Power(x2,6) - \
5936*Power(x1,6)*Power(x2,8) + 2632*Power(x1,4)*Power(x2,10) - \
112*Power(x1,2)*Power(x2,12) + 7*Power(x2,14)) + \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==8) && (k2==6) )
		return (he0*x2*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
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
	else if ( (k1==8) && (k2==7) )
		return (he0*(-(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),15)) + \
2*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(14*Power(x1,4) + 13*Power(x1,2)*Power(x2,2) + \
14*Power(x2,4)) - \
210*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,8) + 2*Power(x1,6)*Power(x2,2) + \
6*Power(x1,4)*Power(x2,4) + 2*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
+ 420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(Power(x1,4) + Power(x2,4))*(Power(x1,8) + \
4*Power(x1,6)*Power(x2,2) + 32*Power(x1,4)*Power(x2,4) + \
4*Power(x1,2)*Power(x2,6) + Power(x2,8)) + 1307674368000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,16) - \
120*Power(x1,14)*Power(x2,2) + 1820*Power(x1,12)*Power(x2,4) - \
8008*Power(x1,10)*Power(x2,6) + 12870*Power(x1,8)*Power(x2,8) - \
8008*Power(x1,6)*Power(x2,10) + 1820*Power(x1,4)*Power(x2,12) - \
120*Power(x1,2)*Power(x2,14) + Power(x2,16)) - \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,16) + 20*Power(x1,14)*Power(x2,2) + \
490*Power(x1,12)*Power(x2,4) - 700*Power(x1,10)*Power(x2,6) + \
2650*Power(x1,8)*Power(x2,8) - 700*Power(x1,6)*Power(x2,10) + \
490*Power(x1,4)*Power(x2,12) + 20*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 630*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,16) + 90*Power(x1,14)*Power(x2,2) - \
910*Power(x1,12)*Power(x2,4) + 4718*Power(x1,10)*Power(x2,6) - \
6870*Power(x1,8)*Power(x2,8) + 4718*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 3150*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,16) - 262*Power(x1,14)*Power(x2,2) + \
4088*Power(x1,12)*Power(x2,4) - 17850*Power(x1,10)*Power(x2,6) + \
28810*Power(x1,8)*Power(x2,8) - 17850*Power(x1,6)*Power(x2,10) + \
4088*Power(x1,4)*Power(x2,12) - 262*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) - 12600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(10*Power(x1,16) - 1207*Power(x1,14)*Power(x2,2) + \
18298*Power(x1,12)*Power(x2,4) - 80521*Power(x1,10)*Power(x2,6) + \
129400*Power(x1,8)*Power(x2,8) - 80521*Power(x1,6)*Power(x2,10) + \
18298*Power(x1,4)*Power(x2,12) - 1207*Power(x1,2)*Power(x2,14) + \
10*Power(x2,16)) - 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==8) && (k2==8) )
		return (he0*x2*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),16) - \
4*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(9*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(63*Power(x1,8) + 114*Power(x1,6)*Power(x2,2) + \
274*Power(x1,4)*Power(x2,4) + 98*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - \
420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,12) + 9*Power(x1,10)*Power(x2,2) + \
57*Power(x1,8)*Power(x2,4) + 24*Power(x1,6)*Power(x2,6) + \
43*Power(x1,4)*Power(x2,8) + 7*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(17*Power(x1,16) - \
680*Power(x1,14)*Power(x2,2) + 6188*Power(x1,12)*Power(x2,4) - \
19448*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
12376*Power(x1,6)*Power(x2,10) + 2380*Power(x1,4)*Power(x2,12) - \
136*Power(x1,2)*Power(x2,14) + Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),2)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),3)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(9*Power(x1,16) + 60*Power(x1,14)*Power(x2,2) + \
1218*Power(x1,12)*Power(x2,4) - 708*Power(x1,10)*Power(x2,6) + \
3970*Power(x1,8)*Power(x2,8) - 244*Power(x1,6)*Power(x2,10) + \
610*Power(x1,4)*Power(x2,12) + 44*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 1260*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,16) + 195*Power(x1,14)*Power(x2,2) - \
1050*Power(x1,12)*Power(x2,4) + 4359*Power(x1,10)*Power(x2,6) - \
4460*Power(x1,8)*Power(x2,8) + 2849*Power(x1,6)*Power(x2,10) - \
350*Power(x1,4)*Power(x2,12) + 53*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 6300*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(303*Power(x1,16) - 12162*Power(x1,14)*Power(x2,2) + \
110586*Power(x1,12)*Power(x2,4) - 347682*Power(x1,10)*Power(x2,6) + \
434480*Power(x1,8)*Power(x2,8) - 221270*Power(x1,6)*Power(x2,10) + \
42518*Power(x1,4)*Power(x2,12) - 2438*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(183*Power(x1,16) - 5850*Power(x1,14)*Power(x2,2) + \
56028*Power(x1,12)*Power(x2,4) - 172518*Power(x1,10)*Power(x2,6) + \
218710*Power(x1,8)*Power(x2,8) - 109606*Power(x1,6)*Power(x2,10) + \
21700*Power(x1,4)*Power(x2,12) - 1114*Power(x1,2)*Power(x2,14) + \
19*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==8) && (k2==9) )
		return (he0*(-(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),17)) + lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),16)*(45*Power(x1,4) + 39*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) - \
2*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(315*Power(x1,8) + 540*Power(x1,6)*Power(x2,2) + \
1056*Power(x1,4)*Power(x2,4) + 392*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) + \
30*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(105*Power(x1,12) + 273*Power(x1,10)*Power(x2,2) + \
1314*Power(x1,8)*Power(x2,4) + 716*Power(x1,6)*Power(x2,6) + \
791*Power(x1,4)*Power(x2,8) + 147*Power(x1,2)*Power(x2,10) + \
14*Power(x2,12)) - \
105*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(45*Power(x1,16) + 180*Power(x1,14)*Power(x2,2) + \
2706*Power(x1,12)*Power(x2,4) - 72*Power(x1,10)*Power(x2,6) + \
5958*Power(x1,8)*Power(x2,8) + 636*Power(x1,6)*Power(x2,10) + \
810*Power(x1,4)*Power(x2,12) + 72*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 1575*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(1212*Power(x1,18) - 186003*Power(x1,16)*Power(x2,2) + \
3719304*Power(x1,14)*Power(x2,4) - 22564836*Power(x1,12)*Power(x2,6) \
+ 53187408*Power(x1,10)*Power(x2,8) - \
53188290*Power(x1,8)*Power(x2,10) + 22564248*Power(x1,6)*Power(x2,12) \
- 3719556*Power(x1,4)*Power(x2,14) + 185940*Power(x1,2)*Power(x2,16) \
- 1219*Power(x2,18)) + 630*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(183*Power(x1,18) - 22959*Power(x1,16)*Power(x2,2) + \
466110*Power(x1,14)*Power(x2,4) - 2817738*Power(x1,12)*Power(x2,6) + \
6652836*Power(x1,10)*Power(x2,8) - 6644016*Power(x1,8)*Power(x2,10) + \
2823618*Power(x1,6)*Power(x2,12) - 463590*Power(x1,4)*Power(x2,14) + \
23589*Power(x1,2)*Power(x2,16) - 113*Power(x2,18)) + \
(355687428096000*(-1 + he0)*Power(lambdasq,17)*(Power(x1,18) - \
153*Power(x1,16)*Power(x2,2) + 3060*Power(x1,14)*Power(x2,4) - \
18564*Power(x1,12)*Power(x2,6) + 43758*Power(x1,10)*Power(x2,8) - \
43758*Power(x1,8)*Power(x2,10) + 18564*Power(x1,6)*Power(x2,12) - \
3060*Power(x1,4)*Power(x2,14) + 153*Power(x1,2)*Power(x2,16) - \
Power(x2,18)))/he0 + \
44460928512000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
7410154752000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
926269344000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
92626934400*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
7718911200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
551350800*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(9*Power(x1,18) + 135*Power(x1,16)*Power(x2,2) + \
8010*Power(x1,14)*Power(x2,4) - 23646*Power(x1,12)*Power(x2,6) + \
90162*Power(x1,10)*Power(x2,8) - 60174*Power(x1,8)*Power(x2,10) + \
41874*Power(x1,6)*Power(x2,12) - 1710*Power(x1,4)*Power(x2,14) + \
873*Power(x1,2)*Power(x2,16) + 19*Power(x2,18)) - \
630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(6*Power(x1,18) + 1287*Power(x1,16)*Power(x2,2) - \
19860*Power(x1,14)*Power(x2,4) + 131754*Power(x1,12)*Power(x2,6) - \
296388*Power(x1,10)*Power(x2,8) + 308148*Power(x1,8)*Power(x2,10) - \
124404*Power(x1,6)*Power(x2,12) + 22590*Power(x1,4)*Power(x2,14) - \
762*Power(x1,2)*Power(x2,16) + 29*Power(x2,18)) + \
177843714048000*Power(lambdasq,16)*(Power(x1,20) - \
152*Power(x1,18)*Power(x2,2) + 2907*Power(x1,16)*Power(x2,4) - \
15504*Power(x1,14)*Power(x2,6) + 25194*Power(x1,12)*Power(x2,8) - \
25194*Power(x1,8)*Power(x2,12) + 15504*Power(x1,6)*Power(x2,14) - \
2907*Power(x1,4)*Power(x2,16) + 152*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==8) && (k2==10) )
		return (he0*x2*(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),18) - lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),17)*(55*Power(x1,4) + 47*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)*(165*Power(x1,8) + 275*Power(x1,6)*Power(x2,2) + \
447*Power(x1,4)*Power(x2,4) + 168*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - \
6*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(1155*Power(x1,12) + 2805*Power(x1,10)*Power(x2,2) + \
10450*Power(x1,8)*Power(x2,4) + 6556*Power(x1,6)*Power(x2,6) + \
5229*Power(x1,4)*Power(x2,8) + 1015*Power(x1,2)*Power(x2,10) + \
70*Power(x2,12)) + 15*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),14)*(1155*Power(x1,16) + 3696*Power(x1,14)*Power(x2,2) \
+ 38874*Power(x1,12)*Power(x2,4) + 13816*Power(x1,10)*Power(x2,6) + \
63554*Power(x1,8)*Power(x2,8) + 14896*Power(x1,6)*Power(x2,10) + \
7882*Power(x1,4)*Power(x2,12) + 728*Power(x1,2)*Power(x2,14) + \
7*Power(x2,16)) - 31185*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(1048*Power(x1,18) - 53553*Power(x1,16)*Power(x2,2) + \
642384*Power(x1,14)*Power(x2,4) - 2784084*Power(x1,12)*Power(x2,6) + \
5103664*Power(x1,10)*Power(x2,8) - 4176126*Power(x1,8)*Power(x2,10) + \
1498896*Power(x1,6)*Power(x2,12) - 214212*Power(x1,4)*Power(x2,14) + \
9432*Power(x1,2)*Power(x2,16) - 57*Power(x2,18)) + \
6402373705728000*(-1 + 1/he0)*Power(lambdasq,18)*(19*Power(x1,18) - \
969*Power(x1,16)*Power(x2,2) + 11628*Power(x1,14)*Power(x2,4) - \
50388*Power(x1,12)*Power(x2,6) + 92378*Power(x1,10)*Power(x2,8) - \
75582*Power(x1,8)*Power(x2,10) + 27132*Power(x1,6)*Power(x2,12) - \
3876*Power(x1,4)*Power(x2,14) + 171*Power(x1,2)*Power(x2,16) - \
Power(x2,18)) - 3201186852864000*Power(lambdasq,17)*(Power(x1,2) + \
Power(x2,2))*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
1890*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(891*Power(x1,18) - 39281*Power(x1,16)*Power(x2,2) + \
485078*Power(x1,14)*Power(x2,4) - 2080958*Power(x1,12)*Power(x2,6) + \
3837548*Power(x1,10)*Power(x2,8) - 3123152*Power(x1,8)*Power(x2,10) + \
1129562*Power(x1,6)*Power(x2,12) - 158594*Power(x1,4)*Power(x2,14) + \
7529*Power(x1,2)*Power(x2,16) + Power(x2,18)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(11*Power(x1,18) + 55*Power(x1,16)*Power(x2,2) + \
2574*Power(x1,14)*Power(x2,4) - 3894*Power(x1,12)*Power(x2,6) + \
16654*Power(x1,10)*Power(x2,8) - 5838*Power(x1,8)*Power(x2,10) + \
6054*Power(x1,6)*Power(x2,12) + 298*Power(x1,4)*Power(x2,14) + \
147*Power(x1,2)*Power(x2,16) + 3*Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(11*Power(x1,18) + 2134*Power(x1,16)*Power(x2,2) - \
18062*Power(x1,14)*Power(x2,4) + 91872*Power(x1,12)*Power(x2,6) - \
152262*Power(x1,10)*Power(x2,8) + 137318*Power(x1,8)*Power(x2,10) - \
42798*Power(x1,6)*Power(x2,12) + 8116*Power(x1,4)*Power(x2,14) - \
41*Power(x1,2)*Power(x2,16) + \
16*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==9) && (k2==0) )
		return -(he0*x1*x2*(Power(x1,8)*Power(Power(x1,2) + Power(x2,2),9) \
- 18*lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,2) + 2*Power(x2,2)) + (725760*(-1 + \
he0)*Power(lambdasq,9)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 + 362880*Power(lambdasq,8)*(Power(x1,2) + \
Power(x2,2))*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 90720*Power(lambdasq,7)*Power(Power(x1,2) + \
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
5*Power(x2,4)) - \
252*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,4) - 6*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==9) && (k2==1) )
		return (he0*x1*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10) - lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(Power(x1,6) + 8*Power(x1,4)*Power(x2,2) + \
8*Power(x1,2)*Power(x2,4) + 21*Power(x2,6)) + \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 80*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) - 90*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,8) + Power(x1,6)*Power(x2,2) + \
29*Power(x1,4)*Power(x2,4) - 21*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) + 3628800*(-1 + \
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
165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==9) && (k2==2) )
		return -(he0*x1*x2*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) + 2*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(24*Power(x1,6) + 127*Power(x1,4)*Power(x2,2) + \
72*Power(x1,2)*Power(x2,4) + 189*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,8) + Power(x1,6)*Power(x2,2) + \
51*Power(x1,4)*Power(x2,4) - 21*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) - 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(12*Power(x1,10) - 325*Power(x1,8)*Power(x2,2) + \
1044*Power(x1,6)*Power(x2,4) - 1134*Power(x1,4)*Power(x2,6) + \
280*Power(x1,2)*Power(x2,8) - 21*Power(x2,10)) + 159667200*(-1 + \
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
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==9) && (k2==3) )
		return (he0*x1*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12) - \
6*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,8) + 30*Power(x1,6)*Power(x2,2) + \
139*Power(x1,4)*Power(x2,4) + 60*Power(x1,2)*Power(x2,6) + \
126*Power(x2,8)) - 12*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(4*Power(x1,10) + 66*Power(x1,8)*Power(x2,2) + \
4*Power(x1,6)*Power(x2,4) + 591*Power(x1,4)*Power(x2,6) - \
126*Power(x1,2)*Power(x2,8) + 105*Power(x2,10)) + \
135*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,12) - 16*Power(x1,10)*Power(x2,2) + \
275*Power(x1,8)*Power(x2,4) - 548*Power(x1,6)*Power(x2,6) + \
476*Power(x1,4)*Power(x2,8) - 84*Power(x1,2)*Power(x2,10) + \
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
13*Power(x2,12)) + 270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(2*Power(x1,12) - 219*Power(x1,10)*Power(x2,2) + \
1955*Power(x1,8)*Power(x2,4) - 4734*Power(x1,6)*Power(x2,6) + \
3528*Power(x1,4)*Power(x2,8) - 791*Power(x1,2)*Power(x2,10) + \
35*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==9) && (k2==4) )
		return -(he0*x1*x2*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13) - \
2*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(5*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
18*Power(x2,4)) + (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(5*Power(x1,8) + 50*Power(x1,6)*Power(x2,2) + \
211*Power(x1,4)*Power(x2,4) + 84*Power(x1,2)*Power(x2,6) + \
126*Power(x2,8)) - 6*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),10)*(35*Power(x1,10) + 350*Power(x1,8)*Power(x2,2) + \
49*Power(x1,6)*Power(x2,4) + 1686*Power(x1,4)*Power(x2,6) - \
126*Power(x1,2)*Power(x2,8) + 210*Power(x2,10)) + \
270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,12) - 833*Power(x1,10)*Power(x2,2) + \
4424*Power(x1,8)*Power(x2,4) - 7734*Power(x1,6)*Power(x2,6) + \
4424*Power(x1,4)*Power(x2,8) - 833*Power(x1,2)*Power(x2,10) + \
28*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(98*Power(x1,12) - 448*Power(x1,10)*Power(x2,2) + \
5299*Power(x1,8)*Power(x2,4) - 6684*Power(x1,6)*Power(x2,6) + \
5124*Power(x1,4)*Power(x2,8) - 588*Power(x1,2)*Power(x2,10) + \
63*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==9) && (k2==5) )
		return (he0*x1*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(15*Power(x1,4) + 23*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) - (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 - \
43589145600*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
10897286400*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1816214400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
227026800*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
22702680*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
1891890*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
135135*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + 5*Power(x2,4))*(Power(x1,8) - \
92*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,8) + 240*Power(x1,6)*Power(x2,2) + \
905*Power(x1,4)*Power(x2,4) + 360*Power(x1,2)*Power(x2,6) + \
378*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),11)*(5*Power(x1,12) + 195*Power(x1,10)*Power(x2,2) + \
1605*Power(x1,8)*Power(x2,4) + 475*Power(x1,6)*Power(x2,6) + \
4590*Power(x1,4)*Power(x2,8) + 126*Power(x1,2)*Power(x2,10) + \
420*Power(x2,12)) - 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(8*Power(x1,14) - 945*Power(x1,12)*Power(x2,2) + \
12180*Power(x1,10)*Power(x2,4) - 44765*Power(x1,8)*Power(x2,6) + \
57480*Power(x1,6)*Power(x2,8) - 26859*Power(x1,4)*Power(x2,10) + \
4060*Power(x1,2)*Power(x2,12) - 135*Power(x2,14)) - \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(14*Power(x1,14) - 390*Power(x1,12)*Power(x2,2) + \
6375*Power(x1,10)*Power(x2,4) - 21845*Power(x1,8)*Power(x2,6) + \
29340*Power(x1,6)*Power(x2,8) - 13032*Power(x1,4)*Power(x2,10) + \
2175*Power(x1,2)*Power(x2,12) - 45*Power(x2,14)) + \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(10*Power(x1,14) + 300*Power(x1,12)*Power(x2,2) - \
750*Power(x1,10)*Power(x2,4) + 7115*Power(x1,8)*Power(x2,6) - \
5580*Power(x1,6)*Power(x2,8) + 4284*Power(x1,4)*Power(x2,10) - \
240*Power(x1,2)*Power(x2,12) + \
45*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==9) && (k2==6) )
		return -(he0*x1*x2*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15) - \
3*lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
12*Power(x2,4)) + \
21*Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(5*Power(x1,8) + 18*Power(x1,6)*Power(x2,2) + \
59*Power(x1,4)*Power(x2,4) + 24*Power(x1,2)*Power(x2,6) + \
18*Power(x2,8)) - 105*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(Power(x1,12) + 13*Power(x1,10)*Power(x2,2) + \
93*Power(x1,8)*Power(x2,4) + 41*Power(x1,6)*Power(x2,6) + \
174*Power(x1,4)*Power(x2,8) + 18*Power(x1,2)*Power(x2,10) + \
12*Power(x2,12)) - 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(68*Power(x1,14) - 2527*Power(x1,12)*Power(x2,2) + \
19446*Power(x1,10)*Power(x2,4) - 51245*Power(x1,8)*Power(x2,6) + \
51000*Power(x1,6)*Power(x2,8) - 19593*Power(x1,4)*Power(x2,10) + \
2478*Power(x1,2)*Power(x2,12) - 75*Power(x2,14)) - \
315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(38*Power(x1,14) - 490*Power(x1,12)*Power(x2,2) + \
5271*Power(x1,10)*Power(x2,4) - 12155*Power(x1,8)*Power(x2,6) + \
13380*Power(x1,6)*Power(x2,8) - 4536*Power(x1,4)*Power(x2,10) + \
735*Power(x1,2)*Power(x2,12) - 3*Power(x2,14)) + 20922789888000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
315*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(4*Power(x1,14) + 70*Power(x1,12)*Power(x2,2) - \
84*Power(x1,10)*Power(x2,4) + 815*Power(x1,8)*Power(x2,6) - \
360*Power(x1,6)*Power(x2,8) + 336*Power(x1,4)*Power(x2,10) + \
3*Power(x2,14)) - 10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) - 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==9) && (k2==7) )
		return (he0*x1*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),16) - \
4*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(7*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
9*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(35*Power(x1,8) + 98*Power(x1,6)*Power(x2,2) + \
274*Power(x1,4)*Power(x2,4) + 114*Power(x1,2)*Power(x2,6) + \
63*Power(x2,8)) - \
420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,12) + 7*Power(x1,10)*Power(x2,2) + \
43*Power(x1,8)*Power(x2,4) + 24*Power(x1,6)*Power(x2,6) + \
57*Power(x1,4)*Power(x2,8) + 9*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) - 1260*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,16) + 53*Power(x1,14)*Power(x2,2) - \
350*Power(x1,12)*Power(x2,4) + 2849*Power(x1,10)*Power(x2,6) - \
4460*Power(x1,8)*Power(x2,8) + 4359*Power(x1,6)*Power(x2,10) - \
1050*Power(x1,4)*Power(x2,12) + 195*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(Power(x1,16) + 44*Power(x1,14)*Power(x2,2) + \
610*Power(x1,12)*Power(x2,4) - 244*Power(x1,10)*Power(x2,6) + \
3970*Power(x1,8)*Power(x2,8) - 708*Power(x1,6)*Power(x2,10) + \
1218*Power(x1,4)*Power(x2,12) + 60*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) \
+ 2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) \
+ Power(x2,2),2)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(19*Power(x1,16) - 1114*Power(x1,14)*Power(x2,2) + \
21700*Power(x1,12)*Power(x2,4) - 109606*Power(x1,10)*Power(x2,6) + \
218710*Power(x1,8)*Power(x2,8) - 172518*Power(x1,6)*Power(x2,10) + \
56028*Power(x1,4)*Power(x2,12) - 5850*Power(x1,2)*Power(x2,14) + \
183*Power(x2,16)) + 6300*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(17*Power(x1,16) - 2438*Power(x1,14)*Power(x2,2) + \
42518*Power(x1,12)*Power(x2,4) - 221270*Power(x1,10)*Power(x2,6) + \
434480*Power(x1,8)*Power(x2,8) - 347682*Power(x1,6)*Power(x2,10) + \
110586*Power(x1,4)*Power(x2,12) - 12162*Power(x1,2)*Power(x2,14) + \
303*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==9) && (k2==8) )
		return -(he0*x1*x2*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17) - \
2*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)*(18*Power(x1,4) + 19*Power(x1,2)*Power(x2,2) + \
18*Power(x2,4)) + (711374856192000*(-1 + \
he0)*Power(lambdasq,17)*(Power(x1,2) - 3*Power(x2,2))*(3*Power(x1,2) \
- Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 355687428096000*Power(lambdasq,16)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 88921857024000*Power(lambdasq,15)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 14820309504000*Power(lambdasq,14)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1852538688000*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 185253868800*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 15437822400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1102701600*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 68918850*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3828825*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + \
2*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(189*Power(x1,8) + 450*Power(x1,6)*Power(x2,2) + \
1066*Power(x1,4)*Power(x2,4) + 450*Power(x1,2)*Power(x2,6) + \
189*Power(x2,8)) - \
12*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(105*Power(x1,12) + 504*Power(x1,10)*Power(x2,2) + \
2601*Power(x1,8)*Power(x2,4) + 1684*Power(x1,6)*Power(x2,6) + \
2601*Power(x1,4)*Power(x2,8) + 504*Power(x1,2)*Power(x2,10) + \
105*Power(x2,12)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(9*Power(x1,16) + 132*Power(x1,14)*Power(x2,2) + \
1530*Power(x1,12)*Power(x2,4) + 180*Power(x1,10)*Power(x2,6) + \
6250*Power(x1,8)*Power(x2,8) + 180*Power(x1,6)*Power(x2,10) + \
1530*Power(x1,4)*Power(x2,12) + 132*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) - 210*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,16) + 1362*Power(x1,14)*Power(x2,2) - \
4662*Power(x1,12)*Power(x2,4) + 33462*Power(x1,10)*Power(x2,6) - \
34270*Power(x1,8)*Power(x2,8) + 33462*Power(x1,6)*Power(x2,10) - \
4662*Power(x1,4)*Power(x2,12) + 1362*Power(x1,2)*Power(x2,14) + \
45*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(207*Power(x1,16) - 4974*Power(x1,14)*Power(x2,2) + \
61488*Power(x1,12)*Power(x2,4) - 215154*Power(x1,10)*Power(x2,6) + \
341570*Power(x1,8)*Power(x2,8) - 215154*Power(x1,6)*Power(x2,10) + \
61488*Power(x1,4)*Power(x2,12) - 4974*Power(x1,2)*Power(x2,14) + \
207*Power(x2,16)) + 2520*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(666*Power(x1,16) - 31137*Power(x1,14)*Power(x2,2) + \
324954*Power(x1,12)*Power(x2,4) - 1209807*Power(x1,10)*Power(x2,6) + \
1845560*Power(x1,8)*Power(x2,8) - 1209807*Power(x1,6)*Power(x2,10) + \
324954*Power(x1,4)*Power(x2,12) - 31137*Power(x1,2)*Power(x2,14) + \
666*Power(x2,16))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==9) && (k2==9) )
		return (he0*x1*(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),18) - \
9*lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17)*(5*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
4*Power(x2,4)) + \
18*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)*(35*Power(x1,8) + 75*Power(x1,6)*Power(x2,2) + \
151*Power(x1,4)*Power(x2,4) + 64*Power(x1,2)*Power(x2,6) + \
21*Power(x2,8)) - \
18*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(175*Power(x1,12) + 665*Power(x1,10)*Power(x2,2) + \
2850*Power(x1,8)*Power(x2,4) + 2028*Power(x1,6)*Power(x2,6) + \
2257*Power(x1,4)*Power(x2,8) + 483*Power(x1,2)*Power(x2,10) + \
70*Power(x2,12)) + \
135*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(35*Power(x1,16) + 280*Power(x1,14)*Power(x2,2) + \
2674*Power(x1,12)*Power(x2,4) + 1256*Power(x1,10)*Power(x2,6) + \
7554*Power(x1,8)*Power(x2,8) + 1416*Power(x1,6)*Power(x2,10) + \
1554*Power(x1,4)*Power(x2,12) + 168*Power(x1,2)*Power(x2,14) + \
7*Power(x2,16)) - 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(592*Power(x1,18) - 104067*Power(x1,16)*Power(x2,2) + \
2355072*Power(x1,14)*Power(x2,4) - 16490796*Power(x1,12)*Power(x2,6) \
+ 45932976*Power(x1,10)*Power(x2,8) - \
56144714*Power(x1,8)*Power(x2,10) + 30621984*Power(x1,6)*Power(x2,12) \
- 7067484*Power(x1,4)*Power(x2,14) + 588768*Power(x1,2)*Power(x2,16) \
- 11563*Power(x2,18)) - 5670*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(23*Power(x1,18) - 2253*Power(x1,16)*Power(x2,2) + \
54078*Power(x1,14)*Power(x2,4) - 373254*Power(x1,12)*Power(x2,6) + \
1046604*Power(x1,10)*Power(x2,8) - 1272976*Power(x1,8)*Power(x2,10) + \
698226*Power(x1,6)*Power(x2,12) - 159546*Power(x1,4)*Power(x2,14) + \
13677*Power(x1,2)*Power(x2,16) - 227*Power(x2,18)) + \
6402373705728000*(-1 + 1/he0)*Power(lambdasq,18)*(Power(x1,18) - \
171*Power(x1,16)*Power(x2,2) + 3876*Power(x1,14)*Power(x2,4) - \
27132*Power(x1,12)*Power(x2,6) + 75582*Power(x1,10)*Power(x2,8) - \
92378*Power(x1,8)*Power(x2,10) + 50388*Power(x1,6)*Power(x2,12) - \
11628*Power(x1,4)*Power(x2,14) + 969*Power(x1,2)*Power(x2,16) - \
19*Power(x2,18)) - 3201186852864000*Power(lambdasq,17)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,18) + 45*Power(x1,16)*Power(x2,2) + \
1050*Power(x1,14)*Power(x2,4) - 1746*Power(x1,12)*Power(x2,6) + \
13626*Power(x1,10)*Power(x2,8) - 8090*Power(x1,8)*Power(x2,10) + \
9666*Power(x1,6)*Power(x2,12) - 306*Power(x1,4)*Power(x2,14) + \
393*Power(x1,2)*Power(x2,16) + 17*Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(5*Power(x1,18) + 468*Power(x1,16)*Power(x2,2) - \
5610*Power(x1,14)*Power(x2,4) + 50148*Power(x1,12)*Power(x2,6) - \
124578*Power(x1,10)*Power(x2,8) + 166178*Power(x1,8)*Power(x2,10) - \
82170*Power(x1,6)*Power(x2,12) + 22248*Power(x1,4)*Power(x2,14) - \
1119*Power(x1,2)*Power(x2,16) + \
94*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==9) && (k2==10) )
		return -(he0*x1*x2*(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),19) - (486580401635328000*(-1 + \
he0)*Power(lambdasq,19)*(x1 - x2)*(x1 + x2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 - 243290200817664000*Power(lambdasq,18)*(x1 - \
x2)*(x1 + x2)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 60822550204416000*Power(lambdasq,17)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 10137091700736000*Power(lambdasq,16)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),3)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 1267136462592000*Power(lambdasq,15)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 126713646259200*Power(lambdasq,14)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),5)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 10559470521600*Power(lambdasq,13)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),6)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 754247894400*Power(lambdasq,12)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),7)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 47140493400*Power(lambdasq,11)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),8)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 2618916300*Power(lambdasq,10)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),9)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - 130945815*Power(lambdasq,9)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),10)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),18)*(55*Power(x1,4) + 53*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) + \
18*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17)*(55*Power(x1,8) + 110*Power(x1,6)*Power(x2,2) + \
190*Power(x1,4)*Power(x2,4) + 80*Power(x1,2)*Power(x2,6) + \
21*Power(x2,8)) - \
18*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(385*Power(x1,12) + 1265*Power(x1,10)*Power(x2,2) + \
4510*Power(x1,8)*Power(x2,4) + 3400*Power(x1,6)*Power(x2,6) + \
2935*Power(x1,4)*Power(x2,8) + 651*Power(x1,2)*Power(x2,10) + \
70*Power(x2,12)) + 9*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),15)*(1925*Power(x1,16) + 10780*Power(x1,14)*Power(x2,2) \
+ 83050*Power(x1,12)*Power(x2,4) + 59400*Power(x1,10)*Power(x2,6) + \
175590*Power(x1,8)*Power(x2,8) + 52820*Power(x1,6)*Power(x2,10) + \
31346*Power(x1,4)*Power(x2,12) + 3640*Power(x1,2)*Power(x2,14) + \
105*Power(x2,16)) - 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(10340*Power(x1,18) - 599775*Power(x1,16)*Power(x2,2) \
+ 8131992*Power(x1,14)*Power(x2,4) - \
40701540*Power(x1,12)*Power(x2,6) + 88138160*Power(x1,10)*Power(x2,8) \
- 88177850*Power(x1,8)*Power(x2,10) + \
40675080*Power(x1,6)*Power(x2,12) - 8143332*Power(x1,4)*Power(x2,14) \
+ 596940*Power(x1,2)*Power(x2,16) - 10655*Power(x2,18)) - \
1890*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),12)*(935*Power(x1,18) - 34815*Power(x1,16)*Power(x2,2) + \
519222*Power(x1,14)*Power(x2,4) - 2517570*Power(x1,12)*Power(x2,6) + \
5549060*Power(x1,10)*Power(x2,8) - 5469680*Power(x1,8)*Power(x2,10) + \
2570490*Power(x1,6)*Power(x2,12) - 496542*Power(x1,4)*Power(x2,14) + \
40485*Power(x1,2)*Power(x2,16) - 305*Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(44*Power(x1,18) + 2343*Power(x1,16)*Power(x2,2) - \
15048*Power(x1,14)*Power(x2,4) + 109890*Power(x1,12)*Power(x2,6) - \
192280*Power(x1,10)*Power(x2,8) + 232264*Power(x1,8)*Power(x2,10) - \
84312*Power(x1,6)*Power(x2,12) + 25086*Power(x1,4)*Power(x2,14) - \
180*Power(x1,2)*Power(x2,16) + 145*Power(x2,18)) - \
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),14)*(77*Power(x1,18) + 1155*Power(x1,16)*Power(x2,2) + \
21714*Power(x1,14)*Power(x2,4) - 12870*Power(x1,12)*Power(x2,6) + \
168410*Power(x1,10)*Power(x2,8) - 39638*Power(x1,8)*Power(x2,10) + \
92250*Power(x1,6)*Power(x2,12) + 6762*Power(x1,4)*Power(x2,14) + \
3885*Power(x1,2)*Power(x2,16) + \
175*Power(x2,18))))/(2.*Power(lambdasq,19)*Pi*Power(Power(x1,2) + \
Power(x2,2),20));
	else if ( (k1==10) && (k2==0) )
		return (he0*x2*(Power(x1,10)*Power(Power(x1,2) + Power(x2,2),10) - \
5*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(5*Power(x1,2) + 9*Power(x2,2)) + \
90*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) - 90*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,6) + 25*Power(x1,4)*Power(x2,2) - \
7*Power(x1,2)*Power(x2,4) + 35*Power(x2,6)) + \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 80*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + (3628800*(-1 + \
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
	else if ( (k1==10) && (k2==1) )
		return (he0*(-(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)) + lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,4) + 24*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - 5*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(5*Power(x1,6) + 48*Power(x1,4)*Power(x2,2) + \
81*Power(x1,2)*Power(x2,4) + 126*Power(x2,6)) + \
90*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 14*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - 45*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(6*Power(x1,10) + 129*Power(x1,8)*Power(x2,2) - \
600*Power(x1,6)*Power(x2,4) + 1358*Power(x1,4)*Power(x2,6) - \
630*Power(x1,2)*Power(x2,8) + 105*Power(x2,10)) + (39916800*(-1 + \
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
Power(x2,12)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,12) - 180*Power(x1,10)*Power(x2,2) + \
1365*Power(x1,8)*Power(x2,4) - 2536*Power(x1,6)*Power(x2,6) + \
1365*Power(x1,4)*Power(x2,8) - 180*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) + 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==10) && (k2==2) )
		return (he0*x2*(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12) - 3*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(23*Power(x1,6) + 124*Power(x1,4)*Power(x2,2) + \
135*Power(x1,2)*Power(x2,4) + 210*Power(x2,6)) - \
30*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(17*Power(x1,8) + 46*Power(x1,6)*Power(x2,2) + \
234*Power(x1,4)*Power(x2,4) - 42*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) + \
135*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,10) + 131*Power(x1,8)*Power(x2,2) - \
376*Power(x1,6)*Power(x2,4) + 658*Power(x1,4)*Power(x2,6) - \
210*Power(x1,2)*Power(x2,8) + 35*Power(x2,10)) + 479001600*(-1 + \
1/he0)*Power(lambdasq,12)*(13*Power(x1,12) - \
286*Power(x1,10)*Power(x2,2) + 1287*Power(x1,8)*Power(x2,4) - \
1716*Power(x1,6)*Power(x2,6) + 715*Power(x1,4)*Power(x2,8) - \
78*Power(x1,2)*Power(x2,10) + Power(x2,12)) - \
239500800*Power(lambdasq,11)*(Power(x1,2) + \
Power(x2,2))*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 59875200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),2)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 9979200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),3)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 1247400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),4)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 124740*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),5)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),6)*(13*Power(x1,12) - 286*Power(x1,10)*Power(x2,2) + \
1287*Power(x1,8)*Power(x2,4) - 1716*Power(x1,6)*Power(x2,6) + \
715*Power(x1,4)*Power(x2,8) - 78*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(73*Power(x1,12) - 1564*Power(x1,10)*Power(x2,2) + \
7101*Power(x1,8)*Power(x2,4) - 9408*Power(x1,6)*Power(x2,6) + \
3955*Power(x1,4)*Power(x2,8) - 420*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==10) && (k2==3) )
		return (he0*(-(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)) + lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(6*Power(x1,4) + 25*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,8) + 44*Power(x1,6)*Power(x2,2) + \
191*Power(x1,4)*Power(x2,4) + 150*Power(x1,2)*Power(x2,6) + \
210*Power(x2,8)) + 3*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(23*Power(x1,10) + 427*Power(x1,8)*Power(x2,2) + \
763*Power(x1,6)*Power(x2,4) + 3675*Power(x1,4)*Power(x2,6) - \
210*Power(x1,2)*Power(x2,8) + 1050*Power(x2,10)) - \
15*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(34*Power(x1,12) + 56*Power(x1,10)*Power(x2,2) + \
3059*Power(x1,8)*Power(x2,4) - 5376*Power(x1,6)*Power(x2,6) + \
7896*Power(x1,4)*Power(x2,8) - 1680*Power(x1,2)*Power(x2,10) + \
315*Power(x2,12)) - 135*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(73*Power(x1,14) - 6496*Power(x1,12)*Power(x2,2) + \
71603*Power(x1,10)*Power(x2,4) - 214662*Power(x1,8)*Power(x2,6) + \
214767*Power(x1,6)*Power(x2,8) - 71540*Power(x1,4)*Power(x2,10) + \
6517*Power(x1,2)*Power(x2,12) - 70*Power(x2,14)) + 6227020800*(-1 + \
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
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,14) + 448*Power(x1,12)*Power(x2,2) - \
4361*Power(x1,10)*Power(x2,4) + 13587*Power(x1,8)*Power(x2,6) - \
13272*Power(x1,6)*Power(x2,8) + 4550*Power(x1,4)*Power(x2,10) - \
385*Power(x1,2)*Power(x2,12) + 7*Power(x2,14)) - \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==10) && (k2==4) )
		return (he0*x2*(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),13)*(10*Power(x1,4) + 27*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) + (87178291200*(-1 + \
he0)*Power(lambdasq,14)*(3*Power(x1,2) - Power(x2,2))*(5*Power(x1,4) \
- 10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)))/he0 + \
43589145600*Power(lambdasq,13)*(3*Power(x1,2) - \
Power(x2,2))*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
10897286400*Power(lambdasq,12)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1816214400*Power(lambdasq,11)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
227026800*Power(lambdasq,10)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
22702680*Power(lambdasq,9)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
1891890*Power(lambdasq,8)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
135135*Power(lambdasq,7)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,8) - \
28*Power(x1,6)*Power(x2,2) + 134*Power(x1,4)*Power(x2,4) - \
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) + \
Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(15*Power(x1,8) + 220*Power(x1,6)*Power(x2,2) + \
843*Power(x1,4)*Power(x2,4) + 540*Power(x1,2)*Power(x2,6) + \
630*Power(x2,8)) - 21*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),11)*(15*Power(x1,10) + 145*Power(x1,8)*Power(x2,2) + \
183*Power(x1,6)*Power(x2,4) + 765*Power(x1,4)*Power(x2,6) + \
30*Power(x1,2)*Power(x2,8) + 150*Power(x2,10)) + \
21*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(120*Power(x1,12) + 110*Power(x1,10)*Power(x2,2) + \
4899*Power(x1,8)*Power(x2,4) - 5130*Power(x1,6)*Power(x2,6) + \
7140*Power(x1,4)*Power(x2,8) - 900*Power(x1,2)*Power(x2,10) + \
225*Power(x2,12)) + 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(30*Power(x1,14) - 2260*Power(x1,12)*Power(x2,2) + \
12837*Power(x1,10)*Power(x2,4) - 29565*Power(x1,8)*Power(x2,6) + \
21720*Power(x1,6)*Power(x2,8) - 6390*Power(x1,4)*Power(x2,10) + \
405*Power(x1,2)*Power(x2,12) - 9*Power(x2,14)) + \
945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(135*Power(x1,14) - 4060*Power(x1,12)*Power(x2,2) + \
26859*Power(x1,10)*Power(x2,4) - 57480*Power(x1,8)*Power(x2,6) + \
44765*Power(x1,6)*Power(x2,8) - 12180*Power(x1,4)*Power(x2,10) + \
945*Power(x1,2)*Power(x2,12) - \
8*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==10) && (k2==5) )
		return (he0*(-(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)) + \
15*lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,8) + 23*Power(x1,6)*Power(x2,2) + \
79*Power(x1,4)*Power(x2,4) + 45*Power(x1,2)*Power(x2,6) + \
42*Power(x2,8)) + 15*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(Power(x1,12) + 60*Power(x1,10)*Power(x2,2) + \
440*Power(x1,8)*Power(x2,4) + 452*Power(x1,6)*Power(x2,6) + \
1485*Power(x1,4)*Power(x2,8) + 168*Power(x1,2)*Power(x2,10) + \
210*Power(x2,12)) - \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,14) + 30*Power(x1,12)*Power(x2,2) + \
20*Power(x1,10)*Power(x2,4) + 647*Power(x1,8)*Power(x2,6) - \
375*Power(x1,6)*Power(x2,8) + 588*Power(x1,4)*Power(x2,10) - \
30*Power(x1,2)*Power(x2,12) + 15*Power(x2,14)) + (1307674368000*(-1 + \
he0)*Power(lambdasq,15)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) \
+ 1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)))/he0 + \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(2*Power(x1,16) - 555*Power(x1,14)*Power(x2,2) + \
8071*Power(x1,12)*Power(x2,4) - 35903*Power(x1,10)*Power(x2,6) + \
57375*Power(x1,8)*Power(x2,8) - 35889*Power(x1,6)*Power(x2,10) + \
8085*Power(x1,4)*Power(x2,12) - 549*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(8*Power(x1,16) - 60*Power(x1,14)*Power(x2,2) + \
2335*Power(x1,12)*Power(x2,4) - 8306*Power(x1,10)*Power(x2,6) + \
15165*Power(x1,8)*Power(x2,8) - 8376*Power(x1,6)*Power(x2,10) + \
2265*Power(x1,4)*Power(x2,12) - 90*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 14175*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(9*Power(x1,16) - 1072*Power(x1,14)*Power(x2,2) + \
16268*Power(x1,12)*Power(x2,4) - 71568*Power(x1,10)*Power(x2,6) + \
115030*Power(x1,8)*Power(x2,8) - 71568*Power(x1,6)*Power(x2,10) + \
16268*Power(x1,4)*Power(x2,12) - 1072*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) + 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==10) && (k2==6) )
		return -(he0*x2*(-(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),15)*(21*Power(x1,4) + 34*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,8) + 35*Power(x1,6)*Power(x2,2) + \
107*Power(x1,4)*Power(x2,4) + 57*Power(x1,2)*Power(x2,6) + \
42*Power(x2,8)) + 105*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(Power(x1,12) + 20*Power(x1,10)*Power(x2,2) + \
124*Power(x1,8)*Power(x2,4) + 116*Power(x1,6)*Power(x2,6) + \
285*Power(x1,4)*Power(x2,8) + 48*Power(x1,2)*Power(x2,10) + \
30*Power(x2,12)) - \
105*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(19*Power(x1,14) + 290*Power(x1,12)*Power(x2,2) + \
196*Power(x1,10)*Power(x2,4) + 3509*Power(x1,8)*Power(x2,6) - \
885*Power(x1,6)*Power(x2,8) + 2172*Power(x1,4)*Power(x2,10) + \
30*Power(x1,2)*Power(x2,12) + 45*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,16)*(17*Power(x1,16) - \
680*Power(x1,14)*Power(x2,2) + 6188*Power(x1,12)*Power(x2,4) - \
19448*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
12376*Power(x1,6)*Power(x2,10) + 2380*Power(x1,4)*Power(x2,12) - \
136*Power(x1,2)*Power(x2,14) + Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),2)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),3)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(226*Power(x1,16) - 12715*Power(x1,14)*Power(x2,2) + \
108871*Power(x1,12)*Power(x2,4) - 350671*Power(x1,10)*Power(x2,6) + \
431295*Power(x1,8)*Power(x2,8) - 223377*Power(x1,6)*Power(x2,10) + \
41685*Power(x1,4)*Power(x2,12) - 2613*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(58*Power(x1,16) - 220*Power(x1,14)*Power(x2,2) + \
6097*Power(x1,12)*Power(x2,4) - 13894*Power(x1,10)*Power(x2,6) + \
21795*Power(x1,8)*Power(x2,8) - 8736*Power(x1,6)*Power(x2,10) + \
2415*Power(x1,4)*Power(x2,12) - 30*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 1575*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(1219*Power(x1,16) - 48592*Power(x1,14)*Power(x2,2) + \
442540*Power(x1,12)*Power(x2,4) - 1390336*Power(x1,10)*Power(x2,6) + \
1738410*Power(x1,8)*Power(x2,8) - 884688*Power(x1,6)*Power(x2,10) + \
170268*Power(x1,4)*Power(x2,12) - 9696*Power(x1,2)*Power(x2,14) + \
75*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==10) && (k2==7) )
		return (he0*(-(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17)) + lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),16)*(28*Power(x1,4) + 39*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - \
2*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(105*Power(x1,8) + 392*Power(x1,6)*Power(x2,2) + \
1056*Power(x1,4)*Power(x2,4) + 540*Power(x1,2)*Power(x2,6) + \
315*Power(x2,8)) + \
30*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(14*Power(x1,12) + 147*Power(x1,10)*Power(x2,2) + \
791*Power(x1,8)*Power(x2,4) + 716*Power(x1,6)*Power(x2,6) + \
1314*Power(x1,4)*Power(x2,8) + 273*Power(x1,2)*Power(x2,10) + \
105*Power(x2,12)) - \
105*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,16) + 72*Power(x1,14)*Power(x2,2) + \
810*Power(x1,12)*Power(x2,4) + 636*Power(x1,10)*Power(x2,6) + \
5958*Power(x1,8)*Power(x2,8) - 72*Power(x1,6)*Power(x2,10) + \
2706*Power(x1,4)*Power(x2,12) + 180*Power(x1,2)*Power(x2,14) + \
45*Power(x2,16)) - 1575*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(1219*Power(x1,18) - 185940*Power(x1,16)*Power(x2,2) + \
3719556*Power(x1,14)*Power(x2,4) - 22564248*Power(x1,12)*Power(x2,6) \
+ 53188290*Power(x1,10)*Power(x2,8) - \
53187408*Power(x1,8)*Power(x2,10) + 22564836*Power(x1,6)*Power(x2,12) \
- 3719304*Power(x1,4)*Power(x2,14) + 186003*Power(x1,2)*Power(x2,16) \
- 1212*Power(x2,18)) - 630*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(113*Power(x1,18) - 23589*Power(x1,16)*Power(x2,2) + \
463590*Power(x1,14)*Power(x2,4) - 2823618*Power(x1,12)*Power(x2,6) + \
6644016*Power(x1,10)*Power(x2,8) - 6652836*Power(x1,8)*Power(x2,10) + \
2817738*Power(x1,6)*Power(x2,12) - 466110*Power(x1,4)*Power(x2,14) + \
22959*Power(x1,2)*Power(x2,16) - 183*Power(x2,18)) + \
355687428096000*(-1 + 1/he0)*Power(lambdasq,17)*(Power(x1,18) - \
153*Power(x1,16)*Power(x2,2) + 3060*Power(x1,14)*Power(x2,4) - \
18564*Power(x1,12)*Power(x2,6) + 43758*Power(x1,10)*Power(x2,8) - \
43758*Power(x1,8)*Power(x2,10) + 18564*Power(x1,6)*Power(x2,12) - \
3060*Power(x1,4)*Power(x2,14) + 153*Power(x1,2)*Power(x2,16) - \
Power(x2,18)) - 44460928512000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
7410154752000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
926269344000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
92626934400*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
7718911200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
551350800*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(29*Power(x1,18) - 762*Power(x1,16)*Power(x2,2) + \
22590*Power(x1,14)*Power(x2,4) - 124404*Power(x1,12)*Power(x2,6) + \
308148*Power(x1,10)*Power(x2,8) - 296388*Power(x1,8)*Power(x2,10) + \
131754*Power(x1,6)*Power(x2,12) - 19860*Power(x1,4)*Power(x2,14) + \
1287*Power(x1,2)*Power(x2,16) + 6*Power(x2,18)) + \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(19*Power(x1,18) + 873*Power(x1,16)*Power(x2,2) - \
1710*Power(x1,14)*Power(x2,4) + 41874*Power(x1,12)*Power(x2,6) - \
60174*Power(x1,10)*Power(x2,8) + 90162*Power(x1,8)*Power(x2,10) - \
23646*Power(x1,6)*Power(x2,12) + 8010*Power(x1,4)*Power(x2,14) + \
135*Power(x1,2)*Power(x2,16) + 9*Power(x2,18)) - \
177843714048000*Power(lambdasq,16)*(Power(x1,20) - \
152*Power(x1,18)*Power(x2,2) + 2907*Power(x1,16)*Power(x2,4) - \
15504*Power(x1,14)*Power(x2,6) + 25194*Power(x1,12)*Power(x2,8) - \
25194*Power(x1,8)*Power(x2,12) + 15504*Power(x1,6)*Power(x2,14) - \
2907*Power(x1,4)*Power(x2,16) + 152*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==10) && (k2==8) )
		return (he0*x2*(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),18) - \
9*lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17)*(4*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + \
18*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(21*Power(x1,8) + 64*Power(x1,6)*Power(x2,2) + \
151*Power(x1,4)*Power(x2,4) + 75*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - \
18*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(70*Power(x1,12) + 483*Power(x1,10)*Power(x2,2) + \
2257*Power(x1,8)*Power(x2,4) + 2028*Power(x1,6)*Power(x2,6) + \
2850*Power(x1,4)*Power(x2,8) + 665*Power(x1,2)*Power(x2,10) + \
175*Power(x2,12)) + \
135*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,16) + 168*Power(x1,14)*Power(x2,2) + \
1554*Power(x1,12)*Power(x2,4) + 1416*Power(x1,10)*Power(x2,6) + \
7554*Power(x1,8)*Power(x2,8) + 1256*Power(x1,6)*Power(x2,10) + \
2674*Power(x1,4)*Power(x2,12) + 280*Power(x1,2)*Power(x2,14) + \
35*Power(x2,16)) + 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(11563*Power(x1,18) - 588768*Power(x1,16)*Power(x2,2) \
+ 7067484*Power(x1,14)*Power(x2,4) - \
30621984*Power(x1,12)*Power(x2,6) + 56144714*Power(x1,10)*Power(x2,8) \
- 45932976*Power(x1,8)*Power(x2,10) + \
16490796*Power(x1,6)*Power(x2,12) - 2355072*Power(x1,4)*Power(x2,14) \
+ 104067*Power(x1,2)*Power(x2,16) - 592*Power(x2,18)) + \
5670*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(227*Power(x1,18) - 13677*Power(x1,16)*Power(x2,2) + \
159546*Power(x1,14)*Power(x2,4) - 698226*Power(x1,12)*Power(x2,6) + \
1272976*Power(x1,10)*Power(x2,8) - 1046604*Power(x1,8)*Power(x2,10) + \
373254*Power(x1,6)*Power(x2,12) - 54078*Power(x1,4)*Power(x2,14) + \
2253*Power(x1,2)*Power(x2,16) - 23*Power(x2,18)) + \
(6402373705728000*(-1 + he0)*Power(lambdasq,18)*(19*Power(x1,18) - \
969*Power(x1,16)*Power(x2,2) + 11628*Power(x1,14)*Power(x2,4) - \
50388*Power(x1,12)*Power(x2,6) + 92378*Power(x1,10)*Power(x2,8) - \
75582*Power(x1,8)*Power(x2,10) + 27132*Power(x1,6)*Power(x2,12) - \
3876*Power(x1,4)*Power(x2,14) + 171*Power(x1,2)*Power(x2,16) - \
Power(x2,18)))/he0 + 3201186852864000*Power(lambdasq,17)*(Power(x1,2) \
+ Power(x2,2))*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(17*Power(x1,18) + 393*Power(x1,16)*Power(x2,2) - \
306*Power(x1,14)*Power(x2,4) + 9666*Power(x1,12)*Power(x2,6) - \
8090*Power(x1,10)*Power(x2,8) + 13626*Power(x1,8)*Power(x2,10) - \
1746*Power(x1,6)*Power(x2,12) + 1050*Power(x1,4)*Power(x2,14) + \
45*Power(x1,2)*Power(x2,16) + Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(94*Power(x1,18) - 1119*Power(x1,16)*Power(x2,2) + \
22248*Power(x1,14)*Power(x2,4) - 82170*Power(x1,12)*Power(x2,6) + \
166178*Power(x1,10)*Power(x2,8) - 124578*Power(x1,8)*Power(x2,10) + \
50148*Power(x1,6)*Power(x2,12) - 5610*Power(x1,4)*Power(x2,14) + \
468*Power(x1,2)*Power(x2,16) + \
5*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==10) && (k2==9) )
		return (he0*(-(Power(x1,10)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),19)) + lambdasq*Power(x1,8)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),18)*(45*Power(x1,4) + 52*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - \
9*Power(lambdasq,2)*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17)*(70*Power(x1,8) + 185*Power(x1,6)*Power(x2,2) + \
382*Power(x1,4)*Power(x2,4) + 185*Power(x1,2)*Power(x2,6) + \
70*Power(x2,8)) + \
18*Power(lambdasq,3)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(175*Power(x1,12) + 910*Power(x1,10)*Power(x2,2) + \
3680*Power(x1,8)*Power(x2,4) + 3306*Power(x1,6)*Power(x2,6) + \
3680*Power(x1,4)*Power(x2,8) + 910*Power(x1,2)*Power(x2,10) + \
175*Power(x2,12)) - \
9*Power(lambdasq,4)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(525*Power(x1,16) + 6650*Power(x1,14)*Power(x2,2) + \
52010*Power(x1,12)*Power(x2,4) + 53230*Power(x1,10)*Power(x2,6) + \
180066*Power(x1,8)*Power(x2,8) + 53230*Power(x1,6)*Power(x2,10) + \
52010*Power(x1,4)*Power(x2,12) + 6650*Power(x1,2)*Power(x2,14) + \
525*Power(x2,16)) + (121645100408832000*(-1 + \
he0)*Power(lambdasq,19)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) \
+ 4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)))/he0 + \
15205637551104000*Power(lambdasq,17)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 2534272925184000*Power(lambdasq,16)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 316784115648000*Power(lambdasq,15)*Power(Power(x1,2) \
+ Power(x2,2),4)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 31678411564800*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 2639867630400*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 188561973600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 11785123350*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 654729075*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) + 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,20) + 560*Power(x1,18)*Power(x2,2) + \
9345*Power(x1,16)*Power(x2,4) + 131400*Power(x1,12)*Power(x2,8) - \
50368*Power(x1,10)*Power(x2,10) + 131400*Power(x1,8)*Power(x2,12) + \
9345*Power(x1,4)*Power(x2,16) + 560*Power(x1,2)*Power(x2,18) + \
7*Power(x2,20)) - 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(17*Power(x1,20) + 1180*Power(x1,18)*Power(x2,2) - \
7305*Power(x1,16)*Power(x2,4) + 114300*Power(x1,14)*Power(x2,6) - \
285480*Power(x1,12)*Power(x2,8) + 506512*Power(x1,10)*Power(x2,10) - \
285480*Power(x1,8)*Power(x2,12) + 114300*Power(x1,6)*Power(x2,14) - \
7305*Power(x1,4)*Power(x2,16) + 1180*Power(x1,2)*Power(x2,18) + \
17*Power(x2,20)) + 3780*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),12)*(47*Power(x1,20) - 2630*Power(x1,18)*Power(x2,2) + \
82290*Power(x1,16)*Power(x2,4) - 627030*Power(x1,14)*Power(x2,6) + \
2082735*Power(x1,12)*Power(x2,8) - 3010088*Power(x1,10)*Power(x2,10) \
+ 2082735*Power(x1,8)*Power(x2,12) - 627030*Power(x1,6)*Power(x2,14) \
+ 82290*Power(x1,4)*Power(x2,16) - 2630*Power(x1,2)*Power(x2,18) + \
47*Power(x2,20)) + 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(454*Power(x1,20) - 100435*Power(x1,18)*Power(x2,2) + \
2539830*Power(x1,16)*Power(x2,4) - 20352660*Power(x1,14)*Power(x2,6) \
+ 66103620*Power(x1,12)*Power(x2,8) - \
96991666*Power(x1,10)*Power(x2,10) + \
66103620*Power(x1,8)*Power(x2,12) - 20352660*Power(x1,6)*Power(x2,14) \
+ 2539830*Power(x1,4)*Power(x2,16) - 100435*Power(x1,2)*Power(x2,18) \
+ 454*Power(x2,20)) + 2835*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),10)*(11563*Power(x1,20) - \
2193820*Power(x1,18)*Power(x2,2) + 55947135*Power(x1,16)*Power(x2,4) \
- 447569520*Power(x1,14)*Power(x2,6) + \
1454610390*Power(x1,12)*Power(x2,8) - \
2133419752*Power(x1,10)*Power(x2,10) + \
1454610390*Power(x1,8)*Power(x2,12) - \
447569520*Power(x1,6)*Power(x2,14) + \
55947135*Power(x1,4)*Power(x2,16) - 2193820*Power(x1,2)*Power(x2,18) \
+ 11563*Power(x2,20)) + \
60822550204416000*Power(lambdasq,18)*(Power(x1,22) - \
189*Power(x1,20)*Power(x2,2) + 4655*Power(x1,18)*Power(x2,4) - \
33915*Power(x1,16)*Power(x2,6) + 87210*Power(x1,14)*Power(x2,8) - \
58786*Power(x1,12)*Power(x2,10) - 58786*Power(x1,10)*Power(x2,12) + \
87210*Power(x1,8)*Power(x2,14) - 33915*Power(x1,6)*Power(x2,16) + \
4655*Power(x1,4)*Power(x2,18) - 189*Power(x1,2)*Power(x2,20) + \
Power(x2,22))))/(2.*Power(lambdasq,19)*Pi*Power(Power(x1,2) + \
Power(x2,2),20));
	else if ( (k1==10) && (k2==10) )
		return -(he0*x2*(-(Power(x1,10)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),20)) + \
5*lambdasq*Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),19)*(11*Power(x1,4) + 12*Power(x1,2)*Power(x2,2) + \
9*Power(x2,4)) - \
5*Power(lambdasq,2)*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),18)*(198*Power(x1,8) + 473*Power(x1,6)*Power(x2,2) + \
858*Power(x1,4)*Power(x2,4) + 405*Power(x1,2)*Power(x2,6) + \
126*Power(x2,8)) + (2432902008176640000*(-1 + \
he0)*Power(lambdasq,20)*(3*Power(x1,2) - Power(x2,2))*(7*Power(x1,6) \
- 35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)))/he0 + \
1216451004088320000*Power(lambdasq,19)*(3*Power(x1,2) - \
Power(x2,2))*(Power(x1,2) + Power(x2,2))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 304112751022080000*Power(lambdasq,18)*(3*Power(x1,2) \
- Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 50685458503680000*Power(lambdasq,17)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 6335682312960000*Power(lambdasq,16)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 633568231296000*Power(lambdasq,15)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 52797352608000*Power(lambdasq,14)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 3771239472000*Power(lambdasq,13)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 235702467000*Power(lambdasq,12)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),8)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 13094581500*Power(lambdasq,11)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),9)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 654729075*Power(lambdasq,10)*(3*Power(x1,2) - \
Power(x2,2))*Power(Power(x1,2) + Power(x2,2),10)*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6))*(Power(x1,12) - 58*Power(x1,10)*Power(x2,2) + \
655*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
1423*Power(x1,4)*Power(x2,8) - 186*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + \
90*Power(lambdasq,3)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),17)*(77*Power(x1,12) + 330*Power(x1,10)*Power(x2,2) + \
1155*Power(x1,8)*Power(x2,4) + 1036*Power(x1,6)*Power(x2,6) + \
945*Power(x1,4)*Power(x2,8) + 238*Power(x1,2)*Power(x2,10) + \
35*Power(x2,12)) - \
45*Power(lambdasq,4)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),16)*(385*Power(x1,16) + 3234*Power(x1,14)*Power(x2,2) + \
21450*Power(x1,12)*Power(x2,4) + 23870*Power(x1,10)*Power(x2,6) + \
56322*Power(x1,8)*Power(x2,8) + 21630*Power(x1,6)*Power(x2,10) + \
13706*Power(x1,4)*Power(x2,12) + 1890*Power(x1,2)*Power(x2,14) + \
105*Power(x2,16)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),15)*(231*Power(x1,20) + 6160*Power(x1,18)*Power(x2,2) + \
82467*Power(x1,16)*Power(x2,4) + 43560*Power(x1,14)*Power(x2,6) + \
733040*Power(x1,12)*Power(x2,8) - 23352*Power(x1,10)*Power(x2,10) + \
561120*Power(x1,8)*Power(x2,12) + 62888*Power(x1,6)*Power(x2,14) + \
37905*Power(x1,4)*Power(x2,16) + 2520*Power(x1,2)*Power(x2,18) + \
21*Power(x2,20)) + 9450*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),13)*(231*Power(x1,20) - 5390*Power(x1,18)*Power(x2,2) + \
108801*Power(x1,16)*Power(x2,4) - 570636*Power(x1,14)*Power(x2,6) + \
1510894*Power(x1,12)*Power(x2,8) - 1749216*Power(x1,10)*Power(x2,10) \
+ 1050210*Power(x1,8)*Power(x2,12) - 262500*Power(x1,6)*Power(x2,14) \
+ 33579*Power(x1,4)*Power(x2,16) - 546*Power(x1,2)*Power(x2,18) + \
29*Power(x2,20)) - 675*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(231*Power(x1,20) + 8008*Power(x1,18)*Power(x2,2) - \
25179*Power(x1,16)*Power(x2,4) + 370260*Power(x1,14)*Power(x2,6) - \
605528*Power(x1,12)*Power(x2,8) + 1044624*Power(x1,10)*Power(x2,10) - \
399672*Power(x1,8)*Power(x2,12) + 189924*Power(x1,6)*Power(x2,14) - \
735*Power(x1,4)*Power(x2,16) + 2352*Power(x1,2)*Power(x2,18) + \
35*Power(x2,20)) + 4725*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),12)*(5082*Power(x1,20) - 353045*Power(x1,18)*Power(x2,2) \
+ 5323626*Power(x1,16)*Power(x2,4) - \
30557340*Power(x1,14)*Power(x2,6) + 77072380*Power(x1,12)*Power(x2,8) \
- 92637678*Power(x1,10)*Power(x2,10) + \
53349660*Power(x1,8)*Power(x2,12) - 14268156*Power(x1,6)*Power(x2,14) \
+ 1561770*Power(x1,4)*Power(x2,16) - 56805*Power(x1,2)*Power(x2,18) + \
122*Power(x2,20)) + 14175*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),11)*(44121*Power(x1,20) - \
2792020*Power(x1,18)*Power(x2,2) + 42724143*Power(x1,16)*Power(x2,4) \
- 244126080*Power(x1,14)*Power(x2,6) + \
617112650*Power(x1,12)*Power(x2,8) - \
740519304*Power(x1,10)*Power(x2,10) + \
427233870*Power(x1,8)*Power(x2,12) - \
113923488*Power(x1,6)*Power(x2,14) + \
12566925*Power(x1,4)*Power(x2,16) - 440580*Power(x1,2)*Power(x2,18) + \
2131*Power(x2,20))))/(2.*Power(lambdasq,20)*Pi*Power(Power(x1,2) + \
Power(x2,2),21));
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
	else if ( (k1==0) && (k2==8) )
		return -(he0*x1*(Power(x2,8)*Power(Power(x1,2) + Power(x2,2),8) - \
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
	else if ( (k1==0) && (k2==9) )
		return (he0*x1*x2*(Power(x2,8)*Power(Power(x1,2) + Power(x2,2),9) - \
18*lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,2) + Power(x2,2)) + (725760*(-1 + \
he0)*Power(lambdasq,9)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 + 362880*Power(lambdasq,8)*(Power(x1,2) + \
Power(x2,2))*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 90720*Power(lambdasq,7)*Power(Power(x1,2) + \
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
5*Power(x2,4)) - \
252*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,4) - 6*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==0) && (k2==10) )
		return (he0*x1*(-(Power(x2,10)*Power(Power(x1,2) + Power(x2,2),10)) \
+ 5*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,2) + 5*Power(x2,2)) - \
90*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
2*Power(x2,4)) + 90*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) - 7*Power(x1,4)*Power(x2,2) + \
25*Power(x1,2)*Power(x2,4) + 3*Power(x2,6)) - \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 80*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + (3628800*(-1 + he0)*Power(lambdasq,10)*(Power(x1,10) \
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
165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
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
	else if ( (k1==1) && (k2==8) )
		return (he0*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),9) - lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 6*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,6) + 15*Power(x1,2)*Power(x2,4) + \
2*Power(x2,6)) - 42*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(10*Power(x1,8) - 35*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 5*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
+ 21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*(5*Power(x1,10) - 200*Power(x1,8)*Power(x2,2) + \
950*Power(x1,6)*Power(x2,4) - 940*Power(x1,4)*Power(x2,6) + \
205*Power(x1,2)*Power(x2,8) - 4*Power(x2,10)) + (362880*(-1 + \
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
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==1) && (k2==9) )
		return (he0*x2*(-(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10)) + lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(36*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 18*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(21*Power(x1,6) + 8*Power(x1,4)*Power(x2,2) + \
8*Power(x1,2)*Power(x2,4) + Power(x2,6)) + \
90*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(14*Power(x1,8) - 21*Power(x1,6)*Power(x2,2) + \
29*Power(x1,4)*Power(x2,4) + Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 80*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 3628800*(-1 + \
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
	else if ( (k1==1) && (k2==10) )
		return (he0*(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(45*Power(x1,4) + 24*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 5*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(126*Power(x1,6) + 81*Power(x1,4)*Power(x2,2) + \
48*Power(x1,2)*Power(x2,4) + 5*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(35*Power(x1,8) - 14*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) + 8*Power(x1,2)*Power(x2,6) + \
2*Power(x2,8)) + 45*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(105*Power(x1,10) - 630*Power(x1,8)*Power(x2,2) + \
1358*Power(x1,6)*Power(x2,4) - 600*Power(x1,4)*Power(x2,6) + \
129*Power(x1,2)*Power(x2,8) + 6*Power(x2,10)) + 39916800*(-1 + \
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
Power(x2,12)) - 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,12) - 180*Power(x1,10)*Power(x2,2) + \
1365*Power(x1,8)*Power(x2,4) - 2536*Power(x1,6)*Power(x2,6) + \
1365*Power(x1,4)*Power(x2,8) - 180*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) - 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) - 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
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
	else if ( (k1==2) && (k2==8) )
		return -(he0*x1*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
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
	else if ( (k1==2) && (k2==9) )
		return (he0*x1*x2*(Power(x1,2)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(36*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + 2*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(189*Power(x1,6) + 72*Power(x1,4)*Power(x2,2) + \
127*Power(x1,2)*Power(x2,4) + 24*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(14*Power(x1,8) - 21*Power(x1,6)*Power(x2,2) + \
51*Power(x1,4)*Power(x2,4) + Power(x1,2)*Power(x2,6) + 3*Power(x2,8)) \
+ 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(21*Power(x1,10) - 280*Power(x1,8)*Power(x2,2) + \
1134*Power(x1,6)*Power(x2,4) - 1044*Power(x1,4)*Power(x2,6) + \
325*Power(x1,2)*Power(x2,8) - 12*Power(x2,10)) + (159667200*(-1 + \
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
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==2) && (k2==10) )
		return (he0*x1*(-(Power(x1,2)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),12)) + 3*lambdasq*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) - 3*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(210*Power(x1,6) + 135*Power(x1,4)*Power(x2,2) + \
124*Power(x1,2)*Power(x2,4) + 23*Power(x2,6)) + \
30*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(105*Power(x1,8) - 42*Power(x1,6)*Power(x2,2) + \
234*Power(x1,4)*Power(x2,4) + 46*Power(x1,2)*Power(x2,6) + \
17*Power(x2,8)) - 135*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),8)*(35*Power(x1,10) - 210*Power(x1,8)*Power(x2,2) + \
658*Power(x1,6)*Power(x2,4) - 376*Power(x1,4)*Power(x2,6) + \
131*Power(x1,2)*Power(x2,8) + 2*Power(x2,10)) + (479001600*(-1 + \
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
13*Power(x2,12)) + 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(7*Power(x1,12) - 420*Power(x1,10)*Power(x2,2) + \
3955*Power(x1,8)*Power(x2,4) - 9408*Power(x1,6)*Power(x2,6) + \
7101*Power(x1,4)*Power(x2,8) - 1564*Power(x1,2)*Power(x2,10) + \
73*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
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
	else if ( (k1==3) && (k2==8) )
		return (he0*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),11) - \
2*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(14*Power(x1,4) + 6*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(210*Power(x1,8) + 28*Power(x1,6)*Power(x2,2) + \
309*Power(x1,4)*Power(x2,4) + 54*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) - 30*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(14*Power(x1,10) - 42*Power(x1,8)*Power(x2,2) + \
154*Power(x1,6)*Power(x2,4) - 37*Power(x1,4)*Power(x2,6) + \
18*Power(x1,2)*Power(x2,8) + Power(x2,10)) + (39916800*(-1 + \
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
Power(x2,2),7)*(7*Power(x1,12) - 252*Power(x1,10)*Power(x2,2) + \
2100*Power(x1,8)*Power(x2,4) - 3724*Power(x1,6)*Power(x2,6) + \
2115*Power(x1,4)*Power(x2,8) - 240*Power(x1,2)*Power(x2,10) + \
10*Power(x2,12)) + 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) + 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==3) && (k2==9) )
		return -(he0*x2*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),12) - \
6*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(6*Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(126*Power(x1,8) + 60*Power(x1,6)*Power(x2,2) + \
139*Power(x1,4)*Power(x2,4) + 30*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 12*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(105*Power(x1,10) - 126*Power(x1,8)*Power(x2,2) + \
591*Power(x1,6)*Power(x2,4) + 4*Power(x1,4)*Power(x2,6) + \
66*Power(x1,2)*Power(x2,8) + 4*Power(x2,10)) + (479001600*(-1 + \
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
Power(x2,12)) + 270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(35*Power(x1,12) - 791*Power(x1,10)*Power(x2,2) + \
3528*Power(x1,8)*Power(x2,4) - 4734*Power(x1,6)*Power(x2,6) + \
1955*Power(x1,4)*Power(x2,8) - 219*Power(x1,2)*Power(x2,10) + \
2*Power(x2,12)) + 135*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 84*Power(x1,10)*Power(x2,2) + \
476*Power(x1,8)*Power(x2,4) - 548*Power(x1,6)*Power(x2,6) + \
275*Power(x1,4)*Power(x2,8) - 16*Power(x1,2)*Power(x2,10) + \
2*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==3) && (k2==10) )
		return (he0*(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),12)*(45*Power(x1,4) + 25*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(210*Power(x1,8) + 150*Power(x1,6)*Power(x2,2) + \
191*Power(x1,4)*Power(x2,4) + 44*Power(x1,2)*Power(x2,6) + \
Power(x2,8)) - 3*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(1050*Power(x1,10) - 210*Power(x1,8)*Power(x2,2) + \
3675*Power(x1,6)*Power(x2,4) + 763*Power(x1,4)*Power(x2,6) + \
427*Power(x1,2)*Power(x2,8) + 23*Power(x2,10)) + \
15*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(315*Power(x1,12) - 1680*Power(x1,10)*Power(x2,2) + \
7896*Power(x1,8)*Power(x2,4) - 5376*Power(x1,6)*Power(x2,6) + \
3059*Power(x1,4)*Power(x2,8) + 56*Power(x1,2)*Power(x2,10) + \
34*Power(x2,12)) - 135*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(70*Power(x1,14) - 6517*Power(x1,12)*Power(x2,2) + \
71540*Power(x1,10)*Power(x2,4) - 214767*Power(x1,8)*Power(x2,6) + \
214662*Power(x1,6)*Power(x2,8) - 71603*Power(x1,4)*Power(x2,10) + \
6496*Power(x1,2)*Power(x2,12) - 73*Power(x2,14)) + 6227020800*(-1 + \
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
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,14) - 385*Power(x1,12)*Power(x2,2) + \
4550*Power(x1,10)*Power(x2,4) - 13272*Power(x1,8)*Power(x2,6) + \
13587*Power(x1,6)*Power(x2,8) - 4361*Power(x1,4)*Power(x2,10) + \
448*Power(x1,2)*Power(x2,12) + 2*Power(x2,14)) - \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
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
	else if ( (k1==4) && (k2==8) )
		return -(he0*x1*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
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
	else if ( (k1==4) && (k2==9) )
		return (he0*x1*x2*(Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),13) - \
2*lambdasq*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(18*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(126*Power(x1,8) + 84*Power(x1,6)*Power(x2,2) + \
211*Power(x1,4)*Power(x2,4) + 50*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 6*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(210*Power(x1,10) - 126*Power(x1,8)*Power(x2,2) + \
1686*Power(x1,6)*Power(x2,4) + 49*Power(x1,4)*Power(x2,6) + \
350*Power(x1,2)*Power(x2,8) + 35*Power(x2,10)) + \
270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,12) - 833*Power(x1,10)*Power(x2,2) + \
4424*Power(x1,8)*Power(x2,4) - 7734*Power(x1,6)*Power(x2,6) + \
4424*Power(x1,4)*Power(x2,8) - 833*Power(x1,2)*Power(x2,10) + \
28*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(63*Power(x1,12) - 588*Power(x1,10)*Power(x2,2) + \
5124*Power(x1,8)*Power(x2,4) - 6684*Power(x1,6)*Power(x2,6) + \
5299*Power(x1,4)*Power(x2,8) - 448*Power(x1,2)*Power(x2,10) + \
98*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==4) && (k2==10) )
		return (he0*x1*(-(Power(x1,4)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),14)) + lambdasq*Power(x1,2)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),13)*(45*Power(x1,4) + 27*Power(x1,2)*Power(x2,2) + \
10*Power(x2,4)) + (87178291200*(-1 + \
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
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
Power(lambdasq,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(630*Power(x1,8) + 540*Power(x1,6)*Power(x2,2) + \
843*Power(x1,4)*Power(x2,4) + 220*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 21*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),11)*(150*Power(x1,10) + 30*Power(x1,8)*Power(x2,2) + \
765*Power(x1,6)*Power(x2,4) + 183*Power(x1,4)*Power(x2,6) + \
145*Power(x1,2)*Power(x2,8) + 15*Power(x2,10)) - \
21*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(225*Power(x1,12) - 900*Power(x1,10)*Power(x2,2) + \
7140*Power(x1,8)*Power(x2,4) - 5130*Power(x1,6)*Power(x2,6) + \
4899*Power(x1,4)*Power(x2,8) + 110*Power(x1,2)*Power(x2,10) + \
120*Power(x2,12)) + 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(8*Power(x1,14) - 945*Power(x1,12)*Power(x2,2) + \
12180*Power(x1,10)*Power(x2,4) - 44765*Power(x1,8)*Power(x2,6) + \
57480*Power(x1,6)*Power(x2,8) - 26859*Power(x1,4)*Power(x2,10) + \
4060*Power(x1,2)*Power(x2,12) - 135*Power(x2,14)) + \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(9*Power(x1,14) - 405*Power(x1,12)*Power(x2,2) + \
6390*Power(x1,10)*Power(x2,4) - 21720*Power(x1,8)*Power(x2,6) + \
29565*Power(x1,6)*Power(x2,8) - 12837*Power(x1,4)*Power(x2,10) + \
2260*Power(x1,2)*Power(x2,12) - \
30*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
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
	else if ( (k1==5) && (k2==8) )
		return (he0*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),12)*(28*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) + \
3*Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(70*Power(x1,8) + 56*Power(x1,6)*Power(x2,2) + \
229*Power(x1,4)*Power(x2,4) + 50*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(140*Power(x1,12) - 70*Power(x1,10)*Power(x2,2) + \
2758*Power(x1,8)*Power(x2,4) - 483*Power(x1,6)*Power(x2,6) + \
1225*Power(x1,4)*Power(x2,8) + 105*Power(x1,2)*Power(x2,10) + \
5*Power(x2,12)) + 45*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(217*Power(x1,14) - 19502*Power(x1,12)*Power(x2,2) + \
214767*Power(x1,10)*Power(x2,4) - 644056*Power(x1,8)*Power(x2,6) + \
644231*Power(x1,6)*Power(x2,8) - 214662*Power(x1,4)*Power(x2,10) + \
19537*Power(x1,2)*Power(x2,12) - 212*Power(x2,14)) + \
45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,14) - 1267*Power(x1,12)*Power(x2,2) + \
13272*Power(x1,10)*Power(x2,4) - 40516*Power(x1,8)*Power(x2,6) + \
39991*Power(x1,6)*Power(x2,8) - 13587*Power(x1,4)*Power(x2,10) + \
1162*Power(x1,2)*Power(x2,12) - 22*Power(x2,14)) + (6227020800*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,14) - 112*Power(x1,12)*Power(x2,2) + \
2632*Power(x1,10)*Power(x2,4) - 5936*Power(x1,8)*Power(x2,6) + \
7511*Power(x1,6)*Power(x2,8) - 1792*Power(x1,4)*Power(x2,10) + \
322*Power(x1,2)*Power(x2,12) + 8*Power(x2,14)) + \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==5) && (k2==9) )
		return (he0*x2*(-(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14)) + lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),13)*(36*Power(x1,4) + 23*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - (87178291200*(-1 + \
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
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(378*Power(x1,8) + 360*Power(x1,6)*Power(x2,2) + \
905*Power(x1,4)*Power(x2,4) + 240*Power(x1,2)*Power(x2,6) + \
45*Power(x2,8)) + 3*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(420*Power(x1,12) + 126*Power(x1,10)*Power(x2,2) + \
4590*Power(x1,8)*Power(x2,4) + 475*Power(x1,6)*Power(x2,6) + \
1605*Power(x1,4)*Power(x2,8) + 195*Power(x1,2)*Power(x2,10) + \
5*Power(x2,12)) - 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(45*Power(x1,14) - 2175*Power(x1,12)*Power(x2,2) + \
13032*Power(x1,10)*Power(x2,4) - 29340*Power(x1,8)*Power(x2,6) + \
21845*Power(x1,6)*Power(x2,8) - 6375*Power(x1,4)*Power(x2,10) + \
390*Power(x1,2)*Power(x2,12) - 14*Power(x2,14)) - \
945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(135*Power(x1,14) - 4060*Power(x1,12)*Power(x2,2) + \
26859*Power(x1,10)*Power(x2,4) - 57480*Power(x1,8)*Power(x2,6) + \
44765*Power(x1,6)*Power(x2,8) - 12180*Power(x1,4)*Power(x2,10) + \
945*Power(x1,2)*Power(x2,12) - 8*Power(x2,14)) - \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(45*Power(x1,14) - 240*Power(x1,12)*Power(x2,2) + \
4284*Power(x1,10)*Power(x2,4) - 5580*Power(x1,8)*Power(x2,6) + \
7115*Power(x1,6)*Power(x2,8) - 750*Power(x1,4)*Power(x2,10) + \
300*Power(x1,2)*Power(x2,12) + \
10*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==5) && (k2==10) )
		return (he0*(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),15) - \
15*lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),14)*(3*Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
Power(x2,4)) + \
15*Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(42*Power(x1,8) + 45*Power(x1,6)*Power(x2,2) + \
79*Power(x1,4)*Power(x2,4) + 23*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) - 15*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(210*Power(x1,12) + 168*Power(x1,10)*Power(x2,2) + \
1485*Power(x1,8)*Power(x2,4) + 452*Power(x1,6)*Power(x2,6) + \
440*Power(x1,4)*Power(x2,8) + 60*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(15*Power(x1,14) - 30*Power(x1,12)*Power(x2,2) + \
588*Power(x1,10)*Power(x2,4) - 375*Power(x1,8)*Power(x2,6) + \
647*Power(x1,6)*Power(x2,8) + 20*Power(x1,4)*Power(x2,10) + \
30*Power(x1,2)*Power(x2,12) + Power(x2,14)) + 1307674368000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,16) - \
120*Power(x1,14)*Power(x2,2) + 1820*Power(x1,12)*Power(x2,4) - \
8008*Power(x1,10)*Power(x2,6) + 12870*Power(x1,8)*Power(x2,8) - \
8008*Power(x1,6)*Power(x2,10) + 1820*Power(x1,4)*Power(x2,12) - \
120*Power(x1,2)*Power(x2,14) + Power(x2,16)) - \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,16) - 549*Power(x1,14)*Power(x2,2) + \
8085*Power(x1,12)*Power(x2,4) - 35889*Power(x1,10)*Power(x2,6) + \
57375*Power(x1,8)*Power(x2,8) - 35903*Power(x1,6)*Power(x2,10) + \
8071*Power(x1,4)*Power(x2,12) - 555*Power(x1,2)*Power(x2,14) + \
2*Power(x2,16)) - 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,16) - 90*Power(x1,14)*Power(x2,2) + \
2265*Power(x1,12)*Power(x2,4) - 8376*Power(x1,10)*Power(x2,6) + \
15165*Power(x1,8)*Power(x2,8) - 8306*Power(x1,6)*Power(x2,10) + \
2335*Power(x1,4)*Power(x2,12) - 60*Power(x1,2)*Power(x2,14) + \
8*Power(x2,16)) - 14175*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(9*Power(x1,16) - 1072*Power(x1,14)*Power(x2,2) + \
16268*Power(x1,12)*Power(x2,4) - 71568*Power(x1,10)*Power(x2,6) + \
115030*Power(x1,8)*Power(x2,8) - 71568*Power(x1,6)*Power(x2,10) + \
16268*Power(x1,4)*Power(x2,12) - 1072*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) - 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
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
	else if ( (k1==6) && (k2==8) )
		return -(he0*x1*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
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
	else if ( (k1==6) && (k2==9) )
		return (he0*x1*x2*(Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),15) - \
3*lambdasq*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(12*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + \
21*Power(lambdasq,2)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(18*Power(x1,8) + 24*Power(x1,6)*Power(x2,2) + \
59*Power(x1,4)*Power(x2,4) + 18*Power(x1,2)*Power(x2,6) + \
5*Power(x2,8)) - 105*Power(lambdasq,3)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(12*Power(x1,12) + 18*Power(x1,10)*Power(x2,2) + \
174*Power(x1,8)*Power(x2,4) + 41*Power(x1,6)*Power(x2,6) + \
93*Power(x1,4)*Power(x2,8) + 13*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(75*Power(x1,14) - 2478*Power(x1,12)*Power(x2,2) + \
19593*Power(x1,10)*Power(x2,4) - 51000*Power(x1,8)*Power(x2,6) + \
51245*Power(x1,6)*Power(x2,8) - 19446*Power(x1,4)*Power(x2,10) + \
2527*Power(x1,2)*Power(x2,12) - 68*Power(x2,14)) + \
315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,14) - 735*Power(x1,12)*Power(x2,2) + \
4536*Power(x1,10)*Power(x2,4) - 13380*Power(x1,8)*Power(x2,6) + \
12155*Power(x1,6)*Power(x2,8) - 5271*Power(x1,4)*Power(x2,10) + \
490*Power(x1,2)*Power(x2,12) - 38*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
315*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,14) + 336*Power(x1,10)*Power(x2,4) - \
360*Power(x1,8)*Power(x2,6) + 815*Power(x1,6)*Power(x2,8) - \
84*Power(x1,4)*Power(x2,10) + 70*Power(x1,2)*Power(x2,12) + \
4*Power(x2,14)) + 10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) + 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==6) && (k2==10) )
		return (he0*x1*(-(Power(x1,6)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,4)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),15)*(45*Power(x1,4) + 34*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,2)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(42*Power(x1,8) + 57*Power(x1,6)*Power(x2,2) + \
107*Power(x1,4)*Power(x2,4) + 35*Power(x1,2)*Power(x2,6) + \
7*Power(x2,8)) + 105*Power(lambdasq,3)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(30*Power(x1,12) + 48*Power(x1,10)*Power(x2,2) + \
285*Power(x1,8)*Power(x2,4) + 116*Power(x1,6)*Power(x2,6) + \
124*Power(x1,4)*Power(x2,8) + 20*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) - 105*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,14) + 30*Power(x1,12)*Power(x2,2) + \
2172*Power(x1,10)*Power(x2,4) - 885*Power(x1,8)*Power(x2,6) + \
3509*Power(x1,6)*Power(x2,8) + 196*Power(x1,4)*Power(x2,10) + \
290*Power(x1,2)*Power(x2,12) + 19*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,16)*(Power(x1,16) - \
136*Power(x1,14)*Power(x2,2) + 2380*Power(x1,12)*Power(x2,4) - \
12376*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
19448*Power(x1,6)*Power(x2,10) + 6188*Power(x1,4)*Power(x2,12) - \
680*Power(x1,2)*Power(x2,14) + 17*Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) \
+ Power(x2,2),2)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,16) - 30*Power(x1,14)*Power(x2,2) + \
2415*Power(x1,12)*Power(x2,4) - 8736*Power(x1,10)*Power(x2,6) + \
21795*Power(x1,8)*Power(x2,8) - 13894*Power(x1,6)*Power(x2,10) + \
6097*Power(x1,4)*Power(x2,12) - 220*Power(x1,2)*Power(x2,14) + \
58*Power(x2,16)) + 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,16) - 2613*Power(x1,14)*Power(x2,2) + \
41685*Power(x1,12)*Power(x2,4) - 223377*Power(x1,10)*Power(x2,6) + \
431295*Power(x1,8)*Power(x2,8) - 350671*Power(x1,6)*Power(x2,10) + \
108871*Power(x1,4)*Power(x2,12) - 12715*Power(x1,2)*Power(x2,14) + \
226*Power(x2,16)) + 1575*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(75*Power(x1,16) - 9696*Power(x1,14)*Power(x2,2) + \
170268*Power(x1,12)*Power(x2,4) - 884688*Power(x1,10)*Power(x2,6) + \
1738410*Power(x1,8)*Power(x2,8) - 1390336*Power(x1,6)*Power(x2,10) + \
442540*Power(x1,4)*Power(x2,12) - 48592*Power(x1,2)*Power(x2,14) + \
1219*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
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
	else if ( (k1==7) && (k2==8) )
		return (he0*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),15) - \
2*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(14*Power(x1,4) + 13*Power(x1,2)*Power(x2,2) + \
14*Power(x2,4)) + \
210*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,8) + 2*Power(x1,6)*Power(x2,2) + \
6*Power(x1,4)*Power(x2,4) + 2*Power(x1,2)*Power(x2,6) + Power(x2,8)) \
- 420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(Power(x1,4) + Power(x2,4))*(Power(x1,8) + \
4*Power(x1,6)*Power(x2,2) + 32*Power(x1,4)*Power(x2,4) + \
4*Power(x1,2)*Power(x2,6) + Power(x2,8)) + (1307674368000*(-1 + \
he0)*Power(lambdasq,15)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) \
+ 1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)))/he0 + \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,16) + 20*Power(x1,14)*Power(x2,2) + \
490*Power(x1,12)*Power(x2,4) - 700*Power(x1,10)*Power(x2,6) + \
2650*Power(x1,8)*Power(x2,8) - 700*Power(x1,6)*Power(x2,10) + \
490*Power(x1,4)*Power(x2,12) + 20*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 630*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,16) + 90*Power(x1,14)*Power(x2,2) - \
910*Power(x1,12)*Power(x2,4) + 4718*Power(x1,10)*Power(x2,6) - \
6870*Power(x1,8)*Power(x2,8) + 4718*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 3150*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,16) - 262*Power(x1,14)*Power(x2,2) + \
4088*Power(x1,12)*Power(x2,4) - 17850*Power(x1,10)*Power(x2,6) + \
28810*Power(x1,8)*Power(x2,8) - 17850*Power(x1,6)*Power(x2,10) + \
4088*Power(x1,4)*Power(x2,12) - 262*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 12600*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(10*Power(x1,16) - 1207*Power(x1,14)*Power(x2,2) + \
18298*Power(x1,12)*Power(x2,4) - 80521*Power(x1,10)*Power(x2,6) + \
129400*Power(x1,8)*Power(x2,8) - 80521*Power(x1,6)*Power(x2,10) + \
18298*Power(x1,4)*Power(x2,12) - 1207*Power(x1,2)*Power(x2,14) + \
10*Power(x2,16)) + 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==7) && (k2==9) )
		return -(he0*x2*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),16) - \
4*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(9*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(63*Power(x1,8) + 114*Power(x1,6)*Power(x2,2) + \
274*Power(x1,4)*Power(x2,4) + 98*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) - \
420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,12) + 9*Power(x1,10)*Power(x2,2) + \
57*Power(x1,8)*Power(x2,4) + 24*Power(x1,6)*Power(x2,6) + \
43*Power(x1,4)*Power(x2,8) + 7*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(17*Power(x1,16) - \
680*Power(x1,14)*Power(x2,2) + 6188*Power(x1,12)*Power(x2,4) - \
19448*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
12376*Power(x1,6)*Power(x2,10) + 2380*Power(x1,4)*Power(x2,12) - \
136*Power(x1,2)*Power(x2,14) + Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),2)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),3)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(9*Power(x1,16) + 60*Power(x1,14)*Power(x2,2) + \
1218*Power(x1,12)*Power(x2,4) - 708*Power(x1,10)*Power(x2,6) + \
3970*Power(x1,8)*Power(x2,8) - 244*Power(x1,6)*Power(x2,10) + \
610*Power(x1,4)*Power(x2,12) + 44*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 1260*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(3*Power(x1,16) + 195*Power(x1,14)*Power(x2,2) - \
1050*Power(x1,12)*Power(x2,4) + 4359*Power(x1,10)*Power(x2,6) - \
4460*Power(x1,8)*Power(x2,8) + 2849*Power(x1,6)*Power(x2,10) - \
350*Power(x1,4)*Power(x2,12) + 53*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 6300*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(303*Power(x1,16) - 12162*Power(x1,14)*Power(x2,2) + \
110586*Power(x1,12)*Power(x2,4) - 347682*Power(x1,10)*Power(x2,6) + \
434480*Power(x1,8)*Power(x2,8) - 221270*Power(x1,6)*Power(x2,10) + \
42518*Power(x1,4)*Power(x2,12) - 2438*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(183*Power(x1,16) - 5850*Power(x1,14)*Power(x2,2) + \
56028*Power(x1,12)*Power(x2,4) - 172518*Power(x1,10)*Power(x2,6) + \
218710*Power(x1,8)*Power(x2,8) - 109606*Power(x1,6)*Power(x2,10) + \
21700*Power(x1,4)*Power(x2,12) - 1114*Power(x1,2)*Power(x2,14) + \
19*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==7) && (k2==10) )
		return (he0*(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),17) - lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),16)*(45*Power(x1,4) + 39*Power(x1,2)*Power(x2,2) + \
28*Power(x2,4)) + \
2*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(315*Power(x1,8) + 540*Power(x1,6)*Power(x2,2) + \
1056*Power(x1,4)*Power(x2,4) + 392*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) - \
30*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(105*Power(x1,12) + 273*Power(x1,10)*Power(x2,2) + \
1314*Power(x1,8)*Power(x2,4) + 716*Power(x1,6)*Power(x2,6) + \
791*Power(x1,4)*Power(x2,8) + 147*Power(x1,2)*Power(x2,10) + \
14*Power(x2,12)) + \
105*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(45*Power(x1,16) + 180*Power(x1,14)*Power(x2,2) + \
2706*Power(x1,12)*Power(x2,4) - 72*Power(x1,10)*Power(x2,6) + \
5958*Power(x1,8)*Power(x2,8) + 636*Power(x1,6)*Power(x2,10) + \
810*Power(x1,4)*Power(x2,12) + 72*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 1575*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(1212*Power(x1,18) - 186003*Power(x1,16)*Power(x2,2) + \
3719304*Power(x1,14)*Power(x2,4) - 22564836*Power(x1,12)*Power(x2,6) \
+ 53187408*Power(x1,10)*Power(x2,8) - \
53188290*Power(x1,8)*Power(x2,10) + 22564248*Power(x1,6)*Power(x2,12) \
- 3719556*Power(x1,4)*Power(x2,14) + 185940*Power(x1,2)*Power(x2,16) \
- 1219*Power(x2,18)) - 630*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(183*Power(x1,18) - 22959*Power(x1,16)*Power(x2,2) + \
466110*Power(x1,14)*Power(x2,4) - 2817738*Power(x1,12)*Power(x2,6) + \
6652836*Power(x1,10)*Power(x2,8) - 6644016*Power(x1,8)*Power(x2,10) + \
2823618*Power(x1,6)*Power(x2,12) - 463590*Power(x1,4)*Power(x2,14) + \
23589*Power(x1,2)*Power(x2,16) - 113*Power(x2,18)) + \
355687428096000*(-1 + 1/he0)*Power(lambdasq,17)*(Power(x1,18) - \
153*Power(x1,16)*Power(x2,2) + 3060*Power(x1,14)*Power(x2,4) - \
18564*Power(x1,12)*Power(x2,6) + 43758*Power(x1,10)*Power(x2,8) - \
43758*Power(x1,8)*Power(x2,10) + 18564*Power(x1,6)*Power(x2,12) - \
3060*Power(x1,4)*Power(x2,14) + 153*Power(x1,2)*Power(x2,16) - \
Power(x2,18)) - 44460928512000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
7410154752000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
926269344000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
92626934400*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
7718911200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
551350800*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(9*Power(x1,18) + 135*Power(x1,16)*Power(x2,2) + \
8010*Power(x1,14)*Power(x2,4) - 23646*Power(x1,12)*Power(x2,6) + \
90162*Power(x1,10)*Power(x2,8) - 60174*Power(x1,8)*Power(x2,10) + \
41874*Power(x1,6)*Power(x2,12) - 1710*Power(x1,4)*Power(x2,14) + \
873*Power(x1,2)*Power(x2,16) + 19*Power(x2,18)) + \
630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(6*Power(x1,18) + 1287*Power(x1,16)*Power(x2,2) - \
19860*Power(x1,14)*Power(x2,4) + 131754*Power(x1,12)*Power(x2,6) - \
296388*Power(x1,10)*Power(x2,8) + 308148*Power(x1,8)*Power(x2,10) - \
124404*Power(x1,6)*Power(x2,12) + 22590*Power(x1,4)*Power(x2,14) - \
762*Power(x1,2)*Power(x2,16) + 29*Power(x2,18)) - \
177843714048000*Power(lambdasq,16)*(Power(x1,20) - \
152*Power(x1,18)*Power(x2,2) + 2907*Power(x1,16)*Power(x2,4) - \
15504*Power(x1,14)*Power(x2,6) + 25194*Power(x1,12)*Power(x2,8) - \
25194*Power(x1,8)*Power(x2,12) + 15504*Power(x1,6)*Power(x2,14) - \
2907*Power(x1,4)*Power(x2,16) + 152*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==8) && (k2==0) )
		return -(he0*x1*(Power(x1,8)*Power(Power(x1,2) + Power(x2,2),8) - \
4*lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,2) + 9*Power(x2,2)) + \
14*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(7*Power(x1,4) + 18*Power(x1,2)*Power(x2,2) + \
27*Power(x2,4)) + 20160*Power(lambdasq,7)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
33*Power(x1,4)*Power(x2,2) + 27*Power(x1,2)*Power(x2,4) - \
3*Power(x2,6)) - 84*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) + 9*Power(x1,4)*Power(x2,2) - \
9*Power(x1,2)*Power(x2,4) + 15*Power(x2,6)) + (40320*(-1 + \
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
	else if ( (k1==8) && (k2==1) )
		return (he0*x1*x2*(Power(x1,8)*Power(Power(x1,2) + Power(x2,2),9) - \
18*lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,2) + 2*Power(x2,2)) + (725760*(-1 + \
he0)*Power(lambdasq,9)*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 + 362880*Power(lambdasq,8)*(Power(x1,2) + \
Power(x2,2))*(5*Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 90720*Power(lambdasq,7)*Power(Power(x1,2) + \
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
5*Power(x2,4)) - \
252*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(5*Power(x1,4) - 6*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(5*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
21*Power(x2,4))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==8) && (k2==2) )
		return (he0*x1*(-(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),10)) + lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) - 18*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(Power(x1,6) + 8*Power(x1,4)*Power(x2,2) + \
8*Power(x1,2)*Power(x2,4) + 21*Power(x2,6)) - \
315*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(15*Power(x1,8) - 80*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 40*Power(x1,2)*Power(x2,6) + \
3*Power(x2,8)) + 90*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,8) + Power(x1,6)*Power(x2,2) + \
29*Power(x1,4)*Power(x2,4) - 21*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) + (3628800*(-1 + \
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
165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==8) && (k2==3) )
		return (he0*x1*x2*(Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,4) + 17*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) + 2*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(24*Power(x1,6) + 127*Power(x1,4)*Power(x2,2) + \
72*Power(x1,2)*Power(x2,4) + 189*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,8) + Power(x1,6)*Power(x2,2) + \
51*Power(x1,4)*Power(x2,4) - 21*Power(x1,2)*Power(x2,6) + \
14*Power(x2,8)) - 45*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(12*Power(x1,10) - 325*Power(x1,8)*Power(x2,2) + \
1044*Power(x1,6)*Power(x2,4) - 1134*Power(x1,4)*Power(x2,6) + \
280*Power(x1,2)*Power(x2,8) - 21*Power(x2,10)) + 159667200*(-1 + \
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
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==8) && (k2==4) )
		return -(he0*x1*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),12) - \
6*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,4) + 3*Power(x1,2)*Power(x2,2) + \
6*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,8) + 30*Power(x1,6)*Power(x2,2) + \
139*Power(x1,4)*Power(x2,4) + 60*Power(x1,2)*Power(x2,6) + \
126*Power(x2,8)) - 12*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),9)*(4*Power(x1,10) + 66*Power(x1,8)*Power(x2,2) + \
4*Power(x1,6)*Power(x2,4) + 591*Power(x1,4)*Power(x2,6) - \
126*Power(x1,2)*Power(x2,8) + 105*Power(x2,10)) + \
135*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,12) - 16*Power(x1,10)*Power(x2,2) + \
275*Power(x1,8)*Power(x2,4) - 548*Power(x1,6)*Power(x2,6) + \
476*Power(x1,4)*Power(x2,8) - 84*Power(x1,2)*Power(x2,10) + \
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
13*Power(x2,12)) + 270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(2*Power(x1,12) - 219*Power(x1,10)*Power(x2,2) + \
1955*Power(x1,8)*Power(x2,4) - 4734*Power(x1,6)*Power(x2,6) + \
3528*Power(x1,4)*Power(x2,8) - 791*Power(x1,2)*Power(x2,10) + \
35*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==8) && (k2==5) )
		return (he0*x1*x2*(Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13) - \
2*lambdasq*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(5*Power(x1,4) + 10*Power(x1,2)*Power(x2,2) + \
18*Power(x2,4)) + (12454041600*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 6227020800*Power(lambdasq,12)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1556755200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 259459200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3243240*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 270270*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 19305*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6))*(7*Power(x1,6) - \
35*Power(x1,4)*Power(x2,2) + 21*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(5*Power(x1,8) + 50*Power(x1,6)*Power(x2,2) + \
211*Power(x1,4)*Power(x2,4) + 84*Power(x1,2)*Power(x2,6) + \
126*Power(x2,8)) - 6*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),10)*(35*Power(x1,10) + 350*Power(x1,8)*Power(x2,2) + \
49*Power(x1,6)*Power(x2,4) + 1686*Power(x1,4)*Power(x2,6) - \
126*Power(x1,2)*Power(x2,8) + 210*Power(x2,10)) + \
270*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(28*Power(x1,12) - 833*Power(x1,10)*Power(x2,2) + \
4424*Power(x1,8)*Power(x2,4) - 7734*Power(x1,6)*Power(x2,6) + \
4424*Power(x1,4)*Power(x2,8) - 833*Power(x1,2)*Power(x2,10) + \
28*Power(x2,12)) + 15*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(98*Power(x1,12) - 448*Power(x1,10)*Power(x2,2) + \
5299*Power(x1,8)*Power(x2,4) - 6684*Power(x1,6)*Power(x2,6) + \
5124*Power(x1,4)*Power(x2,8) - 588*Power(x1,2)*Power(x2,10) + \
63*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==8) && (k2==6) )
		return (he0*x1*(-(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),14)) + lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(15*Power(x1,4) + 23*Power(x1,2)*Power(x2,2) + \
36*Power(x2,4)) + (87178291200*(-1 + \
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
28*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,8) + 240*Power(x1,6)*Power(x2,2) + \
905*Power(x1,4)*Power(x2,4) + 360*Power(x1,2)*Power(x2,6) + \
378*Power(x2,8)) + 3*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),11)*(5*Power(x1,12) + 195*Power(x1,10)*Power(x2,2) + \
1605*Power(x1,8)*Power(x2,4) + 475*Power(x1,6)*Power(x2,6) + \
4590*Power(x1,4)*Power(x2,8) + 126*Power(x1,2)*Power(x2,10) + \
420*Power(x2,12)) + 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(8*Power(x1,14) - 945*Power(x1,12)*Power(x2,2) + \
12180*Power(x1,10)*Power(x2,4) - 44765*Power(x1,8)*Power(x2,6) + \
57480*Power(x1,6)*Power(x2,8) - 26859*Power(x1,4)*Power(x2,10) + \
4060*Power(x1,2)*Power(x2,12) - 135*Power(x2,14)) + \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(14*Power(x1,14) - 390*Power(x1,12)*Power(x2,2) + \
6375*Power(x1,10)*Power(x2,4) - 21845*Power(x1,8)*Power(x2,6) + \
29340*Power(x1,6)*Power(x2,8) - 13032*Power(x1,4)*Power(x2,10) + \
2175*Power(x1,2)*Power(x2,12) - 45*Power(x2,14)) - \
21*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(10*Power(x1,14) + 300*Power(x1,12)*Power(x2,2) - \
750*Power(x1,10)*Power(x2,4) + 7115*Power(x1,8)*Power(x2,6) - \
5580*Power(x1,6)*Power(x2,8) + 4284*Power(x1,4)*Power(x2,10) - \
240*Power(x1,2)*Power(x2,12) + \
45*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==8) && (k2==7) )
		return (he0*x1*x2*(Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15) - \
3*lambdasq*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,4) + 9*Power(x1,2)*Power(x2,2) + \
12*Power(x2,4)) + \
21*Power(lambdasq,2)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(5*Power(x1,8) + 18*Power(x1,6)*Power(x2,2) + \
59*Power(x1,4)*Power(x2,4) + 24*Power(x1,2)*Power(x2,6) + \
18*Power(x2,8)) - 105*Power(lambdasq,3)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(Power(x1,12) + 13*Power(x1,10)*Power(x2,2) + \
93*Power(x1,8)*Power(x2,4) + 41*Power(x1,6)*Power(x2,6) + \
174*Power(x1,4)*Power(x2,8) + 18*Power(x1,2)*Power(x2,10) + \
12*Power(x2,12)) - 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(68*Power(x1,14) - 2527*Power(x1,12)*Power(x2,2) + \
19446*Power(x1,10)*Power(x2,4) - 51245*Power(x1,8)*Power(x2,6) + \
51000*Power(x1,6)*Power(x2,8) - 19593*Power(x1,4)*Power(x2,10) + \
2478*Power(x1,2)*Power(x2,12) - 75*Power(x2,14)) - \
315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(38*Power(x1,14) - 490*Power(x1,12)*Power(x2,2) + \
5271*Power(x1,10)*Power(x2,4) - 12155*Power(x1,8)*Power(x2,6) + \
13380*Power(x1,6)*Power(x2,8) - 4536*Power(x1,4)*Power(x2,10) + \
735*Power(x1,2)*Power(x2,12) - 3*Power(x2,14)) + 20922789888000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
315*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),11)*(4*Power(x1,14) + 70*Power(x1,12)*Power(x2,2) - \
84*Power(x1,10)*Power(x2,4) + 815*Power(x1,8)*Power(x2,6) - \
360*Power(x1,6)*Power(x2,8) + 336*Power(x1,4)*Power(x2,10) + \
3*Power(x2,14)) - 10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) - 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==8) && (k2==8) )
		return -(he0*x1*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),16) - \
4*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15)*(7*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
9*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(35*Power(x1,8) + 98*Power(x1,6)*Power(x2,2) + \
274*Power(x1,4)*Power(x2,4) + 114*Power(x1,2)*Power(x2,6) + \
63*Power(x2,8)) - \
420*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,12) + 7*Power(x1,10)*Power(x2,2) + \
43*Power(x1,8)*Power(x2,4) + 24*Power(x1,6)*Power(x2,6) + \
57*Power(x1,4)*Power(x2,8) + 9*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) - 1260*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,16) + 53*Power(x1,14)*Power(x2,2) - \
350*Power(x1,12)*Power(x2,4) + 2849*Power(x1,10)*Power(x2,6) - \
4460*Power(x1,8)*Power(x2,8) + 4359*Power(x1,6)*Power(x2,10) - \
1050*Power(x1,4)*Power(x2,12) + 195*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),12)*(Power(x1,16) + 44*Power(x1,14)*Power(x2,2) + \
610*Power(x1,12)*Power(x2,4) - 244*Power(x1,10)*Power(x2,6) + \
3970*Power(x1,8)*Power(x2,8) - 708*Power(x1,6)*Power(x2,10) + \
1218*Power(x1,4)*Power(x2,12) + 60*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) \
+ 2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) \
+ Power(x2,2),2)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(19*Power(x1,16) - 1114*Power(x1,14)*Power(x2,2) + \
21700*Power(x1,12)*Power(x2,4) - 109606*Power(x1,10)*Power(x2,6) + \
218710*Power(x1,8)*Power(x2,8) - 172518*Power(x1,6)*Power(x2,10) + \
56028*Power(x1,4)*Power(x2,12) - 5850*Power(x1,2)*Power(x2,14) + \
183*Power(x2,16)) + 6300*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(17*Power(x1,16) - 2438*Power(x1,14)*Power(x2,2) + \
42518*Power(x1,12)*Power(x2,4) - 221270*Power(x1,10)*Power(x2,6) + \
434480*Power(x1,8)*Power(x2,8) - 347682*Power(x1,6)*Power(x2,10) + \
110586*Power(x1,4)*Power(x2,12) - 12162*Power(x1,2)*Power(x2,14) + \
303*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==8) && (k2==9) )
		return (he0*x1*x2*(Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17) - \
2*lambdasq*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)*(18*Power(x1,4) + 19*Power(x1,2)*Power(x2,2) + \
18*Power(x2,4)) + (711374856192000*(-1 + \
he0)*Power(lambdasq,17)*(Power(x1,2) - 3*Power(x2,2))*(3*Power(x1,2) \
- Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 + 355687428096000*Power(lambdasq,16)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 88921857024000*Power(lambdasq,15)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 14820309504000*Power(lambdasq,14)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1852538688000*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 185253868800*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 15437822400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 1102701600*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 68918850*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + 3828825*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + \
2*Power(lambdasq,2)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(189*Power(x1,8) + 450*Power(x1,6)*Power(x2,2) + \
1066*Power(x1,4)*Power(x2,4) + 450*Power(x1,2)*Power(x2,6) + \
189*Power(x2,8)) - \
12*Power(lambdasq,3)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(105*Power(x1,12) + 504*Power(x1,10)*Power(x2,2) + \
2601*Power(x1,8)*Power(x2,4) + 1684*Power(x1,6)*Power(x2,6) + \
2601*Power(x1,4)*Power(x2,8) + 504*Power(x1,2)*Power(x2,10) + \
105*Power(x2,12)) + 105*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),13)*(9*Power(x1,16) + 132*Power(x1,14)*Power(x2,2) + \
1530*Power(x1,12)*Power(x2,4) + 180*Power(x1,10)*Power(x2,6) + \
6250*Power(x1,8)*Power(x2,8) + 180*Power(x1,6)*Power(x2,10) + \
1530*Power(x1,4)*Power(x2,12) + 132*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) - 210*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(45*Power(x1,16) + 1362*Power(x1,14)*Power(x2,2) - \
4662*Power(x1,12)*Power(x2,4) + 33462*Power(x1,10)*Power(x2,6) - \
34270*Power(x1,8)*Power(x2,8) + 33462*Power(x1,6)*Power(x2,10) - \
4662*Power(x1,4)*Power(x2,12) + 1362*Power(x1,2)*Power(x2,14) + \
45*Power(x2,16)) + 630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(207*Power(x1,16) - 4974*Power(x1,14)*Power(x2,2) + \
61488*Power(x1,12)*Power(x2,4) - 215154*Power(x1,10)*Power(x2,6) + \
341570*Power(x1,8)*Power(x2,8) - 215154*Power(x1,6)*Power(x2,10) + \
61488*Power(x1,4)*Power(x2,12) - 4974*Power(x1,2)*Power(x2,14) + \
207*Power(x2,16)) + 2520*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(666*Power(x1,16) - 31137*Power(x1,14)*Power(x2,2) + \
324954*Power(x1,12)*Power(x2,4) - 1209807*Power(x1,10)*Power(x2,6) + \
1845560*Power(x1,8)*Power(x2,8) - 1209807*Power(x1,6)*Power(x2,10) + \
324954*Power(x1,4)*Power(x2,12) - 31137*Power(x1,2)*Power(x2,14) + \
666*Power(x2,16))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==8) && (k2==10) )
		return (he0*x1*(-(Power(x1,8)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),18)) + \
9*lambdasq*Power(x1,6)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17)*(5*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
4*Power(x2,4)) - \
18*Power(lambdasq,2)*Power(x1,4)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)*(35*Power(x1,8) + 75*Power(x1,6)*Power(x2,2) + \
151*Power(x1,4)*Power(x2,4) + 64*Power(x1,2)*Power(x2,6) + \
21*Power(x2,8)) + \
18*Power(lambdasq,3)*Power(x1,2)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(175*Power(x1,12) + 665*Power(x1,10)*Power(x2,2) + \
2850*Power(x1,8)*Power(x2,4) + 2028*Power(x1,6)*Power(x2,6) + \
2257*Power(x1,4)*Power(x2,8) + 483*Power(x1,2)*Power(x2,10) + \
70*Power(x2,12)) - \
135*Power(lambdasq,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(35*Power(x1,16) + 280*Power(x1,14)*Power(x2,2) + \
2674*Power(x1,12)*Power(x2,4) + 1256*Power(x1,10)*Power(x2,6) + \
7554*Power(x1,8)*Power(x2,8) + 1416*Power(x1,6)*Power(x2,10) + \
1554*Power(x1,4)*Power(x2,12) + 168*Power(x1,2)*Power(x2,14) + \
7*Power(x2,16)) + 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(592*Power(x1,18) - 104067*Power(x1,16)*Power(x2,2) + \
2355072*Power(x1,14)*Power(x2,4) - 16490796*Power(x1,12)*Power(x2,6) \
+ 45932976*Power(x1,10)*Power(x2,8) - \
56144714*Power(x1,8)*Power(x2,10) + 30621984*Power(x1,6)*Power(x2,12) \
- 7067484*Power(x1,4)*Power(x2,14) + 588768*Power(x1,2)*Power(x2,16) \
- 11563*Power(x2,18)) + 5670*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(23*Power(x1,18) - 2253*Power(x1,16)*Power(x2,2) + \
54078*Power(x1,14)*Power(x2,4) - 373254*Power(x1,12)*Power(x2,6) + \
1046604*Power(x1,10)*Power(x2,8) - 1272976*Power(x1,8)*Power(x2,10) + \
698226*Power(x1,6)*Power(x2,12) - 159546*Power(x1,4)*Power(x2,14) + \
13677*Power(x1,2)*Power(x2,16) - 227*Power(x2,18)) + \
(6402373705728000*(-1 + he0)*Power(lambdasq,18)*(Power(x1,18) - \
171*Power(x1,16)*Power(x2,2) + 3876*Power(x1,14)*Power(x2,4) - \
27132*Power(x1,12)*Power(x2,6) + 75582*Power(x1,10)*Power(x2,8) - \
92378*Power(x1,8)*Power(x2,10) + 50388*Power(x1,6)*Power(x2,12) - \
11628*Power(x1,4)*Power(x2,14) + 969*Power(x1,2)*Power(x2,16) - \
19*Power(x2,18)))/he0 + \
3201186852864000*Power(lambdasq,17)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,18) + 45*Power(x1,16)*Power(x2,2) + \
1050*Power(x1,14)*Power(x2,4) - 1746*Power(x1,12)*Power(x2,6) + \
13626*Power(x1,10)*Power(x2,8) - 8090*Power(x1,8)*Power(x2,10) + \
9666*Power(x1,6)*Power(x2,12) - 306*Power(x1,4)*Power(x2,14) + \
393*Power(x1,2)*Power(x2,16) + 17*Power(x2,18)) - \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(5*Power(x1,18) + 468*Power(x1,16)*Power(x2,2) - \
5610*Power(x1,14)*Power(x2,4) + 50148*Power(x1,12)*Power(x2,6) - \
124578*Power(x1,10)*Power(x2,8) + 166178*Power(x1,8)*Power(x2,10) - \
82170*Power(x1,6)*Power(x2,12) + 22248*Power(x1,4)*Power(x2,14) - \
1119*Power(x1,2)*Power(x2,16) + \
94*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==9) && (k2==0) )
		return (he0*(Power(x1,10)*Power(Power(x1,2) + Power(x2,2),9) - \
9*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(3*Power(x1,2) + 5*Power(x2,2)) - \
126*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,2) + 5*Power(x2,2))*(3*Power(x1,4) + \
5*Power(x2,4)) + 18*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(11*Power(x1,4) + 30*Power(x1,2)*Power(x2,2) + \
35*Power(x2,4)) + 189*Power(lambdasq,4)*Power(Power(x1,2) + \
Power(x2,2),5)*Power(Power(x1,5) - 10*Power(x1,3)*Power(x2,2) + \
5*x1*Power(x2,4),2) + (362880*(-1 + \
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
181440*Power(lambdasq,8)*(Power(x1,12) - 44*Power(x1,10)*Power(x2,2) \
+ 165*Power(x1,8)*Power(x2,4) - 165*Power(x1,4)*Power(x2,8) + \
44*Power(x1,2)*Power(x2,10) - \
Power(x2,12))))/(2.*Power(lambdasq,9)*Pi*Power(Power(x1,2) + \
Power(x2,2),10));
	else if ( (k1==9) && (k2==1) )
		return (he0*x2*(-(Power(x1,10)*Power(Power(x1,2) + Power(x2,2),10)) \
+ 5*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(5*Power(x1,2) + 9*Power(x2,2)) - \
90*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
7*Power(x2,4)) + 90*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,6) + 25*Power(x1,4)*Power(x2,2) - \
7*Power(x1,2)*Power(x2,4) + 35*Power(x2,6)) - \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) - 40*Power(x1,6)*Power(x2,2) + \
118*Power(x1,4)*Power(x2,4) - 80*Power(x1,2)*Power(x2,6) + \
15*Power(x2,8)) + 3628800*(-1 + \
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
	else if ( (k1==9) && (k2==2) )
		return (he0*(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),11) - lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,4) + 24*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) + 5*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(5*Power(x1,6) + 48*Power(x1,4)*Power(x2,2) + \
81*Power(x1,2)*Power(x2,4) + 126*Power(x2,6)) - \
90*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,8) + 8*Power(x1,6)*Power(x2,2) + \
45*Power(x1,4)*Power(x2,4) - 14*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + 45*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),7)*(6*Power(x1,10) + 129*Power(x1,8)*Power(x2,2) - \
600*Power(x1,6)*Power(x2,4) + 1358*Power(x1,4)*Power(x2,6) - \
630*Power(x1,2)*Power(x2,8) + 105*Power(x2,10)) + 39916800*(-1 + \
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
Power(x2,12)) - 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,12) - 180*Power(x1,10)*Power(x2,2) + \
1365*Power(x1,8)*Power(x2,4) - 2536*Power(x1,6)*Power(x2,6) + \
1365*Power(x1,4)*Power(x2,8) - 180*Power(x1,2)*Power(x2,10) + \
3*Power(x2,12)) - 19958400*Power(lambdasq,10)*(Power(x1,14) - \
65*Power(x1,12)*Power(x2,2) + 429*Power(x1,10)*Power(x2,4) - \
429*Power(x1,8)*Power(x2,6) - 429*Power(x1,6)*Power(x2,8) + \
429*Power(x1,4)*Power(x2,10) - 65*Power(x1,2)*Power(x2,12) + \
Power(x2,14)) - 4989600*Power(lambdasq,9)*(Power(x1,16) - \
64*Power(x1,14)*Power(x2,2) + 364*Power(x1,12)*Power(x2,4) - \
858*Power(x1,8)*Power(x2,8) + 364*Power(x1,4)*Power(x2,12) - \
64*Power(x1,2)*Power(x2,14) + \
Power(x2,16))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==9) && (k2==3) )
		return (he0*x2*(-(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)) + 3*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
15*Power(x2,4)) - 3*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(23*Power(x1,6) + 124*Power(x1,4)*Power(x2,2) + \
135*Power(x1,2)*Power(x2,4) + 210*Power(x2,6)) + \
30*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(17*Power(x1,8) + 46*Power(x1,6)*Power(x2,2) + \
234*Power(x1,4)*Power(x2,4) - 42*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) - \
135*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,10) + 131*Power(x1,8)*Power(x2,2) - \
376*Power(x1,6)*Power(x2,4) + 658*Power(x1,4)*Power(x2,6) - \
210*Power(x1,2)*Power(x2,8) + 35*Power(x2,10)) + (479001600*(-1 + \
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
Power(x2,12)) + 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),7)*(73*Power(x1,12) - 1564*Power(x1,10)*Power(x2,2) + \
7101*Power(x1,8)*Power(x2,4) - 9408*Power(x1,6)*Power(x2,6) + \
3955*Power(x1,4)*Power(x2,8) - 420*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==9) && (k2==4) )
		return (he0*(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),12)*(6*Power(x1,4) + 25*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) + 3*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,8) + 44*Power(x1,6)*Power(x2,2) + \
191*Power(x1,4)*Power(x2,4) + 150*Power(x1,2)*Power(x2,6) + \
210*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),10)*(23*Power(x1,10) + 427*Power(x1,8)*Power(x2,2) + \
763*Power(x1,6)*Power(x2,4) + 3675*Power(x1,4)*Power(x2,6) - \
210*Power(x1,2)*Power(x2,8) + 1050*Power(x2,10)) + \
15*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(34*Power(x1,12) + 56*Power(x1,10)*Power(x2,2) + \
3059*Power(x1,8)*Power(x2,4) - 5376*Power(x1,6)*Power(x2,6) + \
7896*Power(x1,4)*Power(x2,8) - 1680*Power(x1,2)*Power(x2,10) + \
315*Power(x2,12)) + 135*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),7)*(73*Power(x1,14) - 6496*Power(x1,12)*Power(x2,2) + \
71603*Power(x1,10)*Power(x2,4) - 214662*Power(x1,8)*Power(x2,6) + \
214767*Power(x1,6)*Power(x2,8) - 71540*Power(x1,4)*Power(x2,10) + \
6517*Power(x1,2)*Power(x2,12) - 70*Power(x2,14)) + (6227020800*(-1 + \
he0)*Power(lambdasq,13)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
778377600*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
129729600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
16216200*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
1621620*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
135135*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 91*Power(x1,12)*Power(x2,2) + \
1001*Power(x1,10)*Power(x2,4) - 3003*Power(x1,8)*Power(x2,6) + \
3003*Power(x1,6)*Power(x2,8) - 1001*Power(x1,4)*Power(x2,10) + \
91*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,14) + 448*Power(x1,12)*Power(x2,2) - \
4361*Power(x1,10)*Power(x2,4) + 13587*Power(x1,8)*Power(x2,6) - \
13272*Power(x1,6)*Power(x2,8) + 4550*Power(x1,4)*Power(x2,10) - \
385*Power(x1,2)*Power(x2,12) + 7*Power(x2,14)) + \
3113510400*Power(lambdasq,12)*(Power(x1,16) - \
90*Power(x1,14)*Power(x2,2) + 910*Power(x1,12)*Power(x2,4) - \
2002*Power(x1,10)*Power(x2,6) + 2002*Power(x1,6)*Power(x2,10) - \
910*Power(x1,4)*Power(x2,12) + 90*Power(x1,2)*Power(x2,14) - \
Power(x2,16))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==9) && (k2==5) )
		return (he0*x2*(-(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)) + lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),13)*(10*Power(x1,4) + 27*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - (87178291200*(-1 + \
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
92*Power(x1,2)*Power(x2,6) + Power(x2,8)) - \
Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(15*Power(x1,8) + 220*Power(x1,6)*Power(x2,2) + \
843*Power(x1,4)*Power(x2,4) + 540*Power(x1,2)*Power(x2,6) + \
630*Power(x2,8)) + 21*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),11)*(15*Power(x1,10) + 145*Power(x1,8)*Power(x2,2) + \
183*Power(x1,6)*Power(x2,4) + 765*Power(x1,4)*Power(x2,6) + \
30*Power(x1,2)*Power(x2,8) + 150*Power(x2,10)) - \
21*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(120*Power(x1,12) + 110*Power(x1,10)*Power(x2,2) + \
4899*Power(x1,8)*Power(x2,4) - 5130*Power(x1,6)*Power(x2,6) + \
7140*Power(x1,4)*Power(x2,8) - 900*Power(x1,2)*Power(x2,10) + \
225*Power(x2,12)) - 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(30*Power(x1,14) - 2260*Power(x1,12)*Power(x2,2) + \
12837*Power(x1,10)*Power(x2,4) - 29565*Power(x1,8)*Power(x2,6) + \
21720*Power(x1,6)*Power(x2,8) - 6390*Power(x1,4)*Power(x2,10) + \
405*Power(x1,2)*Power(x2,12) - 9*Power(x2,14)) - \
945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(135*Power(x1,14) - 4060*Power(x1,12)*Power(x2,2) + \
26859*Power(x1,10)*Power(x2,4) - 57480*Power(x1,8)*Power(x2,6) + \
44765*Power(x1,6)*Power(x2,8) - 12180*Power(x1,4)*Power(x2,10) + \
945*Power(x1,2)*Power(x2,12) - \
8*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==9) && (k2==6) )
		return (he0*(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),15) - \
15*lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14)*(Power(x1,4) + 2*Power(x1,2)*Power(x2,2) + \
3*Power(x2,4)) + \
15*Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,8) + 23*Power(x1,6)*Power(x2,2) + \
79*Power(x1,4)*Power(x2,4) + 45*Power(x1,2)*Power(x2,6) + \
42*Power(x2,8)) - 15*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(Power(x1,12) + 60*Power(x1,10)*Power(x2,2) + \
440*Power(x1,8)*Power(x2,4) + 452*Power(x1,6)*Power(x2,6) + \
1485*Power(x1,4)*Power(x2,8) + 168*Power(x1,2)*Power(x2,10) + \
210*Power(x2,12)) + \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,14) + 30*Power(x1,12)*Power(x2,2) + \
20*Power(x1,10)*Power(x2,4) + 647*Power(x1,8)*Power(x2,6) - \
375*Power(x1,6)*Power(x2,8) + 588*Power(x1,4)*Power(x2,10) - \
30*Power(x1,2)*Power(x2,12) + 15*Power(x2,14)) + 1307674368000*(-1 + \
1/he0)*Power(lambdasq,15)*(Power(x1,16) - \
120*Power(x1,14)*Power(x2,2) + 1820*Power(x1,12)*Power(x2,4) - \
8008*Power(x1,10)*Power(x2,6) + 12870*Power(x1,8)*Power(x2,8) - \
8008*Power(x1,6)*Power(x2,10) + 1820*Power(x1,4)*Power(x2,12) - \
120*Power(x1,2)*Power(x2,14) + Power(x2,16)) - \
163459296000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 27243216000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 3405402000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 340540200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 28378350*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 120*Power(x1,14)*Power(x2,2) + \
1820*Power(x1,12)*Power(x2,4) - 8008*Power(x1,10)*Power(x2,6) + \
12870*Power(x1,8)*Power(x2,8) - 8008*Power(x1,6)*Power(x2,10) + \
1820*Power(x1,4)*Power(x2,12) - 120*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) - 1575*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(2*Power(x1,16) - 555*Power(x1,14)*Power(x2,2) + \
8071*Power(x1,12)*Power(x2,4) - 35903*Power(x1,10)*Power(x2,6) + \
57375*Power(x1,8)*Power(x2,8) - 35889*Power(x1,6)*Power(x2,10) + \
8085*Power(x1,4)*Power(x2,12) - 549*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) - 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(8*Power(x1,16) - 60*Power(x1,14)*Power(x2,2) + \
2335*Power(x1,12)*Power(x2,4) - 8306*Power(x1,10)*Power(x2,6) + \
15165*Power(x1,8)*Power(x2,8) - 8376*Power(x1,6)*Power(x2,10) + \
2265*Power(x1,4)*Power(x2,12) - 90*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) - 14175*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(9*Power(x1,16) - 1072*Power(x1,14)*Power(x2,2) + \
16268*Power(x1,12)*Power(x2,4) - 71568*Power(x1,10)*Power(x2,6) + \
115030*Power(x1,8)*Power(x2,8) - 71568*Power(x1,6)*Power(x2,10) + \
16268*Power(x1,4)*Power(x2,12) - 1072*Power(x1,2)*Power(x2,14) + \
9*Power(x2,16)) - 653837184000*Power(lambdasq,14)*(Power(x1,18) - \
119*Power(x1,16)*Power(x2,2) + 1700*Power(x1,14)*Power(x2,4) - \
6188*Power(x1,12)*Power(x2,6) + 4862*Power(x1,10)*Power(x2,8) + \
4862*Power(x1,8)*Power(x2,10) - 6188*Power(x1,6)*Power(x2,12) + \
1700*Power(x1,4)*Power(x2,14) - 119*Power(x1,2)*Power(x2,16) + \
Power(x2,18))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==9) && (k2==7) )
		return (he0*x2*(-(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),15)*(21*Power(x1,4) + 34*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,8) + 35*Power(x1,6)*Power(x2,2) + \
107*Power(x1,4)*Power(x2,4) + 57*Power(x1,2)*Power(x2,6) + \
42*Power(x2,8)) + 105*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(Power(x1,12) + 20*Power(x1,10)*Power(x2,2) + \
124*Power(x1,8)*Power(x2,4) + 116*Power(x1,6)*Power(x2,6) + \
285*Power(x1,4)*Power(x2,8) + 48*Power(x1,2)*Power(x2,10) + \
30*Power(x2,12)) - \
105*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(19*Power(x1,14) + 290*Power(x1,12)*Power(x2,2) + \
196*Power(x1,10)*Power(x2,4) + 3509*Power(x1,8)*Power(x2,6) - \
885*Power(x1,6)*Power(x2,8) + 2172*Power(x1,4)*Power(x2,10) + \
30*Power(x1,2)*Power(x2,12) + 45*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,16)*(17*Power(x1,16) - \
680*Power(x1,14)*Power(x2,2) + 6188*Power(x1,12)*Power(x2,4) - \
19448*Power(x1,10)*Power(x2,6) + 24310*Power(x1,8)*Power(x2,8) - \
12376*Power(x1,6)*Power(x2,10) + 2380*Power(x1,4)*Power(x2,12) - \
136*Power(x1,2)*Power(x2,14) + Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),2)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),3)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(17*Power(x1,16) - 680*Power(x1,14)*Power(x2,2) + \
6188*Power(x1,12)*Power(x2,4) - 19448*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 12376*Power(x1,6)*Power(x2,10) + \
2380*Power(x1,4)*Power(x2,12) - 136*Power(x1,2)*Power(x2,14) + \
Power(x2,16)) + 315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(226*Power(x1,16) - 12715*Power(x1,14)*Power(x2,2) + \
108871*Power(x1,12)*Power(x2,4) - 350671*Power(x1,10)*Power(x2,6) + \
431295*Power(x1,8)*Power(x2,8) - 223377*Power(x1,6)*Power(x2,10) + \
41685*Power(x1,4)*Power(x2,12) - 2613*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(58*Power(x1,16) - 220*Power(x1,14)*Power(x2,2) + \
6097*Power(x1,12)*Power(x2,4) - 13894*Power(x1,10)*Power(x2,6) + \
21795*Power(x1,8)*Power(x2,8) - 8736*Power(x1,6)*Power(x2,10) + \
2415*Power(x1,4)*Power(x2,12) - 30*Power(x1,2)*Power(x2,14) + \
3*Power(x2,16)) + 1575*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(1219*Power(x1,16) - 48592*Power(x1,14)*Power(x2,2) + \
442540*Power(x1,12)*Power(x2,4) - 1390336*Power(x1,10)*Power(x2,6) + \
1738410*Power(x1,8)*Power(x2,8) - 884688*Power(x1,6)*Power(x2,10) + \
170268*Power(x1,4)*Power(x2,12) - 9696*Power(x1,2)*Power(x2,14) + \
75*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==9) && (k2==8) )
		return (he0*(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),17) - lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),16)*(28*Power(x1,4) + 39*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) + \
2*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15)*(105*Power(x1,8) + 392*Power(x1,6)*Power(x2,2) + \
1056*Power(x1,4)*Power(x2,4) + 540*Power(x1,2)*Power(x2,6) + \
315*Power(x2,8)) - \
30*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(14*Power(x1,12) + 147*Power(x1,10)*Power(x2,2) + \
791*Power(x1,8)*Power(x2,4) + 716*Power(x1,6)*Power(x2,6) + \
1314*Power(x1,4)*Power(x2,8) + 273*Power(x1,2)*Power(x2,10) + \
105*Power(x2,12)) + \
105*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(Power(x1,16) + 72*Power(x1,14)*Power(x2,2) + \
810*Power(x1,12)*Power(x2,4) + 636*Power(x1,10)*Power(x2,6) + \
5958*Power(x1,8)*Power(x2,8) - 72*Power(x1,6)*Power(x2,10) + \
2706*Power(x1,4)*Power(x2,12) + 180*Power(x1,2)*Power(x2,14) + \
45*Power(x2,16)) + 1575*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(1219*Power(x1,18) - 185940*Power(x1,16)*Power(x2,2) + \
3719556*Power(x1,14)*Power(x2,4) - 22564248*Power(x1,12)*Power(x2,6) \
+ 53188290*Power(x1,10)*Power(x2,8) - \
53187408*Power(x1,8)*Power(x2,10) + 22564836*Power(x1,6)*Power(x2,12) \
- 3719304*Power(x1,4)*Power(x2,14) + 186003*Power(x1,2)*Power(x2,16) \
- 1212*Power(x2,18)) + 630*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(113*Power(x1,18) - 23589*Power(x1,16)*Power(x2,2) + \
463590*Power(x1,14)*Power(x2,4) - 2823618*Power(x1,12)*Power(x2,6) + \
6644016*Power(x1,10)*Power(x2,8) - 6652836*Power(x1,8)*Power(x2,10) + \
2817738*Power(x1,6)*Power(x2,12) - 466110*Power(x1,4)*Power(x2,14) + \
22959*Power(x1,2)*Power(x2,16) - 183*Power(x2,18)) + \
(355687428096000*(-1 + he0)*Power(lambdasq,17)*(Power(x1,18) - \
153*Power(x1,16)*Power(x2,2) + 3060*Power(x1,14)*Power(x2,4) - \
18564*Power(x1,12)*Power(x2,6) + 43758*Power(x1,10)*Power(x2,8) - \
43758*Power(x1,8)*Power(x2,10) + 18564*Power(x1,6)*Power(x2,12) - \
3060*Power(x1,4)*Power(x2,14) + 153*Power(x1,2)*Power(x2,16) - \
Power(x2,18)))/he0 + \
44460928512000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
7410154752000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
926269344000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
92626934400*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
7718911200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
551350800*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 153*Power(x1,16)*Power(x2,2) + \
3060*Power(x1,14)*Power(x2,4) - 18564*Power(x1,12)*Power(x2,6) + \
43758*Power(x1,10)*Power(x2,8) - 43758*Power(x1,8)*Power(x2,10) + \
18564*Power(x1,6)*Power(x2,12) - 3060*Power(x1,4)*Power(x2,14) + \
153*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
630*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(29*Power(x1,18) - 762*Power(x1,16)*Power(x2,2) + \
22590*Power(x1,14)*Power(x2,4) - 124404*Power(x1,12)*Power(x2,6) + \
308148*Power(x1,10)*Power(x2,8) - 296388*Power(x1,8)*Power(x2,10) + \
131754*Power(x1,6)*Power(x2,12) - 19860*Power(x1,4)*Power(x2,14) + \
1287*Power(x1,2)*Power(x2,16) + 6*Power(x2,18)) - \
105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(19*Power(x1,18) + 873*Power(x1,16)*Power(x2,2) - \
1710*Power(x1,14)*Power(x2,4) + 41874*Power(x1,12)*Power(x2,6) - \
60174*Power(x1,10)*Power(x2,8) + 90162*Power(x1,8)*Power(x2,10) - \
23646*Power(x1,6)*Power(x2,12) + 8010*Power(x1,4)*Power(x2,14) + \
135*Power(x1,2)*Power(x2,16) + 9*Power(x2,18)) + \
177843714048000*Power(lambdasq,16)*(Power(x1,20) - \
152*Power(x1,18)*Power(x2,2) + 2907*Power(x1,16)*Power(x2,4) - \
15504*Power(x1,14)*Power(x2,6) + 25194*Power(x1,12)*Power(x2,8) - \
25194*Power(x1,8)*Power(x2,12) + 15504*Power(x1,6)*Power(x2,14) - \
2907*Power(x1,4)*Power(x2,16) + 152*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==9) && (k2==9) )
		return (he0*x2*(-(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),18)) + \
9*lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17)*(4*Power(x1,4) + 5*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - \
18*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(21*Power(x1,8) + 64*Power(x1,6)*Power(x2,2) + \
151*Power(x1,4)*Power(x2,4) + 75*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + \
18*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(70*Power(x1,12) + 483*Power(x1,10)*Power(x2,2) + \
2257*Power(x1,8)*Power(x2,4) + 2028*Power(x1,6)*Power(x2,6) + \
2850*Power(x1,4)*Power(x2,8) + 665*Power(x1,2)*Power(x2,10) + \
175*Power(x2,12)) - \
135*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,16) + 168*Power(x1,14)*Power(x2,2) + \
1554*Power(x1,12)*Power(x2,4) + 1416*Power(x1,10)*Power(x2,6) + \
7554*Power(x1,8)*Power(x2,8) + 1256*Power(x1,6)*Power(x2,10) + \
2674*Power(x1,4)*Power(x2,12) + 280*Power(x1,2)*Power(x2,14) + \
35*Power(x2,16)) - 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(11563*Power(x1,18) - 588768*Power(x1,16)*Power(x2,2) \
+ 7067484*Power(x1,14)*Power(x2,4) - \
30621984*Power(x1,12)*Power(x2,6) + 56144714*Power(x1,10)*Power(x2,8) \
- 45932976*Power(x1,8)*Power(x2,10) + \
16490796*Power(x1,6)*Power(x2,12) - 2355072*Power(x1,4)*Power(x2,14) \
+ 104067*Power(x1,2)*Power(x2,16) - 592*Power(x2,18)) - \
5670*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(227*Power(x1,18) - 13677*Power(x1,16)*Power(x2,2) + \
159546*Power(x1,14)*Power(x2,4) - 698226*Power(x1,12)*Power(x2,6) + \
1272976*Power(x1,10)*Power(x2,8) - 1046604*Power(x1,8)*Power(x2,10) + \
373254*Power(x1,6)*Power(x2,12) - 54078*Power(x1,4)*Power(x2,14) + \
2253*Power(x1,2)*Power(x2,16) - 23*Power(x2,18)) + \
6402373705728000*(-1 + 1/he0)*Power(lambdasq,18)*(19*Power(x1,18) - \
969*Power(x1,16)*Power(x2,2) + 11628*Power(x1,14)*Power(x2,4) - \
50388*Power(x1,12)*Power(x2,6) + 92378*Power(x1,10)*Power(x2,8) - \
75582*Power(x1,8)*Power(x2,10) + 27132*Power(x1,6)*Power(x2,12) - \
3876*Power(x1,4)*Power(x2,14) + 171*Power(x1,2)*Power(x2,16) - \
Power(x2,18)) - 3201186852864000*Power(lambdasq,17)*(Power(x1,2) + \
Power(x2,2))*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) - \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(19*Power(x1,18) - 969*Power(x1,16)*Power(x2,2) + \
11628*Power(x1,14)*Power(x2,4) - 50388*Power(x1,12)*Power(x2,6) + \
92378*Power(x1,10)*Power(x2,8) - 75582*Power(x1,8)*Power(x2,10) + \
27132*Power(x1,6)*Power(x2,12) - 3876*Power(x1,4)*Power(x2,14) + \
171*Power(x1,2)*Power(x2,16) - Power(x2,18)) + \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(17*Power(x1,18) + 393*Power(x1,16)*Power(x2,2) - \
306*Power(x1,14)*Power(x2,4) + 9666*Power(x1,12)*Power(x2,6) - \
8090*Power(x1,10)*Power(x2,8) + 13626*Power(x1,8)*Power(x2,10) - \
1746*Power(x1,6)*Power(x2,12) + 1050*Power(x1,4)*Power(x2,14) + \
45*Power(x1,2)*Power(x2,16) + Power(x2,18)) - \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(94*Power(x1,18) - 1119*Power(x1,16)*Power(x2,2) + \
22248*Power(x1,14)*Power(x2,4) - 82170*Power(x1,12)*Power(x2,6) + \
166178*Power(x1,10)*Power(x2,8) - 124578*Power(x1,8)*Power(x2,10) + \
50148*Power(x1,6)*Power(x2,12) - 5610*Power(x1,4)*Power(x2,14) + \
468*Power(x1,2)*Power(x2,16) + \
5*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==9) && (k2==10) )
		return (he0*(Power(x1,10)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),19) - lambdasq*Power(x1,8)*Power(x2,8)*Power(Power(x1,2) \
+ Power(x2,2),18)*(45*Power(x1,4) + 52*Power(x1,2)*Power(x2,2) + \
45*Power(x2,4)) + \
9*Power(lambdasq,2)*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17)*(70*Power(x1,8) + 185*Power(x1,6)*Power(x2,2) + \
382*Power(x1,4)*Power(x2,4) + 185*Power(x1,2)*Power(x2,6) + \
70*Power(x2,8)) - \
18*Power(lambdasq,3)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(175*Power(x1,12) + 910*Power(x1,10)*Power(x2,2) + \
3680*Power(x1,8)*Power(x2,4) + 3306*Power(x1,6)*Power(x2,6) + \
3680*Power(x1,4)*Power(x2,8) + 910*Power(x1,2)*Power(x2,10) + \
175*Power(x2,12)) + \
9*Power(lambdasq,4)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(525*Power(x1,16) + 6650*Power(x1,14)*Power(x2,2) + \
52010*Power(x1,12)*Power(x2,4) + 53230*Power(x1,10)*Power(x2,6) + \
180066*Power(x1,8)*Power(x2,8) + 53230*Power(x1,6)*Power(x2,10) + \
52010*Power(x1,4)*Power(x2,12) + 6650*Power(x1,2)*Power(x2,14) + \
525*Power(x2,16)) + 121645100408832000*(-1 + \
1/he0)*Power(lambdasq,19)*(Power(x1,20) - \
190*Power(x1,18)*Power(x2,2) + 4845*Power(x1,16)*Power(x2,4) - \
38760*Power(x1,14)*Power(x2,6) + 125970*Power(x1,12)*Power(x2,8) - \
184756*Power(x1,10)*Power(x2,10) + 125970*Power(x1,8)*Power(x2,12) - \
38760*Power(x1,6)*Power(x2,14) + 4845*Power(x1,4)*Power(x2,16) - \
190*Power(x1,2)*Power(x2,18) + Power(x2,20)) - \
15205637551104000*Power(lambdasq,17)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 2534272925184000*Power(lambdasq,16)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 316784115648000*Power(lambdasq,15)*Power(Power(x1,2) \
+ Power(x2,2),4)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 31678411564800*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 2639867630400*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 188561973600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 11785123350*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 654729075*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,20) - 190*Power(x1,18)*Power(x2,2) + \
4845*Power(x1,16)*Power(x2,4) - 38760*Power(x1,14)*Power(x2,6) + \
125970*Power(x1,12)*Power(x2,8) - 184756*Power(x1,10)*Power(x2,10) + \
125970*Power(x1,8)*Power(x2,12) - 38760*Power(x1,6)*Power(x2,14) + \
4845*Power(x1,4)*Power(x2,16) - 190*Power(x1,2)*Power(x2,18) + \
Power(x2,20)) - 135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,20) + 560*Power(x1,18)*Power(x2,2) + \
9345*Power(x1,16)*Power(x2,4) + 131400*Power(x1,12)*Power(x2,8) - \
50368*Power(x1,10)*Power(x2,10) + 131400*Power(x1,8)*Power(x2,12) + \
9345*Power(x1,4)*Power(x2,16) + 560*Power(x1,2)*Power(x2,18) + \
7*Power(x2,20)) + 945*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(17*Power(x1,20) + 1180*Power(x1,18)*Power(x2,2) - \
7305*Power(x1,16)*Power(x2,4) + 114300*Power(x1,14)*Power(x2,6) - \
285480*Power(x1,12)*Power(x2,8) + 506512*Power(x1,10)*Power(x2,10) - \
285480*Power(x1,8)*Power(x2,12) + 114300*Power(x1,6)*Power(x2,14) - \
7305*Power(x1,4)*Power(x2,16) + 1180*Power(x1,2)*Power(x2,18) + \
17*Power(x2,20)) - 3780*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),12)*(47*Power(x1,20) - 2630*Power(x1,18)*Power(x2,2) + \
82290*Power(x1,16)*Power(x2,4) - 627030*Power(x1,14)*Power(x2,6) + \
2082735*Power(x1,12)*Power(x2,8) - 3010088*Power(x1,10)*Power(x2,10) \
+ 2082735*Power(x1,8)*Power(x2,12) - 627030*Power(x1,6)*Power(x2,14) \
+ 82290*Power(x1,4)*Power(x2,16) - 2630*Power(x1,2)*Power(x2,18) + \
47*Power(x2,20)) - 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(454*Power(x1,20) - 100435*Power(x1,18)*Power(x2,2) + \
2539830*Power(x1,16)*Power(x2,4) - 20352660*Power(x1,14)*Power(x2,6) \
+ 66103620*Power(x1,12)*Power(x2,8) - \
96991666*Power(x1,10)*Power(x2,10) + \
66103620*Power(x1,8)*Power(x2,12) - 20352660*Power(x1,6)*Power(x2,14) \
+ 2539830*Power(x1,4)*Power(x2,16) - 100435*Power(x1,2)*Power(x2,18) \
+ 454*Power(x2,20)) - 2835*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),10)*(11563*Power(x1,20) - \
2193820*Power(x1,18)*Power(x2,2) + 55947135*Power(x1,16)*Power(x2,4) \
- 447569520*Power(x1,14)*Power(x2,6) + \
1454610390*Power(x1,12)*Power(x2,8) - \
2133419752*Power(x1,10)*Power(x2,10) + \
1454610390*Power(x1,8)*Power(x2,12) - \
447569520*Power(x1,6)*Power(x2,14) + \
55947135*Power(x1,4)*Power(x2,16) - 2193820*Power(x1,2)*Power(x2,18) \
+ 11563*Power(x2,20)) - \
60822550204416000*Power(lambdasq,18)*(Power(x1,22) - \
189*Power(x1,20)*Power(x2,2) + 4655*Power(x1,18)*Power(x2,4) - \
33915*Power(x1,16)*Power(x2,6) + 87210*Power(x1,14)*Power(x2,8) - \
58786*Power(x1,12)*Power(x2,10) - 58786*Power(x1,10)*Power(x2,12) + \
87210*Power(x1,8)*Power(x2,14) - 33915*Power(x1,6)*Power(x2,16) + \
4655*Power(x1,4)*Power(x2,18) - 189*Power(x1,2)*Power(x2,20) + \
Power(x2,22))))/(2.*Power(lambdasq,19)*Pi*Power(Power(x1,2) + \
Power(x2,2),20));
	else if ( (k1==10) && (k2==0) )
		return -(he0*x1*(Power(x1,10)*Power(Power(x1,2) + Power(x2,2),10) - \
5*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,2) + 11*Power(x2,2)) + \
90*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(4*Power(x1,4) + 11*Power(x1,2)*Power(x2,2) + \
11*Power(x2,4)) - 90*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),7)*(13*Power(x1,6) + 55*Power(x1,4)*Power(x2,2) + \
55*Power(x1,2)*Power(x2,4) + 77*Power(x2,6)) + \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),6)*(3*Power(x1,8) + 110*Power(x1,4)*Power(x2,4) - \
88*Power(x1,2)*Power(x2,6) + 55*Power(x2,8)) + (3628800*(-1 + \
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
165*Power(x1,2)*Power(x2,8) - \
11*Power(x2,10))))/(2.*Power(lambdasq,10)*Pi*Power(Power(x1,2) + \
Power(x2,2),11));
	else if ( (k1==10) && (k2==1) )
		return (he0*x1*x2*(Power(x1,10)*Power(Power(x1,2) + Power(x2,2),11) \
- 11*lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(3*Power(x1,2) + 5*Power(x2,2)) + \
110*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(3*Power(x1,4) + 8*Power(x1,2)*Power(x2,2) + \
9*Power(x2,4)) - 990*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),8)*(Power(x1,6) + 5*Power(x1,4)*Power(x2,2) + \
3*Power(x1,2)*Power(x2,4) + 7*Power(x2,6)) + \
495*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),7)*(3*Power(x1,8) - 20*Power(x1,6)*Power(x2,2) + \
114*Power(x1,4)*Power(x2,4) - 84*Power(x1,2)*Power(x2,6) + \
35*Power(x2,8)) + (159667200*(-1 + \
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
55*Power(x1,2)*Power(x2,8) - 3*Power(x2,10)) + \
79833600*Power(lambdasq,10)*(3*Power(x1,12) - \
52*Power(x1,10)*Power(x2,2) + 143*Power(x1,8)*Power(x2,4) - \
143*Power(x1,4)*Power(x2,8) + 52*Power(x1,2)*Power(x2,10) - \
3*Power(x2,12))))/(2.*Power(lambdasq,11)*Pi*Power(Power(x1,2) + \
Power(x2,2),12));
	else if ( (k1==10) && (k2==2) )
		return (he0*x1*(-(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),12)) + \
330*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,4) + Power(x2,4))*(Power(x1,4) + \
6*Power(x1,2)*Power(x2,2) + 21*Power(x2,4)) + \
lambdasq*Power(x1,8)*Power(Power(x1,2) + Power(x2,2),11)*(Power(x1,4) \
+ 32*Power(x1,2)*Power(x2,2) + 55*Power(x2,4)) + \
1485*Power(lambdasq,5)*(x1 - x2)*(x1 + x2)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,2) - 4*x1*x2 + Power(x2,2))*(Power(x1,2) + \
4*x1*x2 + Power(x2,2))*(Power(x1,6) - 21*Power(x1,4)*Power(x2,2) + \
35*Power(x1,2)*Power(x2,4) - 7*Power(x2,6)) - \
33*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(Power(x1,6) + 12*Power(x1,4)*Power(x2,2) + \
25*Power(x1,2)*Power(x2,4) + 30*Power(x2,6)) - \
495*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),8)*(2*Power(x1,10) + 19*Power(x1,8)*Power(x2,2) - \
40*Power(x1,6)*Power(x2,4) + 194*Power(x1,4)*Power(x2,6) - \
98*Power(x1,2)*Power(x2,8) + 35*Power(x2,10)) + (479001600*(-1 + \
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
13*Power(x2,12))))/(2.*Power(lambdasq,12)*Pi*Power(Power(x1,2) + \
Power(x2,2),13));
	else if ( (k1==10) && (k2==3) )
		return (he0*x1*x2*(Power(x1,10)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),13) - lambdasq*Power(x1,8)*Power(Power(x1,2) + \
Power(x2,2),12)*(3*Power(x1,4) + 32*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) - (12454041600*(-1 + \
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
Power(x2,6)) + 3*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(31*Power(x1,6) + 184*Power(x1,4)*Power(x2,2) + \
275*Power(x1,2)*Power(x2,4) + 330*Power(x2,6)) - \
66*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) + \
Power(x2,2),10)*(14*Power(x1,8) + 56*Power(x1,6)*Power(x2,2) + \
175*Power(x1,4)*Power(x2,4) + 30*Power(x1,2)*Power(x2,6) + \
105*Power(x2,8)) + \
165*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),9)*(14*Power(x1,10) + 161*Power(x1,8)*Power(x2,2) - \
224*Power(x1,6)*Power(x2,4) + 894*Power(x1,4)*Power(x2,6) - \
294*Power(x1,2)*Power(x2,8) + 105*Power(x2,10)) - \
1485*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),8)*(7*Power(x1,12) - 140*Power(x1,10)*Power(x2,2) + \
833*Power(x1,8)*Power(x2,4) - 1368*Power(x1,6)*Power(x2,6) + \
833*Power(x1,4)*Power(x2,8) - 140*Power(x1,2)*Power(x2,10) + \
7*Power(x2,12))))/(2.*Power(lambdasq,13)*Pi*Power(Power(x1,2) + \
Power(x2,2),14));
	else if ( (k1==10) && (k2==4) )
		return -(he0*x1*(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),14) - lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) \
+ Power(x2,2),13)*(6*Power(x1,4) + 33*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) + (87178291200*(-1 + \
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
Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(3*Power(x1,8) + 180*Power(x1,6)*Power(x2,2) + \
795*Power(x1,4)*Power(x2,4) + 880*Power(x1,2)*Power(x2,6) + \
990*Power(x2,8)) - 3*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),11)*(31*Power(x1,10) + 705*Power(x1,8)*Power(x2,2) + \
2055*Power(x1,6)*Power(x2,4) + 5885*Power(x1,4)*Power(x2,6) + \
990*Power(x1,2)*Power(x2,8) + 2310*Power(x2,10)) + \
231*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),10)*(4*Power(x1,12) + 30*Power(x1,10)*Power(x2,2) + \
285*Power(x1,8)*Power(x2,4) - 250*Power(x1,6)*Power(x2,6) + \
900*Power(x1,4)*Power(x2,8) - 180*Power(x1,2)*Power(x2,10) + \
75*Power(x2,12)) + 10395*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 84*Power(x1,12)*Power(x2,2) + \
1113*Power(x1,10)*Power(x2,4) - 4060*Power(x1,8)*Power(x2,6) + \
5235*Power(x1,6)*Power(x2,8) - 2436*Power(x1,4)*Power(x2,10) + \
371*Power(x1,2)*Power(x2,12) - 12*Power(x2,14)) - \
1155*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),9)*(2*Power(x1,14) + 60*Power(x1,12)*Power(x2,2) - \
501*Power(x1,10)*Power(x2,4) + 2125*Power(x1,8)*Power(x2,6) - \
2520*Power(x1,6)*Power(x2,8) + 1278*Power(x1,4)*Power(x2,10) - \
165*Power(x1,2)*Power(x2,12) + \
9*Power(x2,14))))/(2.*Power(lambdasq,14)*Pi*Power(Power(x1,2) + \
Power(x2,2),15));
	else if ( (k1==10) && (k2==5) )
		return (he0*x1*x2*(Power(x1,10)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),15) - \
5*lambdasq*Power(x1,8)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(2*Power(x1,4) + 7*Power(x1,2)*Power(x2,2) + \
11*Power(x2,4)) + 15*Power(lambdasq,2)*Power(x1,6)*Power(Power(x1,2) \
+ Power(x2,2),13)*(Power(x1,8) + 20*Power(x1,6)*Power(x2,2) + \
75*Power(x1,4)*Power(x2,4) + 66*Power(x1,2)*Power(x2,6) + \
66*Power(x2,8)) - 15*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),12)*(29*Power(x1,10) + 305*Power(x1,8)*Power(x2,2) + \
657*Power(x1,6)*Power(x2,4) + 1705*Power(x1,4)*Power(x2,6) + \
330*Power(x1,2)*Power(x2,8) + 462*Power(x2,10)) + \
315*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),11)*(14*Power(x1,12) + 60*Power(x1,10)*Power(x2,2) + \
467*Power(x1,8)*Power(x2,4) - 220*Power(x1,6)*Power(x2,6) + \
880*Power(x1,4)*Power(x2,8) - 88*Power(x1,2)*Power(x2,10) + \
55*Power(x2,12)) + 17325*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,14) - 224*Power(x1,12)*Power(x2,2) + \
1785*Power(x1,10)*Power(x2,4) - 4630*Power(x1,8)*Power(x2,6) + \
4665*Power(x1,6)*Power(x2,8) - 1764*Power(x1,4)*Power(x2,10) + \
231*Power(x1,2)*Power(x2,12) - 6*Power(x2,14)) + (20922789888000*(-1 \
+ he0)*Power(lambdasq,15)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) \
+ 273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)))/he0 + \
2615348736000*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
54486432000*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
5448643200*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
454053600*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
32432400*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) + \
2027025*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,14) - 35*Power(x1,12)*Power(x2,2) + \
273*Power(x1,10)*Power(x2,4) - 715*Power(x1,8)*Power(x2,6) + \
715*Power(x1,6)*Power(x2,8) - 273*Power(x1,4)*Power(x2,10) + \
35*Power(x1,2)*Power(x2,12) - Power(x2,14)) - \
3465*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),10)*(2*Power(x1,14) + 80*Power(x1,12)*Power(x2,2) - \
381*Power(x1,10)*Power(x2,4) + 1255*Power(x1,8)*Power(x2,6) - \
1080*Power(x1,6)*Power(x2,8) + 486*Power(x1,4)*Power(x2,10) - \
45*Power(x1,2)*Power(x2,12) + 3*Power(x2,14)) + \
10461394944000*Power(lambdasq,14)*(Power(x1,16) - \
34*Power(x1,14)*Power(x2,2) + 238*Power(x1,12)*Power(x2,4) - \
442*Power(x1,10)*Power(x2,6) + 442*Power(x1,6)*Power(x2,10) - \
238*Power(x1,4)*Power(x2,12) + 34*Power(x1,2)*Power(x2,14) - \
Power(x2,16)) + 435891456000*Power(lambdasq,12)*(Power(x1,20) - \
32*Power(x1,18)*Power(x2,2) + 171*Power(x1,16)*Power(x2,4) - \
646*Power(x1,12)*Power(x2,8) + 646*Power(x1,8)*Power(x2,12) - \
171*Power(x1,4)*Power(x2,16) + 32*Power(x1,2)*Power(x2,18) - \
Power(x2,20))))/(2.*Power(lambdasq,15)*Pi*Power(Power(x1,2) + \
Power(x2,2),16));
	else if ( (k1==10) && (k2==6) )
		return (he0*x1*(-(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),16)) + lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),15)*(15*Power(x1,4) + 38*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) - \
15*Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(3*Power(x1,8) + 31*Power(x1,6)*Power(x2,2) + \
103*Power(x1,4)*Power(x2,4) + 77*Power(x1,2)*Power(x2,6) + \
66*Power(x2,8)) + 15*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),13)*(Power(x1,12) + 84*Power(x1,10)*Power(x2,2) + \
620*Power(x1,8)*Power(x2,4) + 1044*Power(x1,6)*Power(x2,6) + \
2365*Power(x1,4)*Power(x2,8) + 528*Power(x1,2)*Power(x2,10) + \
462*Power(x2,12)) - \
15*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),12)*(29*Power(x1,14) + 1006*Power(x1,12)*Power(x2,2) + \
3020*Power(x1,10)*Power(x2,4) + 19771*Power(x1,8)*Power(x2,6) - \
3355*Power(x1,6)*Power(x2,8) + 23892*Power(x1,4)*Power(x2,10) - \
462*Power(x1,2)*Power(x2,12) + 1155*Power(x2,14)) - \
3465*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),10)*(2*Power(x1,16) + 253*Power(x1,14)*Power(x2,2) - \
3745*Power(x1,12)*Power(x2,4) + 20377*Power(x1,10)*Power(x2,6) - \
39145*Power(x1,8)*Power(x2,8) + 31911*Power(x1,6)*Power(x2,10) - \
9891*Power(x1,4)*Power(x2,12) + 1155*Power(x1,2)*Power(x2,14) - \
21*Power(x2,16)) + (20922789888000*(-1 + \
he0)*Power(lambdasq,16)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) \
+ 2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)))/he0 + \
10461394944000*Power(lambdasq,15)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2615348736000*Power(lambdasq,14)*Power(Power(x1,2) \
+ Power(x2,2),2)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 435891456000*Power(lambdasq,13)*Power(Power(x1,2) \
+ Power(x2,2),3)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 54486432000*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 5448643200*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 454053600*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 32432400*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 2027025*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,16) - 136*Power(x1,14)*Power(x2,2) + \
2380*Power(x1,12)*Power(x2,4) - 12376*Power(x1,10)*Power(x2,6) + \
24310*Power(x1,8)*Power(x2,8) - 19448*Power(x1,6)*Power(x2,10) + \
6188*Power(x1,4)*Power(x2,12) - 680*Power(x1,2)*Power(x2,14) + \
17*Power(x2,16)) + 315*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),11)*(14*Power(x1,16) + 76*Power(x1,14)*Power(x2,2) + \
2795*Power(x1,12)*Power(x2,4) - 8066*Power(x1,10)*Power(x2,6) + \
22385*Power(x1,8)*Power(x2,8) - 13728*Power(x1,6)*Power(x2,10) + \
5973*Power(x1,4)*Power(x2,12) - 330*Power(x1,2)*Power(x2,14) + \
33*Power(x2,16)) + 17325*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),9)*(7*Power(x1,16) - 880*Power(x1,14)*Power(x2,2) + \
15484*Power(x1,12)*Power(x2,4) - 80416*Power(x1,10)*Power(x2,6) + \
158050*Power(x1,8)*Power(x2,8) - 126384*Power(x1,6)*Power(x2,10) + \
40236*Power(x1,4)*Power(x2,12) - 4416*Power(x1,2)*Power(x2,14) + \
111*Power(x2,16))))/(2.*Power(lambdasq,16)*Pi*Power(Power(x1,2) + \
Power(x2,2),17));
	else if ( (k1==10) && (k2==7) )
		return (he0*x1*x2*(Power(x1,10)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),17) - lambdasq*Power(x1,8)*Power(x2,4)*Power(Power(x1,2) \
+ Power(x2,2),16)*(21*Power(x1,4) + 42*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) - (711374856192000*(-1 + \
he0)*Power(lambdasq,17)*(Power(x1,2) - 3*Power(x2,2))*(3*Power(x1,2) \
- Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)))/he0 - 355687428096000*Power(lambdasq,16)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*(Power(x1,2) + \
Power(x2,2))*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 88921857024000*Power(lambdasq,15)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 14820309504000*Power(lambdasq,14)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 1852538688000*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 185253868800*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 15437822400*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 1102701600*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 68918850*Power(lambdasq,9)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) - 3828825*Power(lambdasq,8)*(Power(x1,2) - \
3*Power(x2,2))*(3*Power(x1,2) - Power(x2,2))*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,6) - 33*Power(x1,4)*Power(x2,2) + \
27*Power(x1,2)*Power(x2,4) - 3*Power(x2,6))*(3*Power(x1,6) - \
27*Power(x1,4)*Power(x2,2) + 33*Power(x1,2)*Power(x2,4) - \
Power(x2,6)) + \
Power(lambdasq,2)*Power(x1,6)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(105*Power(x1,8) + 693*Power(x1,6)*Power(x2,2) + \
2061*Power(x1,4)*Power(x2,4) + 1375*Power(x1,2)*Power(x2,6) + \
990*Power(x2,8)) - 15*Power(lambdasq,3)*Power(x1,4)*Power(Power(x1,2) \
+ Power(x2,2),14)*(7*Power(x1,12) + 196*Power(x1,10)*Power(x2,2) + \
1176*Power(x1,8)*Power(x2,4) + 1660*Power(x1,6)*Power(x2,6) + \
3179*Power(x1,4)*Power(x2,8) + 792*Power(x1,2)*Power(x2,10) + \
462*Power(x2,12)) + \
105*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),13)*(27*Power(x1,14) + 426*Power(x1,12)*Power(x2,2) + \
972*Power(x1,10)*Power(x2,4) + 5229*Power(x1,8)*Power(x2,6) + \
275*Power(x1,6)*Power(x2,8) + 4356*Power(x1,4)*Power(x2,10) + \
198*Power(x1,2)*Power(x2,12) + 165*Power(x2,14)) + \
315*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),11)*(6*Power(x1,16) + 13203*Power(x1,14)*Power(x2,2) - \
111951*Power(x1,12)*Power(x2,4) + 451623*Power(x1,10)*Power(x2,6) - \
657415*Power(x1,8)*Power(x2,8) + 450153*Power(x1,6)*Power(x2,10) - \
113421*Power(x1,4)*Power(x2,12) + 12573*Power(x1,2)*Power(x2,14) - \
99*Power(x2,16)) - 105*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),12)*(288*Power(x1,16) + 804*Power(x1,14)*Power(x2,2) + \
22743*Power(x1,12)*Power(x2,4) - 39330*Power(x1,10)*Power(x2,6) + \
102245*Power(x1,8)*Power(x2,8) - 41976*Power(x1,6)*Power(x2,10) + \
20097*Power(x1,4)*Power(x2,12) - 330*Power(x1,2)*Power(x2,14) + \
99*Power(x2,16)) - 3465*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),10)*(513*Power(x1,16) - 22416*Power(x1,14)*Power(x2,2) + \
237132*Power(x1,12)*Power(x2,4) - 878256*Power(x1,10)*Power(x2,6) + \
1344230*Power(x1,8)*Power(x2,8) - 878256*Power(x1,6)*Power(x2,10) + \
237132*Power(x1,4)*Power(x2,12) - 22416*Power(x1,2)*Power(x2,14) + \
513*Power(x2,16))))/(2.*Power(lambdasq,17)*Pi*Power(Power(x1,2) + \
Power(x2,2),18));
	else if ( (k1==10) && (k2==8) )
		return -(he0*x1*(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),18) - lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) \
+ Power(x2,2),17)*(28*Power(x1,4) + 47*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) + \
6*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),16)*(35*Power(x1,8) + 168*Power(x1,6)*Power(x2,2) + \
447*Power(x1,4)*Power(x2,4) + 275*Power(x1,2)*Power(x2,6) + \
165*Power(x2,8)) - \
6*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),15)*(70*Power(x1,12) + 1015*Power(x1,10)*Power(x2,2) + \
5229*Power(x1,8)*Power(x2,4) + 6556*Power(x1,6)*Power(x2,6) + \
10450*Power(x1,4)*Power(x2,8) + 2805*Power(x1,2)*Power(x2,10) + \
1155*Power(x2,12)) + \
15*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) + \
Power(x2,2),14)*(7*Power(x1,16) + 728*Power(x1,14)*Power(x2,2) + \
7882*Power(x1,12)*Power(x2,4) + 14896*Power(x1,10)*Power(x2,6) + \
63554*Power(x1,8)*Power(x2,8) + 13816*Power(x1,6)*Power(x2,10) + \
38874*Power(x1,4)*Power(x2,12) + 3696*Power(x1,2)*Power(x2,14) + \
1155*Power(x2,16)) + 31185*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),10)*(57*Power(x1,18) - 9432*Power(x1,16)*Power(x2,2) + \
214212*Power(x1,14)*Power(x2,4) - 1498896*Power(x1,12)*Power(x2,6) + \
4176126*Power(x1,10)*Power(x2,8) - 5103664*Power(x1,8)*Power(x2,10) + \
2784084*Power(x1,6)*Power(x2,12) - 642384*Power(x1,4)*Power(x2,14) + \
53553*Power(x1,2)*Power(x2,16) - 1048*Power(x2,18)) + \
(6402373705728000*(-1 + he0)*Power(lambdasq,18)*(Power(x1,18) - \
171*Power(x1,16)*Power(x2,2) + 3876*Power(x1,14)*Power(x2,4) - \
27132*Power(x1,12)*Power(x2,6) + 75582*Power(x1,10)*Power(x2,8) - \
92378*Power(x1,8)*Power(x2,10) + 50388*Power(x1,6)*Power(x2,12) - \
11628*Power(x1,4)*Power(x2,14) + 969*Power(x1,2)*Power(x2,16) - \
19*Power(x2,18)))/he0 + \
3201186852864000*Power(lambdasq,17)*(Power(x1,2) + \
Power(x2,2))*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
800296713216000*Power(lambdasq,16)*Power(Power(x1,2) + \
Power(x2,2),2)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
133382785536000*Power(lambdasq,15)*Power(Power(x1,2) + \
Power(x2,2),3)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
16672848192000*Power(lambdasq,14)*Power(Power(x1,2) + \
Power(x2,2),4)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
1667284819200*Power(lambdasq,13)*Power(Power(x1,2) + \
Power(x2,2),5)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
138940401600*Power(lambdasq,12)*Power(Power(x1,2) + \
Power(x2,2),6)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
9924314400*Power(lambdasq,11)*Power(Power(x1,2) + \
Power(x2,2),7)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
620269650*Power(lambdasq,10)*Power(Power(x1,2) + \
Power(x2,2),8)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) + \
34459425*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),9)*(Power(x1,18) - 171*Power(x1,16)*Power(x2,2) + \
3876*Power(x1,14)*Power(x2,4) - 27132*Power(x1,12)*Power(x2,6) + \
75582*Power(x1,10)*Power(x2,8) - 92378*Power(x1,8)*Power(x2,10) + \
50388*Power(x1,6)*Power(x2,12) - 11628*Power(x1,4)*Power(x2,14) + \
969*Power(x1,2)*Power(x2,16) - 19*Power(x2,18)) - \
945*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),13)*(3*Power(x1,18) + 147*Power(x1,16)*Power(x2,2) + \
298*Power(x1,14)*Power(x2,4) + 6054*Power(x1,12)*Power(x2,6) - \
5838*Power(x1,10)*Power(x2,8) + 16654*Power(x1,8)*Power(x2,10) - \
3894*Power(x1,6)*Power(x2,12) + 2574*Power(x1,4)*Power(x2,14) + \
55*Power(x1,2)*Power(x2,16) + 11*Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),12)*(16*Power(x1,18) - 41*Power(x1,16)*Power(x2,2) + \
8116*Power(x1,14)*Power(x2,4) - 42798*Power(x1,12)*Power(x2,6) + \
137318*Power(x1,10)*Power(x2,8) - 152262*Power(x1,8)*Power(x2,10) + \
91872*Power(x1,6)*Power(x2,12) - 18062*Power(x1,4)*Power(x2,14) + \
2134*Power(x1,2)*Power(x2,16) + 11*Power(x2,18)) - \
1890*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),11)*(Power(x1,18) + 7529*Power(x1,16)*Power(x2,2) - \
158594*Power(x1,14)*Power(x2,4) + 1129562*Power(x1,12)*Power(x2,6) - \
3123152*Power(x1,10)*Power(x2,8) + 3837548*Power(x1,8)*Power(x2,10) - \
2080958*Power(x1,6)*Power(x2,12) + 485078*Power(x1,4)*Power(x2,14) - \
39281*Power(x1,2)*Power(x2,16) + \
891*Power(x2,18))))/(2.*Power(lambdasq,18)*Pi*Power(Power(x1,2) + \
Power(x2,2),19));
	else if ( (k1==10) && (k2==9) )
		return (he0*x1*x2*(Power(x1,10)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),19) + (486580401635328000*(-1 + \
he0)*Power(lambdasq,19)*(x1 - x2)*(x1 + x2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)))/he0 + 243290200817664000*Power(lambdasq,18)*(x1 - \
x2)*(x1 + x2)*(Power(x1,2) + Power(x2,2))*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 60822550204416000*Power(lambdasq,17)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),2)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 10137091700736000*Power(lambdasq,16)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),3)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 1267136462592000*Power(lambdasq,15)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),4)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 126713646259200*Power(lambdasq,14)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),5)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 10559470521600*Power(lambdasq,13)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),6)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 754247894400*Power(lambdasq,12)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),7)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 47140493400*Power(lambdasq,11)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),8)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 2618916300*Power(lambdasq,10)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),9)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) + 130945815*Power(lambdasq,9)*(x1 - x2)*(x1 + \
x2)*Power(Power(x1,2) + Power(x2,2),10)*(5*Power(x1,4) - \
10*Power(x1,2)*Power(x2,2) + Power(x2,4))*(Power(x1,4) - \
4*Power(x1,3)*x2 - 14*Power(x1,2)*Power(x2,2) - 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) + 4*Power(x1,3)*x2 - \
14*Power(x1,2)*Power(x2,2) + 4*x1*Power(x2,3) + \
Power(x2,4))*(Power(x1,4) - 10*Power(x1,2)*Power(x2,2) + \
5*Power(x2,4)) - lambdasq*Power(x1,8)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),18)*(36*Power(x1,4) + 53*Power(x1,2)*Power(x2,2) + \
55*Power(x2,4)) + \
18*Power(lambdasq,2)*Power(x1,6)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),17)*(21*Power(x1,8) + 80*Power(x1,6)*Power(x2,2) + \
190*Power(x1,4)*Power(x2,4) + 110*Power(x1,2)*Power(x2,6) + \
55*Power(x2,8)) - \
18*Power(lambdasq,3)*Power(x1,4)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),16)*(70*Power(x1,12) + 651*Power(x1,10)*Power(x2,2) + \
2935*Power(x1,8)*Power(x2,4) + 3400*Power(x1,6)*Power(x2,6) + \
4510*Power(x1,4)*Power(x2,8) + 1265*Power(x1,2)*Power(x2,10) + \
385*Power(x2,12)) + 9*Power(lambdasq,4)*Power(x1,2)*Power(Power(x1,2) \
+ Power(x2,2),15)*(105*Power(x1,16) + 3640*Power(x1,14)*Power(x2,2) + \
31346*Power(x1,12)*Power(x2,4) + 52820*Power(x1,10)*Power(x2,6) + \
175590*Power(x1,8)*Power(x2,8) + 59400*Power(x1,6)*Power(x2,10) + \
83050*Power(x1,4)*Power(x2,12) + 10780*Power(x1,2)*Power(x2,14) + \
1925*Power(x2,16)) + 2835*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),11)*(10655*Power(x1,18) - 596940*Power(x1,16)*Power(x2,2) \
+ 8143332*Power(x1,14)*Power(x2,4) - \
40675080*Power(x1,12)*Power(x2,6) + 88177850*Power(x1,10)*Power(x2,8) \
- 88138160*Power(x1,8)*Power(x2,10) + \
40701540*Power(x1,6)*Power(x2,12) - 8131992*Power(x1,4)*Power(x2,14) \
+ 599775*Power(x1,2)*Power(x2,16) - 10340*Power(x2,18)) + \
1890*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),12)*(305*Power(x1,18) - 40485*Power(x1,16)*Power(x2,2) + \
496542*Power(x1,14)*Power(x2,4) - 2570490*Power(x1,12)*Power(x2,6) + \
5469680*Power(x1,10)*Power(x2,8) - 5549060*Power(x1,8)*Power(x2,10) + \
2517570*Power(x1,6)*Power(x2,12) - 519222*Power(x1,4)*Power(x2,14) + \
34815*Power(x1,2)*Power(x2,16) - 935*Power(x2,18)) + \
1890*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),13)*(145*Power(x1,18) - 180*Power(x1,16)*Power(x2,2) + \
25086*Power(x1,14)*Power(x2,4) - 84312*Power(x1,12)*Power(x2,6) + \
232264*Power(x1,10)*Power(x2,8) - 192280*Power(x1,8)*Power(x2,10) + \
109890*Power(x1,6)*Power(x2,12) - 15048*Power(x1,4)*Power(x2,14) + \
2343*Power(x1,2)*Power(x2,16) + 44*Power(x2,18)) - \
135*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),14)*(175*Power(x1,18) + 3885*Power(x1,16)*Power(x2,2) + \
6762*Power(x1,14)*Power(x2,4) + 92250*Power(x1,12)*Power(x2,6) - \
39638*Power(x1,10)*Power(x2,8) + 168410*Power(x1,8)*Power(x2,10) - \
12870*Power(x1,6)*Power(x2,12) + 21714*Power(x1,4)*Power(x2,14) + \
1155*Power(x1,2)*Power(x2,16) + \
77*Power(x2,18))))/(2.*Power(lambdasq,19)*Pi*Power(Power(x1,2) + \
Power(x2,2),20));
	else if ( (k1==10) && (k2==10) )
		return (he0*x1*(-(Power(x1,10)*Power(x2,10)*Power(Power(x1,2) + \
Power(x2,2),20)) + \
5*lambdasq*Power(x1,8)*Power(x2,8)*Power(Power(x1,2) + \
Power(x2,2),19)*(9*Power(x1,4) + 12*Power(x1,2)*Power(x2,2) + \
11*Power(x2,4)) - \
5*Power(lambdasq,2)*Power(x1,6)*Power(x2,6)*Power(Power(x1,2) + \
Power(x2,2),18)*(126*Power(x1,8) + 405*Power(x1,6)*Power(x2,2) + \
858*Power(x1,4)*Power(x2,4) + 473*Power(x1,2)*Power(x2,6) + \
198*Power(x2,8)) + (2432902008176640000*(-1 + \
he0)*Power(lambdasq,20)*(Power(x1,2) - 3*Power(x2,2))*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)))/he0 + \
1216451004088320000*Power(lambdasq,19)*(Power(x1,2) - \
3*Power(x2,2))*(Power(x1,2) + Power(x2,2))*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 304112751022080000*Power(lambdasq,18)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),2)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 50685458503680000*Power(lambdasq,17)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),3)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 6335682312960000*Power(lambdasq,16)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),4)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 633568231296000*Power(lambdasq,15)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),5)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 52797352608000*Power(lambdasq,14)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),6)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 3771239472000*Power(lambdasq,13)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),7)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 235702467000*Power(lambdasq,12)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),8)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 13094581500*Power(lambdasq,11)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),9)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + 654729075*Power(lambdasq,10)*(Power(x1,2) - \
3*Power(x2,2))*Power(Power(x1,2) + Power(x2,2),10)*(Power(x1,6) - \
21*Power(x1,4)*Power(x2,2) + 35*Power(x1,2)*Power(x2,4) - \
7*Power(x2,6))*(Power(x1,12) - 186*Power(x1,10)*Power(x2,2) + \
1423*Power(x1,8)*Power(x2,4) - 1772*Power(x1,6)*Power(x2,6) + \
655*Power(x1,4)*Power(x2,8) - 58*Power(x1,2)*Power(x2,10) + \
Power(x2,12)) + \
90*Power(lambdasq,3)*Power(x1,4)*Power(x2,4)*Power(Power(x1,2) + \
Power(x2,2),17)*(35*Power(x1,12) + 238*Power(x1,10)*Power(x2,2) + \
945*Power(x1,8)*Power(x2,4) + 1036*Power(x1,6)*Power(x2,6) + \
1155*Power(x1,4)*Power(x2,8) + 330*Power(x1,2)*Power(x2,10) + \
77*Power(x2,12)) - \
45*Power(lambdasq,4)*Power(x1,2)*Power(x2,2)*Power(Power(x1,2) + \
Power(x2,2),16)*(105*Power(x1,16) + 1890*Power(x1,14)*Power(x2,2) + \
13706*Power(x1,12)*Power(x2,4) + 21630*Power(x1,10)*Power(x2,6) + \
56322*Power(x1,8)*Power(x2,8) + 23870*Power(x1,6)*Power(x2,10) + \
21450*Power(x1,4)*Power(x2,12) + 3234*Power(x1,2)*Power(x2,14) + \
385*Power(x2,16)) + 9450*Power(lambdasq,7)*Power(Power(x1,2) + \
Power(x2,2),13)*(29*Power(x1,20) - 546*Power(x1,18)*Power(x2,2) + \
33579*Power(x1,16)*Power(x2,4) - 262500*Power(x1,14)*Power(x2,6) + \
1050210*Power(x1,12)*Power(x2,8) - 1749216*Power(x1,10)*Power(x2,10) \
+ 1510894*Power(x1,8)*Power(x2,12) - 570636*Power(x1,6)*Power(x2,14) \
+ 108801*Power(x1,4)*Power(x2,16) - 5390*Power(x1,2)*Power(x2,18) + \
231*Power(x2,20)) + 45*Power(lambdasq,5)*Power(Power(x1,2) + \
Power(x2,2),15)*(21*Power(x1,20) + 2520*Power(x1,18)*Power(x2,2) + \
37905*Power(x1,16)*Power(x2,4) + 62888*Power(x1,14)*Power(x2,6) + \
561120*Power(x1,12)*Power(x2,8) - 23352*Power(x1,10)*Power(x2,10) + \
733040*Power(x1,8)*Power(x2,12) + 43560*Power(x1,6)*Power(x2,14) + \
82467*Power(x1,4)*Power(x2,16) + 6160*Power(x1,2)*Power(x2,18) + \
231*Power(x2,20)) - 675*Power(lambdasq,6)*Power(Power(x1,2) + \
Power(x2,2),14)*(35*Power(x1,20) + 2352*Power(x1,18)*Power(x2,2) - \
735*Power(x1,16)*Power(x2,4) + 189924*Power(x1,14)*Power(x2,6) - \
399672*Power(x1,12)*Power(x2,8) + 1044624*Power(x1,10)*Power(x2,10) - \
605528*Power(x1,8)*Power(x2,12) + 370260*Power(x1,6)*Power(x2,14) - \
25179*Power(x1,4)*Power(x2,16) + 8008*Power(x1,2)*Power(x2,18) + \
231*Power(x2,20)) + 4725*Power(lambdasq,8)*Power(Power(x1,2) + \
Power(x2,2),12)*(122*Power(x1,20) - 56805*Power(x1,18)*Power(x2,2) + \
1561770*Power(x1,16)*Power(x2,4) - 14268156*Power(x1,14)*Power(x2,6) \
+ 53349660*Power(x1,12)*Power(x2,8) - \
92637678*Power(x1,10)*Power(x2,10) + \
77072380*Power(x1,8)*Power(x2,12) - 30557340*Power(x1,6)*Power(x2,14) \
+ 5323626*Power(x1,4)*Power(x2,16) - 353045*Power(x1,2)*Power(x2,18) \
+ 5082*Power(x2,20)) + 14175*Power(lambdasq,9)*Power(Power(x1,2) + \
Power(x2,2),11)*(2131*Power(x1,20) - 440580*Power(x1,18)*Power(x2,2) \
+ 12566925*Power(x1,16)*Power(x2,4) - \
113923488*Power(x1,14)*Power(x2,6) + \
427233870*Power(x1,12)*Power(x2,8) - \
740519304*Power(x1,10)*Power(x2,10) + \
617112650*Power(x1,8)*Power(x2,12) - \
244126080*Power(x1,6)*Power(x2,14) + \
42724143*Power(x1,4)*Power(x2,16) - 2792020*Power(x1,2)*Power(x2,18) \
+ 44121*Power(x2,20))))/(2.*Power(lambdasq,20)*Pi*Power(Power(x1,2) + \
Power(x2,2),21));
	else {
#if HVS_DEBUG
		printf("Err:%d,%d\n",k1,k2);
		hvsdie("Error. Function hb2.\n");
#endif
		return 0.0;
	}
}

