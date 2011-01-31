#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/errorutil.h"
#include "libhvs/matrixutil.h"

#define LF(x) (long double)((x))

int main() {
	int status,i;
	FLOAT_TYPE x,y,l;
	
	FLOAT_TYPE m1[] = {1.0,2.0,3.0,0.0,1.0,2.0,0.0,2.0,3.0};
	FLOAT_TYPE m2[] = {0.0,2.0,3.0,1.0,1.0,2.0,1.0,2.0,3.0};
	FLOAT_TYPE m3[9];
	
	if ((status = matrix_multiply(m1,3,3,m2,3,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),
		LF(m3[3]),LF(m3[4]),LF(m3[5]),
		LF(m3[6]),LF(m3[7]),LF(m3[8]));
	if ((status = matrix_add(m1,m2,3,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),LF(m3[3]),LF(m3[4]),LF(m3[5]),LF(m3[6]),LF(m3[7]),LF(m3[8]));
	if ((status = matrix_scalarmult((FLOAT_TYPE)3.0,m1,3,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),LF(m3[3]),LF(m3[4]),LF(m3[5]),LF(m3[6]),LF(m3[7]),LF(m3[8]));
	if ((status = matrix_inv(m1,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),LF(m3[3]),LF(m3[4]),LF(m3[5]),LF(m3[6]),LF(m3[7]),LF(m3[8]));
	if ((status = matrix_multiply(m1,3,3,m3,3,3,m2))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m2[0]),LF(m2[1]),LF(m2[2]),LF(m2[3]),LF(m2[4]),LF(m2[5]),LF(m2[6]),LF(m2[7]),LF(m2[8]));
	if ((status = matrix_inv(m2,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),LF(m3[3]),LF(m3[4]),LF(m3[5]),LF(m3[6]),LF(m3[7]),LF(m3[8]));
	if ((status = matrix_multiply(m2,3,3,m3,3,3,m1))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m1[0]),LF(m1[1]),LF(m1[2]),LF(m1[3]),LF(m1[4]),LF(m1[5]),LF(m1[6]),LF(m1[7]),LF(m1[8]));
//	gmres(matrix_test,vect_y,vect_x,3,(FLOAT_TYPE)0.01,vector_test);
//	printf("Res=(%Lf,%Lf,%Lf)\n",vector_test[0],vector_test[1],vector_test[2]);
	
	return 0;
}