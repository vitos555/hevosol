#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/errorutil.h"
#include "libhvs/matrixutil.h"

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
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
	if ((status = matrix_add(m1,m2,3,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
	if ((status = matrix_scalarmult((FLOAT_TYPE)3.0,m1,3,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
	if ((status = matrix_inv(m1,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
	if ((status = matrix_multiply(m1,3,3,m3,3,3,m2))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m2[0],m2[1],m2[2],m2[3],m2[4],m2[5],m2[6],m2[7],m2[8]);
	if ((status = matrix_inv(m2,3,m3))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
	if ((status = matrix_multiply(m2,3,3,m3,3,3,m1))!=HVS_OK) {
		hvserror(status, "Error");
		return 1;
	}
	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8]);
//	gmres(matrix_test,vect_y,vect_x,3,(FLOAT_TYPE)0.01,vector_test);
//	printf("Res=(%Lf,%Lf,%Lf)\n",vector_test[0],vector_test[1],vector_test[2]);
	
	return 0;
}