#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/errorutil.h"
#include "libhvs/matrixutil.h"

int main() {
	int status,i;
	FLOAT_TYPE x,y,l;
	
	FLOAT_TYPE m3[] = {1.0,2.0,3.0,0.0,1.0,2.0,0.0,2.0,3.0};
	FLOAT_TYPE m2[] = {0.0,2.0,3.0,1.0,1.0,2.0,1.0,2.0,3.0};
	FLOAT_TYPE m[] = {1.0,1.0,1.0,2.0,1.0,3.0};
	FLOAT_TYPE vy[] = {7.0,7.0,8.0};
	FLOAT_TYPE vect_x0[] = {0.0,0.0,0.0};
	FLOAT_TYPE vect_y[] = {1.0,1.0,1.0};
	FLOAT_TYPE vector_test[3];

	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		m3[0],m3[1],m3[2],m3[3],m3[4],m3[5],m3[6],m3[7],m3[8]);
//	printf("%Lf\t%Lf\n%Lf%Lf\n%Lf\t%Lf\n\n",
//		m[0],m[1],m[2],m[3],m[4],m[5]);
//	if((status = ols(m,vy,3,2,vector_test))!=HVS_OK) {
//		hvserror(status,"Error");
//		return 1;
//	}
//	printf("Res=(%Lf,%Lf,%Lf)\n",vector_test[0],vector_test[1],vector_test[2]);
	if((status = gmres(m3,vect_x0,vect_y,3,(FLOAT_TYPE)0.01,10,vector_test))!=HVS_OK) {
		hvserror(status,"Error");
		return 1;
	}
	printf("Res=(%Lf,%Lf,%Lf)\n",vector_test[0],vector_test[1],vector_test[2]);
	
	return 0;
}