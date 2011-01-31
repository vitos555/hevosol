#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/errorutil.h"
#include "libhvs/matrixutil.h"

#define LF(x) (long double)((x))

int main() {
	int status,i;
	FLOAT_TYPE x,y,l;
	
	FLOAT_TYPE m3[] = {1.0,2.0,3.0,4.0,0.0,1.0,2.0,2.0,0.0,2.0,3.0,-1.0,1.0,-1.0,3.0,1.0};
	FLOAT_TYPE m2[] = {0.0,2.0,3.0,1.0,1.0,2.0,1.0,2.0,3.0};
	FLOAT_TYPE m[] = {1.0,1.0,1.0,2.0,1.0,3.0};
	FLOAT_TYPE vy[] = {7.0,7.0,8.0};
	FLOAT_TYPE vect_x0[] = {0.0,0.0,0.0,0.0};
	FLOAT_TYPE vect_y[] = {1.0,1.0,1.0,1.0};
	FLOAT_TYPE vector_test[3];

	printf("%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n%Lf\t%Lf\t%Lf\n\n",
		LF(m3[0]),LF(m3[1]),LF(m3[2]),
		LF(m3[3]),LF(m3[4]),LF(m3[5]),
		LF(m3[6]),LF(m3[7]),LF(m3[8]));
//	printf("%Lf\t%Lf\n%Lf%Lf\n%Lf\t%Lf\n\n",
//		m[0],m[1],m[2],m[3],m[4],m[5]);
//	if((status = ols(m,vy,3,2,vector_test))!=HVS_OK) {
//		hvserror(status,"Error");
//		return 1;
//	}
//	printf("Res=(%Lf,%Lf,%Lf)\n",vector_test[0],vector_test[1],vector_test[2]);
	if((status = gmres(m3,vect_x0,vect_y,4,(FLOAT_TYPE)0.001,10,vector_test))!=HVS_OK) {
		hvserror(status,"Error");
		return 1;
	}
	printf("Res=(%f,%f,%f,%f)\n",vector_test[0],vector_test[1],vector_test[2],vector_test[3]);
	
	return 0;
}