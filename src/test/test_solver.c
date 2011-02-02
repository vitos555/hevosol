#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/hevosol.h"
#include "libhvs/fileutil.h"
#include "libhvs/hermiteutil.h"
#include "libhvs/vectorutil.h"

#define OUTFILE "test.out"
#define LF(x) (long double)((x))


int main() {
	hvs_state* state;
	hvs_params params;
	int status,i;
	FLOAT_TYPE x,y,l;
	
	// Initialize float type
#if HVS_FLOAT_TYPE==HVS_FLOAT
printf("Float\n");
#elif HVS_FLOAT_TYPE==HVS_DOUBLE
printf("Double\n");
#elif HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
printf("Long double");
#else
printf("Else\n");
#endif


	FLOAT_TYPE matrix_test[] = {1.0,2.0,3.0,0.0,1.0,2.0,0.0,2.0,3.0};
	FLOAT_TYPE vector_test[3],vect_x[]={-1.0,2.0,1.0},vect_y[]={1.0,0.0,0.0};
	FLOAT_TYPE res;
	
	params.initvortfile = "test";
	params.timestep = 0.001;
	params.lambda0 = 1.0;
	params.nu = 0.1;
	params.xmin = -5.0;
	params.xmax = 5.0;
	params.xstep = 0.1;
	params.ymin = -5.0;
	params.ymax = 5.0;
	params.ystep = 0.1;
	hvs_center centers[] = {{1.0,0.0},{0.0,0.0},{-1.0,0.0}};
	hvs_moment moment1 = {-2.0,0.0,0.0,0.0,0.0,0.0};
	hvs_moment moment2 = {1.0,0.0,0.0,0.0,0.0,0.0};
	hvs_moment moments[3];
	memcpy(moments[0],moment1,sizeof(hvs_moment));
	memcpy(moments[1],moment2,sizeof(hvs_moment));
	memcpy(moments[2],moment1,sizeof(hvs_moment));
	if ((status = init_solver_by_moments(&params, 3, centers, moments,
				&state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 0;
	}
	write_params(&params,OUTFILE);
	append_vorticity(state,OUTFILE);
	for (i=0;i<100;i++) {
		params.t0 = 0.1*i;
		params.t1 = 0.1*(i+1);
		if ((status = run_solver(&params, state)) == HVS_OK) {
			printf("Centers: (%.4Lf,%.4Lf),(%.4Lf,%.4Lf)\n", 
				LF(state->centers[0].x),LF(state->centers[0].y),
				LF(state->centers[1].x),LF(state->centers[1].y));
			printf("Moments: (%Lf,%Lf,%Lf), (%Lf,%Lf,%Lf)\n",
				LF(state->moments[0][MOM_INDEX(0,0)]),
				LF(state->moments[0][MOM_INDEX(0,2)]),
				LF(state->moments[0][MOM_INDEX(2,0)]),
				LF(state->moments[1][MOM_INDEX(0,0)]),
				LF(state->moments[1][MOM_INDEX(0,2)]),
				LF(state->moments[1][MOM_INDEX(2,0)]));
			if ((i+1)%10==0) {
				append_vorticity(state,OUTFILE);
//				append_centers(state,OUTFILE);
//				append_moments(state,OUTFILE);
			}
		} else {
			hvserror(status,"Run error");
			break;
		}
	}
	free_solver(&state);
	x=0.7;
	y=2.301;
	l=15.5;
	printf("hb1(%Lf,%Lf,%Lf,3,2)=%.16Lg\n",LF(x),LF(y),LF(l*l),LF(hb1(x,y,l*l,3,2)));
	printf("hb1(%Lf,%Lf,%Lf,2,3)=%.16Lg\n",LF(x),LF(y),LF(l*l),LF(hb1(x,y,l*l,2,3)));
	printf("hb1(%Lf,%Lf,%Lf,3,3)=%.16Lg\n",LF(x),LF(y),LF(l*l),LF(hb1(x,y,l*l,3,3)));
	return 0;
}