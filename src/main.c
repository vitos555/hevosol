#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/hevosol.h"
#include "libhvs/hermiteutil.h"

int main() {
	hvs_state* state;
	hvs_params params;
	int status,i;
	FLOAT_TYPE x,y,l;
	params.initvortfile = "test";
	params.timestep = 0.001;
	params.lambda0 = 1.0;
	params.nu = 0.1;
	hvs_center centers[] = {{-1.0,0.0},{1.0,0.0}};
	hvs_moment moment1;
	moment1[0] = 1.0;
	moment1[1] = 0.0;
	moment1[2] = 0.0;
	moment1[3] = 0.0;
	moment1[4] = 0.0;
	moment1[5] = 0.0;
	hvs_moment moment2 = {1.0,0.0,0.0,0.0,0.0,0.0};
	hvs_moment moments[2];
	memcpy(moments[0],moment1,sizeof(hvs_moment));
	memcpy(moments[1],moment2,sizeof(hvs_moment));
	if ((status = init_solver_by_moments(&params, 2, centers, moments,
				-5.0, 5.0, 0.1, 
				-5.0, 5.0, 0.1, &state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 0;
	}
	for (i=0;i<20;i++) {
		params.t0 = 0.1*i;
		params.t1 = 0.1*(i+1);
		if ((status = run_solver(&params, state)) == HVS_OK) {
			printf("Centers: (%.4Lf,%.4Lf),(%.4Lf,%.4Lf)\n", 
				state->centers[0].x,state->centers[0].y,
				state->centers[1].x,state->centers[1].y);
			write_output(state,"test.out");
		} else {
			hvserror(status,"Run error");
			break;
		}
	}
	free_solver(&state);
	x=0.7;
	y=2.301;
	l=15.5;
	printf("hb1(%Lf,%Lf,%Lf,3,2)=%.16Lg\n",x,y,l*l,hb1(x,y,l*l,3,2));
	printf("hb1(%Lf,%Lf,%Lf,2,3)=%.16Lg\n",x,y,l*l,hb1(x,y,l*l,2,3));
	printf("hb1(%Lf,%Lf,%Lf,3,3)=%.16Lg\n",x,y,l*l,hb1(x,y,l*l,3,3));
	return 0;
}