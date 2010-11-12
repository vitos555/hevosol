#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include "libhvs/hevosol.h"
#include "libhvs/fileutil.h"
#include "libhvs/hermiteutil.h"

#define SQRT_NCENTERS 7
#define NCENTERS SQRT_NCENTERS*SQRT_NCENTERS

// Output the result to file
#define OUTFILE "main.out"

int main() {
	hvs_state* state;
	hvs_params params;
	int status,i;
	time_t starttime,endtime;
	FLOAT_TYPE x,y,l;
	params.initvortfile = "test"; // File with initial vorticity field, not used at the moment
	params.timestep = 0.001;
	params.lambda0 = 1.0;
	params.nu = 0.1;
	hvs_center centers[NCENTERS];
	hvs_moment moment1 = {1.0/NCENTERS,0.0,0.0,0.0,0.0,0.0}; // Initial moments
	hvs_moment moments[NCENTERS];
	// In the loop copy moment to the array moments
	for(i=0;i<NCENTERS;i++) {
		memcpy(moments[i],moment1,sizeof(hvs_moment));
		centers[i].x=-2.0+4.0/(SQRT_NCENTERS-1)*(i%(SQRT_NCENTERS));
		centers[i].y=-2.0+4.0/(SQRT_NCENTERS-1)*(i/(SQRT_NCENTERS));
		printf("Center[%d] = (%Lf,%Lf)\n",i,centers[i].x,centers[i].y);
	}
	
	// Get the time when the solver starts working
	starttime = time(NULL);
	
	// Initialize solver
	if ((status = init_solver_by_moments(&params, NCENTERS, centers, moments,
				-5.0, 5.0, 0.1, 
				-5.0, 5.0, 0.1, &state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 0;
	}
	write_params(&params,OUTFILE);
	
	// Integrate
	for (i=0;i<10;i++) {
		params.t0 = 0.1*i;
		params.t1 = 0.1*(i+1);
		if ((status = run_solver(&params, state)) == HVS_OK) {
			printf("Centers: (%.4Lf,%.4Lf),(%.4Lf,%.4Lf)\n", 
				state->centers[0].x,state->centers[0].y,
				state->centers[1].x,state->centers[1].y);
			// Write vorticity to output file
			if (i%10==0) append_vorticity(state,OUTFILE);
		} else {
			hvserror(status,"Run error");
			break;
		}
	}
	
	// Deinitialize solver
	free_solver(&state);
	endtime = time(NULL);
	printf("Time: %d sec\n",endtime-starttime);
	return 0;
}