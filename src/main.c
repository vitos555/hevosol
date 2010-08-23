#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "libhvs/hevosol.h"

int main() {
	hvs_state* state;
	hvs_params params;
	int status;
	params.initvortfile = "test";
	params.gridwidth = 100;
	params.gridheight = 100;
	params.timestep = 0.1;
	if ((status = init_solver(&params, &state)) != HVS_OK) {
		hvserror(status, "Init error");
		return 0;
	}
	if ((status = run_solver(&params, state)) == HVS_OK) {
		printf("Size: %i\n", state->size);
		write_output(state,"test.out");
	} else {
		hvserror(status,"Run error");
	}
	free_solver(&state);
	return 0;
}