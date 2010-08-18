#include <stdio.h>
#include "libhvs/hevosol.h"

int main() {
	hvs_state* state;
	hvs_params params;
	params.gridwidth = 100;
	params.gridheight = 100;
	params.timestep = 0.1;
	init_solver(&params, &state);
	if (run_solver(&params, state) == HVS_OK) {
		char string[100];
		printf("%u and %u\n", state->width, state->height);
		format_output(state,100,&string);
		printf("%s", string);
	}
	free_solver(&state);
	return 0;
}