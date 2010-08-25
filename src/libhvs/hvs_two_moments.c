//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

int step_solver(hvs_state *state, FLOAT_TYPE timestep) {
	int i;
	for (i=0; i<state->ncenters; i++) {
		state->moments[i] = rk4_solve(&eval_moments_eq, state->centers[i], timestep);
		state->centers[i] = rk4_solve(&eval_centers_eq, state->centers[i], timestep);
	}
	return HVS_OK;
}

int update_vorticity_field(hvs_state *state, hvs_params *params) {
	int i,j;
	FLOAT_TYPE sum;
	FLOAT_TYPE lambda_sq = params->lambda0*params->lambda0;
	for (i=0; i<state->size; i++) {
		sum = 0.0;
		for (j=0; j<state->ncenters; j++) {
			sum += 	state->moments[j][0]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,0)+
				state->moments[j][1]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,2)+
				state->moments[j][2]*he(state->grid[i].x,state->grid[i].y,lambda_sq,2,0);
		}
		state->vorticity_field[i] = sum;
	}
	return HVS_OK;
}

hvs_moment eval_moments_eq(hvs_position pos, FLOAT_TYPE timestep) {
	return 0.0;
}

hvs_vector eval_centers_eq(hvs_position pos, FLOAT_TYPE timestep) {
	return 0.0;
}