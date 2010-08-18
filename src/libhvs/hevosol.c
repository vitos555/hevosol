//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hevosol.h"

int init_solver(const hvs_params *params, hvs_state **sstate) {
	hvs_state* state = (hvs_state *) malloc(sizeof(hvs_state));
	if (state == NULL) 
		return HVS_ERR;
	// Init grid
	state->width = params->gridwidth;
	state->height = params->gridheight;
	UINT gridsize = params->gridwidth*params->gridheight;
	state->grid = (hvs_position *) malloc(gridsize*sizeof(hvs_position));
	if (state->grid == NULL) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	state->moments = NULL;
	state->velocity_field = (hvs_vector *) malloc(gridsize*sizeof(hvs_vector));
	if (state->velocity_field == NULL) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	memset(state->velocity_field, 0, gridsize*sizeof(hvs_vector));
	state->vorticity_field = NULL;
	*sstate = state;
	return HVS_OK;
}

int advect_points(hvs_state *state, FLOAT_TYPE timestep) {
	UINT i,j;
	for (i=0; i<state->width; i++)
		for (j=0; j<state->height; j++) {
			(state->grid[i+j*state->height]).x = timestep*state->velocity_field[i+j*state->width].xval;
			(state->grid[i+j*state->height]).y = timestep*state->velocity_field[i+j*state->width].yval;
		}
	return HVS_OK;
}

int update_moments(hvs_state *state, FLOAT_TYPE timestep) {
	return HVS_OK;
}

int update_vorticity_field(hvs_state *state) {
	return HVS_OK;
}

int update_velocity_field(hvs_state *state) {
	return HVS_OK;
}

void free_solver(hvs_state **sstate) {
	hvs_state *state = *sstate;
	if (state->grid != NULL) {
		free(state->grid);
		state->grid = NULL;
	}
	if (state->moments != NULL) {
		free(state->moments);
		state->moments = NULL;
	}
	if (state->velocity_field != NULL) {
		free(state->velocity_field);
		state->velocity_field = NULL;
	}
	if (state->vorticity_field != NULL) {
		free(state->vorticity_field);
		state->vorticity_field = NULL;
	}
}

int run_solver(const hvs_params *params, hvs_state *state) {
	return advect_points(state, params->timestep)
			& update_moments(state, params->timestep)
			& update_vorticity_field(state)
			& update_velocity_field(state);
}

int format_output(hvs_state *state, size_t size, char *output) {
	return snprintf(output, size, "Width: %u, height: %u.\n", state->width, state->height);
}
