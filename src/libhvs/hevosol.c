//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hevosol.h"
#include "fileutil.h"

#if HVS_MOMENTS == 2
#include "hvs_two_moments.c"
#elif HVS_MOMENTS == 1
#include "hvs_one_moment.c"
#else
#include "hvs_simple.c"
#endif

int init_solver(const hvs_params *params, hvs_state **sstate) {
	hvs_state* state = (hvs_state *) malloc(sizeof(hvs_state));
	hvs_file* file;
	hvs_position *pos = NULL;
	hvs_vorticity *vort = NULL;
	int status = 0;
	UINT cursize = 0;
	if (state == NULL) 
		return HVS_ERR;
	// Init file handler
	if (initfile(params->initvortfile, &file)) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	
	// Read data
	do {
		cursize += status;
		if ((pos = (hvs_position *)realloc(pos, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_position)))==NULL) {
			free(state);
			state = NULL;
			return HVS_ERR;
		}
		if ((vort = (hvs_vorticity *)realloc(vort, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_vorticity)))==NULL) {
			free(state);
			state = NULL;
			return HVS_ERR;
		}
	} while ((status = readdata(file,HVS_READ_BLOCK_SIZE,
			pos+cursize*sizeof(hvs_position),
			vort+cursize*sizeof(hvs_vorticity))) == HVS_READ_BLOCK_SIZE);
	if (status < 0) {
		free(state);
		state = NULL;
		return status;
	}
	cursize += status;
	
	// Check if we have regular grid
	if (checkgrid(pos) != HVS_OK) {
		return HVS_ERR_IRREGULAR_GRID;
	}
	
	// Copy values to the state structure
	state->size = cursize;
	state->grid = pos;
	state->vorticity_field = vort;
	state->moments = NULL;
	
	// Initialize velocity field
	state->velocity_field = (hvs_vector *) malloc(cursize*sizeof(hvs_vector));
	if (state->velocity_field == NULL) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	memset(state->velocity_field, 0, cursize*sizeof(hvs_vector));
	*sstate = state;
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
	unsigned int i;
	char status;
	unsigned int stepsnum = floor((params->t1-params->t0)/params->timestep);
	for(i=0;i<stepsnum;i++) {
		if (!(status = step_solver(state, params->timestep))) {
			return status;
		}
	}
	return update_vorticity_field(state);
}

int write_output(const hvs_state *state, const char *filename) {
	return writedata(state, filename);
}

int checkgrid(hvs_position *pos, hvs_vorticity *vort) {
	return HVS_OK;
}