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
#include "hermiteutil.h"

#if NMOMENTS == 2
#include "hvs_two_moments.c"
#else
#include "hvs_simple.c"
#endif

int init_solver_by_moments(hvs_params *params, UINT ncenters, const hvs_center *centers, const hvs_moment *moments,
				FLOAT_TYPE xmin, FLOAT_TYPE xmax, FLOAT_TYPE xstep, 
				FLOAT_TYPE ymin, FLOAT_TYPE ymax, FLOAT_TYPE ystep, hvs_state **sstate) {
	hvs_state* state = (hvs_state *) malloc(sizeof(hvs_state));
	int i, j;
	if (state == NULL)
		return HVS_ERR;
	state->ncenters = ncenters;
	
	state->centers = (hvs_center *) malloc(sizeof(hvs_center)*ncenters);
	if(state->centers == NULL) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	if (memcpy(state->centers, centers, sizeof(hvs_center)*ncenters) == NULL) {
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	
	state->moments = (hvs_moment *) malloc(sizeof(hvs_moment)*ncenters);
	if (state->moments == NULL) {
		free(state->centers);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	
	state->xmin = xmin;
	state->xmax = xmax;
	state->xstep = xstep;
	state->ymin = ymin;
	state->ymax = ymax;
	state->ystep = ystep;
	state->sizex = floor(abs(xmax-xmin)/xstep);
	state->sizey = floor(abs(ymax-ymin)/ystep);
	state->size = state->sizex*state->sizey;
	
	state->grid = (hvs_position *) malloc(sizeof(hvs_position)*state->size);
	if(state->grid == NULL) {
		free(state->moments);
		free(state->centers);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	state->vorticity_field = (hvs_vorticity *) malloc(sizeof(hvs_vorticity)*state->size);
	if (state->vorticity_field == NULL) {
		free(state->moments);
		free(state->centers);
		free(state->grid);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	memset(state->vorticity_field, 0, sizeof(hvs_vorticity)*state->size);
	for (i=0;i<state->sizey;i++)
		for(j=0;j<state->sizex;j++) {
			state->grid[i*state->sizex+j].x = xmin+xstep*j;
			state->grid[i*state->sizex+j].y = ymin+ystep*i;
		}

	// Initialize velocity field
	state->velocity_field = (hvs_vector *) malloc(state->size*sizeof(hvs_vector));
	if (state->velocity_field == NULL) {
		free(state->moments);
		free(state->centers);
		free(state->grid);
		free(state->vorticity_field);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	memset(state->velocity_field, 0, state->size*sizeof(hvs_vector));

	// We have moments and thus we can update vorticity field
	update_vorticity_field(params, state);
	(*sstate) = state;
	return HVS_OK;
}

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
	
	state->grid = pos;
	state->vorticity_field = vort;
	state->moments = NULL;
	state->centers = NULL;
	
	// Check if we have regular grid
	if (checkgrid(pos) != HVS_OK) {
		free(pos);
		free(vort);
		free(state);
		return HVS_ERR_IRREGULAR_GRID;
	}
	
	// Copy values to the state structure
	state->size = cursize;
	
	state->xmin = pos[0].x;
	state->xmax = pos[cursize].x;
	state->ymin = pos[0].y;
	state->ymax = pos[cursize].y;
	state->ystep = pos[1].y-pos[0].y;
	state->sizey = floor(abs(state->ymax-state->ymin)/state->ystep);
	state->xstep = pos[state->sizey].x-pos[0].x;
	state->sizex = floor(abs(state->xmax-state->xmin)/state->xstep);
	
	if(cursize != (state->sizex*state->sizey)) {
		
	}

	// Initialize velocity field
	state->velocity_field = (hvs_vector *) malloc(cursize*sizeof(hvs_vector));
	if (state->velocity_field == NULL) {
		free(pos);
		free(vort);
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
	if (state != NULL) {
		if (state->grid != NULL) {
			free(state->grid);
			state->grid = NULL;
		}
		if (state->moments != NULL) {
			free(state->moments);
			state->moments = NULL;
		}
		if (state->centers != NULL) {
			free(state->centers);
			state->centers = NULL;
		}
		if (state->velocity_field != NULL) {
			free(state->velocity_field);
			state->velocity_field = NULL;
		}
		if (state->vorticity_field != NULL) {
			free(state->vorticity_field);
			state->vorticity_field = NULL;
		}
		free(state);
		(*sstate) = NULL;
	}
}

int run_solver(const hvs_params *params, hvs_state *state) {
	unsigned int i;
	char status;
	unsigned int stepsnum = floor((params->t1-params->t0)/params->timestep);
	for(i=0;i<stepsnum;i++) {
		if ((status = step_solver(state, params->timestep)) != HVS_OK) {
			return status;
		}
	}
	return update_vorticity_field(params, state);
}

int write_output(const hvs_state *state, const char *filename) {
	return writedata(state, filename);
}

int checkgrid(hvs_position *pos, hvs_vorticity *vort) {
	return HVS_OK;
}