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
#include "factorialutil.h"

#if NMOMENTS <= 4
#include "hvs_multi_moments.c"
#else
#error "Currently implemented only 4 moments."
#endif

#ifdef HVS_PROFILE
#include <time.h>

Timings timings = {.init = 0, .integration=0, .init_moments=0, .init_coefs=0, .rk_step=0, .eval_equation=0, .vorticity_update=0 };
#endif

int init_solver_by_moments(const hvs_params *params, UINT ncenters, const hvs_centers centers, const hvs_moments moments,
				hvs_state **sstate) {
	hvs_state* state = (hvs_state *) malloc(sizeof(hvs_state));
	int i, j;
	if (state == NULL)
		return HVS_ERR;
	state->ncenters = ncenters;
	state->lambdasq = params->lambda0*params->lambda0;
	
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
	
	state->moments = (hvs_moments) malloc(sizeof(hvs_moment)*ncenters);
	if (state->moments == NULL) {
		free(state->centers);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	if(memcpy(state->moments, moments, sizeof(hvs_moment)*ncenters) == NULL) {
		free(state->centers);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	
	state->xmin = params->xmin;
	state->xmax = params->xmax;
	state->xstep = params->xstep;
	state->ymin = params->ymin;
	state->ymax = params->ymax;
	state->ystep = params->ystep;
	state->sizex = floor(abs(state->xmax-state->xmin)/state->xstep)+1;
	state->sizey = floor(abs(state->ymax-state->ymin)/state->ystep)+1;
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
			state->grid[i*state->sizex+j].x = state->xmin+state->xstep*j;
			state->grid[i*state->sizex+j].y = state->ymin+state->ystep*i;
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

	// Initialize coefficients
	state->coefs = (hvs_coefs *) malloc(sizeof(hvs_coefs));
	if (state->coefs == NULL) {
		free(state->moments);
		free(state->centers);
		free(state->grid);
		free(state->vorticity_field);
		free(state->velocity_field);
		free(state);
		state = NULL;
	}
	memset(state->coefs, 0, sizeof(hvs_coefs));
	init_coefs(state->coefs);

	// We have moments and thus we can update vorticity field
	update_vorticity_field(state);
	(*sstate) = state;
	return HVS_OK;
}

int init_solver(hvs_params *params, hvs_state **sstate) {
#ifdef HVS_PROFILE
	time_t starttime, endtime;
	starttime = time(NULL);
#endif
	hvs_file* file;
	hvs_position *pos = NULL;
	hvs_vorticity *vort = NULL;
	int status = 0;
	UINT cursize = 0;
	hvs_state* state = NULL;
	hvs_moments moments = NULL;
	hvs_centers centers = NULL;
	if (params->initvortfile==NULL)
		if (params->initmomentsfile==NULL)
			return HVS_ERR_WRONG_USAGE;
	if ((params->initmomentsfile!=NULL) || (params->tmpmomentsfile!=NULL)) {
		if (params->tmpmomentsfile!=NULL) {
			// First try to init file handler of tmp moments file
			if (initfile(params->tmpmomentsfile, &file)) {
				// On failure use file handler of init moments file
				if (initfile(params->initmomentsfile, &file)) {
					return HVS_ERR;
				}
			} else {
				if (read_time(file, &params->tmpmomentstime)<1) {
					return HVS_ERR;
				}
			}
		} else {
			// Init file handler of init moments file
			if (initfile(params->initmomentsfile, &file)) {
				return HVS_ERR;
			}
		}
		// Read data
		do {
			cursize += status;
			if ((centers = (hvs_position *)realloc(centers, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_center)))==NULL) {
				return HVS_ERR;
			}
			if ((moments = realloc(moments, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_moment)))==NULL) {
				return HVS_ERR;
			}
		} while ((status = read_moments(file,HVS_READ_BLOCK_SIZE,
			&centers[cursize],
			&moments[cursize])) == HVS_READ_BLOCK_SIZE);
		if (status < 0) {
			return status;
		}
		cursize += status;
		closefile(&file);
		return init_solver_by_moments(params, cursize, centers, moments, sstate);
	}

	// Init file handler
	if (initfile(params->initvortfile, &file)) {
		return HVS_ERR;
	}
	
	// Read data
	do {
		cursize += status;
		if ((pos = (hvs_position *)realloc(pos, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_position)))==NULL) {
			return HVS_ERR;
		}
		if ((vort = (hvs_vorticity *)realloc(vort, (cursize+HVS_READ_BLOCK_SIZE)*sizeof(hvs_vorticity)))==NULL) {
			return HVS_ERR;
		}
	} while ((status = read_vorticity(file,HVS_READ_BLOCK_SIZE,
			&pos[cursize],
			&vort[cursize])) == HVS_READ_BLOCK_SIZE);
	if (status < 0) {
		return status;
	}
	cursize += status;
	closefile(&file);
	
	if ((state = (hvs_state *) malloc(sizeof(hvs_state))) == NULL) 
		return HVS_ERR;
	state->grid = pos;
	state->vorticity_field = vort;
	state->ncenters = 0;
	state->lambdasq = params->lambda0*params->lambda0;
	
	// Read initial centers positions
	if(params->initcentersfile!=NULL) {
		if((status=initfile(params->initcentersfile, &file))==HVS_OK) {
			do {
				state->ncenters += status;
				if ((centers = (hvs_position *)realloc(centers, (state->ncenters+HVS_READ_BLOCK_SIZE)*sizeof(hvs_center)))==NULL) {
					free(state->grid);
					free(state);
					state = NULL;
					return HVS_ERR;
				}
			} while ((status = read_centers(file,HVS_READ_BLOCK_SIZE,
					&centers[cursize])) == HVS_READ_BLOCK_SIZE);
			if (status < 0) {
				free(state->grid);
				free(state->vorticity_field);
				free(state);
				state = NULL;
				return status;
			}
			closefile(&file);
			state->centers = centers;
		} else {
			free(state->grid);
			free(state->vorticity_field);
			free(state);
			state = NULL;
			return status;
		}
	}
	// If there is no initcentersfile then centers coinside with grid
	if (state->ncenters==0) {
		if ((state->centers = malloc(sizeof(hvs_position)*cursize)) == NULL) {
			free(state->grid);
			free(state->vorticity_field);
			free(state);
			return HVS_ERR;
		}
		memcpy(state->centers,state->grid,sizeof(hvs_position)*cursize);
		state->ncenters = cursize;
	}
	
	// Allocate memory for moments
	if ((state->moments = malloc(sizeof(hvs_moment)*state->ncenters))==NULL) {
		free(state->centers);
		free(state->vorticity_field);
		free(state->grid);
		free(state);
		return HVS_ERR;
	}
	memset(state->moments,0,sizeof(hvs_moment)*state->ncenters);
	
	// Check if we have regular grid
	if (checkgrid(pos) != HVS_OK) {
		free(state->centers);
		free(state->moments);
		free(state->vorticity_field);
		free(state->grid);
		free(state);
		return HVS_ERR_IRREGULAR_GRID;
	}
	
	// Copy values to the state structure
	state->size = cursize;
	
	state->xmin = MIN(state->grid[0].x,state->grid[cursize-1].x);
	state->xmax = MAX(state->grid[0].x,state->grid[cursize-1].x);
	state->ymin = MIN(state->grid[0].y,state->grid[cursize-1].y);
	state->ymax = MAX(state->grid[0].y,state->grid[cursize-1].y);
	state->ystep = M_ABS(state->grid[1].y-state->grid[0].y);
	if (state->ystep>HVS_EPS) {
		state->sizey = round((state->ymax-state->ymin)/state->ystep)+1;
		state->xstep = M_ABS(pos[state->sizey].x-pos[0].x);
		state->sizex = round((state->xmax-state->xmin)/state->xstep)+1;
	} else {
		state->xstep = M_ABS(state->grid[1].x-state->grid[0].x);
		state->sizex = round((state->xmax-state->xmin)/state->xstep)+1;
		state->ystep = M_ABS(pos[state->sizex].y-pos[0].y);
		state->sizey = round((state->ymax-state->ymin)/state->ystep)+1;
	}
	
	if(cursize != (state->sizex*state->sizey)) {
		free(state->centers);
		free(state->moments);
		free(state->vorticity_field);
		free(state->grid);
		free(state);
		return HVS_ERR_WRONG_SIZE;
	}
	if ((status = init_moments(state))!=HVS_OK) {
		free(state->centers);
		free(state->moments);
		free(state->vorticity_field);
		free(state->grid);
		free(state);
		return status;
	}

	// Initialize velocity field
	state->velocity_field = (hvs_vector *) malloc(cursize*sizeof(hvs_vector));
	if (state->velocity_field == NULL) {
		free(state->centers);
		free(state->moments);
		free(state->vorticity_field);
		free(state->grid);
		free(state);
		state = NULL;
		return HVS_ERR;
	}
	memset(state->velocity_field, 0, cursize*sizeof(hvs_vector));

	// Initialize coefficients
	state->coefs = (hvs_coefs *) malloc(sizeof(hvs_coefs));
	if (state->coefs == NULL) {
		free(state->moments);
		free(state->centers);
		free(state->grid);
		free(state->vorticity_field);
		free(state->velocity_field);
		free(state);
		state = NULL;
	}
	memset(state->coefs, 0, sizeof(hvs_coefs));
	init_coefs(state->coefs);

	*sstate = state;
#ifdef HVS_PROFILE
	endtime = time(NULL);
	timings.init = (int)(endtime-starttime);
#endif
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
		if (state->coefs != NULL) {
			free(state->coefs);
			state->coefs = NULL;
		}
		free(state);
		(*sstate) = NULL;
	}
}

int run_solver(const hvs_params *params, hvs_state *state) {
#ifdef HVS_PROFILE
    time_t starttime, endtime;
    time_t startstep, endstep;
    starttime = time(NULL);
#endif
	unsigned int i;
	char status;
	unsigned int stepsnum = round((params->t1-params->t0)/params->timestep);
	FLOAT_TYPE tn;
	tn = params->t0;
	if (params->tmpmomentstime > 0)
		if (params->tmpmomentstime < params->t1) {
			if (params->tmpmomentstime > params->t0) {
				tn = params->tmpmomentstime;
				stepsnum = round((params->t1-tn)/params->timestep);
			}
		} else {
			// Skip steps until we get into the interval of previous stop
			return HVS_SKIP_STEP;
		}
	for(i=0;i<stepsnum;i++) {
#ifdef HVS_PROFILE
		int duration;
		startstep = time(NULL);
#endif
		if ((status = step_solver(state, &tn, params)) != HVS_OK) {
			return status;
		}
		if (params->tmpmomentsfile!=NULL) {
			write_tmp_moments(state, tn, params->tmpmomentsfile);
		}
#ifdef HVS_PROFILE
		endstep = time(NULL);
		duration = (int)(endstep-startstep);
		if (duration>5) {
			printf("Step %d: %d sec\n", i, duration);
		}
#endif
	}
	state->curtime = tn;
#ifdef HVS_PROFILE
	endtime = time(NULL);
	timings.integration = (int)(endtime-starttime);
#endif
	return update_vorticity_field(state);
}

int checkgrid(hvs_position *pos, hvs_vorticity *vort) {
	return HVS_OK;
}

