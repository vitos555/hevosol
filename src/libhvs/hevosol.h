//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#define HEVOSOL_H 1

// Custom types. To make further scaling easier.
#define FLOAT_TYPE float
#define UINT unsigned int

// Number of moments to include
#define NMOMENTS 2

#include <unistd.h>

#ifndef HVS_ERRORUTIL_H
#include "errorutil.h"
#endif

typedef struct s_hvs_position {
	FLOAT_TYPE x,y;
} hvs_position;

typedef struct s_hvs_vector {
	FLOAT_TYPE	xval, yval;
} hvs_vector;

typedef FLOAT_TYPE hvs_moment[NMOMENTS*(NMOMENTS-1)/2];
typedef struct s_hvs_position hvs_center;
typedef struct s_hvs_vector hvs_velocity;
typedef FLOAT_TYPE hvs_vorticity;

typedef struct {
	char	*initvortfile; /* Initial vorticity field file name */
	FLOAT_TYPE	t0;
	FLOAT_TYPE	t1;
	FLOAT_TYPE	timestep;
	FLOAT_TYPE	gridwidth, gridheight;
} hvs_params;

typedef struct {
	UINT	size, sizex, sizey, ncenters;
	FLOAT_TYPE	xmin, xmax, xstep;
	FLOAT_TYPE	ymin, ymax, ystep;
	hvs_position*	grid;
	hvs_velocity*	velocity_field;
	hvs_vorticity*	vorticity_field;
	hvs_moment*	moments;
	hvs_center*	centers;
} hvs_state;

int init_solver(const hvs_params *params, hvs_state **sstate);
int init_solver_by_moments(UINT ncenters, const hvs_center *centers, const hvs_moment *moments,
				FLOAT_TYPE xmin, FLOAT_TYPE xmax, FLOAT_TYPE xstep, 
				FLOAT_TYPE ymin, FLOAT_TYPE ymax, FLOAT_TYPE ystep, hvs_state **sstate);
void free_solver(hvs_state **sstate);
int run_solver(const hvs_params *params, hvs_state *state);
ssize_t write_data(const hvs_state *state, const char *filename);
