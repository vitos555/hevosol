//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

// Custom types. To make further scaling easier.
#define FLOAT_TYPE float
#define UINT unsigned int

// Status values
#define HVS_OK 0
#define HVS_ERR 1

typedef struct {
	FLOAT_TYPE x,y;
} hvs_position;

typedef struct {
	FLOAT_TYPE	xval, yval;
} hvs_vector;

typedef FLOAT_TYPE hvs_moment;

typedef struct {
	char	*initvortfile; /* Initial vorticity field file name */
	FLOAT_TYPE	t0;
	FLOAT_TYPE	t1;
	FLOAT_TYPE	timestep;
	FLOAT_TYPE	gridwidth, gridheight;
} hvs_params;

typedef struct {
	UINT	width, height;
	hvs_position*	grid;
	hvs_vector*	velocity_field;
	hvs_vector*	vorticity_field;
	hvs_moment*	moments;
} hvs_state;

int init_solver(const hvs_params *params, hvs_state **sstate);
void free_solver(hvs_state **sstate);
int run_solver(const hvs_params *params, hvs_state *state);
int format_output(hvs_state *state, size_t size, char *output);
