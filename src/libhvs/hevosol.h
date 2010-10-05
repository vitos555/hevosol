//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#define HEVOSOL_H 1

// Number of moments to include
#define NMOMENTS 2
#define NMOMENTS_2 4
#define NMOMENTS_3 8
#define NMOMENTS_4 16
#define NMOMENTS_5 32
#define NMOMENTS_6 64
#define NMOMENTS_7 128
#define NMOMENTS_8 256
#define NMOMENTS_9 512
#define NMOMENTS_10 1024

#define MOMENTS_LEN 4

#define MIN(a,b) ((a)>(b)?(b):(a))
#define POWN1(a) ((a)%2==0?1:-1)

#define COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j) ((k1)*NMOMENTS_7+(k2)*NMOMENTS_6+\
	(l1)*NMOMENTS_5+(l2)*NMOMENTS_4+\
	(m1)*NMOMENTS_3+(m2)*NMOMENTS_2+\
	(i)*NMOMENTS+(j))
#define MOM_INDEX(i,j) ((i)*((i)-1)/2+(j))

#include <unistd.h>

#include "gentypes.h"

#ifndef HVS_ERRORUTIL_H
#include "errorutil.h"
#endif

typedef struct s_hvs_position {
	FLOAT_TYPE x, y;
} hvs_position;

typedef struct s_hvs_vector {
	FLOAT_TYPE	xval, yval;
} hvs_vector;

typedef FLOAT_TYPE hvs_moment[MOMENTS_LEN];
typedef struct s_hvs_position hvs_center;
typedef struct s_hvs_vector hvs_velocity;
typedef FLOAT_TYPE hvs_vorticity;

typedef struct {
	char	*initvortfile; /* Initial vorticity field file name */
	FLOAT_TYPE	t0;
	FLOAT_TYPE	t1;
	FLOAT_TYPE	timestep;
	FLOAT_TYPE	lambda0;
	FLOAT_TYPE	nu;
} hvs_params;

typedef struct {
	FLOAT_TYPE gamma1a[NMOMENTS_8];
	FLOAT_TYPE gamma2a[NMOMENTS_8];
	FLOAT_TYPE gamma1b[NMOMENTS_8];
	FLOAT_TYPE gamma2b[NMOMENTS_8];
} hvs_coefs;

typedef struct {
	UINT	size, sizex, sizey, ncenters;
	FLOAT_TYPE	xmin, xmax, xstep;
	FLOAT_TYPE	ymin, ymax, ystep;
	hvs_position*	grid;
	hvs_velocity*	velocity_field;
	hvs_vorticity*	vorticity_field;
	hvs_moment*	moments;
	hvs_center*	centers;
	hvs_coefs*	coefs;
} hvs_state;

typedef struct {
	hvs_moment	*moments;
	hvs_position	*centers;
	FLOAT_TYPE	lambdasq;
	UINT		ncenters;
} hvs_ode_data;

int init_solver(const hvs_params *params, hvs_state **sstate);
int init_solver_by_moments(hvs_params *params, UINT ncenters, const hvs_center *centers, const hvs_moment *moments,
				FLOAT_TYPE xmin, FLOAT_TYPE xmax, FLOAT_TYPE xstep, 
				FLOAT_TYPE ymin, FLOAT_TYPE ymax, FLOAT_TYPE ystep, hvs_state **sstate);
void free_solver(hvs_state **sstate);
int run_solver(const hvs_params *params, hvs_state *state);
ssize_t write_data(const hvs_state *state, const char *filename);
