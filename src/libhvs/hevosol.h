//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HEVOSOL_H
#define HEVOSOL_H 1

// Number of moments to include
#define NMOMENTS 2
#define NMOMENTS_1 (NMOMENTS+1)
#define NMOMENTS_2 (NMOMENTS_1*NMOMENTS_1)
#define NMOMENTS_3 (NMOMENTS_2*NMOMENTS_1)
#define NMOMENTS_4 (NMOMENTS_3*NMOMENTS_1)
#define NMOMENTS_5 (NMOMENTS_4*NMOMENTS_1)
#define NMOMENTS_6 (NMOMENTS_5*NMOMENTS_1)
#define NMOMENTS_7 (NMOMENTS_6*NMOMENTS_1)
#define NMOMENTS_8 (NMOMENTS_7*NMOMENTS_1)
#define NMOMENTS_9 (NMOMENTS_8*NMOMENTS_1)
#define NMOMENTS_10 (NMOMENTS_9*NMOMENTS_1)
// Number of combinations of moments
#define NCOMBS (NMOMENTS+1)*(NMOMENTS+2)/2
// Combination index:
// 0 -> (0,0)
// 1 -> (1,0)
// 2 -> (0,1)
// 3 -> (2,0)
// 4 -> (1,1)
// 5 -> (0,2)
// 6 -> (3,0)
// 7 -> (2,1)
// 8 -> (1,2)
// 9 -> (0,3)
#define COMBS_IND1(i) (i)==0?0:(i)==1?1:(i)==2?0:(i)==3?2:(i)==4?1:(i)==5?0:(i)==6?3:(i)==7?2:(i)==8?1:(i)==9?0:0
#define COMBS_IND2(i) (i)==0?0:(i)==1?0:(i)==2?1:(i)==3?0:(i)==4?1:(i)==5?2:(i)==6?0:(i)==7?1:(i)==8?2:(i)==9?3:0

#define MOMENTS_LEN NCOMBS

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define POWN1(a) ((a)%2==0?1:-1)
#define POW2(a) ((a)==0?1:2<<((a)-1))
#define POW(a,b) ((b)==0?1:((b)==1?(a):((b)==2?(a)*(a):((b)==3?(a)*(a)*(a):M_POW((a),(b))))))

#define COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j) ((k1)*NMOMENTS_7+(k2)*NMOMENTS_6+\
	(l1)*NMOMENTS_5+(l2)*NMOMENTS_4+\
	(m1)*NMOMENTS_3+(m2)*NMOMENTS_2+\
	(i)*NMOMENTS_1+(j))
#define MOM_INDEX(i,j) (((i)+(j))*((i)+(j)+1)/2+(j))

#define HVS_GMRES_PRECISION (FLOAT_TYPE)1e-10
#define HVS_GMRES_ITERATIONS 250

#include <unistd.h>

#include "gentypes.h"
#include "errorutil.h"

typedef struct s_hvs_position {
	FLOAT_TYPE x, y;
} hvs_position;

typedef struct s_hvs_vector {
	FLOAT_TYPE	xval, yval;
} hvs_vector;

typedef FLOAT_TYPE hvs_moment[MOMENTS_LEN];
typedef hvs_moment* hvs_moments;
typedef struct s_hvs_position hvs_center;
typedef hvs_center* hvs_centers;
typedef struct s_hvs_vector hvs_velocity;
typedef FLOAT_TYPE hvs_vorticity;

typedef struct {
	char	*initvortfile; /* Initial vorticity field file name */
	char	*initcentersfile; /* Initial centers positions file name */
	char	*initmomentsfile; /* Initial moments file name */
	FLOAT_TYPE	t0;
	FLOAT_TYPE	t1;
	FLOAT_TYPE	timestep;
	FLOAT_TYPE	lambda0;
	FLOAT_TYPE	nu;
	FLOAT_TYPE	xmin,xmax,xstep;
	FLOAT_TYPE	ymin,ymax,ystep;
} hvs_params;

typedef struct {
	FLOAT_TYPE gamma1[NMOMENTS_8];
	FLOAT_TYPE gamma2[NMOMENTS_8];
} hvs_coefs;

typedef struct {
	UINT	size, sizex, sizey, ncenters;
	FLOAT_TYPE	xmin, xmax, xstep;
	FLOAT_TYPE	ymin, ymax, ystep;
	FLOAT_TYPE	curtime;
	FLOAT_TYPE	lambdasq;
	hvs_position*	grid;
	hvs_velocity*	velocity_field;
	hvs_vorticity*	vorticity_field;
	hvs_moments	moments;
	hvs_centers	centers;
	hvs_coefs*	coefs;
} hvs_state;

typedef struct {
	hvs_moments	moments;
	hvs_centers	centers;
	FLOAT_TYPE	lambdasq;
	UINT		ncenters;
} hvs_ode_data;

int init_solver(const hvs_params *params, hvs_state **sstate);
int init_solver_by_moments(const hvs_params *params, UINT ncenters, const hvs_centers centers, const hvs_moments moments,
				hvs_state **sstate);
void free_solver(hvs_state **sstate);
int run_solver(const hvs_params *params, hvs_state *state);

#endif