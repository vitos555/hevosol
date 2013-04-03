//
// Copyright (C) 2010-2011, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//
#include <omp.h>
#if HVS_DEBUG>2
#include <stdio.h>
#endif

#ifdef HVS_PROFILE
#include <time.h>
#endif

int init_ode_data(hvs_ode_data *data, hvs_state *input) {
	data->ncenters = input->ncenters;
	data->lambdasq = input->lambdasq;
	data->moments = (hvs_moment *) malloc(input->ncenters*sizeof(hvs_moment));
	if (data->moments != NULL) {
		memset(data->moments, 0, input->ncenters*sizeof(hvs_moment));
	} else {
		return HVS_ERR;
	}
	data->centers = (hvs_center *) malloc(input->ncenters*sizeof(hvs_center));
	if (data->centers != NULL) {
		memset(data->centers, 0, input->ncenters*sizeof(hvs_center));
	} else {
		free(data->moments);
		return HVS_ERR;
	}
	return HVS_OK;
}

int free_ode_data(hvs_ode_data *data) {
	if (data->moments != NULL)
	    free(data->moments);
	data->moments = NULL;
	if (data->centers != NULL)
	    free(data->centers);
	data->centers = NULL;
	return HVS_OK;
}

int update_vorticity_field(hvs_state *state) {
	int i,j;
	FLOAT_TYPE sum;
#ifdef HVS_PROFILE
    time_t starttime, endtime;
    starttime = time(NULL);
#endif

#if HVS_DEBUG
	printf("update_vorticity_field\n");
#endif

	for (i=0; i<state->size; i++) {
		sum = 0.0;
		#pragma omp parallel for reduction(+:sum) private(j) shared(state)
		for (j=0; j<state->ncenters; j++) {
			sum += 	state->moments[j][MOM_INDEX(0,0)]*he(state->grid[i].x-state->centers[j].x,state->grid[i].y-state->centers[j].y,state->lambdasq,0,0)+
				state->moments[j][MOM_INDEX(1,1)]*he(state->grid[i].x-state->centers[j].x,state->grid[i].y-state->centers[j].y,state->lambdasq,1,1)+
				state->moments[j][MOM_INDEX(0,2)]*he(state->grid[i].x-state->centers[j].x,state->grid[i].y-state->centers[j].y,state->lambdasq,0,2)+
				state->moments[j][MOM_INDEX(2,0)]*he(state->grid[i].x-state->centers[j].x,state->grid[i].y-state->centers[j].y,state->lambdasq,2,0);
		}
		state->vorticity_field[i] = sum;
	}
#ifdef HVS_PROFILE
    endtime = time(NULL);
    timings.vorticity_update += (int)(endtime-starttime);
    timings.vorticity_update = (int)(timings.vorticity_update/2);
#endif
	return HVS_OK;
}

int init_moments(hvs_state *state) {
	FLOAT_TYPE *A,*x;
	UINT i,j,k;
	int status;

#ifdef HVS_PROFILE
    time_t starttime, endtime;
    starttime = time(NULL);
#endif

#if HVS_DEBUG
	printf("init_moments\n");
#endif

	if (state->size!=state->ncenters) {
		return HVS_ERR;
	}
	if ((A=malloc(sizeof(FLOAT_TYPE)*state->size*state->size))==NULL) {
		return HVS_ERR;
	}
	if ((x=malloc(sizeof(FLOAT_TYPE)*state->size))==NULL) {
		free(A);
		return HVS_ERR;
	}
	memset(x,0,sizeof(FLOAT_TYPE)*state->size);
	#pragma omp parallel for collapse(2) private(i,j) shared(A,state)
	for(i=0;i<state->size;i++)
		for(j=0;j<state->size;j++)
				A[state->size*i+j] = 
					he(state->grid[i].x-state->centers[j].x,
					   state->grid[i].y-state->centers[j].y,
					   state->lambdasq,
					   0,0);
	if ((status=gmres(A,x,state->vorticity_field,state->size,HVS_GMRES_PRECISION,HVS_GMRES_MAX_INNER_MATRIX,x))!=HVS_OK) {
		free(A);
		free(x);
		return status;
	}
	// Get moments
	for(i=0;i<state->ncenters;i++) {
		state->moments[i][0]=x[i];
	}
	free(A);
	free(x);
#ifdef HVS_PROFILE
    endtime = time(NULL);
    timings.init_moments = (int)(endtime-starttime);
#endif
	return HVS_OK;
}

int eval_eq(hvs_ode_data *input, hvs_ode_data *output, hvs_coefs *coefs, FLOAT_TYPE time_val) {
	FLOAT_TYPE lambdasq = input->lambdasq;
	FLOAT_TYPE gamma1,gamma2;
	int i0,j0,k,k1,k2,l,l1,l2,m,m1,m2,i,j;
#ifdef HVS_PROFILE
    time_t starttime, endtime;
    starttime = time(NULL);
#endif
	output->lambdasq = input->lambdasq;

#if HVS_DEBUG
	printf("eval_eq\n");
#endif

	// Parallelize execution
	#pragma omp parallel default(none) shared(input,output,coefs,time_val,lambdasq)
	{
	// First evaluate centers equation
	#pragma omp for private(i0,j0,l,m,l1,l2,m1,m2)
	for (i0=0; i0<input->ncenters; i0++) {
		output->centers[i0].x=(FLOAT_TYPE)0.0;
		output->centers[i0].y=(FLOAT_TYPE)0.0;
		for(j0=0; j0<input->ncenters; j0++) {
		if (i0==j0) continue;
		for(l=0;l<NCOMBS;l++) {
			l1=COMBS_IND1(l);
			l2=COMBS_IND2(l);
			for(m=0;m<NCOMBS;m++) {
				m1=COMBS_IND1(m);
				m2=COMBS_IND2(m);
				output->centers[i0].x += 
				input->moments[j0][MOM_INDEX(l1,l2)]*
				input->moments[i0][MOM_INDEX(m1,m2)]*
				POWN1(m1+m2)*
				hb1(
				input->centers[i0].x-input->centers[j0].x,
				input->centers[i0].y-input->centers[j0].y,
				input->lambdasq,
				m1+l1,m2+l2);

				output->centers[i0].y += 
				input->moments[j0][MOM_INDEX(l1,l2)]*
				input->moments[i0][MOM_INDEX(m1,m2)]*
				POWN1(m1+m2)*
				hb2(
				input->centers[i0].x-input->centers[j0].x,
				input->centers[i0].y-input->centers[j0].y,
				input->lambdasq,
				m1+l1,m2+l2);
			}
		}
		}
		output->centers[i0].x/=input->moments[i0][MOM_INDEX(0,0)];
		output->centers[i0].y/=input->moments[i0][MOM_INDEX(0,0)];
	}
	
	// Now get the moments equations
	#pragma omp for collapse(2)  private(i0,j0,l,m,k,k1,k2,l1,l2,m1,m2,gamma1,gamma2,i,j) 
	for (i0=0; i0<input->ncenters; i0++) {
		for(j0=0; j0<input->ncenters; j0++)
			for(k=0;k<NCOMBS;k++) {
				k1=COMBS_IND1(k);
				k2=COMBS_IND2(k);
				for(l=0;l<NCOMBS;l++) {
					l1=COMBS_IND1(l);
					l2=COMBS_IND2(l);
					for(m=0;m<NCOMBS;m++) {
						m1=COMBS_IND1(m);
						m2=COMBS_IND2(m);
						gamma1=0.0;
						gamma2=0.0;
						for (i=0;i<=NMOMENTS;i++)
						for (j=0;j<=NMOMENTS;j++) {
							if (j0==i0) {
// If j'==j gammas will come from A
							if ( (i<=MIN(l1,k1-1)) && (j<=MIN(l2,k2)) )
gamma1+=coefs->gamma1[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	POW(lambdasq,i+j+1)*h1(m1+k1-i-1+l1-i,m2+k2-j+l2-j,lambdasq);
							if ( (i<=MIN(l1,k1)) && (j<=MIN(l2,k2-1)) )
gamma2+=coefs->gamma2[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	POW(lambdasq,i+j+1)*h2(m1+k1-i+l1-i,m2+k2-j-1+l2-j,lambdasq);
							} else {
// If j'<>j gammas will come from B 
								if ((i<=MIN(l1,k1-1)) && (j<=MIN(l2,k2)))
gamma1+=coefs->gamma1[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	POW(lambdasq,i+j+1)*hb1(input->centers[i0].x-input->centers[j0].x,
				input->centers[i0].y-input->centers[j0].y,
				lambdasq,m1+k1-i-1+l1-i,m2+k2-j+l2-j);
								if ((i<=MIN(l1,k1)) && (j<=MIN(l2,k2-1)))
gamma2+=coefs->gamma2[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	POW(lambdasq,i+j+1)*hb2(input->centers[i0].x-input->centers[j0].x,
				input->centers[i0].y-input->centers[j0].y,
				lambdasq,m1+k1-i+l1-i,m2+k2-j-1+l2-j);
							}
						}
					// Add A+B
		output->moments[i0][MOM_INDEX(k1,k2)] += (gamma1+gamma2)*
			input->moments[i0][MOM_INDEX(l1,l2)]*
			input->moments[j0][MOM_INDEX(m1,m2)]*
			POWN1(k1+k2)*POW(input->lambdasq,k1+k2)/
			(POW2(k1+k2)*factorial(k1)*factorial(k2));
					}
				}
				// Add C once
				if (j0==0) {
					if (k1>0)
						output->moments[i0][MOM_INDEX(k1,k2)]+=
							output->centers[i0].x*
							input->moments[i0][MOM_INDEX(k1-1,k2)];
					if (k2>0)
						output->moments[i0][MOM_INDEX(k1,k2)]+=
							output->centers[i0].y*
							input->moments[i0][MOM_INDEX(k1,k2-1)];
				}
			}
	} // end for
	} // end pragma omp parallel
#ifdef HVS_PROFILE
    endtime = time(NULL);
    timings.eval_equation += (int)(endtime-starttime);
    timings.eval_equation = (int)(timings.eval_equation/2);
#endif
	return HVS_OK;
}

int rk4_hvs_solve(hvs_state *curdata, FLOAT_TYPE tn, FLOAT_TYPE timestep, FLOAT_TYPE nu) {
	hvs_ode_data k1, k2, k3, k4, kt;
	int status = HVS_OK;
	int i0,i,i1,i2;

#ifdef HVS_PROFILE
    time_t starttime, endtime;
    starttime = time(NULL);
#endif

#if HVS_DEBUG
	printf("rk4_hvs_solve\n");
#endif

	// Initialize temporary variable for moments and centers
	if ((status=init_ode_data(&kt,curdata))!=HVS_OK) {
		return status;
	}
	memcpy(kt.moments,curdata->moments,sizeof(hvs_moment)*curdata->ncenters);
	memcpy(kt.centers,curdata->centers,sizeof(hvs_center)*curdata->ncenters);

	// Initialize k1,k2,k3,k4 - rk4 function evaluations
	if ((status=init_ode_data(&k1,curdata))!=HVS_OK) {
		return status;
	}
	if ((status=init_ode_data(&k2,curdata))!=HVS_OK) {
		return status;
	}
	if ((status=init_ode_data(&k3,curdata))!=HVS_OK) {
		return status;
	}
	if ((status=init_ode_data(&k4,curdata))!=HVS_OK) {
		return status;
	}
	
	// Get k1
	if ((status=eval_eq(&kt, &k1, curdata->coefs, tn))!=HVS_OK) {
		return status;
	}
	
#if HVS_DEBUG>2
	printf("m00=%f,m20=%f,m11=%f,m02=%f\n",k1.moments[0][0],k1.moments[0][3],k1.moments[0][4],k1.moments[0][5]);
	printf("m00=%f,m20=%f,m11=%f,m02=%f\n",k1.moments[1][0],k1.moments[1][3],k1.moments[1][4],k1.moments[1][5]);
	printf("y1=%.12f,y2=%.1f\n", k1.centers[0].y, k1.centers[1].y);
#endif

	// Update kt
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i=0; i<NCOMBS; i++) {
			i1=COMBS_IND1(i);
			i2=COMBS_IND2(i);
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k1.moments[i0][MOM_INDEX(i1,i2)]*0.5*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k1.centers[i0].x*0.5*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k1.centers[i0].y*0.5*timestep;
	}
	kt.lambdasq += 4*0.5*timestep*nu;
	
	// Get k2
	if ((status=eval_eq(&kt, &k2, curdata->coefs, tn+0.5*timestep))!=HVS_OK) {
		return status;
	}
	
	// Update kt
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i=0; i<NCOMBS; i++) {
			i1=COMBS_IND1(i);
			i2=COMBS_IND2(i);
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k2.moments[i0][MOM_INDEX(i1,i2)]*0.5*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k2.centers[i0].x*0.5*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k2.centers[i0].y*0.5*timestep;
	}
	
	// Get k3
	if ((status=eval_eq(&kt, &k3, curdata->coefs, tn+0.5*timestep))!=HVS_OK) {
		return status;
	}
	
	// Update kt
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i=0; i<NCOMBS; i++) {
			i1=COMBS_IND1(i);
			i2=COMBS_IND2(i);
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k3.moments[i0][MOM_INDEX(i1,i2)]*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k3.centers[i0].x*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k3.centers[i0].y*timestep;
	}
	kt.lambdasq += 4*0.5*timestep*nu;
	
	// Get k4
	if ((status=eval_eq(&kt, &k4, curdata->coefs, tn+timestep))!=HVS_OK) {
		return status;
	}
	
	
	// Use k1,k2,k3,k4 to find rk4 value
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i=0; i<NCOMBS; i++) {
			i1=COMBS_IND1(i);
			i2=COMBS_IND2(i);
			curdata->moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+1.0/6*timestep*
				(k1.moments[i0][MOM_INDEX(i1,i2)]+2*k2.moments[i0][MOM_INDEX(i1,i2)]+2*k3.moments[i0][MOM_INDEX(i1,i2)]+k4.moments[i0][MOM_INDEX(i1,i2)]);
		}
		curdata->centers[i0].x = curdata->centers[i0].x+1.0/6*timestep*(k1.centers[i0].x+2*k2.centers[i0].x+2*k3.centers[i0].x+k4.centers[i0].x);
		curdata->centers[i0].y = curdata->centers[i0].y+1.0/6*timestep*(k1.centers[i0].y+2*k2.centers[i0].y+2*k3.centers[i0].y+k4.centers[i0].y);
	}
	curdata->lambdasq += 4.0*timestep*nu;

	// Free all the Intermediate values
	free_ode_data(&kt);
	free_ode_data(&k1);
	free_ode_data(&k2);
	free_ode_data(&k3);
	free_ode_data(&k4);
	
#if HVS_DEBUG>2
	printf("m00=%f,m20=%f,m11=%f,m02=%f\n",curdata->moments[0][0],curdata->moments[0][3],curdata->moments[0][4],curdata->moments[0][5]);
	printf("m00=%f,m20=%f,m11=%f,m02=%f\n",curdata->moments[1][0],curdata->moments[1][3],curdata->moments[1][4],curdata->moments[1][5]);
	printf("y1=%f,y2=%f\n", curdata->centers[0].y, curdata->centers[1].y);
#endif

#ifdef HVS_PROFILE
    endtime = time(NULL);
    timings.rk_step += (int)(endtime-starttime);
    timings.rk_step = (int)(timings.rk_step/2);
#endif

	return HVS_OK;
}

int init_coefs(hvs_coefs *coefs) {
	int k,k1,k2,m,m1,m2,l,l1,l2,i,j;
#ifdef HVS_PROFILE
    time_t starttime, endtime;
    starttime = time(NULL);
#endif

#if HVS_DEBUG
	printf("init_coefs\n");
#endif
	#pragma omp parallel for default(none) shared(coefs) private(k,k1,k2,l,l1,l2,m,m1,m2,i,j)
	for(k=0;k<NCOMBS;k++) {
		k1=COMBS_IND1(k);
		k2=COMBS_IND2(k);
		for(l=0;l<NCOMBS;l++) {
			l1=COMBS_IND1(l);
			l2=COMBS_IND2(l);
			for(m=0;m<NCOMBS;m++) {
				m1=COMBS_IND1(m);
				m2=COMBS_IND2(m);
				for(i=0;i<=NMOMENTS;i++)
				for(j=0;j<=NMOMENTS;j++) {
					if ((i<=MIN(l1,k1-1))&&(j<=MIN(l2,k2))) {
coefs->gamma1[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=(FLOAT_TYPE)POWN1(l1+l2)*
	binomial(l1,i)*binomial(l2,j)*(POW2(i+1)*factorial(k1))/factorial(k1-i-1)*
	(POW2(j)*factorial(k2))/factorial(k2-j);
					}
					if ((i<=MIN(l1,k1))&&(j<=MIN(l2,k2-1))) {
coefs->gamma2[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=(FLOAT_TYPE)POWN1(l1+l2)*
	binomial(l1,i)*binomial(l2,j)*(POW2(i)*factorial(k1))/factorial(k1-i)*
	(POW2(j+1)*factorial(k2))/factorial(k2-j-1);
					}
				}
			}
		}
	}
#ifdef HVS_PROFILE
    endtime = time(NULL);
    timings.init_coefs = (int)(endtime-starttime);
#endif
}


int step_solver(hvs_state *state, FLOAT_TYPE *tn, const hvs_params *params) {
	int i,status;

#if HVS_DEBUG
	printf("step_solver, t=%f\n",(float)(*tn));
#endif

	if ((status = rk4_hvs_solve(state, (*tn), params->timestep, params->nu))!=HVS_OK) {
		return status;
	}
	(*tn) = (*tn)+params->timestep;

#if HVS_DEBUG>2
	printf("%f=%f\n",state->lambdasq,(*tn)*4.0*params->nu+params->lambda0*params->lambda0);
#endif

	return HVS_OK;
}

