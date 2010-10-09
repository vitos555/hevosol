//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

int init_ode_data(hvs_ode_data *data, hvs_ode_data *input) {
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
	free(data->moments);
	free(data->centers);
	return HVS_OK;
}

int update_vorticity_field(hvs_state *state, const hvs_params *params) {
	int i,j;
	FLOAT_TYPE sum;
	FLOAT_TYPE lambda_sq = params->lambda0*params->lambda0;
	for (i=0; i<state->size; i++) {
		sum = 0.0;
		for (j=0; j<state->ncenters; j++) {
			sum += 	state->moments[j][MOM_INDEX(0,0)]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,0)+
				state->moments[j][MOM_INDEX(0,2)]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,2)+
				state->moments[j][MOM_INDEX(2,0)]*he(state->grid[i].x,state->grid[i].y,lambda_sq,2,0);
		}
		state->vorticity_field[i] = sum;
	}
	return HVS_OK;
}

int eval_eq(hvs_ode_data *input, hvs_ode_data *output, hvs_coefs *coefs, FLOAT_TYPE time) {
	FLOAT_TYPE lambdasq = input->lambdasq;
	FLOAT_TYPE gamma1,gamma2;
	int i0,j0,k1,k2,l1,l2,m1,m2,i,j;
// Debug code
/*	FLOAT_TYPE b1,b2;
	k1=1;k2=1;l1=1;l2=1;m1=1;m2=1;
	i0=0;j0=1;
	gamma1 = 0.0;
	for(i=0;i<NMOMENTS;i++)
	for(j=0;j<NMOMENTS;j++)
		if ((i<=MIN(l1,k1)) && (j<=MIN(l2,k2-1))) {
		FLOAT_TYPE test;
		lambdasq=1.2*1.2;
		b1=0.6;
		b2=-0.39;
		printf("i=%d,j=%d\n",i,j);
		printf("COEF_INDEX=%d\n",COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j));
		printf("lambdasq=%f\n",lambdasq);
		printf("gamma1=%f\n",coefs->gamma1[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]);
		printf("gamma2=%f\n",coefs->gamma2[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]);
		printf("hb1=%f\n",hb2(b1,b2,lambdasq,m1+k1-i+l1-i,m2+k2-j-1+l2-j));
		printf("POW=%f\n",1/POW(lambdasq,i+j+1));
		
	    test=coefs->gamma2[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
		POW(lambdasq,i+j+1)*hb2(b1,b2,lambdasq,m1+k1-i+l1-i,m2+k2-j-1+l2-j);
		printf("test=%.14f\n",test);
		gamma1+=test;
	    printf("gamma_test=%.10f\n",gamma1);
	}
	printf("gamma_test=%.10f\n",gamma1);
*/
	for (i0=0; i0<input->ncenters; i0++) {
		output->centers[i0].x=0;
		output->centers[i0].y=0;
		for(j0=0; j0<input->ncenters; j0++)
		for(l1=0;l1<NMOMENTS;l1++)
		for(l2=0;l2<NMOMENTS;l2++)
			for(m1=0;m1<NMOMENTS;m1++)
			for(m2=0;m2<NMOMENTS;m2++) {
				output->centers[i0].x += 
				input->moments[j0][MOM_INDEX(l1,l2)]*
				input->moments[j0][MOM_INDEX(m1,m2)]*
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
		output->centers[i0].x/=input->moments[i0][0];
		output->centers[i0].y/=input->moments[i0][0];
	}
	for (i0=0; i0<input->ncenters; i0++) {
		memset(output->moments[i0],0,MOMENTS_LEN*sizeof(hvs_moment));
		for(j0=0; j0<input->ncenters; j0++)
			for(k1=0;k1<NMOMENTS;k1++)
			for(k2=0;k2<NMOMENTS;k2++) {
				for(l1=0;l1<NMOMENTS;l1++)
				for(l2=0;l2<NMOMENTS;l2++)
				for(m1=0;m1<NMOMENTS;m1++)
				for(m2=0;m2<NMOMENTS;m2++) {
					gamma1=0.0;
					gamma2=0.0;
					for (i=0;i<NMOMENTS;i++)
					for (j=0;j<NMOMENTS;j++) {
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
					output->moments[i0][MOM_INDEX(k1,k2)] = (gamma1+gamma2)*
						input->moments[i0][MOM_INDEX(k1,k2)]*
						input->moments[j0][MOM_INDEX(k1,k2)];
				}
				// Add C once
				if (j0==0) {
					output->moments[i0][MOM_INDEX(k1,k2)]*=
					POWN1(k1+k2)*POW(input->lambdasq,k1+k2)/
					((2<<(k1+k2-1))*factorial(k1)*factorial(k2));
					output->moments[i0][MOM_INDEX(k1,k2)]+=
					output->centers[i0].x*output->moments[i0][MOM_INDEX(k1-1,k2)]+
					output->centers[i0].y*output->moments[i0][MOM_INDEX(k1,k2-1)];
				}
			}
	}
}

int rk4_hvs_solve(hvs_ode_data *curdata, hvs_coefs *coefs, FLOAT_TYPE tn, FLOAT_TYPE timestep, FLOAT_TYPE nu) {
	hvs_ode_data k1, k2, k3, k4, kt;
	int status;
	int i0,i1,i2;

	if ((status=init_ode_data(&kt,curdata))!=HVS_OK) {
		return status;
	}
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

	if ((status=eval_eq(curdata, &k1, coefs, tn))!=HVS_OK) {
		return status;
	}
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i1=0; i1<NMOMENTS; i1++) 
		for (i2=0; i2<NMOMENTS; i2++) {
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k1.moments[i0][MOM_INDEX(i1,i2)]*0.5*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k1.centers[i0].x*0.5*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k1.centers[i0].y*0.5*timestep;
		kt.lambdasq += 4*0.5*timestep*nu;
		if ((status=eval_eq(&kt, &k2, coefs, tn))!=HVS_OK) {
			return status;
		}
		for (i1=0; i1<NMOMENTS; i1++) 
		for (i2=0; i2<NMOMENTS; i2++) {
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k2.moments[i0][MOM_INDEX(i1,i2)]*0.5*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k2.centers[i0].x*0.5*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k2.centers[i0].y*0.5*timestep;
		kt.lambdasq += 4*0.5*timestep*nu;
		if ((status=eval_eq(&kt, &k3, coefs, tn))!=HVS_OK) {
			return status;
		}
		for (i1=0; i1<NMOMENTS; i1++) 
		for (i2=0; i2<NMOMENTS; i2++) {
			kt.moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+k3.moments[i0][MOM_INDEX(i1,i2)]*timestep;
		}
		kt.centers[i0].x = curdata->centers[i0].x+k3.centers[i0].x*timestep;
		kt.centers[i0].y = curdata->centers[i0].y+k3.centers[i0].y*timestep;
		kt.lambdasq += 4*timestep*nu;
		if ((status=eval_eq(&kt, &k4, coefs, tn))!=HVS_OK) {
			return status;
		}
		for (i1=0; i1<NMOMENTS; i1++) 
		for (i2=0; i2<NMOMENTS; i2++) {
			curdata->moments[i0][MOM_INDEX(i1,i2)] = curdata->moments[i0][MOM_INDEX(i1,i2)]+1/6*timestep*
				(k1.moments[i0][MOM_INDEX(i1,i2)]+2*k2.moments[i0][MOM_INDEX(i1,i2)]+2*k3.moments[i0][MOM_INDEX(i1,i2)]+k4.moments[i0][MOM_INDEX(i1,i2)]);
		}
		curdata->centers[i0].x = curdata->centers[i0].x+1/6*timestep*(k1.centers[i0].x+2*k2.centers[i0].x+2*k3.centers[i0].x+k4.centers[i0].x);
		curdata->centers[i0].y = curdata->centers[i0].y+1/6*timestep*(k1.centers[i0].y+2*k2.centers[i0].y+2*k3.centers[i0].y+k4.centers[i0].y);
	}
	free_ode_data(&kt);
	free_ode_data(&k1);
	free_ode_data(&k2);
	free_ode_data(&k3);
	free_ode_data(&k4);
	return HVS_OK;
}

int init_coefs(hvs_coefs *coefs) {
	int k1,k2,m1,m2,l1,l2,i,j;
	for(k1=0;k1<NMOMENTS;k1++)
	for(k2=0;k2<NMOMENTS;k2++)
		for(l1=0;l1<NMOMENTS;l1++)
		for(l2=0;l2<NMOMENTS;l2++)
			for(m1=0;m1<NMOMENTS;m1++)
			for(m2=0;m2<NMOMENTS;m2++) {
				for(i=0;i<NMOMENTS;i++)
				for(j=0;j<NMOMENTS;j++) {
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


int step_solver(hvs_state *state, FLOAT_TYPE *tn, const hvs_params *params) {
	int i,status;
	hvs_ode_data oderetval = {state->moments, state->centers, params->lambda0*params->lambda0+4*params->nu*(*tn), state->ncenters};
	if ((status = rk4_hvs_solve(&oderetval, state->coefs, (*tn), params->timestep, params->nu))!=HVS_OK) {
		return status;
	}
	(*tn) = (*tn)+params->timestep;
	return HVS_OK;
}

