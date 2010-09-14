//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

int step_solver(hvs_state *state, FLOAT_TYPE *tn, FLOAT_TYPE timestep) {
	int i;
	hvs_ode_data oderetval = {state->moments, state->centers, state->ncenters, state->lambda0+4*nu*(*tn)};
	oderetval = rk4_hvs_solve(oderetval, state->coefs, (*tn), timestep);
	(*tn) = (*tn)+timestep;
	return HVS_OK;
}

int update_vorticity_field(hvs_state *state, hvs_params *params) {
	int i,j;
	FLOAT_TYPE sum;
	FLOAT_TYPE lambda_sq = params->lambda0*params->lambda0;
	for (i=0; i<state->size; i++) {
		sum = 0.0;
		for (j=0; j<state->ncenters; j++) {
			sum += 	state->moments[j][0]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,0)+
				state->moments[j][1]*he(state->grid[i].x,state->grid[i].y,lambda_sq,0,2)+
				state->moments[j][2]*he(state->grid[i].x,state->grid[i].y,lambda_sq,2,0);
		}
		state->vorticity_field[i] = sum;
	}
	return HVS_OK;
}

hvs_ode_data eval_eq(hvs_ode_data input, hvs_coefs *coefs, FLOAT_TYPE time) {
	hvs_moment_center_combined output;
	int i0,j0,k1,k2,l1,l2,m1,m2,i,j;
	for (i0=0; i0<input->ncenters; i0++) {
		output[i0].center.x=0;
		output[i0].center.y=0;
		for(j0=0; j0<input->ncenters; j0++)
		for(l1=0;l1<NMOMENTS;l1++)
		for(l2=0;l2<NMOMENTS;l2++)
			for(m1=0;m1<NMOMENTS;m1++)
			for(m2=0;m2<NMOMENTS;m2++) {
				output[i0].center.x += 
				input[j0].moment[MOM_INDEX(l1,l2)]*
				input[i0].moment[MOM_INDEX(m1,m2)]*
				POWN1(m1+m2)*
				hb1(
				input[i0].center.x-input[j0].center.x,
				input[i0].center.y-input[j0].center.y,
				input.lambdasq,
				m1+l1,m2+l2);

				output[i0].center.y += 
				input[j0].moment[MOM_INDEX(l1,l2)]*
				input[i0].moment[MOM_INDEX(m1,m2)]*
				POWN1(m1+m2)*
				hb2(
				input[i0].center.x-input[j0].center.x,
				input[i0].center.y-input[j0].center.y,
				input.lambdasq,
				m1+l1,m2+l2);
			}
		output[i0].center.x/=input[i0].moments[0];
		output[i0].center.y/=input[i0].moments[0];
	}
	for (i0=0; i0<input->ncenters; i0++) {
		memset(output[i0].moment,0,MOMENT_LEN*sizeof(hvs_moment));
		for(j0=0; j0<input->ncenters; j0++)
			for(k1=0;k1<NMOMENTS;k1++)
			for(k2=0;k2<NMOMENTS;k2++) {
				for(l1=0;l1<NMOMENTS;l1++)
				for(l2=0;l2<NMOMENTS;l2++)
				for(m1=0;m1<NMOMENTS;m1++)
				for(m2=0;m2<NMOMENTS;m2++) {
					FLOAT_TYPE gamma1=0.0,gamma2=0.0;
					for (i=0;i<NMOMENTS;i++)
					for (j=0;j<NMOMENTS;j++) {
						if (j0==i0) {
// If j'==j gammas will come from A
gamma1+=coefs->gamma1a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	pow(lambdasq,i+j)*h1(m1+k1-i-1+l1-i,m2+k2-j+l2-j,lambdasq);
gamma2+=coefs->gamma2a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	pow(lambdasq,i+j)*h2(m1+k1-i+l1-i,m2+k2-j-1+l2-j,lambdasq);
						} else {
// If j'<>j gammas will come from B 
gamma1+=coefs->gamma1a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	pow(lambdasq,i+j+1)*hb1(input[i0].center.x-input[j0].center.x,
				input[i0].center.y-input[j0].center.y,
				lambdasq,m1+k1-i-1+l1-i,m2+k2-j+l2-j);
gamma2+=coefs->gamma2a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]/
	pow(lambdasq,i+j+1)*hb2(input[i0].center.x-input[j0].center.x,
				input[i0].center.y-input[j0].center.y,
				lambdasq,m1+k1-i+l1-i,m2+k2-j-1+l2-j);
						}
					}
					output[i0].moment[MOM_INDEX(k1,k2)] = (gamma1+gamma2)*
						input[i0].moment[MOM_INDEX(k1,k2)]*
						input[j0].moment[MOM_INDEX(k1,k2)];
				}
				// Add C once
				if (j0==0) {
					output[i0].moment[MOM_INDEX(k1,k2)]*=
					POWN1(k1+k2)*pow(input->lambdasq,k1+k2)/
					((2<<(k1+k2-1))*factorial(k1)*factorial(k2));
					output[i0].moment[MOM_INDEX(k1,k2)]+=
					output[i0].center.x*output[i0].momets[MOM_INDEX(k1-1,k2)]+
					output[i0].center.y*output[i0].momets[MOM_INDEX(k1,k2-1)];
				}
			}
	}
}

hvs_ode_data rk4_hvs_solve(hvs_ode_data curdata, hvs_coefs *coefs, FLOAT_TYPE tn, FLOAT_TYPE timestep) {
	hvs_ode_data retval, k1, k2, k3, k4, kt;
	int i,i0;
	k1=eval_eq(curdata, coefs, tn);
	for (i0=0; i0<curdata->ncenters; i0++) {
		for (i=0; i<MOMENTS_LEN; i++) {
			kt[i0].moment[i] = curdata[i0].moment[i]+k1[i0].moment[i]*0.5*timestep;
		}
		kt[i0].center.x = curdata[i0].center.x+k1[i0].center.x*0.5*timestep;
		kt.center.y = curdata.center.y+k1.center.y*0.5*timestep;
		k2=eval_eq(kt, coefs, tn+0.5*timestep);
		for (i=0; i<MOMENTS_LEN; i++) {
			kt.moment[i] = curdata.moment[i]+k2.moment[i]*0.5*timestep;
		}
		kt.center.x = curdata.center.x+k2.center.x*0.5*timestep;
		kt.center.y = curdata.center.y+k2.center.y*0.5*timestep;
		k3=eval_eq(kt, coefs, tn+0.5*timestep);
		for (i=0; i<MOMENTS_LEN; i++) {
			kt.moment[i] = curdata.moment[i]+k3.moment[i]*timestep;
		}
		kt.center.x = curdata.center.x+k3.center.x*timestep;
		kt.center.y = curdata.center.y+k3.center.y*timestep;
		k4=eval_eq(kt, coefs, tn+timestep);
		for (i=0; i<MOMENTS_LEN; i++) {
			retval.moment[i] = curdata.moment[i]+1/6*timestep*
				(k1.moment[i]+2*k2.moment[i]+2*k3.moment[i]+k4.moment[i]);
		}
		kt.center.x = curdata.center.x+1/6*timestep*(k1.center.x+2*k2.center.x+2*k3.center.x+k4.center.x);
		kt.center.y = curdata.center.y+1/6*timestep*(k1.center.y+2*k2.center.y+2*k3.center.y+k4.center.y);
	}
	return retval;
}


int init_coefs(hvs_coefs *coefs) {
	int k1,k2,m1,m2,l1,l2,i,j;
	for(k1=0;k1<NMOMENTS;k1++)
	for(k2=0;k2<NMOMENTS;k2++)
		for(l1=0;l1<NMOMENTS;l1++)
		for(l2=0;l2<NMOMENTS;l2++)
			for(m1=0;m1<NMOMENTS;m1++)
			for(m2=0;m2<NMOMENTS;m2++)
				for(i=0;i<NMOMENTS;i++)
				for(j=0;j<NMOMENTS;j++) {
					if ((i<MIN(l1,k1-1))&&(j<MIN(l2,k2))) {
coefs->gamma1a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=POWN1(l1+l2)*
	binomial(l1,i)*binomial(l2,j)*((2<<(i+1-1))*factorial(k1))/factorial(k1-i-1)*
	((2<<(j-1))*factorial(k2))/factorial(k2-j);
coefs->gamma1b[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=POWN1(l1+l2)*
	binomial(l1,i)*binomial(l2,j)*((2<<(i+1-1))*factorial(k1))/factorial(k1-i-1)*
	((2<<(j-1))*factorial(k2))/factorial(k2-j);
					}
					if ((i<MIN(l1,k1))&&(j<MIN(l2,k2-1))) {
coefs->gamma2a[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=POWN1(l1+l1)*
	binomial(l1,i)*binomial(l2,j)*((2<<(i-1))*factorial(k1))/factorial(k1-i)*
	((2<<(j+1-1))*factorial(k2))/factorial(k2-j-1);
coefs->gamma2b[COEF_INDEX(k1,k2,l1,l2,m1,m2,i,j)]=POWN1(l1+l1)*
	binomial(l1,i)*binomial(l2,j)*((2<<(i-1))*factorial(k1))/factorial(k1-i)*
	((2<<(j+1-1))*factorial(k2))/factorial(k2-j-1);
					}
				}
}