#include "gmres.h"
#include "vectorutil.h"
#include "hevosol.h"
#include <math.h>
#include <stdlib.h>

int gmres(const FLOAT_TYPE *A, const FLOAT_TYPE *x0, const FLOAT_TYPE *b, UINT size,
		FLOAT_TYPE accuracy, FLOAT_TYPE *result) {
	FLOAT_TYPE beta;
	FLOAT_TYPE residualsq;
	FLOAT_TYPE *betae = NULL,*yn = NULL,*bnorm = NULL;
	FLOAT_TYPE *H = NULL,*Q = NULL;
	
	int status;
	UINT n,i;
	
	if ((bnorm = malloc(sizeof(FLOAT_TYPE)*size))==NULL) {
		return HVS_ERR;
	}

	// Normalize b
	if ((status=vect_scalarmult(1/M_SQRT(vect_normsq(b,size)),b,size,bnorm))!=HVS_OK) {
		return status;
	}
	
	if((status = vect_matrixmultsub(A,x0,bnorm,size,result))!=HVS_OK) {
		return status;
	}
	beta = M_SQRT(vect_normsq(result,size));

	// Initialize result with x0, residualsq to be bigger then accuracy
	memcpy(result,x0,sizeof(FLOAT_TYPE)*size);
	residualsq = 2*accuracy*accuracy;
	n = 0;

	while((residualsq>(accuracy*accuracy)) && (n<MAX_GMRES_ITERATIONS)) {
		n++;
		if ((betae = realloc(betae,sizeof(FLOAT_TYPE)*(n+1)))==NULL) {
			free(bnorm);
			return HVS_ERR;
		}
		if ((yn = realloc(yn,sizeof(FLOAT_TYPE)*n))==NULL) {
			free(betae);
			free(bnorm);
			return HVS_ERR;
		}
		if ((H = realloc(H,sizeof(FLOAT_TYPE)*(n+1)*n))==NULL) {
			free(betae);
			free(bnorm);
			free(yn);
			return HVS_ERR;
		}
		if ((Q = realloc(Q,sizeof(FLOAT_TYPE)*size*n))==NULL) {
			free(betae);
			free(bnorm);
			free(yn);
			free(H);
			return HVS_ERR;
		}
		memset(H,sizeof(FLOAT_TYPE)*(n+1)*n,0);
		memset(betae,sizeof(FLOAT_TYPE)*(n+1),0);
		betae[0]=beta;
		
		// One step Arnoldi
		if ((status = arnoldi(A,b,size,n,H,Q))!=HVS_OK) {
			return status;
		}
		
		// Now OLS
		if ((status = ols(H,betae,n+1,n,yn))!=HVS_OK) {
			return status;
		}
		
		// Get rn into result
		for(i=0;i<n+1;i++) {
			result[i] = vect_dotproduct(&H[i*n],yn,n)-betae[i];
		}
		residualsq = vect_normsq(result,n+1);
		printf("Residual: %Lf\n",residualsq);
	}
	
	// Get xn
	for(i=0;i<size;i++) {
		result[i]=M_SQRT(vect_normsq(b,size))*vect_dotproduct(&Q[i*n],yn,n);
	}
	free(H);
	free(Q);
	free(yn);
	free(betae);
	free(bnorm);
	return HVS_OK;
}

int ols(const FLOAT_TYPE *X, const FLOAT_TYPE *y, UINT m, UINT n, FLOAT_TYPE *result) {
	FLOAT_TYPE inv_factor=0.0;
	UINT i,j;
	int status;
	for(i=0;i<n;i++) {
		for(j=0;j<m;j++) inv_factor += vect_normsq(&X[i*n],n);
		result[i] = 0;
	}
	inv_factor = 1/inv_factor;
	for(i=0;i<m;i++) 
		for(j=0;j<n;j++) result[i] += X[i*n+j]*inv_factor*y[j];
	return HVS_OK;
}

int arnoldi(const FLOAT_TYPE *A, const FLOAT_TYPE *q1, UINT m, UINT n, FLOAT_TYPE *H,
		FLOAT_TYPE *Q) {
	UINT k,j,i;
	FLOAT_TYPE *qk,*qkp;
	int status;
	if ((qk = malloc(sizeof(FLOAT_TYPE)*m)) ==NULL) {
		return HVS_ERR;
	}
	if ((qkp = malloc(sizeof(FLOAT_TYPE)*m)) ==NULL) {
		free(qk);
		return HVS_ERR;
	}
	memcpy(qkp,q1,sizeof(FLOAT_TYPE)*m);
	for (k=1;k<n+1;k++) {
		for(i=0;i<m;i++) Q[i*n+k-1]=qkp[i];
		if ((status=vect_matrixmult(A,qkp,m,qk))!=HVS_OK) {
			return status;
		}
		for(j=0;j<k-1;j++) {
			for(i=0;i<m;i++) qkp[i]=Q[i*n+j];
			H[j*n+k-1]=vect_dotproduct(qk,qkp,m);
			if ((status=vect_scalarmultsub(H[j*n+k-1],qkp,qk,m,qk))!=HVS_OK) {
				return status;
			}
		}
		H[k*n+k-1]=M_SQRT(vect_normsq(qk,m));
		if((status=vect_scalarmult(-1/H[k*n+k-1],qk,m,qkp))!=HVS_OK) {
			return status;
		}
	}
	free(qk);
	free(qkp);
	return HVS_OK;
}