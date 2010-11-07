#include "gmres.h"
#include "vectorutil.h"
#include "hevosol.h"
#include <math.h>
#include <stdlib.h>

int gmres(const FLOAT_TYPE *A, const FLOAT_TYPE *x0, const FLOAT_TYPE *b, UINT size,
		FLOAT_TYPE accuracy, UINT restart, FLOAT_TYPE *result) {
	FLOAT_TYPE beta;
	FLOAT_TYPE residualsq;
	FLOAT_TYPE *betae = NULL,*yn = NULL,*bnorm = NULL, *tmp = NULL;
	FLOAT_TYPE *H = NULL,*Q = NULL;
	
	int status;
	UINT n,i,globaln=0;
	
	if ((bnorm = malloc(sizeof(FLOAT_TYPE)*size))==NULL) {
		return HVS_ERR;
	}
	
	if ((tmp = malloc(sizeof(FLOAT_TYPE)*size))==NULL) {
		free(bnorm);
		return HVS_ERR;
	}

	// Normalize b
	if ((status=vect_scalarmult(1/M_SQRT(vect_normsq(b,size)),b,size,bnorm))!=HVS_OK) {
		return status;
	}
	
	// Initialize result with x0, residualsq to be bigger then accuracy
	memcpy(result,x0,sizeof(FLOAT_TYPE)*size);
	residualsq = 2*accuracy*accuracy;

	while((residualsq>(accuracy*accuracy)) && (globaln<MAX_GMRES_RESTARTS)) {
		globaln++;
		if((status = vect_matrixmultsub(A,result,bnorm,size,tmp))!=HVS_OK) {
			return status;
		}
		beta = M_SQRT(vect_normsq(tmp,size));
		n = 0;
		
		while((residualsq>(accuracy*accuracy)) && (n<restart)) {
			n++;
			if (globaln==1) {
				if ((betae = realloc(betae,sizeof(FLOAT_TYPE)*(n+1)))==NULL) {
					free(bnorm);
					return HVS_ERR;
				}
				if ((yn = realloc(yn,sizeof(FLOAT_TYPE)*n))==NULL) {
					free(betae);
					free(bnorm);
					free(tmp);
					return HVS_ERR;
				}
				if ((H = realloc(H,sizeof(FLOAT_TYPE)*(n+1)*n))==NULL) {
					free(betae);
					free(bnorm);
					free(tmp);
					free(yn);
					return HVS_ERR;
				}
				if ((Q = realloc(Q,sizeof(FLOAT_TYPE)*size*n))==NULL) {
					free(betae);
					free(bnorm);
					free(tmp);
					free(yn);
					free(H);
					return HVS_ERR;
				}
			}
			memset(H,0,sizeof(FLOAT_TYPE)*(n+1)*n);
			memset(betae,0,sizeof(FLOAT_TYPE)*(n+1));
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
	}
	free(H);
	free(Q);
	free(yn);
	free(betae);
	free(tmp);
	free(bnorm);
	return HVS_OK;
}

int ols(const FLOAT_TYPE *X, const FLOAT_TYPE *y, UINT m, UINT n, FLOAT_TYPE *result) {
	FLOAT_TYPE inv_factor=0.0;
	FLOAT_TYPE *xx,*inv;
	UINT i,j,k;
	int status;
	
	for(i=0;i<n;i++) {
		for(j=0;j<m;j++)
			printf("%Lf\t",X[i*n+j]);
		printf("\n");
	}
	printf("\n");

	if ((inv=malloc(sizeof(FLOAT_TYPE)*n*n))==NULL) {
		return HVS_ERR;
	}
	if ((xx=malloc(sizeof(FLOAT_TYPE)*n*n))==NULL) {
		free(inv);
		return HVS_ERR;
	}
	
	// Build X'X
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			xx[i*n+j]=0;
			for (k=0;k<m;k++) 
				xx[i*n+j] = X[i*n+k]*X[j*n+k];
		}
		result[i] = 0;
	}
	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++)
			printf("%Lf\t",xx[i*n+j]);
		printf("\n");
	}
	// Find inverse
	if((status=matrix_inv(xx,n,inv))!=HVS_OK) {
		return status;
	}
	
			printf("Here\n");
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			for (k=0;k<n;k++)
				result[i] += inv[k*n+i]*X[j*n+k]*y[j];
	free(inv);
	free(xx);
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