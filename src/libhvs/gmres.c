#include "gmres.h"
#include "vectorutil.h"
#include "hevosol.h"
#include <math.h>
#include <stdlib.h>

int gmres(const FLOAT_TYPE *A, const FLOAT_TYPE *x0, const FLOAT_TYPE *b, UINT size,
		FLOAT_TYPE accuracy, UINT restart, FLOAT_TYPE *result) {
	FLOAT_TYPE beta,bnormval;
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
	bnormval = M_SQRT(vect_normsq(b,size));
	if ((status=vect_scalarmult((FLOAT_TYPE)1/bnormval,b,size,bnorm))!=HVS_OK) {
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
				if ((H = realloc(H,sizeof(FLOAT_TYPE)*n*(n+1)))==NULL) {
					free(betae);
					free(bnorm);
					free(tmp);
					free(yn);
					return HVS_ERR;
				}
				if ((Q = realloc(Q,sizeof(FLOAT_TYPE)*size*(n+1)))==NULL) {
					free(betae);
					free(bnorm);
					free(tmp);
					free(yn);
					free(H);
					return HVS_ERR;
				}
			}
			if (n>1) {
				memset(&H[n*(n-1)],0,sizeof(FLOAT_TYPE)*n*2);
				betae[n-1]=0;
				for(i=1;i<=n;i++) {
					memcpy(betae,&H[(n-i)*(n-1)],sizeof(FLOAT_TYPE)*(n-1));
					memcpy(&H[(n-i)*n],betae,sizeof(FLOAT_TYPE)*n);
				}
			}
			if (n==1) {
				memcpy(Q,bnorm,sizeof(FLOAT_TYPE)*size);
			}
			memset(betae,0,sizeof(FLOAT_TYPE)*(n+1));
			betae[0]=beta;

			// One step Arnoldi
			if ((status = arnoldi(A,size,n,H,Q))!=HVS_OK) {
				return status;
			}
			// Now OLS
			if ((status = ols(H,betae,n+1,n,yn))!=HVS_OK) {
				return status;
			}
			
			// Get rn into result
			for(i=0;i<n+1;i++) {
				result[i] = vect_dotproduct(&H[i*n],yn,n);
				if (i==0) result[i] -= beta;
			}
			residualsq = vect_normsq(result,n+1);
			printf("Residual: %Lf\n",residualsq);
		}
	
		// Get xn
		for(i=0;i<size;i++) {
			int j;
			result[i] = (FLOAT_TYPE)0.0;
			for (j=0;j<n;j++) {
				result[i]+=Q[j*size+i]*yn[j];
			}
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
				xx[i*n+j] += X[k*n+i]*X[k*n+j];
		}
		result[i] = 0;
	}
	
	// Find inverse
	if((status=matrix_inv(xx,n,inv))!=HVS_OK) {
		return status;
	}
	
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			for (k=0;k<n;k++)
				result[i] += inv[k*n+i]*X[j*n+k]*y[j];
	free(inv);
	free(xx);
	return HVS_OK;
}

int arnoldi(const FLOAT_TYPE *A, UINT m, UINT n, FLOAT_TYPE *H,
		FLOAT_TYPE *Q) {
	UINT k,j,i;
	int status;

	if ((status=vect_matrixmult(A,&Q[(n-1)*m],m,&Q[n*m]))!=HVS_OK) {
		return status;
	}
	for(j=0;j<n;j++) {
		H[j*n+n-1]=vect_dotproduct(&Q[j*m],&Q[n*m],m);
		if ((status=vect_scalarmultsub(H[j*n+n-1],&Q[j*m],&Q[n*m],m,&Q[n*m]))!=HVS_OK) {
			return status;
		}
	}
	H[n*n+n-1]=M_SQRT(vect_normsq(&Q[n*m],m));
	if (H[n*n+n-1]>=0.0000001) {
		if((status=vect_scalarmult(-1/H[n*n+n-1],&Q[n*m],m,&Q[n*m]))!=HVS_OK) {
			return status;
		}
	}
	printf("Q:\n");
	for(i=0;i<n+1;i++) {
		for(k=0;k<m;k++)
		    printf("%Lf\t",Q[i*m+k]);
		printf("\n");
	}
	printf("H:\n");
	for(i=0;i<n+1;i++) {
		for(k=0;k<n;k++)
		    printf("%Lf\t",H[i*n+k]);
		printf("\n");
	}
	return HVS_OK;
}