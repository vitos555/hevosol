#include "matrixutil.h"
#include "errorutil.h"
#include <stdlib.h>
#include <math.h>

int matrix_multiply(const FLOAT_TYPE *A1, UINT m1, UINT n1, 
			const FLOAT_TYPE *A2, UINT m2, UINT n2, FLOAT_TYPE *result) {
	UINT i,j,k;
	if (n1!=m2) return HVS_ERR;
	for(i=0;i<m1;i++)
		for(j=0;j<n2;j++) {
			result[i*n2+j]=0;
			for(k=0;k<n1;k++)
				result[i*n2+j]+=A1[i*n1+k]*A2[k*n2+j];
		}
	return HVS_OK;
}

int matrix_add(const FLOAT_TYPE *A1, const FLOAT_TYPE *A2,
		UINT m, UINT n, FLOAT_TYPE *result) {
	UINT i,j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			result[i*n+j]=A1[i*n+j]+A2[i*n+j];
	return HVS_OK;
}

int matrix_scalarmult(const FLOAT_TYPE scalar, const FLOAT_TYPE *m1,
			UINT m, UINT n, FLOAT_TYPE *result) {
	UINT i,j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			result[i*n+j]=scalar*m1[i*n+j];
	return HVS_OK;
}

int matrix_inv(const FLOAT_TYPE *m, UINT n, FLOAT_TYPE *result) {
	int i,j,k;
	UINT *l1;
	short int *l2;
	FLOAT_TYPE *m1;
	
	if ((l1=malloc(n*sizeof(UINT)))==NULL) {
		return HVS_ERR;
	}
	memset(l1,0,n*sizeof(UINT));
	if ((l2=malloc(n*sizeof(short int)))==NULL) {
		free(l1);
		return HVS_ERR;
	}
	memset(l2,0,n*sizeof(short int));
	if ((m1=malloc(n*n*sizeof(FLOAT_TYPE)))==NULL) {
		free(l1);
		free(l2);
		return HVS_ERR;
	}
	memcpy(m1,m,n*n*sizeof(FLOAT_TYPE));
	
	// Initialize identity
	memset(result,0,sizeof(FLOAT_TYPE)*n*n);
	for(i=0;i<n;i++) result[i*n+i]=1;

	// Do Gaussian elimination with partial scaling
	for(i=0;i<n;i++) {
		FLOAT_TYPE max=(FLOAT_TYPE)0.0;
		for(j=0;j<n;j++) {
			if (l2[j]) continue;
			if (M_ABS(m1[j*n+i])>max) {
				max=m1[j*n+i];
				l1[i]=j;
			}
		}
		if (M_ABS(max)>0) {
			for (j=0;j<n;j++) {
				if (l2[j]) continue;
				if ((j==l1[j]) && (m1[j*n+i]==(FLOAT_TYPE)1.0)) {
					l2[j]=1;
					continue;
				}
				if (m1[j*n+i]==(FLOAT_TYPE)0.0) continue;
				FLOAT_TYPE first = m1[j*n+i];
				for(k=0;k<n;k++) {
					if (j==l1[i]) {
						if (k>=i) m1[j*n+k]=m1[j*n+k]/max;
						result[j*n+k]=result[j*n+k]/max;
					} else {
						if (k>=i) 
							m1[j*n+k]=
							m1[j*n+k]-first*m1[l1[i]*n+k]/max;
						result[j*n+k]=
							result[j*n+k]-first*result[l1[i]*n+k]/max;
					}
				}
				if (j==l1[i]) {
					max=(FLOAT_TYPE)1.0;
					l2[j]=1;
				}
			}
		} else {
			return HVS_ERR;
		}
	}

	// Backward elimination
	for (i=n-1;i>0;i--) {
		for(j=0;j<n;j++) {
			if (j==l1[i]) continue;
			FLOAT_TYPE last = m1[j*n+i];
			for (k=0;k<n;k++) {
			if (k>=i)
				m1[j*n+k]=m1[j*n+k]-m1[l1[i]*n+k]*last;
				result[j*n+k]=result[j*n+k]-result[l1[i]*n+k]*last;
			}
		}
	}
	memcpy(m1,result,n*n*sizeof(FLOAT_TYPE));
	for (i=0;i<n;i++)
		memcpy(&result[i*n],&m1[l1[i]*n],n*sizeof(FLOAT_TYPE));
	free(l1);
	free(l2);
	free(m1);
	return HVS_OK;
}

int matrix_transpose(const FLOAT_TYPE *A, UINT m, UINT n, FLOAT_TYPE *result) {
	UINT i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++) 
			result[j*m+i]=A[i*n+j];
	return HVS_OK;
}
