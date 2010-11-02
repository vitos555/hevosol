#include "vectorutil.h"
#include "errorutil.h"

FLOAT_TYPE vect_dotproduct(const FLOAT_TYPE *x, const FLOAT_TYPE *y, UINT n) {
	FLOAT_TYPE retval=0;
	UINT i;
	for(i=0;i<n;i++) {
		retval += x[i]*y[i];
	}
	return retval;
}

int vect_scalarmult(const FLOAT_TYPE scalar, const FLOAT_TYPE *vector, UINT n, 
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++) {
		result[i] = vector[i]*scalar;
	}
	return HVS_OK;
}

// scalar*x+y
int vect_scalarmultadd(const FLOAT_TYPE scalar, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n,
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++) {
		result[i] = scalar*x[i]+y[i];
	}
	return HVS_OK;
}

// scalar*x-y
int vect_scalarmultsub(const FLOAT_TYPE scalar, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n,
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++) {
		result[i] = scalar*x[i]-y[i];
	}
	return HVS_OK;
}

FLOAT_TYPE vect_normsq(const FLOAT_TYPE *vector,UINT n) {
	return vect_dotproduct(vector,vector,n);
}

int vect_matrixmult(const FLOAT_TYPE *matrix, const FLOAT_TYPE *vector, UINT n, 
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++)
		result[i] = vect_dotproduct(&matrix[i*n],vector,n);
	return HVS_OK;
}

// Ax+y
int vect_matrixmultadd(const FLOAT_TYPE *matrix, const FLOAT_TYPE *x,
			const FLOAT_TYPE *y, UINT n, 
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++)
		result[i] = vect_dotproduct(&matrix[i*n],x,n)+y[i];
	return HVS_OK;
}

// Ax-y
int vect_matrixmultsub(const FLOAT_TYPE *matrix, const FLOAT_TYPE *x,
			const FLOAT_TYPE *y, UINT n, 
			FLOAT_TYPE *result) {
	UINT i;
	for(i=0;i<n;i++)
		result[i] = vect_dotproduct(&matrix[i*n],x,n)-y[i];
	return HVS_OK;
}
