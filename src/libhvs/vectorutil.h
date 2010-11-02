//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_VECTORUTIL_H
#define HVS_VECTORUTIL_H 1

#include "gentypes.h"

FLOAT_TYPE vect_dotproduct(const FLOAT_TYPE *x, const FLOAT_TYPE *y, UINT n);
int vect_scalarmult(const FLOAT_TYPE scalar, const FLOAT_TYPE *vector, UINT n,
			FLOAT_TYPE *result);
int vect_scalarmultadd(const FLOAT_TYPE scalar, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n,
			FLOAT_TYPE *result);
int vect_scalarmultsub(const FLOAT_TYPE scalar, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n,
			FLOAT_TYPE *result);
FLOAT_TYPE vect_normsq(const FLOAT_TYPE *vector,UINT n);
int vect_matrixmult(const FLOAT_TYPE *matrix, const FLOAT_TYPE *vector, UINT n, 
			FLOAT_TYPE *result);
int vect_matrixmultadd(const FLOAT_TYPE *matrix, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n, 
			FLOAT_TYPE *result);
int vect_matrixmultsub(const FLOAT_TYPE *matrix, const FLOAT_TYPE *x, 
			const FLOAT_TYPE *y, UINT n, 
			FLOAT_TYPE *result);

#endif