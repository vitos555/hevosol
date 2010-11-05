//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_MATRIXUTIL_H
#define HVS_MATRIXUTIL_H 1

#include "gentypes.h"

int matrix_multiply(const FLOAT_TYPE *A1, UINT m1, UINT n1, 
			const FLOAT_TYPE *A2, UINT m2, UINT n2, FLOAT_TYPE *result);
int matrix_add(const FLOAT_TYPE *A1, const FLOAT_TYPE *A2,
		UINT m, UINT n, FLOAT_TYPE *result);
int matrix_scalarmult(const FLOAT_TYPE scalar, const FLOAT_TYPE *m1,
			UINT m, UINT n, FLOAT_TYPE *result);
int matrix_inv(const FLOAT_TYPE *m, UINT n, FLOAT_TYPE *result);

#endif