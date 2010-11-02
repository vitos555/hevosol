//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_GMRES_H
#define HVS_GMRES_H 1

#ifndef MAX_GMRES_ITERATIONS
#define MAX_GMRES_ITERATIONS 1000
#endif

#include "gentypes.h"

int gmres(const FLOAT_TYPE *A, const FLOAT_TYPE *x0, const FLOAT_TYPE *b, UINT size, 
		FLOAT_TYPE accuracy, FLOAT_TYPE *result);

#endif