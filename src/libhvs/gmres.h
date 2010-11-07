//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_GMRES_H
#define HVS_GMRES_H 1

#ifndef MAX_GMRES_RESTARTS
#define MAX_GMRES_RESTARTS 10
#endif

#include "gentypes.h"

int gmres(const FLOAT_TYPE *A, const FLOAT_TYPE *x0, const FLOAT_TYPE *b, UINT size, 
		FLOAT_TYPE accuracy, const UINT restart, FLOAT_TYPE *result);

#endif