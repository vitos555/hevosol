//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_HERMITEUTIL_H
#define HVS_HERMITEUTIL_H 1

#include "gentypes.h"

FLOAT_TYPE he(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambda, unsigned short k1, unsigned short k2);
FLOAT_TYPE hb1(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambdasq, UINT k1, UINT k2);
FLOAT_TYPE hb2(FLOAT_TYPE x1, FLOAT_TYPE x2, FLOAT_TYPE lambdasq, UINT k1, UINT k2);
FLOAT_TYPE h1(UINT alpha1,UINT alpha2,FLOAT_TYPE lambdasq);
FLOAT_TYPE h2(UINT alpha1,UINT alpha2,FLOAT_TYPE lambdasq);

#endif