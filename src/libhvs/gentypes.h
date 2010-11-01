//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_GENTYPES_H
#define HVS_GENTYPES_H 1

// Custom types. To make further scaling easier.
#define HVS_FLOAT 1
#define HVS_DOUBLE 2
#define HVS_LONG_DOUBLE 3

// Change this constant to change float type
#ifndef HVS_FLOAT_TYPE
#define HVS_FLOAT_TYPE HVS_LONG_DOUBLE
#endif

// Initialize float type
#if HVS_FLOAT_TYPE==HVS_FLOAT
#define FLOAT_TYPE float
#elif HVS_FLOAT_TYPE==HVS_DOUBLE
#define FLOAT_TYPE double
#elif HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
#define FLOAT_TYPE long double
#else
#define HVS_FLOAT_TYPE HVS_DOUBLE
#define FLOAT_TYPE double
#endif

// Initialize unsigned int
#define UINT unsigned int

// Define math type dependant macros
#if HVS_FLOAT_TYPE==HVS_FLOAT
#define M_EXP(x) expf((x))
#define M_POW(x,y) powf((x),(y))
#define M_SQRT(x) sqrtf((x))
#elif HVS_FLOAT_TYPE==HVS_LONG_DOUBLE
#define M_EXP(x) expl((x))
#define M_POW(x,y) powl((x),(y))
#define M_SQRT(x) sqrtl((x))
#else
#define M_EXP(x) exp((x))
#define M_POW(x,y) pow((x),(y))
#define M_SQRT(x) sqrt((x))
#endif

#endif