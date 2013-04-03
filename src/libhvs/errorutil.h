//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#include <errno.h>
#include <string.h>

#ifndef HVS_ERRORUTIL_H
#define HVS_ERRORUTIL_H 1

// Status values
#define HVS_OK 0
#define HVS_SKIP_STEP 1
#define HVS_ERR -1
#define HVS_ERR_WRONG_FILE_FORMAT -2
#define HVS_ERR_IRREGULAR_GRID -3
#define HVS_ERR_WRONG_USAGE -4
#define HVS_ERR_WRONG_SIZE -5
#define HVS_ERR_GMRES_GOT_NAN -6
#define HVS_ERR_GMRES_MAX_RESTART -7
#define HVS_ERR_NOT_ENOUGH_DATA -8
#define HVS_ERR_FIRST_MOMENT_ZERO -9

void hvserror(int status, const char *string);

#endif