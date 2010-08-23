#include <errno.h>
#include <string.h>

#define HVS_ERRORUTIL_H 1

// Status values
#define HVS_OK 0
#define HVS_ERR -1
#define HVS_ERR_WRONG_FILE_FORMAT -2
#define HVS_ERR_IRREGULAR_GRID -3

void hvserror(int status, const char *string);