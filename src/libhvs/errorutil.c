#include "hevosol.h"
#include <stdio.h>

void hvserror(int status, const char *string) {
	switch(status) {
		case HVS_ERR:
			perror(string);
			break;
		case HVS_ERR_WRONG_FILE_FORMAT:
			fprintf(stderr, "%s: Wrong file format.\n", string);
			break;
		case HVS_ERR_IRREGULAR_GRID:
			fprintf(stderr, "%s: This version doesn't support irregular grids for initial data.\n", string);
			break;
	}
}