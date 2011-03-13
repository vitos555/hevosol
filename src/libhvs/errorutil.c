#include "hevosol.h"
#include <stdio.h>
#include <stdlib.h>

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
		case HVS_ERR_WRONG_SIZE:
			fprintf(stderr, "%s: Given data is not rectangular.\n", string);
			break;
		case HVS_ERR_GOT_NAN:
			fprintf(stderr, "%s: Got NaN.\n", string);
			break;
		case HVS_ERR_GMRES_MAX_RESTART:
			fprintf(stderr, "%s: GMRES: Got max restarts. Try decreasing lambda.\n", string);
			break;
	}
}

void hvsdie(const char *string) {
    perror(string);
    exit(1);
}