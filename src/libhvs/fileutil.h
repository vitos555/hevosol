#define HVS_READ_BLOCK_SIZE 100

#ifndef HEVOSOL_H
#include "hevosol.h"
#endif
#include <stdio.h>

typedef struct {
	FILE *fh;
} hvs_file;

int initfile(const char *filename, hvs_file **file);

ssize_t readdata(const hvs_file *filename, size_t count, hvs_position *pos, hvs_vorticity *vort);

void closefile(hvs_file **file);

ssize_t writedata(const hvs_state *state, const char *filename);