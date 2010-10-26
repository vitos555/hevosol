//
// Copyright (C) 2010, Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
// Author: Vitalii Ostrovskyi <vitalii@ostrovskyi.org.ua>
//

#ifndef HVS_FILEUTIL_H
#define HVS_FILEUTIL_H 1

#define HVS_READ_BLOCK_SIZE 100

#include <stdio.h>
#include "hevosol.h"

typedef struct {
	FILE *fh;
} hvs_file;

int initfile(const char *filename, hvs_file **file);
void closefile(hvs_file **file);

ssize_t read_vorticity(const hvs_file *filename, size_t count, hvs_position *pos, hvs_vorticity *vort);
ssize_t write_vorticity(const hvs_state *state, const char *filename);
ssize_t append_vorticity(const hvs_state *state, const char *filename);
ssize_t write_params(const hvs_params *params, const char *filename);
ssize_t append_params(const hvs_params *params, const char *filename);
ssize_t write_moments(const hvs_state *state, const char *filename);
ssize_t append_moments(const hvs_state *state, const char *filename);
ssize_t write_centers(const hvs_state *state, const char *filename);
ssize_t append_centers(const hvs_state *state, const char *filename);

#endif